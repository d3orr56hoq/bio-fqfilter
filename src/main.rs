use std::env;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::process::ExitCode;

use fqfilter::window::trim_ends;
use fqfilter::{classify, fastq, parse_min_mean, Config, Decision, Reader, Stats};

const USAGE: &str = "\
usage: fqfilter [-q MIN_MEAN] [-l MIN_LENGTH] [-w WINDOW] [INPUT]

Streaming FASTQ quality filter. Reads INPUT (or stdin) and writes passing
records to stdout. Summary counters go to stderr.

  -q MIN_MEAN    minimum mean Phred score to keep a read (default 20)
  -l MIN_LENGTH  minimum sequence length to keep a read (default 1)
  -w WINDOW      trim low-quality ends with a sliding WINDOW before scoring
  -h, --help     show this help

Exit status: 0 on success, 1 on usage or parse error.
";

fn parse_args() -> Result<(Config, Option<usize>, Option<String>), String> {
    let mut cfg = Config::default();
    let mut input = None;
    let mut trim_window = None;
    let mut args = env::args().skip(1);

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-h" | "--help" => return Err(USAGE.to_string()),
            "-q" => {
                let value = args.next().ok_or("-q needs a value")?;
                cfg.min_mean = parse_min_mean(&value)?;
            }
            "-l" => {
                let value = args.next().ok_or("-l needs a value")?;
                cfg.min_length = value
                    .parse()
                    .map_err(|_| format!("invalid -l value: {value}"))?;
            }
            "-w" | "--trim-window" => {
                let value = args.next().ok_or("-w needs a value")?;
                let width: usize = value
                    .parse()
                    .map_err(|_| format!("invalid -w value: {value}"))?;
                if width > 0 {
                    trim_window = Some(width);
                }
            }
            other if other.starts_with('-') => return Err(format!("unknown flag: {other}")),
            other => {
                if input.replace(other.to_string()).is_some() {
                    return Err("more than one input file given".into());
                }
            }
        }
    }

    Ok((cfg, trim_window, input))
}

fn open_input(path: Option<&str>) -> io::Result<Box<dyn Read>> {
    match path {
        Some(p) => Ok(Box::new(File::open(p)?)),
        None => Ok(Box::new(io::stdin())),
    }
}

fn run() -> Result<(), String> {
    let (cfg, trim_window, input) = parse_args()?;
    let source = |e: io::Error| format!("cannot read input: {e}");
    let mut reader = Reader::new(BufReader::with_capacity(
        1 << 16,
        open_input(input.as_deref()).map_err(source)?,
    ));
    let stdout = io::stdout();
    let mut out = BufWriter::with_capacity(1 << 16, stdout.lock());
    let mut stats = Stats::default();

    loop {
        let mut record = match reader.next_record() {
            Ok(Some(record)) => record,
            Ok(None) => break,
            Err(fastq::ParseError::Io(e)) => return Err(format!("read error: {e}")),
            Err(e) => return Err(e.to_string()),
        };
        if let Some(width) = trim_window {
            stats.trimmed_bases += trim_ends(&mut record, width, cfg.min_mean) as u64;
        }
        let decision = classify(&record, &cfg);
        stats.record(decision);
        if decision == Decision::Keep {
            record
                .write_to(&mut out)
                .map_err(|e| format!("write error: {e}"))?;
        }
    }

    out.flush().map_err(|e| format!("write error: {e}"))?;
    eprintln!(
        "seen={} kept={} low_quality={} too_short={} trimmed_bases={}",
        stats.seen, stats.kept, stats.low_quality, stats.too_short, stats.trimmed_bases
    );
    Ok(())
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(msg) => {
            let help_only = msg.starts_with("usage:");
            if help_only {
                print!("{msg}");
                ExitCode::SUCCESS
            } else {
                eprintln!("fqfilter: {msg}");
                eprintln!("try 'fqfilter --help'");
                ExitCode::FAILURE
            }
        }
    }
}
