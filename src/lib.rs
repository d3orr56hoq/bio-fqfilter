//! `fqfilter` — a streaming FASTQ quality filter.
//!
//! Records are parsed one at a time from a [`std::io::BufRead`] and either
//! forwarded unchanged or dropped based on their mean Phred score. Only the
//! current record and the reader/writer buffers are ever resident, so memory
//! use is independent of input size.

pub mod fastq;
pub mod quality;
pub mod window;

pub use fastq::{ParseError, Reader, Record};
pub use quality::mean_phred;

/// Default minimum acceptable mean Phred score (Q20, ~1% error).
pub const DEFAULT_MIN_MEAN: f64 = 20.0;

/// Why a record was rejected.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Reason {
    /// Mean Phred score below the configured threshold.
    LowQuality,
    /// Sequence shorter than the configured minimum after filtering.
    TooShort,
}

/// Decision for one record.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Decision {
    Keep,
    Drop(Reason),
}

/// Filter settings shared by the CLI and library callers.
#[derive(Debug, Clone, Copy)]
pub struct Config {
    pub min_mean: f64,
    pub min_length: usize,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            min_mean: DEFAULT_MIN_MEAN,
            min_length: 1,
        }
    }
}

/// Parse the `-q` threshold.
///
/// `NaN` is rejected because `mean < NaN` is always false, which would keep
/// every read and silently disable filtering. Infinities are accepted: they
/// are valid `f64` values and keep the usual comparison (`-q inf` drops
/// every finite-scoring read).
pub fn parse_min_mean(value: &str) -> Result<f64, String> {
    let min_mean: f64 = value
        .parse()
        .map_err(|_| format!("invalid -q value: {value}"))?;
    if min_mean.is_nan() {
        return Err(format!("invalid -q value: {value}"));
    }
    Ok(min_mean)
}

/// Scores one record against `cfg`.
///
/// Checks are ordered so that the cheapest rejection (length) wins and
/// `mean_phred` runs at most once per record.
pub fn classify(record: &Record, cfg: &Config) -> Decision {
    if record.seq.len() < cfg.min_length {
        return Decision::Drop(Reason::TooShort);
    }
    if mean_phred(&record.qual) < cfg.min_mean {
        return Decision::Drop(Reason::LowQuality);
    }
    Decision::Keep
}

/// Running counters for the CLI summary.
#[derive(Debug, Default)]
pub struct Stats {
    pub seen: u64,
    pub kept: u64,
    pub low_quality: u64,
    pub too_short: u64,
    pub trimmed_bases: u64,
}

impl Stats {
    pub fn record(&mut self, decision: Decision) {
        self.seen += 1;
        match decision {
            Decision::Keep => self.kept += 1,
            Decision::Drop(Reason::LowQuality) => self.low_quality += 1,
            Decision::Drop(Reason::TooShort) => self.too_short += 1,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(seq: &str, qual: &str) -> Record {
        Record {
            header: "r".into(),
            seq: seq.as_bytes().to_vec(),
            plus: "+".into(),
            qual: qual.as_bytes().to_vec(),
        }
    }

    #[test]
    fn keeps_high_quality_record() {
        assert_eq!(
            classify(&record("ACGT", "IIII"), &Config::default()),
            Decision::Keep
        );
    }

    #[test]
    fn drops_low_quality_record() {
        let cfg = Config {
            min_mean: 20.0,
            min_length: 1,
        };
        assert_eq!(
            classify(&record("ACGT", "!!!!"), &cfg),
            Decision::Drop(Reason::LowQuality)
        );
    }

    #[test]
    fn length_check_runs_first() {
        let cfg = Config {
            min_mean: 0.0,
            min_length: 4,
        };
        assert_eq!(
            classify(&record("AC", "!!"), &cfg),
            Decision::Drop(Reason::TooShort)
        );
    }

    #[test]
    fn parse_min_mean_rejects_nan() {
        for value in ["nan", "NaN", "-nan", "+nan"] {
            let err = parse_min_mean(value).unwrap_err();
            assert!(err.contains(value), "{err}");
        }
    }

    #[test]
    fn parse_min_mean_keeps_finite_and_infinite_values() {
        assert_eq!(parse_min_mean("20").unwrap(), 20.0);
        assert_eq!(parse_min_mean("0").unwrap(), 0.0);
        assert_eq!(parse_min_mean("-1.5").unwrap(), -1.5);
        assert!(
            parse_min_mean("inf").unwrap().is_infinite()
                && parse_min_mean("inf").unwrap().is_sign_positive()
        );
        assert!(
            parse_min_mean("+inf").unwrap().is_infinite()
                && parse_min_mean("+inf").unwrap().is_sign_positive()
        );
        assert!(
            parse_min_mean("-inf").unwrap().is_infinite()
                && parse_min_mean("-inf").unwrap().is_sign_negative()
        );
        assert!(parse_min_mean("notanumber").is_err());
    }
}
