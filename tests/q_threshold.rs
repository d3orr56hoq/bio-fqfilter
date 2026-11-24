//! Integration test: `-q` rejects NaN but still accepts infinities.

use std::io::Write;
use std::process::{Command, Stdio};

const TWO_READS: &str = "\
@high
ACGT
+
IIII
@low
ACGT
+
!!!!
";

fn run_q(value: &str, input: &str) -> (i32, String, String) {
    let mut child = Command::new(env!("CARGO_BIN_EXE_fqfilter"))
        .args(["-q", value])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn fqfilter");
    child
        .stdin
        .as_mut()
        .unwrap()
        .write_all(input.as_bytes())
        .unwrap();
    let out = child.wait_with_output().unwrap();
    (
        out.status.code().unwrap_or(-1),
        String::from_utf8_lossy(&out.stdout).into_owned(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

#[test]
fn finite_threshold_keeps_only_high_quality_reads() {
    let (code, stdout, stderr) = run_q("20", TWO_READS);
    assert_eq!(code, 0, "stderr={stderr}");
    assert!(stdout.contains("@high"));
    assert!(!stdout.contains("@low"));
    assert!(
        stderr.contains("seen=2 kept=1 low_quality=1"),
        "stderr={stderr}"
    );
}

#[test]
fn nan_threshold_is_rejected() {
    for value in ["nan", "NaN", "-nan"] {
        let (code, stdout, stderr) = run_q(value, TWO_READS);
        assert_eq!(code, 1, "value={value} stderr={stderr}");
        assert!(stdout.is_empty(), "value={value} wrote {stdout:?}");
        assert!(
            stderr.contains(&format!("invalid -q value: {value}")),
            "value={value} stderr={stderr}"
        );
    }
}

#[test]
fn infinite_threshold_still_filters() {
    let (code, stdout, stderr) = run_q("inf", TWO_READS);
    assert_eq!(code, 0, "stderr={stderr}");
    assert!(
        stdout.is_empty(),
        "inf should drop every finite-scoring read"
    );
    assert!(
        stderr.contains("seen=2 kept=0 low_quality=2"),
        "stderr={stderr}"
    );
}

#[test]
fn unparseable_threshold_is_still_rejected() {
    let (code, stdout, stderr) = run_q("notanumber", TWO_READS);
    assert_eq!(code, 1, "stderr={stderr}");
    assert!(stdout.is_empty());
    assert!(
        stderr.contains("invalid -q value: notanumber"),
        "stderr={stderr}"
    );
}
