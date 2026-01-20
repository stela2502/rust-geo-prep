// tests/example_test.rs

use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};

use flate2::write::GzEncoder;
use flate2::Compression;
use tempfile::TempDir;

fn keep_dir_on_err(tmp: TempDir, err: impl std::fmt::Display) -> ! {
    let path: PathBuf = tmp.into_path(); // prevents cleanup
    panic!("{err}\n\nTest workspace kept at:\n  {}", path.display());
}

/// Write gzipped bytes to `path` (real .gz content, not just a filename suffix).
fn write_gzip_text<P: AsRef<Path>>(path: P, text: &str) -> io::Result<()> {
    let path = path.as_ref();
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    let f = fs::File::create(path)?;
    let mut gz = GzEncoder::new(f, Compression::default());
    gz.write_all(text.as_bytes())?;
    gz.finish()?;
    Ok(())
}

/// Write plain text to a file, creating parents.
fn write_text<P: AsRef<Path>>(path: P, text: &str) -> io::Result<()> {
    let path = path.as_ref();
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    fs::write(path, text)?;
    Ok(())
}

/// Build the example INPUT/ tree described in the README.
/// Returns the root folder (`.../INPUT`).
fn create_example_tree(base_dir: &Path) -> io::Result<PathBuf> {
    let input = base_dir.join("INPUT");

    // folders
    fs::create_dir_all(input.join("experiment_1/data"))?;
    fs::create_dir_all(input.join("experiment_1/sampleA/outs/filtered_features_bc_matrix"))?;
    fs::create_dir_all(input.join("experiment_1/geo_downloaded_data"))?;
    fs::create_dir_all(input.join("experiment_2"))?;
    fs::create_dir_all(input.join("old_runs"))?;

    // FASTQ gz (dummy but structurally plausible)
    write_gzip_text(
        input.join("experiment_1/data/sampleA_R1.fastq.gz"),
        "@SEQ_ID\nACGTACGTACGT\n+\nFFFFFFFFFFFF\n",
    )?;
    write_gzip_text(
        input.join("experiment_1/data/sampleA_R2.fastq.gz"),
        "@SEQ_ID\nTGCATGCATGCA\n+\nFFFFFFFFFFFF\n",
    )?;

    // 10x triplet (gz)
    write_gzip_text(
        input.join("experiment_1/sampleA/outs/filtered_features_bc_matrix/barcodes.tsv.gz"),
        "AAACCTGAGAAACCAT-1\nAAACCTGAGAAACCAA-1\n",
    )?;
    write_gzip_text(
        input.join("experiment_1/sampleA/outs/filtered_features_bc_matrix/features.tsv.gz"),
        "GeneA\tGeneA\tExpression\nGeneB\tGeneB\tExpression\n",
    )?;
    write_gzip_text(
        input.join("experiment_1/sampleA/outs/filtered_features_bc_matrix/matrix.mtx.gz"),
        "%%MatrixMarket matrix coordinate integer general\n2 2 2\n1 1 5\n2 2 3\n",
    )?;

    // 10x H5 placeholder (not real HDF5; enough for filename-based collection tests)
    write_text(
        input.join("experiment_1/sampleA/outs/filtered_feature_bc_matrix.h5"),
        "Dummy 10x HDF5 placeholder\n",
    )?;

    Ok(input)
}

#[test]
fn example_tree_is_created_with_expected_files() -> io::Result<()> {
    let tmp = TempDir::new()?;
    let input = create_example_tree(tmp.path())?;

    // Existence checks (this is the “unit-test” part)
    assert!(input.join("experiment_1/data/sampleA_R1.fastq.gz").is_file());
    assert!(input.join("experiment_1/data/sampleA_R2.fastq.gz").is_file());

    assert!(input
        .join("experiment_1/sampleA/outs/filtered_features_bc_matrix/barcodes.tsv.gz")
        .is_file());
    assert!(input
        .join("experiment_1/sampleA/outs/filtered_features_bc_matrix/features.tsv.gz")
        .is_file());
    assert!(input
        .join("experiment_1/sampleA/outs/filtered_features_bc_matrix/matrix.mtx.gz")
        .is_file());

    assert!(input
        .join("experiment_1/sampleA/outs/filtered_feature_bc_matrix.h5")
        .is_file());

    Ok(())
}

fn must_exist(p: &Path) -> Result<(), String> {
    if p.is_file() {
        Ok(())
    } else {
        Err(format!("Missing expected file: {}", p.display()))
    }
}

/// End-to-end test: runs the CLI on the example tree and checks expected outputs.
/// Requires dev-dep: assert_cmd = "2"
#[test]
fn cli_runs_on_example_tree_and_creates_outputs() {
    let tmp = TempDir::new().expect("TempDir");

    let result: Result<(), String> = (|| {
        // Setup input tree
        let input = create_example_tree(tmp.path()).map_err(|e| e.to_string())?;

        // Output prefix
        let prefix = tmp.path().join("out").join("example");
        fs::create_dir_all(prefix.parent().unwrap()).map_err(|e| e.to_string())?;

        // Run binary
        let mut cmd = assert_cmd::Command::cargo_bin("rust-geo-prep")
            .map_err(|e| format!("binary rust-geo-prep not built: {e}"))?;

        cmd.arg("--input")
            .arg(input.as_os_str())
            .arg("--exclude")
            .arg("geo_downloaded_data")
            .arg("--exclude")
            .arg("old_runs")
            .arg("--suffix")
            .arg(".fastq.gz")
            .arg("--suffix")
            .arg("filtered_feature_bc_matrix.h5")
            .arg("--suffix")
            .arg("matrix.mtx.gz")
            .arg("--prefix")
            .arg(prefix.to_string_lossy().to_string());

        // Run + assert success (assert_cmd will print nice diagnostics on failure)
        cmd.assert().success();

        // Expected artifacts (new output names)
        let sample_collection = PathBuf::from(format!("{}.tsv", prefix.display()));
        let md5_table = PathBuf::from(format!("{}_md5sum.tsv", prefix.display()));
        let pairs_table = PathBuf::from(format!("{}_pairs.tsv", prefix.display()));
        let script_path = if cfg!(windows) {
            PathBuf::from(format!("{}_collection_script.ps1", prefix.display()))
        } else {
            PathBuf::from(format!("{}_collection_script.sh", prefix.display()))
        };

        must_exist(&sample_collection)?;
        must_exist(&md5_table)?;
        must_exist(&pairs_table)?;
        must_exist(&script_path)?;

        Ok(())
    })();

    if let Err(e) = result {
        keep_dir_on_err(tmp, e);
    }

    // On success: TempDir is dropped and cleaned up.
}
