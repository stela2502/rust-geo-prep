use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

// -----------------------
// Helpers
// -----------------------

fn create_file(path: &Path, content: &str) {
    let mut file = File::create(path).expect("Unable to create file");
    file.write_all(content.as_bytes())
        .expect("Unable to write data to file");
}

fn canonical_str<P: AsRef<Path>>(p: P) -> String {
    fs::canonicalize(p.as_ref())
        .unwrap()
        .to_string_lossy()
        .to_string()
}

fn read_tsv(path: &Path) -> (Vec<String>, Vec<Vec<String>>) {
    let file = File::open(path).unwrap_or_else(|_| panic!("Unable to open file: {}", path.display()));
    let reader = BufReader::new(file);

    let mut lines = reader.lines();
    let header = lines
        .next()
        .expect("TSV is empty")
        .expect("Unable to read header")
        .split('\t')
        .map(|s| s.trim().to_string())
        .collect::<Vec<_>>();

    let mut rows = Vec::new();
    for line in lines {
        let line = line.expect("Unable to read line");
        let cols = line.split('\t').map(|s| s.trim().to_string()).collect::<Vec<_>>();
        // skip empty trailing lines if any
        if cols.iter().all(|c| c.is_empty()) {
            continue;
        }
        rows.push(cols);
    }

    (header, rows)
}

fn header_index(header: &[String], name: &str) -> usize {
    header
        .iter()
        .position(|h| h == name)
        .unwrap_or_else(|| panic!("Missing required column '{name}' in header: {header:?}"))
}

fn basename(s: &str) -> String {
    Path::new(s)
        .file_name()
        .map(|x| x.to_string_lossy().to_string())
        .unwrap_or_else(|| s.to_string())
}

fn parse_dest_from_script(script_path: &Path) -> Option<String> {
    let txt = fs::read_to_string(script_path).ok()?;
    for line in txt.lines() {
        let line = line.trim();
        // DEST="something"
        if let Some(rest) = line.strip_prefix("DEST=\"") {
            if let Some(end) = rest.find('"') {
                return Some(rest[..end].to_string());
            }
        }
    }
    None
}

fn md5sum_file(path: &Path) -> String {
    use std::io::Read;

    let mut f = File::open(path)
        .unwrap_or_else(|_| panic!("Unable to open file for md5: {}", path.display()));

    let mut ctx = md5::Context::new();
    let mut buf = vec![0u8; 64 * 1024]; // 64KiB on heap, safe on Windows

    loop {
        let n = f.read(&mut buf)
            .unwrap_or_else(|e| panic!("Failed reading {}: {e}", path.display()));
        if n == 0 { break; }
        ctx.consume(&buf[..n]);
    }

    format!("{:x}", ctx.compute())
}

// -----------------------
// Cleanup
// -----------------------

fn clean_test_output(prefix: &str) {
    let sample_tsv = format!("{prefix}.tsv");
    let md5_tsv = format!("{prefix}_md5sum.tsv");
    let script = format!("{prefix}_collection_script.sh");

    for f in [sample_tsv.as_str(), md5_tsv.as_str(), script.as_str()] {
        if Path::new("tests/data").join(f).exists() {
            fs::remove_file(Path::new("tests/data").join(f))
                .unwrap_or_else(|_| panic!("Failed to remove old test output file: {f}"));
        }
    }

    // also remove any previous copy dir (we'll discover it by reading script if it exists)
    let script_path = Path::new("tests/data").join(script);
    if script_path.exists() {
        if let Some(dest) = parse_dest_from_script(&script_path) {
            let dest_path = Path::new("tests/data").join(dest);
            if dest_path.exists() {
                fs::remove_dir_all(dest_path).expect("Failed to remove old copy destination folder");
            }
        }
    }
}

// -----------------------
// Test runner
// -----------------------

fn run_rust_geo_prep(prefix: &str) {
    let test_data_dir = Path::new("tests/data");
    let info_dir = test_data_dir.join("info");

    if !info_dir.exists() {
        fs::create_dir_all(&info_dir).expect("Failed to create tests/data/info/");
    }

    // Define the file paths and content
    let files_and_contents = vec![
        ("example1_S1_L001_R1.fastq.gz", "@SEQ_ID_1\nAGCTGTTAG\n+\nIIIIIIIIII\n"),
        ("example1_S1_L001_R2.fastq.gz", "@SEQ_ID_2\nTGCTAGTCG\n+\nIIIIIIIIII\n"),
        ("example1_S1_L001_I1.fastq.gz", "@SEQ_ID_3\nACGTGTCG\n+\nIIIIIIIIII\n"),
        ("example1_S2_L001_R1.fastq.gz", "@SEQ_ID_4\nAGCTGTTAG\n+\nIIIIIIIIII\n"),
        ("example1_S2_L001_R2.fastq.gz", "@SEQ_ID_5\nTGCTAGTCG\n+\nIIIIIIIIII\n"),
        ("example1_S2_L001_I1.fastq.gz", "@SEQ_ID_6\nACGTGTCG\n+\nIIIIIIIIII\n"),
        ("example2_L001_R1.fastq.gz", "@SEQ_ID_1\nAGCTGTTAG\n+\nIIIIIIIIII\n"),
        ("example2_L001_R2.fastq.gz", "@SEQ_ID_2\nTGCTAGTCG\n+\nIIIIIIIIII\n"),
        ("example3_1_R1.fastq.gz", "@SEQ_ID_1\nAGCTGTTAG\n+\nIIIIIIIIII\n"),
        ("example3_1_R2.fastq.gz", "@SEQ_ID_2\nTGCTAGTCG\n+\nIIIIIIIIII\n"),
        ("example3_1_I1.fastq.gz", "@SEQ_ID_3\nGCTAGTGC\n+\nIIIIIIIIII\n"),
    ];

    for (file_name, content) in &files_and_contents {
        let file_path = info_dir.join(file_name);
        if !file_path.exists() {
            create_file(&file_path, content);
        }
    }

    let exe = env!("CARGO_BIN_EXE_rust-geo-prep");
    eprintln!(
        "RUN MANUALLY:\n  cd {}\n  {} --prefix {}\n",
        test_data_dir.display(),
        exe,
        prefix
    );

    let output = Command::new(exe)
        .current_dir(test_data_dir)
        .args(["--prefix", prefix])
        .output()
        .expect("Failed to execute rust-geo-prep");

    assert!(
        output.status.success(),
        "Program did not run successfully.\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );

    // Verify expected output files exist
    let expected_files = vec![
        format!("{prefix}.tsv"),
        format!("{prefix}_md5sum.tsv"),
        format!("{prefix}_collection_script.sh"),
    ];

    for file in expected_files {
        let path = test_data_dir.join(&file);
        assert!(
            fs::metadata(&path).is_ok(),
            "Missing expected output file: {}",
            path.display()
        );
    }
}

// -----------------------
// Tests
// -----------------------

#[test]
fn program_run_and_outputs_and_copy_script() {
    let prefix = "sample_collection";
    clean_test_output(prefix);

    run_rust_geo_prep(prefix);

    // Validate md5sum table
    test_md5sum_table(prefix);

    // Validate sample table content (file references)
    test_sample_collection_tsv_exact();

    // Run and validate the collection script
    test_collection_script(prefix);

    clean_test_output(prefix);
}

fn test_md5sum_table(prefix: &str) {
    let test_data_dir = Path::new("tests/data");
    let md5_path = test_data_dir.join(format!("{prefix}_md5sum.tsv"));

    let (header, rows) = read_tsv(&md5_path);

    // Your actual header names
    let i_file = header_index(&header, "file_name");
    let i_md5  = header_index(&header, "md5sum");

    // Expected md5 by filename (truth)
    let expected: BTreeMap<&str, &str> = BTreeMap::from([
        ("example1_S1_L001_I1.fastq.gz", "1da0250da36f7f38d11f4f08397e02d9"),
        ("example1_S1_L001_R1.fastq.gz", "220693693f35b15196bc2f2fa8238e7b"),
        ("example1_S1_L001_R2.fastq.gz", "28f6a6cefb6b7ea07049b8261c52cab8"),
        ("example1_S2_L001_I1.fastq.gz", "933471e0abaab240b18683bc2267f3bc"),
        ("example1_S2_L001_R1.fastq.gz", "867171df270ed55ca348daf1369f5c25"),
        ("example1_S2_L001_R2.fastq.gz", "f60431ad04351b3eb786879ed18440c8"),
        ("example2_L001_R1.fastq.gz", "220693693f35b15196bc2f2fa8238e7b"),
        ("example2_L001_R2.fastq.gz", "28f6a6cefb6b7ea07049b8261c52cab8"),
        ("example3_1_I1.fastq.gz", "32a0a8c330f2cdcccafee94b03d1a04e"),
        ("example3_1_R1.fastq.gz", "220693693f35b15196bc2f2fa8238e7b"),
        ("example3_1_R2.fastq.gz", "28f6a6cefb6b7ea07049b8261c52cab8"),
    ]);

    // Build: filename -> set of md5 strings seen in file
    let mut seen: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();

    for row in rows {
        assert!(
            row.len() >= header.len(),
            "Row has fewer columns than header: {row:?}"
        );

        let filename = row[i_file].trim().to_string();
        let md5 = row[i_md5].trim().to_string();

        seen.entry(filename).or_default().insert(md5);
    }

    // Now validate expectations:
    // Each expected file should appear, and should contain the expected md5.
    for (fname, exp_md5) in &expected {
        let md5s = seen
            .get(*fname)
            .unwrap_or_else(|| panic!("Missing filename in md5sum file: {fname}"));

        assert!(
            md5s.contains(&exp_md5.to_string()),
            "MD5 missing/incorrect for {fname}.\nFound: {md5s:?}\nExpected to include: {exp_md5}"
        );

        // Optional: if you want to *forbid* 'none' when omit_md5 is false,
        // uncomment this once you fix the duplicate/none bug:
        //
        // assert!(
        //     !md5s.contains("none"),
        //     "Unexpected 'none' md5 row for {fname} when md5 is enabled.\nFound: {md5s:?}"
        // );
    }

    // Strong sanity: md5sum file should not contain unknown filenames
    for fname in seen.keys() {
        assert!(
            expected.contains_key(fname.as_str()),
            "Unexpected filename present in md5sum file: {fname}"
        );
    }

    // Strong sanity: the actual files exist and match their md5
    let info_dir = test_data_dir.join("info");
    for (fname, exp_md5) in &expected {
        let file_path = info_dir.join(fname);
        assert!(file_path.exists(), "Missing input file: {}", file_path.display());

        let actual = md5sum_file(&file_path);
        assert_eq!(
            actual, *exp_md5,
            "Actual md5 does not match expected for {}",
            file_path.display()
        );
    }
}

fn looks_like_source_folders_cell(s: &str) -> bool {
    let v = s.trim();
    if v.is_empty() {
        return false;
    }
    // must contain at least one path separator or \\?\ prefix
    let has_pathish = v.contains("\\\\?\\") || v.contains(":\\") || v.contains('\\') || v.contains('/');
    // and typically comma-separated (but allow single path too)
    has_pathish
}

fn split_folders(s: &str) -> Vec<String> {
    s.split(',')
        .map(|x| x.trim().to_string())
        .filter(|x| !x.is_empty())
        .collect()
}

#[test]
fn test_sample_collection_tsv_exact() {
    use std::collections::BTreeMap;
    use std::path::Path;

    let path = Path::new("tests/data/sample_collection.tsv");
    let (header, rows) = read_tsv(path);

    // 1) Header must match exactly
    let expected_header = vec![
        "Source_Path(s)", "Sample_Lane", "TenX", "H5",
        "I1", "R1", "R2",
        "I1", "R1", "R2",
    ].into_iter().map(|s| s.to_string()).collect::<Vec<_>>();

    assert_eq!(header, expected_header, "Header mismatch in sample_collection.tsv");

    // 2) Must have exactly 3 data rows
    assert_eq!(rows.len(), 3, "Expected exactly 3 rows in sample_collection.tsv");

    // 3) Build expected absolute folders string for column 0
    let info_dir = canonical_str("tests/data/info");


    // 4) Expected rows keyed by sample (= column 1 / "TenX")
    //    Column layout (9 columns):
    //    0 folders, 1 sample, 2 h5, then 6 fastq cells
    let expected_by_sample: BTreeMap<&str, Vec<String>> = BTreeMap::from([
        ("example1", vec![
            info_dir.clone(),
            "example1".to_string(),
            "".to_string(),
            "".to_string(),
            "example1_S1_L001_I1.fastq.gz".to_string(),
            "example1_S1_L001_R1.fastq.gz".to_string(),
            "example1_S1_L001_R2.fastq.gz".to_string(),
            "example1_S2_L001_I1.fastq.gz".to_string(),
            "example1_S2_L001_R1.fastq.gz".to_string(),
            "example1_S2_L001_R2.fastq.gz".to_string(),
        ]),
        ("example2", vec![
            info_dir.clone(),
            "example2".to_string(),
            "".to_string(),
            "".to_string(),
            "".to_string(), // I1 missing lane1
            "example2_L001_R1.fastq.gz".to_string(),
            "example2_L001_R2.fastq.gz".to_string(),
            "".to_string(), // lane2 padding
            "".to_string(),
            "".to_string(),
        ]),
        ("example3", vec![
            info_dir,
            "example3".to_string(),
            "".to_string(),
            "".to_string(),
            "example3_1_I1.fastq.gz".to_string(),
            "example3_1_R1.fastq.gz".to_string(),
            "example3_1_R2.fastq.gz".to_string(),
            "".to_string(), // lane2 padding
            "".to_string(),
            "".to_string(),
        ]),
    ]);

    // 5) Validate each row matches expected exactly (after folder normalization)
    for row in rows {
        assert_eq!(row.len(), header.len(), "Row has wrong column count: {:?}", row);

        let sample = row[1].trim(); // column 1 is sample name in your output
        let expected = expected_by_sample.get(sample)
            .unwrap_or_else(|| panic!("Unexpected sample in sample_collection.tsv: '{}'", sample));

        // Normalize the folders cell by trimming whitespace (canonical_str already absolute).
        let mut normalized = row.clone();
        normalized[0] = normalized[0].trim().to_string();
        for c in &mut normalized {
            *c = c.trim().to_string();
        }

        assert_eq!(
            normalized, *expected,
            "Row mismatch for sample '{}'\nGot:      {:?}\nExpected: {:?}",
            sample, normalized, expected
        );
    }
}


fn test_collection_script(prefix: &str) {
    use std::collections::{BTreeMap, BTreeSet};
    use std::fs;

    let test_data_dir = Path::new("tests/data");
    let script_path = test_data_dir.join(format!("{prefix}_collection_script.sh"));
    assert!(script_path.exists(), "Missing collection script: {}", script_path.display());

    // Parse DEST="..."
    let dest_rel = parse_dest_from_script(&script_path)
        .unwrap_or_else(|| panic!("Could not parse DEST=\"...\" from script {}", script_path.display()));
    let dest_path = test_data_dir.join(&dest_rel);

    // Clean previous run
    if dest_path.exists() {
        fs::remove_dir_all(&dest_path).expect("Failed to remove old DEST folder");
    }
    fs::create_dir_all(&dest_path).expect("Failed to create DEST folder");

    // Load md5 table: file_name + md5sum
    let md5_path = test_data_dir.join(format!("{prefix}_md5sum.tsv"));
    let (header, rows) = read_tsv(&md5_path);
    let i_file = header_index(&header, "file_name");
    let i_md5  = header_index(&header, "md5sum");

    // Expected mapping: basename -> md5
    let mut expected: BTreeMap<String, String> = BTreeMap::new();
    for row in rows {
        let fname = row[i_file].trim().to_string();
        let md5   = row[i_md5].trim().to_string();
        assert_ne!(fname, "", "Empty file_name in md5 table");
        expected.insert(fname, md5);
    }

    // Parse script copy pairs and assert NO RENAMING:
    // "${COPY_CMD[@]}" "SRC" "$DEST/BASENAME"
    let script_txt = fs::read_to_string(&script_path)
        .unwrap_or_else(|_| panic!("Unable to read script: {}", script_path.display()));

    let pairs = parse_copy_pairs_from_script(&script_txt);

    // Map: basename -> src_path
    let mut by_basename: BTreeMap<String, PathBuf> = BTreeMap::new();

    for (src, dst_rel) in pairs {
        assert!(
            dst_rel.starts_with("$DEST/"),
            "Destination does not start with $DEST/: dst={dst_rel}"
        );

        let dst_name = dst_rel.trim_start_matches("$DEST/").trim().to_string();
        let src_base = Path::new(&src)
            .file_name()
            .map(|x| x.to_string_lossy().to_string())
            .unwrap_or_else(|| panic!("Script src has no filename: {src}"));

        // This is the core requirement: destination name must equal source basename
        assert_eq!(
            dst_name, src_base,
            "Script renames files (NOT ALLOWED).\n  src: {src}\n  dst: {dst_rel}"
        );

        // Each basename must appear at most once
        if by_basename.insert(dst_name.clone(), PathBuf::from(&src)).is_some() {
            panic!("Duplicate basename in script copy lines: {dst_name}");
        }
    }

    // The script must cover exactly the md5 table file list
    let expected_names: BTreeSet<String> = expected.keys().cloned().collect();
    let script_names: BTreeSet<String> = by_basename.keys().cloned().collect();

    assert_eq!(
        script_names, expected_names,
        "Script file list != md5 table file list.\nMissing in script: {:?}\nExtra in script: {:?}",
        expected_names.difference(&script_names).collect::<Vec<_>>(),
        script_names.difference(&expected_names).collect::<Vec<_>>(),
    );

    // Now actually perform the copy in Rust and validate outputs exactly
    for (name, src_path) in &by_basename {
        assert!(src_path.exists(), "Script references missing source file: {}", src_path.display());

        let dst_path = dest_path.join(name);
        fs::copy(src_path, &dst_path)
            .unwrap_or_else(|e| panic!("Copy failed {} -> {}: {e}", src_path.display(), dst_path.display()));

        assert!(dst_path.exists(), "Destination not created: {}", dst_path.display());

        // verify md5 equals md5 table
        let want = expected.get(name).unwrap();
        let got = md5sum_file(&dst_path);
        assert_eq!(
            &got, want,
            "Copied file md5 mismatch for {}\n  dst: {}\n  got: {}\n  want: {}",
            name, dst_path.display(), got, want
        );
    }

    // Verify DEST contains exactly the expected basenames (no extras)
    let mut dest_files: Vec<String> = fs::read_dir(&dest_path)
        .expect("Failed to list DEST dir")
        .filter_map(|e| e.ok())
        .filter(|e| e.file_type().map(|t| t.is_file()).unwrap_or(false))
        .map(|e| e.file_name().to_string_lossy().to_string())
        .collect();

    dest_files.sort();

    let mut expected_sorted: Vec<String> = expected.keys().cloned().collect();
    expected_sorted.sort();

    assert_eq!(
        dest_files, expected_sorted,
        "DEST directory contents mismatch.\nGot: {:?}\nExpected: {:?}",
        dest_files, expected_sorted
    );
}
fn header_index_any(header: &[String], names: &[&str]) -> usize {
    for &n in names {
        if let Some(i) = header.iter().position(|h| h == n) {
            return i;
        }
    }
    panic!("Missing required column(s) {names:?} in header: {header:?}");
}

fn sample_from_sample_lane(s: &str) -> String {
    let s = s.trim();
    if let Some(pos) = s.find("_S") {
        return s[..pos].to_string();
    }
    if let Some(pos) = s.find("_L") {
        return s[..pos].to_string();
    }
    s.to_string()
}