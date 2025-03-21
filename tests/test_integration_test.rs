use std::process::Command;
use std::fs::{self, File};
use std::path::Path;
use std::io::Write;
use std::io::BufReader;
use std::io::BufRead;


fn create_fastq_file(path: &Path, content: &str) {
    let mut file = File::create(path).expect("Unable to create file");
    file.write_all(content.as_bytes())
        .expect("Unable to write data to file");
}


fn clean_test_output() {
    let files_to_remove = vec![
        "sample_collection_basename_.files_md5sum_lines.tsv",
        "sample_collection_basename_.sample_lines.tsv",
        "sample_collection_.files_md5sum_lines.tsv",
        "sample_collection_.sample_lines.tsv",
    ];

    for file in files_to_remove {
        if Path::new(file).exists() {
            fs::remove_file(file).expect("Failed to remove old test output file");
        }
    }
}

#[test]
fn program_run() {
    clean_test_output();
    rust_geo_prep() ;
    // test the different outfiles
    sample_collection_files_md5sum_lines();
    sample_collection_basename_files_md5sum_lines();

    clean_test_output()
}
fn rust_geo_prep() {
    let test_data_dir = "tests/data";

    // first make sure the files exists and have the content that is expected:
    let data_dir = Path::new("tests/data/info/");

    // Check if the data directory exists, create it if not
    if !data_dir.exists() {
        fs::create_dir_all(data_dir).expect("Failed to create data directory");
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

    // Create files with the provided content if they don't exist
    for (file_name, content) in &files_and_contents {
        let file_path = data_dir.join(file_name);
        if !file_path.exists() {
            //println!("Creating file: {}", file_name);  // Optional: For debugging
            create_fastq_file(&file_path, content);
        } else {
            //println!("File {} already exists", file_name);  // Optional: For debugging
        }
    }
    
    // Run the binary
    let output = Command::new(env!("CARGO_BIN_EXE_rust-geo-prep"))
        .current_dir(test_data_dir) // Run inside test data folder
        .output()
        .expect("Failed to execute rust-geo-prep");

    // Check if execution was successful
    assert!(output.status.success(), "Program did not run successfully");

    // List expected output files
    let expected_files = vec![
        "sample_collection_basename_files_md5sum_lines.tsv",
        "sample_collection_basename_sample_lines.tsv",
        "sample_collection_files_md5sum_lines.tsv",
        "sample_collection_sample_lines.tsv",
    ];

    // Verify output files exist
    for file in expected_files {
        let path = format!("{}/{}", test_data_dir, file);
        assert!(fs::metadata(&path).is_ok(), "Missing expected output file: {}", file);
    }

}


#[test]
fn sample_collection_files_md5sum_lines() {
    // Path to the test file
    let path = "tests/data/sample_collection_files_md5sum_lines.tsv";
    
    // Open the file
    let file = File::open(path).expect("Unable to open file");

    // Create a buffered reader for efficient reading
    let reader = BufReader::new(file);

    // Expected values based on your sample file content
    let expected_contents = vec![
        ("./info/example1_S1_L001_I1.fastq.gz", "1da0250da36f7f38d11f4f08397e02d9"),
        ("./info/example1_S1_L001_R1.fastq.gz", "220693693f35b15196bc2f2fa8238e7b"),
        ("./info/example1_S1_L001_R2.fastq.gz", "28f6a6cefb6b7ea07049b8261c52cab8"),
        ("./info/example1_S2_L001_I1.fastq.gz", "933471e0abaab240b18683bc2267f3bc"),
        ("./info/example1_S2_L001_R1.fastq.gz", "867171df270ed55ca348daf1369f5c25"),
        ("./info/example1_S2_L001_R2.fastq.gz", "f60431ad04351b3eb786879ed18440c8"),
        ("./info/example2_L001_R1.fastq.gz", "220693693f35b15196bc2f2fa8238e7b"),
        ("./info/example2_L001_R2.fastq.gz", "28f6a6cefb6b7ea07049b8261c52cab8"),
        ("./info/example3_1_I1.fastq.gz", "32a0a8c330f2cdcccafee94b03d1a04e"),
        ("./info/example3_1_R1.fastq.gz", "220693693f35b15196bc2f2fa8238e7b"),
        ("./info/example3_1_R2.fastq.gz", "28f6a6cefb6b7ea07049b8261c52cab8"),
    ];

    // Iterate over each line in the file
    for (index, line) in reader.lines().enumerate() {
        let line = line.expect("Unable to read line");
        let parts: Vec<&str> = line.split('\t').collect();

        // Skip header line
        if index == 0 {
            continue;
        }

        // Check that the line contains exactly two parts (file_name, md5sum)
        assert_eq!(parts.len(), 2, "Line does not have exactly two columns");

        // Extract file_name and md5sum
        let file_name = parts[0].trim();
        let md5sum = parts[1].trim();

        // Check if the file name and md5sum match the expected ones
        assert_eq!(file_name, expected_contents[index - 1].0, "File name mismatch at line {}", index);
        assert_eq!(md5sum, expected_contents[index - 1].1, "MD5 sum mismatch at line {}", index);
    }
}


#[test]
fn sample_collection_basename_files_md5sum_lines() {
    // Path to the test file
    let path = "tests/data/sample_collection_basename_files_md5sum_lines.tsv";
    
    // Open the file
    let file = File::open(path).expect("Unable to open file");

    // Create a buffered reader for efficient reading
    let reader = BufReader::new(file);

    // Expected values based on your sample file content
    let expected_contents = vec![
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
    ];

    // Iterate over each line in the file
    for (index, line) in reader.lines().enumerate() {
        let line = line.expect("Unable to read line");
        let parts: Vec<&str> = line.split('\t').collect();

        // Skip header line
        if index == 0 {
            continue;
        }

        // Check that the line contains exactly two parts (file_name, md5sum)
        assert_eq!(parts.len(), 2, "Line does not have exactly two columns");

        // Extract file_name and md5sum
        let file_name = parts[0].trim();
        let md5sum = parts[1].trim();

        // Check if the file name and md5sum match the expected ones
        assert_eq!(file_name, expected_contents[index - 1].0, "BN File name mismatch at line {}", index);
        assert_eq!(md5sum, expected_contents[index - 1].1, "BN MD5 sum mismatch at line {}", index);
    }
}

#[test]
fn test_sample_collection_sample_lines() {
    // Path to the test file
    let path = "tests/data/sample_collection_sample_lines.tsv";

    // Expected values based on your sample file content
    let expected_contents = vec![
        ("example1", vec![
            "./info/example1_S1_L001_I1.fastq.gz", "./info/example1_S1_L001_R1.fastq.gz", "./info/example1_S1_L001_R2.fastq.gz",
            "./info/example1_S2_L001_I1.fastq.gz", "./info/example1_S2_L001_R1.fastq.gz", "./info/example1_S2_L001_R2.fastq.gz"
        ]),
        ("example2", vec![
            "./info/example2_L001_R1.fastq.gz", "./info/example2_L001_R2.fastq.gz"
        ]),
        ("example3_1", vec![
            "./info/example3_1_I1.fastq.gz", "./info/example3_1_R1.fastq.gz", "./info/example3_1_R2.fastq.gz"
        ])
    ];

    // Open the file
    let file = File::open(path).expect("Unable to open file");
    // Create a buffered reader for efficient reading
    let reader = BufReader::new(file);

    // Iterate over each line in the file
    for (index, line) in reader.lines().enumerate() {
        let line = line.expect("Unable to read line");
        let parts: Vec<&str> = line.split('\t').collect();

        // Skip header line
        if index == 0 {
            continue;
        }

        // Check that the line has at least one sample column (adjust this if necessary)
        assert!(parts.len() >= 2, "Line does not have enough columns");

        // Extract sample name
        let sample_name = parts[0].trim();

        // Get the actual filenames from the line (starting from index 1 onward)
        let filenames: Vec<String> = parts[1..]
            .iter()
            .map(|&filename| filename.trim().to_string())
            .collect();

        // Find the expected filenames for the sample
        let expected = expected_contents.iter().find(|(name, _)| name == &sample_name);

        assert!(expected.is_some(), "Sample name {} not found in expected contents", sample_name);

        let expected_files = expected.unwrap().1.clone();

        // Sort both expected and actual filenames for a flexible comparison
        let mut filenames = filenames.clone();
        let mut expected_files = expected_files.clone();

        filenames.sort();
        expected_files.sort();

        // Compare filenames (R1, R2, I1)
        assert_eq!(filenames, expected_files, "File mismatch for sample {}", sample_name);
    }
}

#[test]
fn test_sample_collection_sample_lines_basename() {
    // Path to the test file
    let path = "tests/data/sample_collection_basename_sample_lines.tsv";

    // Expected values based on your sample file content
    let expected_contents = vec![
        ("example1", vec![
            "example1_S1_L001_I1.fastq.gz", "example1_S1_L001_R1.fastq.gz", "example1_S1_L001_R2.fastq.gz",
            "example1_S2_L001_I1.fastq.gz", "example1_S2_L001_R1.fastq.gz", "example1_S2_L001_R2.fastq.gz"
        ]),
        ("example2", vec![
            "example2_L001_R1.fastq.gz", "example2_L001_R2.fastq.gz"
        ]),
        ("example3_1", vec![
            "example3_1_I1.fastq.gz", "example3_1_R1.fastq.gz", "example3_1_R2.fastq.gz"
        ])
    ];

    // Open the file
    let file = File::open(path).expect("Unable to open file");
    // Create a buffered reader for efficient reading
    let reader = BufReader::new(file);

    // Iterate over each line in the file
    for (index, line) in reader.lines().enumerate() {
        let line = line.expect("Unable to read line");
        let parts: Vec<&str> = line.split('\t').collect();

        // Skip header line
        if index == 0 {
            continue;
        }

        // Check that the line has at least one sample column (adjust this if necessary)
        assert!(parts.len() >= 2, "Line does not have enough columns");

        // Extract sample name
        let sample_name = parts[0].trim();

        // Get the actual filenames from the line (starting from index 1 onward)
        let filenames: Vec<String> = parts[1..]
            .iter()
            .map(|&filename| filename.trim().to_string())
            .collect();

        // Find the expected filenames for the sample
        let expected = expected_contents.iter().find(|(name, _)| name == &sample_name);

        assert!(expected.is_some(), "Sample name {} not found in expected contents", sample_name);

        let expected_files = expected.unwrap().1.clone();

        // Sort both expected and actual filenames for a flexible comparison
        let mut filenames = filenames.clone();
        let mut expected_files = expected_files.clone();

        filenames.sort();
        expected_files.sort();

        // Compare filenames (R1, R2, I1)
        assert_eq!(filenames, expected_files, "File mismatch for sample {}", sample_name);
    }
}