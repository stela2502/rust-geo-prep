use std::path::Path;
use std::fs::{self, File, };
use std::io::Write;

use rust_geo_prep::*;

fn create_fastq_file(path: &Path, content: &str) {
    let mut file = File::create(path).expect("Unable to create file");
    file.write_all(content.as_bytes())
        .expect("Unable to write data to file");
}


#[test]
fn test_extract_basename_valid_file() {
    // Test with a valid file path
    let file_path = "/tests/data/info/example1_R1.fastq.gz".to_string();
    let result = extract_basename(Some(&file_path));
    assert_eq!(result, Some("example1_R1.fastq.gz".to_string()));
}

#[test]
fn test_extract_basename_valid_file2() {
    // Test with a valid file path
    let file_path = "tests/data/info/example1_R1.fastq.gz".to_string();
    let result = extract_basename(Some(&file_path));
    assert_eq!(result, Some("example1_R1.fastq.gz".to_string()));
}

#[test]
fn create_and_check_fastq_files() {
    // Define the path to the `data` directory
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

    // Verify that the files were created
    for (file_name, _) in files_and_contents {
        let file_path = data_dir.join(file_name);
        assert!(file_path.exists(), "File {} does not exist", file_name);
    }

}

#[cfg(test)]
mod tests {
	use rust_geo_prep::parse_filename_split;
    #[test]
    fn test_parse_filename_split() {
        let files_and_contents = vec![
            ("example1_S1_L001_R1.fastq.gz", ("example1".to_string(),  "R1".to_string())),
            ("example1_S1_L001_R2.fastq.gz", ("example1".to_string(),  "R2".to_string())),
            ("example1_S1_L001_I1.fastq.gz", ("example1".to_string(),  "I1".to_string())),
            ("example2_L001_R1.fastq.gz", ("example2".to_string(),  "R1".to_string())),
            ("example2_L001_R2.fastq.gz", ("example2".to_string(),  "R2".to_string())),
            ("example3_1_R1.fastq.gz", ("example3_1".to_string(),  "R1".to_string())),
            ("example3_1_R2.fastq.gz", ("example3_1".to_string(),  "R2".to_string())),
            ("example3_1_I1.fastq.gz", ("example3_1".to_string(),  "I1".to_string())),
        ];

        for (file, expected) in files_and_contents {
            let result = parse_filename_split(file);
            assert_eq!(result, Some(expected), "Failed for file: {}", file);
        }
    }
}