use std::path::Path;
use std::fs::{self, File, };
use std::io::Write;

fn create_fastq_file(path: &Path, content: &str) {
    let mut file = File::create(path).expect("Unable to create file");
    file.write_all(content.as_bytes())
        .expect("Unable to write data to file");
}

#[test]
fn create_and_check_fastq_files() {
    // Define the path to the `data` directory
    let data_dir = Path::new("tests/data");

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
            println!("Creating file: {}", file_name);  // Optional: For debugging
            create_fastq_file(&file_path, content);
        } else {
            println!("File {} already exists", file_name);  // Optional: For debugging
        }
    }

    // Verify that the files were created
    for (file_name, _) in files_and_contents {
        let file_path = data_dir.join(file_name);
        assert!(file_path.exists(), "File {} does not exist", file_name);
    }
}