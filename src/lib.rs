use std::process::{Command, exit};
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, Read, BufRead, BufReader, Write};
use std::path::Path;



pub fn parse_filename_split(file_path: &str) -> Option<(String, String)> {
    // Split the path into parts based on '/' and then split by '_'
    let parts: Vec<&str> = file_path.split('/').last()?.split('_').collect();

    // Find the index of "S" and "L" if present
    let mut sample_parts = Vec::new();
    let mut read_type = None;

    for part in &parts {
        if part.starts_with('S') && part.len() > 1 && part[1..].chars().all(|c| c.is_digit(10) ) {
            // Skip over "S" and "L" parts, which are lane-related
            //println!("Skipping {}", part);
            continue;
        } else if part.starts_with('L') && part.len() == 4 && part[1..].chars().all(|c| c.is_digit(10) ) {
            // Skip the "L" part, which is lane-related
            //println!("Skipping {}", part);
            continue;
        } else if part.starts_with("R1") || part.starts_with("R2") || part.starts_with("I1") {
            // This is the read type (R1, R2, I1)
            read_type = Some(part[0..2].to_string());
            //println!("Found an READ part {}", part);
            break
        } else {
            // Otherwise, it's part of the sample name
            //println!("Found an sample part {}", part);
            sample_parts.push(part.to_string());
        }
    }

    if let Some(read) = read_type {
        // Join the sample parts to form the sample name
        let sample_name = sample_parts.join("_");
        Some((sample_name.to_string(), read.to_string()))
    } else {
        None // In case there's no read type found
    }
}


pub fn write_sample_files(path: &str, data: &HashMap<String, HashMap<String, String>>) {
    let mut file = File::create(path).expect("Could not create sample file");
    writeln!(file, "Sample_Lane\tR1\tR2\tI1").unwrap();
    // Sort the keys of the outer HashMap (sample_lane)
    let mut sorted_keys: Vec<String> = data.keys().cloned().collect();
    sorted_keys.sort();

    // Iterate through the sorted keys and write the corresponding data
    for sample_lane in sorted_keys {
        if let Some(reads) = data.get(&sample_lane) {
            writeln!(file, "{}\t{}\t{}\t{}", 
                     sample_lane,
                     reads.get("R1").unwrap_or(&"MISSING_R1".to_string()),
                     reads.get("R2").unwrap_or(&"MISSING_R2".to_string()),
                     reads.get("I1").unwrap_or(&"MISSING_I1".to_string())
            ).unwrap();
        }
    }
}

// Helper function to extract the basename
pub fn extract_basename(file_path: Option<&String>) -> Option<String> {
    file_path
        .and_then(|path| Path::new(path).file_name()) // Extract the file name
        .and_then(|name| name.to_str())               // Convert OsStr to &str
        .map(|s| s.to_string())                       // Convert &str to String
}

pub fn write_sample_files_basename(path: &str, data: &HashMap<String, HashMap<String, String>>) {
    let mut file = File::create(path).expect("Could not create sample file");
    writeln!(file, "Sample_Lane\tR1\tR2\tI1").unwrap();
    
    // Sort the keys of the outer HashMap (sample_lane)
    let mut sorted_keys: Vec<String> = data.keys().cloned().collect();
    sorted_keys.sort();

    // Iterate through the sorted keys and write the corresponding data
    for sample_lane in sorted_keys {
        if let Some(reads) = data.get(&sample_lane) {
            writeln!(file, "{}\t{}\t{}\t{}", 
                sample_lane,
                extract_basename(reads.get("R1")).unwrap_or("MISSING_R1".to_string()),
                extract_basename(reads.get("R2")).unwrap_or("MISSING_R2".to_string()),
                extract_basename(reads.get("I1")).unwrap_or("MISSING_I1".to_string())
            ).unwrap();
        }
    }
}

pub fn generate_md5_file_data(data: &HashMap<String, HashMap<String, String>>) -> Vec<(String, String)> {
    // Collect all (basename, md5sum) tuples in sorted order in one step
    let mut all_files: Vec<(String, String)> = data
        .values()  // Iterating over values (inner HashMap)
        .flat_map(|reads| {
            // Sort file paths directly here
            let mut sorted_file_paths: Vec<String> = reads.values().cloned().collect();
            sorted_file_paths.sort(); // Sort the file paths lexicographically
            sorted_file_paths.into_iter() // Convert the sorted file paths into an iterator
                .filter_map(|file_path| {
                    // For each file path, extract the basename and calculate the MD5sum
                       Some((file_path.clone(), get_md5sum(&file_path)))
                })
        })
        .collect();

    // Sort all the (basename, md5sum) tuples by md5sum
    all_files.sort_by(|a, b| a.0.cmp(&b.0)); // Sort by md5sum (tuple.0)

    all_files // Return the sorted (basename, md5sum) vector
}

pub fn write_md5_files(path: &str, data: &Vec::<(String, String)> )-> io::Result<()> {
    let mut file = File::create(path)?;
    writeln!(file, "file_name\tmd5sum").unwrap();

    // Iterate through the sorted keys and write the corresponding data
    for (file_path, md5sum) in data {
        writeln!(file, "{}\t{}", file_path, md5sum).unwrap();
    }
    Ok(())
}

pub fn write_md5_files_basename(path: &str, data: &Vec::<(String, String)> ) -> io::Result<()> {
    let mut file = File::create(path).expect("Could not create md5 file");
    writeln!(file, "file_name\tmd5sum").unwrap();

    // Iterate through the sorted keys and write the corresponding data
    for (file_path, md5sum) in data {
        writeln!(file, "{}\t{}", extract_basename(Some(&file_path)).unwrap(), md5sum).unwrap();
    }
    Ok(())
}



pub fn compute_file_md5_incremental( file_path:&str ) -> io::Result<String> {
    // Run the md5sum command
    let output = Command::new("md5sum")
        .arg(file_path)
        .output()?;
    // Check if the command was successful
    if !output.status.success() {
        return Err(io::Error::new(io::ErrorKind::Other, "md5sum command failed"));
    }

    let hash = String::from_utf8_lossy(&output.stdout);
    Ok( format!("{}", hash.split_whitespace().next().unwrap() ) )
}


pub fn get_md5sum(file_path: &str) -> String {
    let path = Path::new(file_path);
    let md5_file = path.with_extension("fastq.gz.md5sum");
    if md5_file.exists() {
        if let Ok(file) = File::open(&md5_file) {
            let reader = BufReader::new(file);
            if let Some(Ok(line)) = reader.lines().next() {
                return line;
            }
        }
    }

    if let Ok(md5sum) = compute_file_md5_incremental(file_path) {
        let _ = fs::write(&md5_file, &md5sum);
        return md5sum;
    }
    "none".to_string()
}
