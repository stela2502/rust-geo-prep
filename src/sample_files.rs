use std::collections::HashMap;
use std::process::{Command, exit};
use std::fs::{self, File};
use std::io::{self, Read, BufRead, BufReader, Write};
use std::path::Path;

#[derive(Debug)]
pub struct SampleFiles {
    filenames: Vec<(String, String)>,  // A single vector of (filenames, md5sums)
    filenames_by_sample: HashMap<String, Vec<usize>>,  // Maps sample_name -> list of indices into `filenames`
    filenames_by_sample_tech: HashMap<(String, String), Vec<usize>>,  // Maps (sample_name, technicalities) -> list of indices into `filenames`
}

impl SampleFiles {
    pub fn new() -> Self {
        SampleFiles {
            filenames: Vec::new(),
            filenames_by_sample: HashMap::new(),
            filenames_by_sample_tech: HashMap::new(),
        }
    }

    fn samples(&self) -> Vec<&String> {
        let mut sample_keys: Vec<&String> = self.filenames_by_sample.keys().collect();
        sample_keys.sort();
        sample_keys

    }


    pub fn write_sample_files(&self, path: &str) {
        let mut file = match File::create(path){
            Ok(f) => f,
            Err(e) => {
                panic!( "Could not create sample file {}:\n{}", path, e)
            },
        };
        writeln!(file, "Sample_Lane\tR1\tR2\tI1").unwrap();
        // Sort the keys of the outer HashMap (sample_lane)

        for sample in self.samples() {
            match self.get_files_by_sample( sample ){
                Some(files) => {
                    let entry = files.into_iter()
                        .map(|(fname ,_) | fname )
                        .collect::<Vec<_>>()
                        .join("\t");
                    writeln!(file, "{}\t{}", sample, entry ).unwrap();
                },
                None =>{} ,
            }
        }
    }

    pub fn write_sample_files_basename(&self, path: &str ) {
        let mut file = File::create(path).expect("Could not create sample file");
        writeln!(file, "Sample_Lane\tR1\tR2\tI1").unwrap();

        for sample in self.samples() {
            match self.get_files_by_sample( sample ){
                Some(files) => {
                    let entry = files.into_iter()
                        .map(|(fname ,_) | self.extract_basename( &fname ).unwrap() )
                        .collect::<Vec<_>>()
                        .join("\t");
                    writeln!(file, "{}\t{}", sample, entry ).unwrap();
                },
                None =>{} ,
            }
        }
    }

    pub fn write_md5_files(&self, path: &str )-> io::Result<()> {
        let mut file = File::create(path)?;
        writeln!(file, "file_name\tmd5sum").unwrap();

        // Iterate through the sorted keys and write the corresponding data
        for (file_path, md5sum) in self.get_files(None) {
            writeln!(file, "{}\t{}", file_path, md5sum).unwrap();
        }
        Ok(())
    }

    pub fn write_md5_files_basename(&self, path: &str ) -> io::Result<()> {
        let mut file = File::create(path).expect("Could not create md5 file");
        writeln!(file, "file_name\tmd5sum").unwrap();

        // Iterate through the sorted keys and write the corresponding data
        for (file_path, md5sum) in self.get_files(None) {
            writeln!(file, "{}\t{}", self.extract_basename( &file_path ).unwrap(), md5sum).unwrap();
        }
        Ok(())
    }

    // Add a file with its sample name and technicalities
    pub fn add_file(&mut self, file_path: &str ) {
        // Add the file to the main filenames vector
        if let Some(basename) = self.extract_basename( file_path ) {
            if basename.starts_with("Undetermined") || basename.starts_with("Unmapped") || basename.starts_with("Umapped")  {
                // igonore that crap
                return
            }
        }

        let index = self.filenames.len();
        let md5sum = self.get_md5sum( file_path );
        self.filenames.push((file_path.to_string(), md5sum) );
        if let Some( (sample, technicalities, _read) ) = self.parse_filename_split( &file_path ){
            // Add the index to the sample_name entry
            self.filenames_by_sample.entry(sample.clone())
                .or_insert_with(Vec::new)
                .push(index);

            // Add the index to the (sample_name, technicalities) entry
            self.filenames_by_sample_tech.entry((sample, technicalities))
                .or_insert_with(Vec::new)
                .push(index);
        }

        

    }

    fn parse_filename_split(&self, file_path: &str) -> Option<(String, String, String)> {
        // Split the path into parts based on '/' and then split by '_'

        if let Some(ret) = self.parse_matrix_triplets( file_path ){
            return Some(ret);
        }

        let parts: Vec<&str> = file_path.split('/').last()?.split('_').collect();

        // Find the index of "S" and "L" if present
        let mut sample_parts = Vec::new();
        let mut tech = Vec::new();
        let mut read_type = None;

        for part in &parts {
            if sample_parts.is_empty() {
                // some less clever sample names might be S1 S2 S345223
                sample_parts.push(part.to_string());
            }
            else if part.starts_with('S') && part.len() > 1 && part[1..].chars().all(|c| c.is_digit(10) ) {
                // Skip over "S" and "L" parts, which are lane-related
                tech.push( part.to_string() );
                //println!("Skipping {}", part);
                continue;
            } else if part.starts_with('L') && part.len() == 4 && part[1..].chars().all(|c| c.is_digit(10) ) {
                // Skip the "L" part, which is lane-related
                //println!("Skipping {}", part);
                tech.push( part.to_string() );
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
            let technicalities = tech.join("_");
            Some((sample_name.to_string(), technicalities.to_string(), read.to_string()))
        } else {
            None // In case there's no read type found
        }
    }

    /// fix the barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz files
    fn parse_matrix_triplets(&self, file_path: &str) -> Option<(String, String, String)> {
        let suffixes = [
            ("barcodes.tsv.gz", "barcodes", "tsv"),
            ("features.tsv.gz", "features", "tsv"),
            ("matrix.mtx.gz",   "matrix",   "mtx"),
        ];

        for (suffix, kind, ext) in suffixes {
            if let Some(stripped) = file_path.strip_suffix(suffix) {
                return Some((stripped.to_string(), kind.to_string(), ext.to_string()));
            }
        }

        None
    }

    /// Retrieve filenames by sample name, sorted by filename lexicographically
    fn get_files(&self, ids: Option<&Vec<usize>> ) -> Vec<(String, String)> {
        let mut sorted_files = match ids {
            Some(id_s ) => id_s.into_iter().map(|id| self.filenames[*id].clone()).collect(),
            None => self.filenames.clone(),
        };
        // Sort by the filename (first element of the tuple)
        sorted_files.sort_by(|a, b| a.0.cmp(&b.0)); 
        sorted_files
    }

    // Retrieve filenames by sample name, sorted lexicographically
    fn get_files_by_sample(&self, sample: &str) -> Option<Vec<(String, String)>> {
        self.filenames_by_sample.get(sample).map(|indices| {
            // Sort the indices lexicographically
            self.get_files( Some(indices) )
        })
    }

    // Retrieve filenames by sample name + technicalities, sorted lexicographically
    fn get_files_by_sample_tech(&self, sample: &str, tech: &str) -> Option<Vec<(String,String)>> {
        self.filenames_by_sample_tech.get(&(sample.to_string(), tech.to_string())).map(|indices| {
            // Sort the indices lexicographically
            self.get_files( Some(indices) )
        })
    }

    // Helper function to extract the basename
    fn extract_basename(&self, file_path: &str ) -> Option<String> {
        Path::new(file_path).file_name() // Extract the file name
            .and_then(|name| name.to_str())               // Convert OsStr to &str
            .map(|s| s.to_string())                       // Convert &str to String
    }

    fn compute_file_md5_incremental( &self, file_path:&str ) -> io::Result<String> {
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


    fn get_md5sum(&self, file_path: &str) -> String {
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

        if let Ok(md5sum) = self.compute_file_md5_incremental(file_path) {
            let _ = fs::write(&md5_file, &md5sum);
            return md5sum;
        }
        "none".to_string()
    }





}
