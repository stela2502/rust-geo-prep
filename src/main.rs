
use walkdir::WalkDir;
use clap::Parser;

use std::collections::HashSet;
#[cfg(unix)]
use std::os::unix::fs::MetadataExt;

use rust_geo_prep::sample_files::SampleFiles;

/// Submitting data to GEO is complex. 
/// This tool helps by collecting the different fastq files and grouping them into samples groups.
/// It also calculates the md5sums and reports them for every fastq file.
#[derive(Parser)]
#[clap(version = "1.0.0", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the output prefix
    #[clap(short, long, default_value="sample_collection")]
    prefix: String,

    /// Omit md5sums (default calculate all)
    #[clap(short, long )]
    omit_md5: bool,

    /// File suffixes treated as target files
    ///
    /// Can be specified multiple times:
    ///   --suffix .fastq.gz --suffix .fq.gz
    #[clap(
        short = 's',
        long = "suffix",
        multiple_occurrences = true,
        default_values = &[".fastq.gz", ".fq.gz"]
    )]
    suffixes: Vec<String>,

    /// File suffixes treated as target files
    ///
    /// Can be specified multiple times:
    ///   --suffix .fastq.gz --suffix .fq.gz
    #[clap(
        short = 'e',
        long = "exclude",
        multiple_occurrences = true,
    )]
    exclude: Vec<String>,

    /// Allow collecting files with identical basenames in different directories
    /// (required for 10x matrix / features / barcodes layouts)
    #[clap(short, long )]
    allow_duplicates: bool,

}

fn is_excluded_path(p: &std::path::Path, excludes: &[String]) -> bool {
    if excludes.is_empty() {
        return false;
    }

    // Compare on a stable string form (works even if parts are non-utf8-ish: fallback lossy)
    let p_str = p.to_string_lossy();

    // Also check components to support excluding by folder name (e.g. "outs", "fastq_path")
    for ex in excludes {
        if ex.is_empty() {
            continue;
        }

        // 1) Full substring match on the whole path (most flexible)
        if p_str.contains(ex) {
            return true;
        }

        // 2) Folder/component exact match (exclude by directory name)
        if p.components()
            .any(|c| c.as_os_str().to_string_lossy() == ex.as_str())
        {
            return true;
        }
    }
    false
}

fn main(){
    let opts: Opts = Opts::parse();
    
    let sample_file_path = format!("{}_sample_lines.tsv", opts.prefix);
    let files_file_path = format!("{}_files_md5sum_lines.tsv", opts.prefix);
    
    let sample_file_path_basename = format!("{}_basename_sample_lines.tsv", opts.prefix);
    let files_file_path_basename = format!("{}_basename_files_md5sum_lines.tsv", opts.prefix);

    println!("We are searching for files ending on either of these strings {:?}", opts.suffixes );

    
    let mut data = SampleFiles::new( opts.omit_md5 );
    let mut id= 0;
    
    #[cfg(unix)]
    let mut visited: HashSet<(u64, u64)> = HashSet::new();

    #[cfg(not(unix))]
    let mut visited: std::collections::HashSet<std::path::PathBuf> = HashSet::new();

    // Parse files and group by sample name, technicalities, and read type
    for entry in WalkDir::new( "." ).follow_links(true).into_iter().filter_map(Result::ok) {
        let file_path = entry.path();

        // --- exclude folders/files EARLY (before metadata/canonicalize) ---
        if is_excluded_path(file_path, &opts.exclude) {
            continue;
        }

        // Only directories need loop protection
        if let Ok(md) = file_path.metadata() {

            if md.is_dir() {
                #[cfg(unix)]
                let key = (md.dev(), md.ino());

                #[cfg(not(unix))]
                let key = {
                    // Windows fallback: canonicalized path
                    use std::path::PathBuf;
                    std::fs::canonicalize(&file_path)
                        .unwrap_or_else(|_| PathBuf::from(file_path))
                };


                if !visited.insert(key) {
                    // Already seen → skip this directory entirely
                    // Prevents infinite recursion
                    continue;
                }
            }
        }
        if let Some(file_name) = file_path.file_name().and_then(|n| n.to_str()) {
            id +=1;
            if opts.suffixes.iter().any(|s| file_name.ends_with(s)) {
                let fname = match std::fs::canonicalize(file_path) {
                    Ok(real) => real.to_string_lossy().to_string() ,
                    Err(e) => {
                        // If canonicalize fails, use the original path
                        eprintln!("canonicalize - failed for {} with error {e:?}", file_path.display() );
                        file_path.to_string_lossy().to_string()
                    }
                };
                // this now internally checks if the file is a duplicate!
                // it also handles 10x triplets so that they become unique, too.
                data.add_file(&fname);
            }
        }
    }

    let _ = data.write_sample_files(&sample_file_path);
    let _ = data.write_md5_files(&files_file_path);
    let _ = data.write_sample_files_basename(&sample_file_path_basename);
    let _ = data.write_md5_files_basename(&files_file_path_basename);

    println!("{}/{} files detected - data written to '{}', '{}', '{}' and '{}'", 
        data.len(),
        id,
        &sample_file_path,
        &files_file_path, 
        &sample_file_path_basename,
        &files_file_path_basename
    );
}


