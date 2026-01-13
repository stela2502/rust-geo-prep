
use walkdir::WalkDir;
use clap::Parser;
use std::path::PathBuf;

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
    
    let sample_file_path = format!("{}.tsv", opts.prefix);
    let files_file_path = format!("{}_md5sum.tsv", opts.prefix);
    let collection_script_path = if cfg!(windows) {
        format!("{}_collection_script.ps1", opts.prefix)
    } else {
        format!("{}_collection_script.sh", opts.prefix)
    };
    let collection_dest = format!("{}_all_files_copied", opts.prefix);
    
    //let sample_file_path_basename = format!("{}_basename_sample_lines.tsv", opts.prefix);
    //let files_file_path_basename = format!("{}_basename_files_md5sum_lines.tsv", opts.prefix);

    println!("We are searching for files ending on either of these strings {:?}", opts.suffixes );

    
    let mut data = SampleFiles::new();
    
    let (added, visited) = match data.ingest_dir(".", &opts.suffixes, &opts.exclude) {
        Err(e) => {
            eprintln!("\n❌ Failed while scanning input directories:");
            eprintln!("   {e}\n");
            std::process::exit(1);
        },
        Ok(i) => i,
    };

    let _ = data.write_sample_files_basename(&sample_file_path);
    let _ = data.write_md5_files_basename(&files_file_path);
    let _ = if cfg!(windows) {
        data.write_collect_all_files_script_ps1(&collection_script_path, &collection_dest)
    } else {
        data.write_collect_all_files_script_sh(&collection_script_path, &collection_dest)
    };

    //let _ = data.write_sample_files_basename(&sample_file_path_basename);
    //let _ = data.write_md5_files_basename(&files_file_path_basename);

    let run_cmd = if cfg!(windows) {
        format!(".\\{}", collection_script_path)
    } else {
        format!("bash {}", collection_script_path)
    };

    println!(
        "\n{} target files detected; {} added and {} samples identified.\n\
         \nOutput files:\n\
         - Sample table      : {}\n\
         - MD5 checksum table: {}\n\
         - Collection script : {}\n\
         - Copy destination  : {}\n\
         \nNext steps:\n\
         1) Review the TSV files for correctness.\n\
         2) Run the collection script to gather all referenced files:\n\
            {}\n\
         3) Use the TSV files to fill in the official GEO submission spreadsheets.\n\
         \nNote: These files are intermediate manifests. Experimental metadata must be added manually.\n",
        visited, 
        added,
        data.len(),
        sample_file_path,
        files_file_path,
        collection_script_path,
        collection_dest,
        run_cmd
    );
}


