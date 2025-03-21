# rust-geo-prep

The tool is designed to create two types of reference tables for a GEO submission. The initial project for this was a multi-year, multi-sample 10x-based single-cell experiment with a large number of FASTQ files.

GEO expects two types of information:

1. A list of all files with their corresponding MD5 checksums.
2. A table mapping each sample and sequencing lane to its associated files.

This is a direct re-implementation of the [original Perl script](./get_sample_fastq_files.pl).\
I decided to implement it in Rust because it is more convenient at the moment.

In terms of speed, rust might be a little faster but as the main bottleneck is reading the files and data from disk this difference will not be too impressive. 
However, I prefer Rust’s testing framework over Perl’s, and distributing a compiled binary is simpler than setting up a Perl package.

## Installation

You need to have the [Rust compiler](https://www.rust-lang.org/tools/install) installed.

To install `rust-geo-prep`, use Cargo:

```sh
cargo install --git https://github.com/stela2502/rust-geo-prep.git
```

## Usage

When started, the tool scans all subfolders of the current folder and collects all "fastq.gz" files. It processes each file name by splitting it by '\_', dropping any `S\d`, `L\d\d\d`, and `[RI]` entries, and merging the rest into a sample name.

The tool generates two output files for each type of reference table: one containing only the file basenames and another with the full file paths. This ensures that even if sample names are duplicated across sequencing runs, the absolute file path provides an additional level of distinction. Both output files maintain the same structure, meaning that corresponding samples are on the same line in each output, allowing easy cross-referencing.

Run the script in the folder containing your data:
```
rust-geo-prep
```

This will scan all subdirectories for fastq.gz files and process their filenames to extract sample and sequencing lane information.

The script will generate two sets of output files, each consisting of:

 - One version with only basenames (for GEO table submission).
 - One version with full paths (to ensure unique identification when needed).

Output Files
Files containing only basenames:

 - sample_collection_basename_.files_md5sum_lines.tsv  
    (List of MD5 checksums for each file, using only basenames.)
 - sample_collection_basename_.sample_lines.tsv  
    (List of sequencing files grouped by sample and lane, using only basenames.)

Files containing full paths:

 - sample_collection_.files_md5sum_lines.tsv  
    (List of MD5 checksums for each file, using absolute paths.)
 - sample_collection_.sample_lines.tsv  
    (List of sequencing files grouped by sample and lane, using absolute paths.)

This dual-output structure helps mitigate sample name collisions across sequencing runs while still providing a clean format for GEO submissions.
