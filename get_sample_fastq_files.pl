#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;
use Digest::MD5;


# Hash to store organized data
my %data;

# Find and process .fastq.gz files
my $path_to_search = shift || "."; # Default to current directory
my $sampleid = shift || ""; # Default to none ""

my $path = $path_to_search;
$path =~ s/[\/\\\s]+/_/g; 

open my $SAMPLE ,"> $path/sample_collection_$sampleid.sample_lines.tsv" or die "couln not create '$path./sample_collection_$sampleid.sample_lines.tsv'\n";
open my $FILES ,"> $path/sample_collection_$sampleid.files_md5sum_lines.tsv" or die "couln not create '$path./sample_collection_$sampleid.files_md5sum_lines.tsv'\n";

open my $SAMPLE_BN ,"> $path/sample_collection_basename_$sampleid.sample_lines.tsv" or die "couln not create '$path./sample_collection_basename_$sampleid.sample_lines.tsv'\n";
open my $FILES_BN ,"> $path/sample_collection_basename_$sampleid.files_md5sum_lines.tsv" or die "couln not create '$path./sample_collection_basename_$sampleid.files_md5sum_lines.tsv'\n";


open my $find, "-|", "find $path_to_search -name '$sampleid*.fastq.gz'" or die "Failed to execute find: $!\n";

while (<$find>) {
    chomp;
    # Parse filename, e.g., sample_L001_R1_001.fastq.gz
    print( "." );
    #NextSeq7/HVHNYBGXB/outs/fastq_path/Chromium_20191011/SI-GA-B1_3/Sample_16_S7_L001_R1_001.fastq.gz
    #data/files/2025_001/0_fastq/AGCYWLM5/negative_1_S8_R2_001.fastq.gz
    if (m!.*/(.*)_(S.*L\d{3}).*(R[12]|I1).*\.fastq\.gz!) {
        my ($sample, $lane, $read_type) = ($1, $2, $3);
	#print("${sample}_$lane\n");
        $data{"${sample}_$lane"}{$read_type} = $_;
    }elsif (m!.*/(.*)_(S\d+).*(R[12]|I1).*\.fastq\.gz!) {
	my ($sample, $lane, $read_type) = ($1, $2, $3);
	#print("${sample}_$lane\n");
        $data{"${sample}_$lane"}{$read_type} = $_;
    }
    $| = 1;  # Force output flushin
}
close $find;

# Print organized data
print $SAMPLE join("\t", "Sample_Lane", "R1", "R2", "I1"), "\n";
print $SAMPLE_BN join("\t", "Sample_Lane", "R1", "R2", "I1"), "\n";

foreach my $sample_lane (sort keys %data) {
    #for my $file in %{$data{$sample_lane}}{
    #    check_and_add_md5sum( $file );
    #}    
    print $SAMPLE join("\t", 
        $sample_lane,
        $data{$sample_lane}{R1} // "MISSING_R1",
        $data{$sample_lane}{R2} // "MISSING_R2",
        $data{$sample_lane}{I1} // "MISSING_I1"
    ), "\n";
    print $SAMPLE_BN join("\t",
        $sample_lane,
        basename($data{$sample_lane}{R1} // "MISSING_R1"),
        basename($data{$sample_lane}{R2} // "MISSING_R2"),
        basename($data{$sample_lane}{I1} // "MISSING_I1")
    ), "\n";


}


# Print organized data
print $FILES join("\t", "file_name", "md5sum"), "\n";
print $FILES_BN join("\t", "file_name", "md5sum"), "\n";

foreach my $sample_lane (sort keys %data) {
    #for my $file in %{$data{$sample_lane}}{
    #    check_and_add_md5sum( $file );
    #}
    print $FILES join("\n",
        $data{$sample_lane}{R1}."\t".get_md5sum($data{$sample_lane}{R1}) // "MISSING_R1\tnone",
        $data{$sample_lane}{R2}."\t".get_md5sum($data{$sample_lane}{R2}) // "MISSING_R2\tnone",
        $data{$sample_lane}{I1}."\t".get_md5sum($data{$sample_lane}{I1}) // "MISSING_I1\tnone"
    ), "\n";
    print $FILES_BN join("\n",
        basename($data{$sample_lane}{R1}// "MISSING_R1")."\t".get_md5sum($data{$sample_lane}{R1}) // "MISSING_R1\tnone",
        basename($data{$sample_lane}{R2}// "MISSING_R1")."\t".get_md5sum($data{$sample_lane}{R2}) // "MISSING_R2\tnone",
        basename($data{$sample_lane}{I1}// "MISSING_R1")."\t".get_md5sum($data{$sample_lane}{I1}) // "MISSING_I1\tnone"
    ), "\n";

}

close $SAMPLE ;
close $FILES ;

print("Data written to '$path.$sampleid.sample_lines.tsv' and '$path.$sampleid.files_md5sum_lines.tsv'\n");

# Function to check and add md5sum calculation command to the script
sub check_and_add_md5sum {
    my $file = shift;
    
    # Check if the file and its corresponding .md5sum file exist
    my $md5file = $file . ".md5sum";
    
    if (-e $md5file) {
        print STDERR "MD5sum already exists for: $file\n";
    } else {
        # If the .md5sum file doesn't exist, add the command to the calculate.md5sums.sh script
        open my $sh_script, '>>', 'calculate.md5sums.sh' or die "Cannot open calculate.md5sums.sh: $!\n";
        print $sh_script "md5sum $file > $md5file\n";
        close $sh_script;
        print STDERR "Added md5sum calculation command for $file to calculate.md5sums.sh\n";
    }
}

# Function to check if the md5sum file exists and return it, or calculate it if not
sub get_md5sum {
    my $file = shift;

    if (undef $file) {
	return "";
    }
    unless ( -f $file ) {
        return "" ;
    }
    # Create the expected .md5sum file name
    my $md5file = $file . ".md5sum";

    # If the .md5sum file exists, read and return its content
    if (-e $md5file) {
        open my $fh, '<', $md5file or die "Cannot open $md5file: $!\n";
        my $md5sum = <$fh>;
        close $fh;
        chomp $md5sum;  # Remove the newline
        return $md5sum;
    }
    else {
        # If no .md5sum file, calculate and return the MD5 sum
        open my $fh, '<', $file or die "Cannot open $file: $!\n";
        my $md5 = Digest::MD5->new;

        # Read the file and calculate its MD5 sum
        $md5->addfile($fh);
        close $fh;

        # Get the MD5 sum as a hexadecimal string
        my $md5sum = $md5->hexdigest;

        # Optionally, save this MD5 sum to a file
        open my $out_fh, '>', $md5file or die "Cannot write to $md5file: $!\n";
        print $out_fh "$md5sum\n";
        close $out_fh;

        # Return the calculated MD5 sum
        return $md5sum;
    }
}

