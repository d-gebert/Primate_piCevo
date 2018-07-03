#!/usr/bin/perl
# Call of BLASTN on files from a given list

use strict;
use warnings;
use lib '/Users/dgebert/Dropbox/Perlmodules';
use FileIO;

# Global variables
my @blast_outs = ();
# Global constants
my $outfile = 'list_of_files_blastout.txt';
# Options
$|=1; #Autoflush

# Genome identifiers
my @genomes = ('GRCh38','rheMac8','macFas5','calJac3','micMur2','otoGar3');

# Program name
print("\n--- $0 ---\n");

# Collect command line arguments
my $USAGE = "perl $0 <list_of_files_cds.txt>\n";
unless ($ARGV[0]) {
	die("\nUsage: $USAGE\n");
}
my $files_list = $ARGV[0];
# Extract filenames from list
my @files = FileIO::get_file_data_array($files_list);

# Compare each p1 file to each p5 file
print("Processing");
foreach my $x (0..$#files) {
	foreach my $y (0..$#files) {
		# Extract species names with regex (match all from start up to digit)
		my $spec_1  = $genomes[$x];
		my $spec_2  = $genomes[$y];
		if ($spec_1 eq $spec_2) { next }
		# Name output file
		my $blast_out = $spec_1.'_vs_'.$spec_2.'.cds.blast';
		push(@blast_outs,$blast_out);
		print(".");
		# System call to start blastn
		if (-e $blast_out) { next }
		system("blastn -query $files[$x] -subject $files[$y] -out $blast_out -outfmt '7 qseqid sseqid pident evalue qcovs qlen qstart qend slen sstart send sstrand length bitscore' -task dc-megablast");
	}
}
print("done.\n");

my $out = FileIO::open_outfile($outfile);
foreach my $blast_out (@blast_outs) {
	print($out "$blast_out\n");
}

exit;