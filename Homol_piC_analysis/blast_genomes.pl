#!/usr/bin/perl

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

# Species identifiers
my @species = ('Hsap','Mmul','Mfas','Cjac','Mmur','Ogar');
my %spec_id = ();
foreach my $i (0..$#species) { $spec_id{$species[$i]} = $i }

# Program name
print("\n--- $0 ---\n");
# Collect command line arguments
my $USAGE = "perl $0 <list_of_files_genomes.txt>\n";
unless ($ARGV[0]) {
	die("\nUsage: $USAGE\n");
}
my $list_of_files = $ARGV[0];
# Extract filenames from list
my @files = FileIO::get_file_data_array($list_of_files);

# Compare each file to each file
foreach my $file_x (@files) {
	foreach my $file_y (@files) {
		# Extract species names with regex (match all from start up to digit)
		my($spec_1) = ($file_x =~ /(^.+?)\./);
		my($spec_2) = ($file_y =~ /(^.+?)\./);
		# Skip if same species
		if ($spec_1 eq $spec_2) { next }
		# Name output file
		my $blast_out = $spec_1.'_VS_'.$spec_2.'.dcm.out';
		push(@blast_outs,$blast_out);		
		# Call system to start blastn
		print("$blast_out ... ");
		unless (-e $blast_out) {
			system("blastn -query $file_x -subject $file_y -out $blast_out -outfmt '7 qseqid sseqid pident evalue qcovs qlen qstart qend slen sstart send sstrand length bitscore' -task dc-megablast");
		}
		print("ok\n");
	}
}

my $out = FileIO::open_outfile($outfile);
foreach my $blast_out (@blast_outs) {
	print($out "$blast_out\n");
}
print("Done.\n");

exit;