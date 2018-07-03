#!/usr/bin/perl
use strict;
use warnings;
use lib '/Users/dgebert/Dropbox/Perlmodules/';
use FileIO;
use BioStat;

# Constants
my $outfile = 'Genome_identities.txt';
$|=1; #Autoflush
# Variables
my @identity_matrix = ();

# Species identifiers
my @species = ('Hsap','Mmul','Mfas','Cjac','Mmur','Ogar');
my %spec_id = ();
foreach my $i (0..$#species) { $spec_id{$species[$i]} = $i }

# Program name
print("\n--- $0 ---\n");
# Collect command line arguments
my $USAGE = "perl $0 <blast_output_file_list.txt>\n";
unless ($ARGV[0]) {
	die("\nUsage: $USAGE\n");
}
# Input files
my $blast_file_list = $ARGV[0];
# Get blast file list
my @blast_files = FileIO::get_file_data_array($blast_file_list);

# Extract and save total identity for each species combination
foreach my $blast_file (@blast_files) {
	# Get species names and ids
	my($spec_x) = ($blast_file =~ /^(\w+)\_VS/);
	my($spec_y) = ($blast_file =~ /VS_(\w+)\./);
	my $x = $spec_id{$spec_x};
	my $y = $spec_id{$spec_y};
	# Get genomic blast hits
	my $blast_hits = get_blast_hits($blast_file);
	# Get total identity
	my $identity = get_blast_total_identity($blast_hits);
	# Save identity in matrix
	$identity_matrix[$x][$y] = $identity;
}

# Output matrix
my $out = FileIO::open_outfile($outfile);
foreach my $y (0..$#species) { print($out "\t$species[$y]") }
print($out "\n");
foreach my $x (0..$#species) {
	print($out "$species[$x]\t");
	foreach my $y (0..$#species) {
		if ($x == $y) { $identity_matrix[$x][$y] = 100 }
		printf($out "%.1f\t", $identity_matrix[$x][$y]);
	}
	print($out "\n");
}

exit;

################################# subroutines #################################

sub get_blast_hits {
	# Take infile name
	my($file) = @_;
	# Get file data
	my @file_data = FileIO::get_file_data_array($file);
	# Initialize variables
	my @blast_hits = ();
	# Parse blast output
	foreach my $line (@file_data) {
		# Results line
		if ($line =~ /^[^#]/) {
			# Get info
			my @d = split(/\t/, $line);
			my $query_id = $d[0];
			# Save hits
			push(@blast_hits,[@d]);
		}
	}
	return \@blast_hits;
}

sub get_blast_total_identity {
	# Take blast hits array ref
	my($blast_hits) = @_;
	# Initialize variables
	my $n_match_total = 0;
	my $n_mismatch_total = 0;
	# Extract total counts of matches and mismatches from each hit
	foreach my $hit (sort {$a->[6] <=> $b->[6]} @{$blast_hits}) {
		my $ident = $hit->[2];
		my $a_len = $hit->[12];
		my $n_match = int(($a_len*$ident/100)+0.5);
		my $n_mismatch = $a_len-$n_match;
		$n_match_total += $n_match;
		$n_mismatch_total += $n_mismatch;
	}
	# Calculate and return total identity
	my $n_bases_total = $n_match_total+$n_mismatch_total;
	my $identity_total = $n_match_total/$n_bases_total*100;
	return $identity_total;
}