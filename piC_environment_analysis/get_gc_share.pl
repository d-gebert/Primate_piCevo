#!/usr/bin/perl
use strict;
use warnings;
use lib '/Users/dgebert/Dropbox/Perlmodules/';
use FileIO;
use FastaIO;

# Constants
$|=1; #Autoflush
# Variables
my %count = ();

# Program name
print("\n--- $0 ---\n");
# Collect command line arguments
my $USAGE = "perl $0 <genome.fas>\n";
unless ($ARGV[0]) {
	die("\nUsage: $USAGE\n");
}
# Input files
my $genome_file = $ARGV[0];

print("Processing...");
# Get genome and chromosome names
my $genome = FastaIO::get_fasta_seqs($genome_file);
my $chroms = FastaIO::get_fasta_headers($genome_file);

# Go through each chromosome
foreach my $chr (@{$chroms}) {
	while ($genome->{$chr} =~ /[AT]/ig){$count{'AT'}++};
	while ($genome->{$chr} =~ /[GC]/ig){$count{'GC'}++};
}
print("done\n");

my $gc_share = $count{'GC'}/($count{'AT'}+$count{'GC'})*100;
my $at_share = $count{'AT'}/($count{'AT'}+$count{'GC'})*100;

print("piC GC share: $gc_share\npiC AT share: $at_share\n");