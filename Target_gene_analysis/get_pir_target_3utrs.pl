#!/usr/bin/perl
use strict;
use warnings;
use lib '/Users/dgebert/Dropbox/Perlmodules/';
use FileIO;
use FastaIO;

# Constants
$|=1; #Autoflush
# Variables

# Program name
print("\n--- $0 ---\n");
# Collect command line arguments
my $USAGE = "perl $0 <pi_target_genes.txt> <3UTRs.fas> <expression.fas>\n";
unless ($ARGV[0] && $ARGV[1]) {
	die("\nUsage: $USAGE\n");
}
# Input files
my $target_file = $ARGV[0];
my $utrs_3_file = $ARGV[1];
my $expres_file = $ARGV[2];

# Get piRNA target genes
my $targets = get_target_genes($target_file);

# Get 3'-UTR lengths
my $utrlens = get_3utr_lengths($utrs_3_file);

# Get gene expressions
my $express = get_gene_expression($expres_file);

# Output for target genes
my $out = FileIO::open_outfile('pi_tar_utr_lengths.txt');
# Go through each gene
foreach my $gid (sort {$targets->{$b}->[1] <=> $targets->{$a}->[1]} keys %{$targets}) {
	print($out "$targets->{$gid}->[-1]\t$targets->{$gid}->[1]\t$utrlens->{$gid}\n");
}
# Output for non-target genes
my $out2 = FileIO::open_outfile('pi_nontar_utr_lengths.txt');
# Go through each gene
foreach my $gid (sort {$utrlens->{$b} <=> $utrlens->{$a}} keys %{$utrlens}) {
	# If testis-expressed and non-target
	if ($express->{$gid}->[15]) {
		print($out2 "$gid\t$utrlens->{$gid}\n") unless $targets->{$gid};
	}
}

exit;

################################# subroutines #################################

sub get_3utr_lengths {
	# Take name of pic pseudogenes file
	my($infile) = @_;
	# Variables
	my %utr_lens = ();
	# Extract file data
	my @file_data = FileIO::get_file_data_array($infile);
	# Parse file data
	foreach my $line (@file_data) {
		# Fasta header
		if ($line =~ /^>/) {
			# Delete arrow and split entries
			$line =~ s/>//;
			my @d = split(/\|/,$line);
			# Get data
			my $gid = $d[0];
			my $tid = $d[1];
			my $beg = $d[2] ? $d[2] : 0;
			my $end = $d[3] ? $d[3] : 0;
			$beg =~ s/\;.*//;
			$end =~ s/\;.*//;
			my $len = $end-$beg;
			# Save 3’-UTR length
			if ($utr_lens{$gid}) {
				$utr_lens{$gid} = $len if $len > $utr_lens{$gid};
			} else {
				$utr_lens{$gid} = $len;
			}
		}
	}
	return \%utr_lens;
}

sub get_target_genes {
	# Take name of pic pseudogenes file
	my($infile, $rpkm_th) = @_;
	# RPKM threshold
	unless ($rpkm_th) { $rpkm_th = 0 }
	# Variables
	my %tab_data = ();
	# Extract file data
	my @file_data = FileIO::get_file_data_array($infile);
	# Parse file data
	foreach my $line (@file_data) {
		# Skip empty line
		if ($line =~ /^$/) { next }
		if ($line =~ /^Gene_id/) { next }
		# Split tab-separated entries
		my @d = split(/\t/,$line);
		# Get gene id
		my $g_id = $d[0];
		# Get RPKM
		my $rpkm = $d[1];
		# Save data
		$tab_data{$g_id} = [@d] if $rpkm >= $rpkm_th;
	}
	return \%tab_data;
}

sub get_gene_expression {
	# Take name of pic pseudogenes file
	my($infile) = @_;
	# Variables
	my %expression = ();
	# Extract file data
	my @file_data = FileIO::get_file_data_array($infile);
	# Parse file data
	foreach my $line (@file_data) {
		# Skip comment and title lines
		if ($line =~ /^#/ || $line =~ /^Gene/) { next }
		# Split entries
		my @d = split(/\t/,$line);
		# Get gene id and delete from array
		my $gid = $d[0];
		shift(@d);
		# Save gene expression
		$expression{$gid} = [@d];
	}
	return \%expression;
}