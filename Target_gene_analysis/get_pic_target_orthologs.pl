#!/usr/bin/perl
use strict;
use warnings;
use lib '/Users/dgebert/Dropbox/Perlmodules/';
use FileIO;
use FastaIO;

# Constants
$|=1; #Autoflush
# Variables
my @pp_infos = ();

# Species identifiers
my @species = ('Hsap','Mmul','Mfas','Cjac','Mmur','Ltar');
my %spec_id = ();
foreach my $i (0..$#species) {
	$spec_id{$species[$i]} = $i;
}

# Program name
print("\n--- $0 ---\n");
# Collect command line arguments
my $USAGE = "perl $0 <files_list_pp_all_pic.txt> <orthologs.txt>\n";
unless ($ARGV[0] && $ARGV[1]) {
	die("\nUsage: $USAGE\n");
}
# Input files
my $files_list = $ARGV[0];
my $ortho_file = $ARGV[1];

# Get pp genes files
my @files = FileIO::get_file_data_array($files_list);
# Go through each homloc file
foreach my $file (@files) {
	# Get species and genome IDs
	my($spec) = ($file =~ /^(\D+)\d*\_/);
	my $i = $spec_id{$spec};
	# Get pp info
	$pp_infos[$i] = get_target_genes($file);
}

# Get orthologs
my $orthologs = get_othologs($ortho_file);

# Get pp patterns
my %pp_patterns = ();
my %pi_patterns = ();
my %pp_gid_patterns = ();
my %pi_gid_patterns = ();
my @proc_genes = ();
# Go through each gene
foreach my $oid (sort keys %{$orthologs}) {
	my $pp_pattern = '';
	my $pi_pattern = '';
	my $pp_gid_pattern = '';
	my $pi_gid_pattern = '';
	my @pp_gid_pattern = ();
	my @pi_gid_pattern = ();
	# Go through each species
	foreach my $sp_i (0..$#species) {
		if ($sp_i == 2 || $sp_i == 5) { next }
		# Get gene id
		my $g_id = $orthologs->{$oid}->[$sp_i];
		unless ($g_id) { $g_id = '' }
		# Get pp / pi status
		my $pp = 0;
		my $pi = 0;
		$pp = $pp_infos[$sp_i]->{$g_id}->[3] ? $pp_infos[$sp_i]->{$g_id}->[10] : 0;
		$pi = $pp_infos[$sp_i]->{$g_id}->[3] ? 1 : 0;
		$pp_pattern .= $pp;
		$pi_pattern .= $pi;
		my $pp_gid = $pp ? $g_id : '-';
		my $pi_gid = $pi ? $g_id : '-';
		$pp_gid_pattern .= $pp_gid.' ';
		$pi_gid_pattern .= $pi_gid.' ';
		push(@pp_gid_pattern, $pp_gid);
		push(@pi_gid_pattern, $pi_gid);
		# Save processed genes
		$proc_genes[$sp_i]{$g_id}++ if $g_id;
	}
	$pp_gid_patterns{$pp_gid_pattern}++;
	$pi_gid_patterns{$pi_gid_pattern}++;
	$pp_patterns{$pp_pattern}++ if $pp_gid_patterns{$pp_gid_pattern} == 1;
	$pi_patterns{$pi_pattern}++ if $pi_gid_patterns{$pi_gid_pattern} == 1;
}

# Add orphan target genes
my $pat_pos = 0;
my @orph_genes_pp = ();
my @orph_genes_pi = ();
# Go through each species
foreach my $sp_i (0..$#species) {
	if ($sp_i == 2 || $sp_i == 5) { next }
	# Go through each target gene
	foreach my $g_id (sort keys %{$pp_infos[$sp_i]}) {
		# Skip title line
		if ($g_id eq 'Gene_id') { next }
		# Get pp / pi status
		my $pp = 0;
		my $pi = 0;
		$pp = $pp_infos[$sp_i]->{$g_id}->[3] ? $pp_infos[$sp_i]->{$g_id}->[10] : 0;
		$pi = $pp_infos[$sp_i]->{$g_id}->[3] ? 1 : 0;
		# Check if gene already processed in previous step
		unless ($proc_genes[$sp_i]{$g_id}) {
			# Count orphan genes
			$orph_genes_pp[$sp_i] += $pp;
			$orph_genes_pi[$sp_i] += $pi;
		}
	}
	# Get pattern
	my $pattern = '0000';
	substr($pattern, $pat_pos, 1, '1');
	# Save to pattern count
	$pp_patterns{$pattern} += $orph_genes_pp[$sp_i];
	$pi_patterns{$pattern} += $orph_genes_pi[$sp_i];
	# Pattern pos
	$pat_pos++;
}

foreach my $pat (sort {$pp_patterns{$b} <=> $pp_patterns{$a}} keys %pp_patterns) {
	print("$pat\t$pp_patterns{$pat}\n") unless $pat eq '0000';
}
print("\n");
foreach my $pat (sort {$pi_patterns{$b} <=> $pi_patterns{$a}} keys %pi_patterns) {
	print("$pat\t$pi_patterns{$pat}\n") unless $pat eq '0000';
}
print("\n");

exit;

################################# subroutines #################################

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
		my @d = split(/\s/,$line);
		# Get gene id
		my $g_id = $d[0];
		$g_id =~ s/\..*//;
		# Get RPKM
		my $rpkm = $d[2];
		# Save data
		$tab_data{$g_id} = [@d] if $rpkm >= $rpkm_th;
	}
	return \%tab_data;
}

sub get_othologs {
	# Take name of pic pseudogenes file
	my($infile,$key_i) = @_;
	# Key index
	unless ($key_i) { $key_i = 0 }
	# Variables
	my %tab_data = ();
	# Extract file data
	my @file_data = FileIO::get_file_data_array($infile);
	# Parse file data
	foreach my $line (@file_data) {
		# Skip empty line
		if ($line =~ /^$/) { next }
		if ($line =~ /^Gene/) { next }
		# Split tab-separated entries
		my @d = split(/\t/,$line);
		my $key = $d[$key_i];
		shift(@d);
		$tab_data{$key} = [@d];
	}
	return \%tab_data;
}