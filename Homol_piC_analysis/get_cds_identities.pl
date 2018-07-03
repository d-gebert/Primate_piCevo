#!/usr/bin/perl

use strict;
use warnings;
use lib '/Users/dgebert/Dropbox/Perlmodules';
use FileIO;

# Global constants
my $outfile = 'cds_identities.txt';
# Global variables
my @identities = ();
# Options
$|=1; #Autoflush

# Species identifiers
my @species = ('Hsap','Mmul','Mfas','Cjac','Mmur','Ltar');
my %spec_id = ();
foreach my $i (0..$#species) {
	$spec_id{$species[$i]} = $i;
}
# Genome identifiers
my @genomes = ('GRCh38','rheMac8','macFas5','calJac3','micMur2','otoGar3');
my %gnom_id = ();
foreach my $i (0..$#genomes) {
	$gnom_id{$genomes[$i]} = $i;
}

# Program name
print("\n--- $0 ---\n");

# Collect command line arguments
my $USAGE = "perl $0 <list_of_files_blastout>\n";
unless ($ARGV[0]) {
	die("\nUsage: $USAGE\n");
}
my $list_of_files = $ARGV[0];
## Get blast files
my @files = FileIO::get_file_data_array($list_of_files);

## Extract homologous loci from homloc files
# Go through each homloc file
foreach my $file (@files) {
	# Get species and genome IDs
	my($gnom_x) = ($file =~ /^(\w+)_vs_/);
	my($gnom_y) = ($file =~ /_vs_(\w+)\./);
	my $x = $gnom_id{$gnom_x};
	my $y = $gnom_id{$gnom_y};
	# Get total identities
	$identities[$x][$y] = get_total_identity($file);
}

output_identity_matrix(\@identities,$outfile);

exit;

################################# subroutines #################################

sub get_total_identity {
	# File name
	my($blast_file) = @_;
	# Get blast hits per gene and query descriptions
	my($gene_hits,$qgene_infos,$qtran_infos,$qtran_lens) = get_blast_hits_per_gene($blast_file);
	# (Sub-)Global variables
	my $out_suffix = 'gene_hits.txt';
	my $global_ident_x_a_len_sum = 0;
	my $global_a_len_sum = 0;
	# Open output file
	my($prefix) = ($blast_file =~ /(\w+)\./);
	my $outfile = $prefix.'.'.$out_suffix;
	my $out = FileIO::open_outfile($outfile);
	# Go through each gene
	foreach my $qgene_id (sort keys %{$gene_hits}) {
		# Ignore non-protein coding genes
		#if ($qgene_infos->{$qgene_id} !~ /gene_biotype:protein_coding/) { next }
		print($out "$qgene_infos->{$qgene_id}\t");
		# Choose best transcript
		my $max_len = 0;
		my $qtran_best = '';
		foreach my $q_id (keys %{$gene_hits->{$qgene_id}}) {
			# Ignore non-protein coding transcripts
			#if ($qtran_infos->{$q_id} !~ /transcript_biotype:protein_coding/) { next }
			# Find largest transcript
			if ($qtran_lens->{$q_id} > $max_len) {
				$max_len = $qtran_lens->{$q_id};
				$qtran_best = $q_id;
			}
		}
		# Choose best hit transcript
		my $max_score = 0;
		my $stran_best = '';
		foreach my $s_id (sort keys %{$gene_hits->{$qgene_id}->{$qtran_best}}) {
			# Go through each hit and sum subject transcript score
			my $stran_score = 0;
			foreach my $hit (@{$gene_hits->{$qgene_id}->{$qtran_best}->{$s_id}}) {
				# Get hit score
				my @d = split(/\t/,$hit);
				$stran_score += $d[13];
			}
			# Find highest scoring hit transcript
			if ($stran_score > $max_score) {
				$max_score = $stran_score;
				$stran_best = $s_id;
			}
		}
		# Calculate total identity between two homologs
		my $q_cov = 0;
		my %pos_identities = ();
		foreach my $hit (@{$gene_hits->{$qgene_id}->{$qtran_best}->{$stran_best}}) {
			# Get hit data
			my @d = split(/\t/,$hit);
			my $ident = $d[2];
			$q_cov = $d[4];
			# Get identities for each position
			my $q_beg = $d[6];
			my $q_end = $d[7];
			foreach my $pos ($q_beg..$q_end) {
				if ($pos_identities{$pos}) {
					$pos_identities{$pos} = $ident if $ident>$pos_identities{$pos};
				} else {
					$pos_identities{$pos} = $ident;
				}
			}
		}
		# Get alignment length
		my $a_len = scalar(keys %pos_identities);
		# Get total identity for the current gene
		my $homol_identity = 0;
		foreach my $pos (keys %pos_identities) {
			$homol_identity += $pos_identities{$pos}/$a_len;
		}
		# Add to global counts for global total identity
		$global_a_len_sum += $a_len;
		$global_ident_x_a_len_sum += $homol_identity*$a_len;
	
		# Print output
		if ($homol_identity) {
			printf($out "%s\t%s\t%.2f\t%d\t%d\n", $qtran_best,$stran_best,$homol_identity,$q_cov,$a_len);
		} else {
			print($out "N/A\n");
		}
	}
	# Calculate global identity
	my $global_identity = $global_ident_x_a_len_sum/$global_a_len_sum;
	return $global_identity;
}

sub get_blast_hits_per_gene {
	# Take blast file name
	my($blast_file) = @_;
	# Get file data
	my @blast_data = FileIO::get_file_data_array($blast_file);
	# Initialize variables
	my $qtran_id = '';
	my $qgene_info = '';
	my $qtran_info = '';
	my %gene_hits = ();
	my %qgene_infos = ();
	my %qtran_infos = ();
	my %qtran_lens = ();
	# Parse file data
	foreach my $line (@blast_data) {
		# Save query name and id
		if ($line =~ /^#\sQuery/) {
			# Get line data
			my @d = split(/\s/, $line);
			# Save gene info
			$qgene_info = "$d[5] $d[6] $d[8]";
			# Save transcript info
			$qtran_id = $d[2];
			$qtran_info = "$d[2] $d[4] $d[5] $d[6] $d[7] $d[8]";
		}
		# Collect data from match report lines
		elsif ($line !~ /^#/ and $line !~ /^\s*$/) {
			# Get line data
			my @d = split(/\t/, $line);
			my $s_id = $d[1];
			my $ident = $d[2];
			my $q_cov = $d[4];
			my $q_len = $d[5];
			my $s_len = $d[8];
			my $s_start = $d[9];
			my $s_end = $d[10];
			my $s_strand = $d[11];
			my $score = $d[13];
			# Save hit and link to gene id
			if ($s_strand eq 'plus') {
				push(@{$gene_hits{$qgene_info}{$qtran_id}{$s_id}},$line);
				$qgene_infos{$qgene_info} = $qgene_info;
				$qtran_infos{$qtran_id} = $qtran_info;
				$qtran_lens{$qtran_id} = $q_len;
			}
		}
	}
	return \(%gene_hits,%qgene_infos,%qtran_infos,%qtran_lens);
}

sub output_identity_matrix {
	# Take identity matrix and outfile name
	my($idents,$outfile) = @_;
	# Open output file
	$outfile = FileIO::find_unique_filename($outfile);
	my $out = FileIO::open_outfile($outfile);
	# Print table
	foreach my $y (0..$#species) {
		print($out "\t$species[$y]");
	}
	print($out "\n");
	foreach my $x (0..$#species) {
		print($out "$species[$x]\t");
		foreach my $y (0..$#species) {
			if ($x eq $y) { $idents->[$x]->[$y] = 100 }
			printf($out "%.2f\t", $idents->[$x]->[$y]);
		}
		print($out "\n");
	}
}