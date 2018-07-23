#!/usr/bin/perl

use strict;
use warnings;
use lib '/Users/dgebert/Dropbox/Perlmodules/';
use FileIO;
use FastaIO;

# Constants
my $min_len = 24;
my $max_len = 32;
$|=1; #Autoflush
# Variables
my %tar_gene_hits = ();
my %tar_gene_syms = ();
my %hits_per_target = ();
my %tar_gene_syms_all = ();

# Collect command line arguments
my $USAGE = "perl $0 <piCpseudogenes.txt> <pic_reads_folder> <cdna.fas>\n";
unless ($ARGV[0] && $ARGV[1] && $ARGV[2]) {
	die("\nUsage: $USAGE\n");
}
# Input files
my $pic_pseudos_file = $ARGV[0];
my $pic_reads_folder = $ARGV[1];
my $coding_cdna_file = $ARGV[2];

## Get pic pseudogenes
my($pic_info,$pic_pgenes) = get_pic_pseudogenes($pic_pseudos_file);
## Get pic reads
my $pic_reads = get_pic_reads($pic_reads_folder);

# Get species name
my($sample_id) = ($pic_pseudos_file =~ /([^_]+)_/);
# Name output file
my $outfile_tg = $sample_id.'_pic_target_genes.txt';
my $out_tg = FileIO::open_outfile($outfile_tg);
# Name folder for map files
my $map_files_dir = $sample_id.'_mapfiles';
mkdir($map_files_dir) unless (-d $map_files_dir);

## Go through each pseudogene-containing piC
foreach my $pic (sort keys %{$pic_pgenes}) {
	## Get piRNA reads from piC
	# Get pic ids (if merged pics)
	my @pic_ids = split(/_/,$pic);
	shift(@pic_ids);
	# Go through each pseudogene
	foreach my $pg_id (sort keys %{$pic_pgenes->{$pic}}) {
		# Name reads output file
		my $reads_outfile = $pic.'.'.$pg_id.'_PG.hits.fas';
		my $out = FileIO::open_outfile($reads_outfile);
		# Bool: Pseudogene is in parallel orientation
		my $pg_parallel = 1;
		# Go through each pseudogene exon
		my $pgp_sym = '';
		foreach my $pg_exon (@{$pic_pgenes->{$pic}->{$pg_id}}) {
			# Get pseudogene exon coordinates
			my $pge_beg = $pg_exon->[7];
			my $pge_end = $pg_exon->[8];
			my $pge_str = $pg_exon->[9];
			my $pge_dir = $pg_exon->[11];
			$pgp_sym = $pg_exon->[1];
			# Skip if pseudogene orientation is not reverse to pic expression
			if ($pge_dir =~ /par/) { next }
			if ($pge_dir =~ /rev/) { $pg_parallel = 0 }
			# Go through each pic id (if merged pics)
			foreach my $pic_id (@pic_ids) {
				# Go through each pic hit
				foreach my $hit_id (sort {$a <=> $b} keys %{$pic_reads->{$pic_id}}) {
					my $pos   = $pic_reads->{$pic_id}->{$hit_id}->[1];
					my $reads = $pic_reads->{$pic_id}->{$hit_id}->[4];
					my $seq   = $pic_reads->{$pic_id}->{$hit_id}->[6];
					# Print hit to file if (fully) located in exon
					if ($pos >= $pge_beg && ($pos+length($seq)-1) <= $pge_end) {
						print($out ">$reads\n$seq\n");
					}
				}
			}
		}
		# Skip pseudogene if entirely in parallel orientation
		if ($pg_parallel) { unlink($reads_outfile); next; }
		## Map reads to gene set
		my $map_outfile = $map_files_dir.'/'.$pic.'.'.$pg_id.'_PG.hits.cdna.map';
		unless (-e $map_outfile) {
			system("./seqmap 2 $reads_outfile $coding_cdna_file $map_outfile /output_all_matches >.stdout");
		}
		unlink($reads_outfile);
		#unlink('.stdout');
		## Get target gene names
		# Get mapfile data if map file exists
		unless (-e $map_outfile) { next }
		my @map_data = FileIO::get_file_data_array($map_outfile);
		# Parse file data
		foreach my $line (@map_data) {
			# Extract line data
			my @d = split(/\t/,$line);
			my $gen_inf = $d[0];
			my $hit_pos = $d[1];
			my $hit_seq = $d[2];
			my $hit_cnt = $d[3];
			my $pir_seq = $d[4];
			my $hit_mms = $d[5];
			my $hit_str = $d[6];
			my $hit_len = length($hit_seq);
			my($gen_eid) = ($gen_inf =~ /gene:(\S+)/);
			my($gen_sym) = ($gen_inf =~ /gene_symbol:(\S+)/);
			unless ($gen_sym) { $gen_sym = $gen_eid }
			# Skip if hit is not on minus strand
			unless ($hit_str eq '-') { next }
			unless ($hit_len >= $min_len && $hit_len <= $max_len) { next }
			# Save target
			$tar_gene_hits{$pic}{$pg_id}{$gen_eid} += $hit_cnt;
			$tar_gene_syms{$pic}{$pg_id}{$gen_eid} = $gen_eid."\t".$gen_sym;
			$tar_gene_syms_all{$gen_eid} = $gen_eid."\t".$gen_sym;
		}
		print($out_tg "#$pic\t$pg_id\t$pgp_sym\n");
		# Go through each target gene
		foreach my $tar_id (sort {$tar_gene_hits{$pic}{$pg_id}{$b} <=> $tar_gene_hits{$pic}{$pg_id}{$a} } keys %{$tar_gene_hits{$pic}{$pg_id}}) {
			# Get sum hits per target
			$hits_per_target{$tar_id} += $tar_gene_hits{$pic}{$pg_id}{$tar_id};
			
			unless ($tar_gene_hits{$pic}{$pg_id}{$tar_id} >= 1) { next }
			printf($out_tg "%s\t%.2f\n", $tar_gene_syms{$pic}{$pg_id}{$tar_id},$tar_gene_hits{$pic}{$pg_id}{$tar_id});
		}
		print($out_tg "\n");
	}
	## Check if piC target genes have ping-pong signature
}

my $outfile_tg2 = $sample_id.'_all_target_genes.txt';
my $out_tg2 = FileIO::open_outfile($outfile_tg2);
# Go through each target gene
foreach my $tar_id (sort {$hits_per_target{$b} <=> $hits_per_target{$a} } keys %hits_per_target) {
	printf($out_tg2 "%s\t%.2f\n", $tar_gene_syms_all{$tar_id},$hits_per_target{$tar_id});
}

exit;

################################# subroutines #################################

sub get_pic_pseudogenes {
	# Take name of pic pseudogenes file
	my($infile) = @_;
	# Variables
	my %pic_info = ();
	my %pic_pgenes = ();
	my $pic_id = '';
	# Extract file data
	my @file_data = FileIO::get_file_data_array($infile);
	# Parse file data
	foreach my $line (@file_data) {
		# Get pic info
		if ($line =~ /^#/) {
			my @d = split(/\t/,$line);
			$pic_id = $d[0];
			$pic_id =~ s/#//;
			my $chr = $d[1];
			my @loc = split(/-/,$d[2]);
			my $beg = $loc[0];
			my $end = $loc[1];
			my $dir = $d[3];
			$pic_info{$pic_id} = [$chr,$beg,$end,$dir];
		} 
		# Get pseudogene
		elsif ($line !~ /^$/) {
			my @d = split(/\t/,$line);
			my $parent_id = $d[0];
			my $parent_gs = $d[1];
			my $pgene_ide = $d[4];
			my $pgene_beg = $d[7];
			my $pgene_end = $d[8];
			my $pgene_str = $d[9];
			my $pgene_dir = $d[11];
			push(@{$pic_pgenes{$pic_id}{$parent_id}},[@d]);
		}
	}
	return \(%pic_info,%pic_pgenes);
}

sub get_pic_reads {
	# Take name of pic reads folder
	my($folder) = @_;
	# Storage variable
	my %pic_reads = ();
	# Get file names
	my @files = FileIO::get_folder_file_names($folder);
	# Extract reads from each file and save with pic id
	foreach my $file (@files) {
		# Get file data
		my @file_data = FileIO::get_file_data_array($folder.'/'.$file);
		# Get pic id
		my($pic_id) = ($file =~ /(\d+)\./);
		# Initialize seq id
		my $hit_id = 0;
		# Parse file data
		foreach my $line (@file_data) {
			# Get read info
			if ($line =~ /^>/) {
				my @d = split(/\t/,$line);
				my($chr)   = ($d[0] =~ /:(.*)/);
				my($pos)   = ($d[1] =~ /:(.*)/);
				my($reads) = ($d[2] =~ /:(.*)/);
				my($hits)  = ($d[3] =~ /:(.*)/);
				my($alloc) = ($d[4] =~ /:(.*)/);
				my($str)   = ($d[5] =~ /:(.*)/);
				# Save data
				@{$pic_reads{$pic_id}{$hit_id}} = ($chr,$pos,$reads,$hits,$alloc,$str);
			}
			# Get read sequence
			elsif ($line !~ /^$/) {
				my $seq = $line;
				# Save data
				push(@{$pic_reads{$pic_id}{$hit_id}},$seq);
				$hit_id++;
			}
		}
	}
	return \%pic_reads;
}

sub HammingDistance {
	my($p, $q) = @_;
	my $hamming_dist = ($p ^ $q) =~ tr/\0//c;
	return $hamming_dist;
}
