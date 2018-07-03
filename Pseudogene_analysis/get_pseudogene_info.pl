#!/usr/bin/perl
use strict;
use warnings;
use lib '/Users/dgebert/Dropbox/Perlmodules/';
use FileIO;
use FastaIO;
use BioStat;

# Constants
$|=1; #Autoflush
# Variables
my %parents_gnm = ();
my %pgtypes_gnm = ();
my @identys_gnm = ();
my %parents_pic = ();
my %pgtypes_pic = ();
my @identys_pic = ();
my %ide_pic_dir = ();
my $pgene_bps_gnm = 0;
my $pgene_bps_pic = 0;

# Collect command line arguments
my $USAGE = "perl $0 <genome.gff> <piCpseudogenes.txt>\n";
unless ($ARGV[0] && $ARGV[1]) {
	die("\nUsage: $USAGE\n");
}
# Input files
my $gene_gff_file = $ARGV[0];
my $pic_pgen_file = $ARGV[1];

# Get gff gene and cds locations
my($pgene_locs,$pexon_locs) = get_gff_gene_locs($gene_gff_file);
# Get pic pseudogenes
my($pic_info,$pic_pgenes) = get_pic_pseudogenes($pic_pgen_file);

# Go through each chromosome
foreach my $chr (sort keys %{$pgene_locs}) {
	# Go through each pgene on chromosome
	foreach my $pgene (sort { $pgene_locs->{$chr}->{$a}->[3] <=> $pgene_locs->{$chr}->{$b}->[3]} keys %{$pgene_locs->{$chr}}) {
		# Gene pgene info
		my $gen_sym = $pgene_locs->{$chr}->{$pgene}->[-1];
		my $gen_ide = $pgene_locs->{$chr}->{$pgene}->[9];
		my $gen_inf = $pgene_locs->{$chr}->{$pgene}->[8];
		my($gen_typ) = ($gen_inf =~ /gene_biotype=([^;,]+)/);
		# Count parent gene offspring
		$parents_gnm{$gen_sym}++;
		$pgtypes_gnm{$gen_typ}++;
		# Save identities
		push(@identys_gnm,$gen_ide);
		# Count pseudogenic base pairs
		foreach my $pexon (@{$pexon_locs->{$chr}->{$pgene}}) {
			my $pex_beg = $pexon->[3];
			my $pex_end = $pexon->[4];
			my $pex_len = $pex_end-$pex_beg+1;
			$pgene_bps_gnm += $pex_len;
		}
	}
}

# Go through each pseudogene-containing piC
foreach my $pic (sort keys %{$pic_pgenes}) {
	# Go through each pseudogene
	foreach my $pgene (sort keys %{$pic_pgenes->{$pic}}) {
		# Get pgene info
		my $gen_sym = $pic_pgenes->{$pic}->{$pgene}->[0]->[1];
		my $gen_ide = $pic_pgenes->{$pic}->{$pgene}->[0]->[4];
		my $gen_typ = $pic_pgenes->{$pic}->{$pgene}->[0]->[10];
		my $gen_dir = $pic_pgenes->{$pic}->{$pgene}->[0]->[11];
		# Count parent gene offspring
		$parents_pic{$gen_sym}++;
		$pgtypes_pic{$gen_typ}++;
		# Save identities
		push(@identys_pic,$gen_ide);
		# Save identities per directionality
		push(@{$ide_pic_dir{$gen_dir}},$gen_ide);
		# Count pseudogenic base pairs
		foreach my $pexon (@{$pic_pgenes->{$pic}->{$pgene}}) {
			my $pex_beg = $pexon->[7];
			my $pex_end = $pexon->[8];
			my $pex_len = $pex_end-$pex_beg+1;
			$pgene_bps_pic += $pex_len;
		}
	}
}

## Output
my($genome_name) = ($gene_gff_file =~ /(\w+)\./);
my $outfile_sum = $genome_name.'.piCpseudo_info.txt';
my $out_sum = FileIO::open_outfile($outfile_sum);

# Output pseudogenic base pair count
print($out_sum "\tGenome\tpiClusters\n");
print($out_sum "Pgene bps\t$pgene_bps_gnm\t$pgene_bps_pic\n\n");

print($out_sum "\tGenome\tpiClusters\tpar_dir\trev_dir\n");
# Output median identities
my $med_ide_gnm = BioStat::get_median(\@identys_gnm);
my $med_ide_pic = BioStat::get_median(\@identys_pic);
my $med_ide_par = BioStat::get_median($ide_pic_dir{'par'});
my $med_ide_rev = BioStat::get_median($ide_pic_dir{'rev'});
printf($out_sum "Median identity\t%.2f\t%.2f\t%.2f\t%.2f\n", $med_ide_gnm,$med_ide_pic,$med_ide_par,$med_ide_rev);

# Output mean identities
my $mea_ide_gnm = BioStat::get_mean(\@identys_gnm);
my $mea_ide_pic = BioStat::get_mean(\@identys_pic);
my $mea_ide_par = BioStat::get_mean($ide_pic_dir{'par'});
my $mea_ide_rev = BioStat::get_mean($ide_pic_dir{'rev'});
printf($out_sum "Mean identity\t%.2f\t%.2f\t%.2f\t%.2f\n", $mea_ide_gnm,$mea_ide_pic,$mea_ide_par,$mea_ide_rev);

# Output pseudogene counts
my $pgs_cnt_par = @{$ide_pic_dir{'par'}};
my $pgs_cnt_rev = @{$ide_pic_dir{'rev'}};
my $pgs_cnt_tot = @identys_pic;
my $gnm_cnt_tot = @identys_gnm;
print($out_sum "Pseudog. count\t$gnm_cnt_tot\t$pgs_cnt_tot\t$pgs_cnt_par\t$pgs_cnt_rev\n\n");

# Output pseudogene type counts
print($out_sum "\tGenome\tpiClusters\n");
foreach my $type (sort { $pgtypes_gnm{$b} <=> $pgtypes_gnm{$a} } keys %pgtypes_gnm) {
	print($out_sum "$type\t$pgtypes_gnm{$type}\t$pgtypes_pic{$type}\n");
}
print($out_sum "\n");

# Output parent gene counts
print($out_sum "\tGenome\tpiClusters\n");
foreach my $parent (sort { $parents_gnm{$b} <=> $parents_gnm{$a} } keys %parents_gnm) {
	unless ($parents_pic{$parent}) { $parents_pic{$parent} = 0 }
	print($out_sum "$parent\t$parents_gnm{$parent}\t$parents_pic{$parent}\n");
}
print($out_sum "\n");

## Output identity array
my $outfile_ide_gnm = $genome_name.'.GnmPseudo_idents_all.txt';
my $out_ide_gnm = FileIO::open_outfile($outfile_ide_gnm);
foreach my $ide (@identys_gnm) {
	print($out_ide_gnm "$ide\n");
}

my $outfile_ide_pic = $genome_name.'.piCpseudo_idents_all.txt';
my $out_ide_pic = FileIO::open_outfile($outfile_ide_pic);
foreach my $ide (@identys_pic) {
	print($out_ide_pic "$ide\n");
}

my $outfile_ide_par = $genome_name.'.piCpseudo_idents_par.txt';
my $out_ide_par = FileIO::open_outfile($outfile_ide_par);
foreach my $ide (@{$ide_pic_dir{'par'}}) {
	print($out_ide_par "$ide\n");
}

my $outfile_ide_rev = $genome_name.'.piCpseudo_idents_rev.txt';
my $out_ide_rev = FileIO::open_outfile($outfile_ide_rev);
foreach my $ide (@{$ide_pic_dir{'rev'}}) {
	print($out_ide_rev "$ide\n");
}

## Output genomic pseudogene occurence per parentgene for whole genome
my $outfile_pgp_gnm = $genome_name.'.GnmPseudo_parents.txt';
my $out_pgp_gnm = FileIO::open_outfile($outfile_pgp_gnm);
foreach my $parent (sort { $parents_gnm{$b} <=> $parents_gnm{$a} } keys %parents_gnm) {
	print($out_pgp_gnm "$parents_gnm{$parent}\n");
}

## Output genomic pseudogene occurence per parentgene for pics
my $outfile_pgp_pic = $genome_name.'.piCpseudo_parents.txt';
my $out_pgp_pic = FileIO::open_outfile($outfile_pgp_pic);
foreach my $parent (sort { $parents_pic{$b} <=> $parents_pic{$a} } keys %parents_pic) {
	print($out_pgp_pic "$parents_gnm{$parent}\n") if $parents_pic{$parent}; # x $parents_pic{$parent}
}

exit;

################################# subroutines #################################

sub get_gff_gene_locs {
	# Take infile name and keyword
	my($gff_file) = @_;
	# Get file data
	my $gff = FileIO::open_infile($gff_file);
	# Initialize variables
	my %ggene_locs = ();
	my %pexon_locs = ();
	my $gene = '';
	my $name = '';
	# Parse gff file
	while (my $line = <$gff>) {
		$line =~ s/\s+$//; #better chomp
		# Results line
		if ($line =~ /^[^#]\w+/) {
			# Get info
			my @d = split(/\t/, $line);
			my $chr  = $d[0];
			my $type = $d[2];
			my $beg  = $d[3];
			my $end  = $d[4];
			my $info = $d[8];
			my $iden = $d[9] ? $d[9] : 0;
			$d[9] = $iden;
			# Gene line
			if ($type =~ /gene/i) {
				# Skip non pseudogenic
				if ($info !~ /pseudogene/i) { next }
				# Get gene id
				# GFF
				if ($info =~ /GeneID:([^;,]+)/) {
					($gene) = ($info =~ /GeneID:([^;,]+)/);
				}
				# GTF
				if ($info =~ /gene_id\s\"([^\"]+)\";/) {
					($gene) = ($info =~ /gene_id\s\"([^\"]+)\";/);
				}
				# Get gene name 
				# GFF
				if ($info =~ /Name=([^;,]+)/) {
					($name) = ($info =~ /Name=([^;,]+)/);
				}
				# GTF
				if ($info =~ /gene_name\s"([^\"]+)"/) {
					($name) = ($info =~ /gene_name\s"([^\"]+)"/);
				}
				# If pseudogenes are predicted by Pseudon, add unit id
				if ($d[1] eq 'Pseudon') {
					my($id) = ($info =~ /ID=(\d+);/);
					$gene = $gene.'_'.$id;
				}
				# Save entries
				push(@d,$name);
				$ggene_locs{$chr}{$gene} = \@d;
			}
			# Exon line
			elsif ($type =~ /exon/i) {
				# Save entries
				push(@d,$name);
				push(@{$pexon_locs{$chr}{$gene}},\@d);
			}
		}
	}
	return \(%ggene_locs,%pexon_locs);
}

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
			my($beg,$end) = split(/-/,$d[2]);
			my $dir = $d[3];
			$pic_info{$pic_id} = [$chr,$beg,$end,$dir];
		} 
		# Get pseudogene
		elsif ($line !~ /^$/) {
			my @d = split(/\t/,$line);
			my $parent_id = $d[0];
			push(@{$pic_pgenes{$pic_id}{$parent_id}},\@d);
		}
	}
	return \(%pic_info,%pic_pgenes);
}