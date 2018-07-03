#!/usr/bin/perl

use strict;
use warnings;
use lib '/Users/dgebert/Dropbox/Perlmodules';
use FileIO;

# Global constants
my $outfile_iden = 'homloc_identities.txt';
my $outfile_idex = 'homloc_identities_ex.txt';
my $outfile_synt = 'homloc_syntenies.txt';
my $outfile_homl = 'homloc_homologies.txt';
my $outfile_expr = 'homloc_expressed_homs.txt';
# Global variables
my @identities = ();
my @identities_expr = ();
my @syntenies = ();
my @homologies = ();
my @homexpress = ();
my $pic_loci;
# Options
$|=1; #Autoflush

# Species identifiers
my @species = ('Hsap','Mmul','Mfas','Cjac','Mmur','Ltar');
my %spec_id = ();
foreach my $i (0..$#species) {
	$spec_id{$species[$i]} = $i;
}
# Genome identifiers
my @genomes = ('GRCh38','rheMac8','macFas5','calJac3','micMur3','otoGar3');
my %gnom_id = ();
foreach my $i (0..$#genomes) {
	$gnom_id{$genomes[$i]} = $i;
}

# Program name
print("\n--- $0 ---\n");

# Collect command line arguments
my $USAGE = "perl $0 <list_of_files_homlocs> <piC_loci_list.txt>\n";
unless ($ARGV[0]) {
	die("\nUsage: $USAGE\n");
}
my $list_of_files = $ARGV[0];
my $picloci_file = $ARGV[1];
## Get blast files
my @files = FileIO::get_file_data_array($list_of_files);
# Get piC loci if given
if ($ARGV[1]) {
	$pic_loci = get_pic_loci($picloci_file);
}

## Extract homologous loci from homloc files
# Go through each homloc file
foreach my $file (@files) {
	# Get species and genome IDs
	my($spec_x) = ($file =~ /^(\D+)\_/);
	my($gnom_y) = ($file =~ /\_(.+)\_/);
	my $x = $spec_id{$spec_x};
	my $y = $gnom_id{$gnom_y};
	# Get total identities
	$identities[$x][$y] = get_total_identity($file);
	if ($ARGV[1]) {
		$identities_expr[$x][$y] = get_total_identity_expr($file,$pic_loci,$y);
	}
	$syntenies[$x][$y] = get_total_synteny($file);
	$homologies[$x][$y] = get_total_homology($file);
	$homexpress[$x][$y] = get_total_homexpr($file);
}

output_matrix(\@identities,$outfile_iden);
output_matrix(\@identities_expr,$outfile_idex) if ($ARGV[1]);
output_matrix(\@syntenies,$outfile_synt);
output_matrix(\@homologies,$outfile_homl);
output_matrix(\@homexpress,$outfile_expr);

exit;

################################# subroutines #################################

sub get_total_identity {
	# File name
	my($homlocs_file) = @_;
	# Get homlocs data
	my @homlocs_data = FileIO::get_file_data_array($homlocs_file);
	# Initialize variables
	my $len_align_sum = 0;
	my $len_align_x_ident_sum = 0;
	# Parse homlocs file
	foreach my $line (@homlocs_data) {
		# Result line
		if ($line !~ /^pic/) {
			my @d = split(/\t/,$line);
			# Get fields
			my $pic_x = $d[1];
			my $loc_x = $d[2];
			my $len_x = $d[3];
			my $loc_y = $d[4];
			my $len_y = $d[5];
			my $relen = $d[6];
			my $qcovg = $d[7];
			my $ident = $d[8];
			my $score = $d[9];
			my $legit = $d[10];
			# Include homologous loci in calculation
			if ($legit) {
				# Get alignment length
				my $len_align = $len_x*($qcovg/100);
				# Multiply identity with alignment length
				my $len_align_x_ident = $ident*$len_align;
				# Sum alignment length and product
				$len_align_sum += $len_align;
				$len_align_x_ident_sum += $len_align_x_ident;
			}
		}
	}
	my $total_identity = $len_align_x_ident_sum/$len_align_sum;
	return $total_identity;
}

sub get_total_synteny {
	# File name
	my($homlocs_file) = @_;
	# Get homlocs data
	my @homlocs_data = FileIO::get_file_data_array($homlocs_file);
	# Initialize variables
	my $n_pics = 0;
	my $n_synt = 0;
	# Parse homlocs file
	foreach my $line (@homlocs_data) {
		# Result line
		if ($line !~ /^pic/) {
			my @d = split(/\t/,$line);
			# Get fields
			my $pic_x = $d[1];
			my $loc_x = $d[2];
			my $len_x = $d[3];
			my $loc_y = $d[4];
			my $len_y = $d[5];
			my $relen = $d[6];
			my $qcovg = $d[7];
			my $ident = $d[8];
			my $score = $d[9];
			my $legit = $d[10];
			my $s_loc = $d[11];
			my $synty = $d[12];
			# Ignore unreliable loci
			if ($loc_x =~ /chrUn_/) { next }
			if ($loc_x =~ /random/) { next }
			# Include homologous loci in calculation
			$n_pics++;
			if ($synty) { $n_synt++ }
		}
	}
	my $total_synteny = $n_synt/$n_pics*100;
	return $total_synteny;
}

sub get_total_homology {
	# File name
	my($homlocs_file) = @_;
	# Get homlocs data
	my @homlocs_data = FileIO::get_file_data_array($homlocs_file);
	# Initialize variables
	my $n_pics = 0;
	my $n_homl = 0;
	# Parse homlocs file
	foreach my $line (@homlocs_data) {
		# Result line
		if ($line !~ /^pic/) {
			my @d = split(/\t/,$line);
			# Get fields
			my $pic_x = $d[1];
			my $loc_x = $d[2];
			my $len_x = $d[3];
			my $loc_y = $d[4];
			my $len_y = $d[5];
			my $relen = $d[6];
			my $qcovg = $d[7];
			my $ident = $d[8];
			my $score = $d[9];
			my $legit = $d[10];
			my $s_loc = $d[11];
			my $synty = $d[12];
			# Ignore unreliable loci
			if ($loc_x =~ /chrUn_/) { next }
			if ($loc_x =~ /random/) { next }
			# Include homologous loci in calculation
			$n_pics++;
			if ($legit) { $n_homl++ }
		}
	}
	my $total_homology = $n_homl/$n_pics*100;
	return $total_homology;
}

sub get_total_homexpr {
	# File name
	my($homlocs_file) = @_;
	# Get homlocs data
	my @homlocs_data = FileIO::get_file_data_array($homlocs_file);
	# Initialize variables
	my $n_pics = 0;
	my $n_expr = 0;
	# Parse homlocs file
	foreach my $line (@homlocs_data) {
		# Result line
		if ($line !~ /^pic/) {
			my @d = split(/\t/,$line);
			# Get fields
			my $pic_x = $d[1];
			my $loc_x = $d[2];
			my $len_x = $d[3];
			my $loc_y = $d[4];
			my $len_y = $d[5];
			my $relen = $d[6];
			my $qcovg = $d[7];
			my $ident = $d[8];
			my $score = $d[9];
			my $legit = $d[10];
			my $s_loc = $d[11];
			my $synty = $d[12];
			my $hompi = $d[13];
			# Ignore unreliable loci
			if ($loc_x =~ /chrUn_/) { next }
			if ($loc_x =~ /random/) { next }
			# Include homologous loci in calculation
			$n_pics++;
			if ($hompi) { $n_expr++ }
		}
	}
	my $total_homexpr = $n_expr/$n_pics*100;
	return $total_homexpr;
}

sub get_total_identity_expr {
	# File name
	my($homlocs_file,$pic_loci,$sp) = @_;
	# Get homlocs data
	my @homlocs_data = FileIO::get_file_data_array($homlocs_file);
	# Initialize variables
	my $len_align_sum = 0;
	my $len_align_x_ident_sum = 0;
	# Parse homlocs file
	foreach my $line (@homlocs_data) {
		# Result line
		if ($line !~ /^$/) {
			my @d = split(/\t/,$line);
			# Get fields
			my $pic_x = $d[0];
			my $loc_x = $d[1];
			my $len_x = $d[2];
			my $loc_y = $d[3];
			my $len_y = $d[4];
			my $relen = $d[5];
			my $qcovg = $d[6];
			my $ident = $d[7];
			my $score = $d[8];
			my $legit = $d[9];
			# Include homologous loci in calculation
			if ($legit) {
				# Go though expressed pics list
				foreach my $pic (keys %{$pic_loci}) {
					# Get species
					my($spec) = ($pic =~ /^([^-_]+)/);
					#print("$spec_id{$spec}\t$sp\n");
					unless ($spec_id{$spec} eq $sp) { next }
					# Get loc_pic coordinates
					my $loc_pic = $pic_loci->{$pic}->[0];
					my @loc_pic = split(/-|:/,$loc_pic);
					my $chr_pic = $loc_pic[0];
					my $beg_pic = $loc_pic[1];
					my $end_pic = $loc_pic[2];
					# Get loc_y coordinates
					my @loc_y = split(/-|:/,$loc_y);
					my $chr_y = $loc_y[0];
					my $beg_y = $loc_y[1];
					my $end_y = $loc_y[2];
					# Check if loci overlap
					if ($chr_y eq $chr_pic && $beg_pic <= $end_y && $end_pic >= $beg_y) {
						# Get alignment length
						my $len_align = $len_x*($qcovg/100);
						# Multiply identity with alignment length
						my $len_align_x_ident = $ident*$len_align;
						# Sum alignment length and product
						$len_align_sum += $len_align;
						$len_align_x_ident_sum += $len_align_x_ident;
					}
				}
			}
		}
	}
	my $total_identity = $len_align_x_ident_sum/$len_align_sum if $len_align_sum;
	return $total_identity;
}

sub get_pic_loci {
	# Take pic loci file name
	my($infile) = @_;
	my %pic_locs = ();
	# Get file data
	my @in_data = FileIO::get_file_data_array($infile);
	# Parse file data
	foreach my $line (@in_data) {
		if ($line !~ /^$/) {
			# Get piC locus data
			my @d = split(/\t/,$line);
			my $pic  = $d[0];
			my $loc  = $d[1];
			my $size = $d[2];
			my $hits = $d[3];
			my $dens = $d[4];
			my $dir  = $d[5];
			@{$pic_locs{$pic}} = ($loc,$size,$hits,$dens,$dir);
		}
	}
	return \%pic_locs;
}

sub output_matrix {
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
			if ($x eq $y) { 
				$idents->[$x]->[$y] = '';
				printf($out "%s\t", $idents->[$x]->[$y]);
			} else {
				printf($out "%.2f\t", $idents->[$x]->[$y]);
			}
		}
		print($out "\n");
	}
}