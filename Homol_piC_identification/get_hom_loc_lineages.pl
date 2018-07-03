#!/usr/bin/perl

use strict;
use warnings;
use lib '/Users/dgebert/Dropbox/Perlmodules';
use FileIO;
use BioStat;

# Global constants
my $outfile = 'homolocs_tree.txt';
my $iterations_max = 10;
# Global variables
my @homlocs = ();
my @tree = ();
my @new_locs = ();
my $br = 0;
my @branches = ();
my %deleted_br = ();
my $n_pics_all_prev = 0;
my $n_pics_uni_prev = 0;
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
my $USAGE = "perl $0 <list_of_files_homlocs>\n";
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
	my($spec_x) = ($file =~ /^(\D+)\_/);
	my($gnom_y) = ($file =~ /\_(\w+)\_homoloc/);
	my $x = $spec_id{$spec_x};
	my $y = $gnom_id{$gnom_y};
	# Get homologous loci
	$homlocs[$x][$y] = get_homlocs($file);
}

## Build initial piC tree
# Go through each x species
foreach my $x (0..$#species) {
	# Go through each y species
	foreach my $y (0..$#species) {
		# Go through each x species loc
		foreach my $loc_x (sort keys %{$homlocs[$x][$y]}) {
			# Increment branch id number
			$br++;
			push(@branches, $br);
			# Save x species loc in current tree branch
			push(@{$tree[$x]{$br}}, $loc_x);
			# Save y species homloc in current tree branch
			my $loc_y = $homlocs[$x][$y]->{$loc_x};
			push(@{$tree[$y]{$br}}, $loc_y) if $loc_y;
		}
	}
}

## Collapse tree
# Get overall and non-redundant piC counts
my($n_pics_all, $n_pics_uni) = count_tree_pics(\@tree, \@branches);
# Repeat until overall and non-redundant counts are equal
while ($n_pics_all != $n_pics_uni) {
	# Go through each tree branch (p)
	foreach my $br_p (@branches) {
		# Skip if branch is deleted
		if ($deleted_br{$br_p}) { next }
		# Go through each following tree branch (q)
		foreach my $br_q (($br_p+1)..scalar(@branches)) {
			# Skip if branch is deleted
			if ($deleted_br{$br_q}) { next }
			# Check if loci in branches p and q overlap, set default: false
			my $loci_overlap = 0;
			my $p_len = 0;
			my $q_len = 0;
			# Go through each species
			foreach my $i (0..$#species) {
				# Get each loc of branch p
				foreach my $loc_p (@{$tree[$i]{$br_p}}) {
					unless ($loc_p) { next }
					# Get loc p coordinates 
					my($chr_p,$beg_p,$end_p) = split(/-|:/,$loc_p);
					# Get each loc of branch q
					foreach my $loc_q (@{$tree[$i]{$br_q}}) {
						# Get loc q coordinates 
						my($chr_q,$beg_q,$end_q) = split(/-|:/,$loc_q);
						# Check if loci overlap
						if ($chr_p eq $chr_q && $beg_q <= $end_p && $end_q >= $beg_p) {
							# Save positions
							my %p_positions = ();
							foreach my $pos ($beg_p..$end_p) {
								$p_positions{$chr_p}{$pos} = 1;
							}
							# Compare positions
							foreach my $pos ($beg_q..$end_q) {
								$loci_overlap++ if $p_positions{$chr_q}{$pos};
							}
							$p_len = $end_p-$beg_p+1;
							$q_len = $end_q-$beg_q+1;
						}
					}
				}
			}
			# If loci overlap add loci in branch q to branch p
			if ($loci_overlap) {
				if ($loci_overlap/$p_len > 0.66 || $loci_overlap/$q_len > 0.66) {
					# Go through each species
					foreach my $i (0..$#species) {
						# Get each loc of branch q
						foreach my $loc_q (@{$tree[$i]{$br_q}}) {
							# Get loc q coordinates 
							my($chr_q,$beg_q,$end_q) = split(/-|:/,$loc_q);
							# Check if q branch loc is included in p branch, set default: false
							my $loc_included = 0;
							# Get each loc of branch p
							foreach my $loc_p (@{$tree[$i]{$br_p}}) {
								# Get loc p coordinates 
								my($chr_p,$beg_p,$end_p) = split(/-|:/,$loc_p);
								# Check if loci overlap
								if ($chr_p eq $chr_q && $beg_q <= $end_p && $end_q >= $beg_p) {
									# Locus included: true
									$loc_included = 1;
									# Locus might need update if loci not identical
									my $loc_needs_update = 0;
									# Loci overlap but are not identical
									if ($loc_p ne $loc_q) {
										# Get new start and end coordinates
										my $beg_new = $beg_p;
										my $end_new = $end_p;
										if ($beg_q < $beg_p) {
											$beg_new = $beg_q;
											# Locus needs update: true
											$loc_needs_update = 1;
										}
										if ($end_q > $end_p) {
											$end_new = $end_q;
											# Locus needs update: true
											$loc_needs_update = 1;
										}
										if ($loc_needs_update) {
											# Define new locus
											my $loc_new = $chr_p.':'.$beg_new.'-'.$end_new;
											# Get index of old locus to replace
											my $n_locs_p = scalar(@{$tree[$i]{$br_p}})-1;
											my($i_old) = grep {$tree[$i]{$br_p}[$_] eq $loc_p} 0..$n_locs_p;
											# Update locus
											$new_locs[$i]{$br_p}{$i_old} = $loc_new;
										}
									}
								}
							}
							# If not already included, add branch q locus to branch p
							unless ($loc_included) {
								push(@{$tree[$i]{$br_p}}, $loc_q) if $tree[$i]{$br_p};
							}
						}
						# Delete branch q
						@{$tree[$i]{$br_q}} = ();
						$deleted_br{$br_q} = 1;
					}
				}
			}
			
		}
		foreach my $i (0..$#species) {
			foreach my $i_old (keys %{$new_locs[$i]{$br_p}}) {
				$tree[$i]{$br_p}[$i_old] = $new_locs[$i]{$br_p}{$i_old};
			}
		}
	}
	# Get new overall and non-redundant piC counts
	$n_pics_all_prev = $n_pics_all;
	$n_pics_uni_prev = $n_pics_uni;
	($n_pics_all, $n_pics_uni) = count_tree_pics(\@tree, \@branches);
	# Abort if overall and non-redundant piC counts stay the same
	print("$n_pics_all\t$n_pics_uni\n");
	if ($n_pics_all == $n_pics_all_prev && $n_pics_uni == $n_pics_uni_prev) { last }
}

## Collapse overlapping loci within a branch and species
# Go through each tree branch
foreach my $br_p (@branches) {
	# Skip if branch is deleted
	if ($deleted_br{$br_p}) { next }
	# Go through each species
	foreach my $i (0..$#species) {
		# Skip if branch is not defined for species
		unless ($tree[$i]{$br_p}) { next }
		# Get locs
		my @update_locs = @{$tree[$i]{$br_p}};
		# Compare all locs
		for (my $a = 0; $a < (scalar(@{$tree[$i]{$br_p}})-1); $a++) {
			for (my $b = ($a+1); $b <= (scalar(@{$tree[$i]{$br_p}})-1); $b++) {
				# Skip obsolete locs
				unless ($update_locs[$a] && $update_locs[$b]) { next }
				# Get each loc p coordinates (a and b)
				my($chr_a,$beg_a,$end_a) = split(/-|:/,$update_locs[$a]);
				my($chr_b,$beg_b,$end_b) = split(/-|:/,$update_locs[$b]);
				# Check overlap
				if ($chr_a eq $chr_b && $beg_b <= $end_a && $end_b >= $beg_a) {
					# Get new start and end coordinates
					my $beg_new = $beg_a;
					my $end_new = $end_a;
					if ($beg_b < $beg_a) { $beg_new = $beg_b }
					if ($end_b > $end_a) { $end_new = $end_b }
					# Update loc (a) and delete loc (b)
					$update_locs[$a] = $chr_a.' '.$beg_new.'-'.$end_new;
					splice(@update_locs, $b, 1);
				}
			}
		}
		@{$tree[$i]{$br_p}} = @update_locs;
	}
}
($n_pics_all, $n_pics_uni) = count_tree_pics(\@tree, \@branches);
print("$n_pics_all\t$n_pics_uni\n");

## Output piC tree
# Sort tree branches by number of homologs
my @sorted_branches = sort_tree_by_homologs(\@tree, \@branches);
# Create plain text representation of piC tree
draw_global_tree_txt(\@tree, \@sorted_branches, 'homlocs_tree.txt');
# Create html heatmap of piC tree
draw_global_tree_html_heatmap(\@tree, \@sorted_branches, 'homlocs_tree_hm.html');

exit;

################################# subroutines #################################

sub get_homlocs {
	# File name
	my($homlocs_file) = @_;
	# Get homlocs data
	my @homlocs_data = FileIO::get_file_data_array($homlocs_file);
	# Initialize variables
	my %homlocs = ();
	# Parse homlocs file
	foreach my $line (@homlocs_data) {
		# Result line
		if ($line !~ /^$/) {
			if ($line =~ /^pic/) { next }
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
			# Ignore unreliable loci
			if ($loc_x =~ /chrUn_/) { next }
			if ($loc_x =~ /random/) { next }
			# Link homologous loci
			if ($legit) {
				$homlocs{$loc_x} = $loc_y;
			} else {
				$homlocs{$loc_x} = '';
			}
		}
	}
	return \%homlocs;
}

sub count_tree_pics {
	# Take piC tree
	my($pictree, $treebranches) = @_;
	# Count variables
	my $n_pics_all = 0;
	my $n_pics_uni = 0;
	my %uni_tree_pics = ();
	# Go through each tree branch
	foreach my $br (@{$treebranches}) {
		# Go through each species
		foreach my $i (0..(scalar(@{$pictree})-1)) {
			# Skip if no piCs
			if ($pictree->[$i]->{$br}) {
				# Add piC counts
				$n_pics_all += scalar(@{$pictree->[$i]->{$br}});
				foreach my $pic (@{$pictree->[$i]->{$br}}) {
					$uni_tree_pics{$pic} = 1;
				}
			}
		}
		# Get non-redundant piC counts
		my @uni_tree_pics = keys %uni_tree_pics;
		$n_pics_uni = @uni_tree_pics;
	}
	return ($n_pics_all, $n_pics_uni);
}

sub sort_tree_by_homologs {
	# Take piC tree and its branches
	my($tree, $branches) = @_;
	# Variables initialization
	my %n_homs = ();
	my @sorted_brs = ();
	# Go through each tree branch
	foreach my $br (@{$branches}) {
		# Get number of branch homologs
		my $homs = 0;
		foreach my $i (0..(scalar(@{$tree})-1)) {
			if ($tree->[$i]->{$br}) {
				if (scalar(@{$tree->[$i]->{$br}})>0) { $homs++ }
			}
		}
		$n_homs{$br} = $homs;
	}
	# Get branch list sorted by number of homologs
	foreach my $br (sort { $n_homs{$b} <=> $n_homs{$a} } keys %n_homs) {
		push(@sorted_brs,$br) if $n_homs{$br};
	}
	return @sorted_brs;
}

sub draw_global_tree_txt {
	# Take piC tree and its branches
	my($tree, $branches, $outfile) = @_;
	# Open output file
	$outfile = FileIO::find_unique_filename($outfile);
	my $out = FileIO::open_outfile($outfile);
	# Go through each global tree branch
	foreach my $br (@{$branches}) {
		# Get max array size for the groups in this branch
		my @sizes = ();
		foreach my $i (0..(scalar(@{$tree})-1)) {
			push(@sizes,scalar(@{$tree->[$i]->{$br}})) if $tree->[$i]->{$br};
		}
		my $size_max = BioStat::get_maximum(\@sizes);
		# Skip empty branch
		if ($size_max == 0) { next }
		# Go through each putative array position
		foreach my $j (0..($size_max-1)) {
			# If this position in the group array exists print the element
			foreach my $i (0..(scalar(@{$tree})-1)) {
				if (${$tree->[$i]->{$br}}[$j]) {
					print($out "${$tree->[$i]->{$br}}[$j]\t");
				} else {
					print($out "\t");
				}
			}
			print($out "\n");
		}
		print($out "\n");
	}
}

sub draw_global_tree_html_heatmap {
	# Take piC tree and its branches
	my($tree, $branches, $outfile) = @_;
	# Option for skipping single piCs
	my $skip_single_clusters = 0;
	# Open output file
	$outfile = FileIO::find_unique_filename($outfile);
	my $out = FileIO::open_outfile($outfile);
	# Parameters
	my $width = 200;
	my $height = 500;
	my $cell_spacing = 0;
	my $cell_border = 0;
	# RGB colors
	my @no_pics = (0,0,0);
	my @one_pic = (3,101,192);
	my @many_pics = (81,167,249);
	# Print html header
	print($out "\n<table width=\"${width}px\" height=\"${height}px\" border=\"${cell_border}px\" cellspacing=\"${cell_spacing}px\">");
	
	# Go through each global tree branch
	foreach my $br (@{$branches}) {
		# Skip single clusters if option true
		if ($skip_single_clusters == 1) {
			# Get number of branch homologs
			my $homs = 0;
			foreach my $i (0..(scalar(@{$tree})-1)) {
				if ($tree->[$i]->{$br}) {
					if (scalar(@{$tree->[$i]->{$br}})>0) { $homs++ }
				}
			}
			if ($homs <= 1) { next }
		}
		# Initialize RGB
		my @rgb = (0,0,0);
		# Draw heatmap fields for current branch
		print($out "\n<tr>");
		foreach my $i (0..(scalar(@{$tree})-1)) {
			# Choose field color according to number of grouped intra-species pics
			if (scalar(@{$tree->[$i]->{$br}}) == 1) { @rgb = @one_pic }
			if (scalar(@{$tree->[$i]->{$br}}) >= 2) { @rgb = @many_pics }
			if (scalar(@{$tree->[$i]->{$br}}) == 0) { @rgb = @no_pics }
			print($out "<td style=\"background-color:rgb($rgb[0],$rgb[1],$rgb[2]);\"></td>");
		}
		print($out "</tr>");
	}
	print($out "</table>");
}