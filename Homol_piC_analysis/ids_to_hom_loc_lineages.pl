#!/usr/bin/perl

use strict;
use warnings;
use lib '/Users/dgebert/Dropbox/Perlmodules';
use FileIO;
use BioStat;
use GradientRGB;

# Global constants
my $outfile_expr = "homlocs_tree_piCexpr.txt";
# Global variables
my @syn_tree = ();
my @hits_min = (1000000,1000000,1000000,1000000,1000000,1000000);
my @hits_max = (0,0,0,0,0,0);
# Options
$|=1; #Autoflush
my $limit_max_hits = 0;

# Species identifiers
my @species = ('Hsap','Mmul','Mfas','Cjac','Mmur','Ltar');
my %spec_id = ();
foreach my $i (0..$#species) {
	$spec_id{$species[$i]} = $i;
}

# Program name
print("\n--- $0 ---\n");

# Collect command line arguments
my $USAGE = "perl $0 <homlocs_tree.txt> <piC_loci_list.txt>\n";
unless ($ARGV[0]&&$ARGV[1]) {
	die("\nUsage: $USAGE\n");
}
# Input files
my $homtree_file = $ARGV[0];
my $picloci_file = $ARGV[1];

## Extract date
# Get homlocs tree
my($hom_tree, $hom_branches) = get_tree_data($homtree_file);
# Get piC loci
my $pic_loci = get_pic_loci($picloci_file);

## Compare loci of piC tree to homlocs tree
print("Processing");
# Go though homlocs tree branches
foreach my $br_hom (@{$hom_branches}) {
	# Progress
	print(".");
	# Go through each species
	foreach my $sp (0..$#species) {
		# Save loci in new tree
		@{$syn_tree[$sp]{$br_hom}} = @{$hom_tree->[$sp]->{$br_hom}};
		# Get loci from homlocs tree
		foreach my $loc_hom (@{$hom_tree->[$sp]->{$br_hom}}) {
			# Get loc_hom coordinates 
			my @loc_hom = split(/-|:|\s/,$loc_hom);
			my $chr_hom = $loc_hom[0];
			my $beg_hom = $loc_hom[1];
			my $end_hom = $loc_hom[2];
			# Go though piC tree branches
			foreach my $pic (keys %{$pic_loci}) {
				# Get species
				my($spec) = ($pic =~ /^([^-_]+)/);
				unless ($spec_id{$spec} eq $sp) { next }
				# Get loc_pic coordinates
				my $loc_pic = $pic_loci->{$pic}->[0];
				my @loc_pic = split(/-|:|\s/,$loc_pic);
				my $chr_pic = $loc_pic[0];
				my $beg_pic = $loc_pic[1];
				my $end_pic = $loc_pic[2];
				# Check if loci overlap
				if ($chr_hom eq $chr_pic && $beg_pic <= $end_hom && $end_pic >= $beg_hom) {					
					# Get index of current locus
					my $n_inds = scalar(@{$hom_tree->[$sp]->{$br_hom}})-1;
					my($index) = grep {$hom_tree->[$sp]->{$br_hom}->[$_] eq $loc_hom} 0..$n_inds;
					# Add piC id to locus of homlocs tree
					if ($syn_tree[$sp]{$br_hom}[$index] =~ /\(.*\)/) {
						if ($syn_tree[$sp]{$br_hom}[$index] =~ /\(.*-p5.*\)/) {
							$syn_tree[$sp]{$br_hom}[$index] =~ s/\(.*\)/\($pic\)/;
						} else {
							if ($pic !~ /-p5/) {
								$syn_tree[$sp]{$br_hom}[$index] =~ s/\)/\;$pic\)/;
							}
						}
					} else {
						$syn_tree[$sp]{$br_hom}[$index] = $syn_tree[$sp]{$br_hom}[$index].' ('.$pic.')';
					}
				}
			}
		}
	}
}
# Process finished
print("done\n");

# Get global max and min hits
foreach my $sp (0..$#species) {
	foreach my $br (@{$hom_branches}) {
		foreach my $loc (@{$syn_tree[$sp]{$br}}) {
			# Get piC and its hits
			my($pic) = ($loc =~ /\((.*)\)/);
			my $hits = $pic_loci->{$pic}->[2] if $pic;
			# Get max and min hits
			if ($hits) {
				if ($hits > $hits_max[$sp]) { $hits_max[$sp] = $hits }
				if ($hits < $hits_min[$sp]) { $hits_min[$sp] = $hits }
			}
		}
	}
}
# Limit global max
if ($limit_max_hits) {
	foreach my $sp (0..$#species) {
		$hits_max[$sp] = 5000;
		$hits_min[$sp] = 10;
	}
}

## Output piC tree
# Sort tree branches by number of homologs
my($sorted_branches,$sorted_brinfo) = sort_tree_by_homologs_expr(\@syn_tree, $hom_branches);
# Create plain text representation of piC tree
draw_global_tree_txt(\@syn_tree, $sorted_branches, 'homlocs_tree_piCids.txt');
# Create html heatmap of piC tree
draw_global_tree_html_heatmap(\@syn_tree, $sorted_branches, 'homlocs_tree_piCids_hm.html');

draw_global_tree_expression(\@syn_tree, $sorted_branches, 'homlocs_tree_expression.txt');

# Open output file
$outfile_expr = FileIO::find_unique_filename($outfile_expr);
my $out = FileIO::open_outfile($outfile_expr);
foreach my $br (@{$sorted_brinfo}) {
	print($out "$br->[0]\t$br->[1]\t$br->[2]\n");
}

exit;

################################# subroutines #################################

sub get_tree_data {
	# Take name of pic tree file
	my($infile) = @_;
	# Get file data
	my @in_data = FileIO::get_file_data_array($infile);
	# Global tree hashes for each species
	my @tree = ();
	my $branch = 1;
	my @branches = (1);
	# Get number of species
	my $n_spec = scalar(split(/\t/,$in_data[0]));
	# Go through file data
	foreach my $line (@in_data) {
		# Increment branch number if line is empty
		if ($line =~ /^\s*$/) {
			$branch++;
			push(@branches, $branch);
		} 
		# Save clusters of each species to branch in tree hashes
		else {
			# Get clusters from branch for each species
			my @clusters = split(/\t/, $line);
			# Save cluster to branch in corresponding tree hash
			for (my $sp_i=0; $sp_i<$n_spec; $sp_i++) {
				push(@{$tree[$sp_i]{$branch}}, $clusters[$sp_i]) if $clusters[$sp_i];
			}
		}
	}
	# Fill empty branch entries
	foreach my $br (@branches) {
		for (my $sp_i=0; $sp_i<$n_spec; $sp_i++) {
			@{$tree[$sp_i]{$br}} = () unless $tree[$sp_i]{$br};
		}
	}
	
	return \(@tree, @branches);
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

sub sort_tree_by_homologs_expr {
	# Take piC tree and its branches
	my($tree, $branches) = @_;
	# Variables initialization
	my %grouped_brs = ();
	my @sorted_brs = ();
	my @sorted_brs_info = ();
	# Go through each tree branch
	foreach my $br (@{$branches}) {
		# Get number of branch homologs and expressed piCs
		my $homs = 0;
		my $expr = 0;
		my $hits = 0;
		# Go through each species
		foreach my $i (0..(scalar(@{$tree})-1)) {
			# Skip empty branches
			unless ($tree->[$i]->{$br}) { next }
			unless (scalar(@{$tree->[$i]->{$br}})>0) { next }
			# Increment number of homologs
			$homs++;
			# Check if loc is expressed piC
			if (grep {$_ =~ /\(.*\)/} @{$tree->[$i]->{$br}}) {
				# Increment number of expressed homologous piCs
				$expr++;
				# Get total length of homloc group
				my $group_len = 0;
				foreach my $homloc (@{$tree->[$i]->{$br}}) {
					my($chr_hom,$beg_hom,$end_hom) = split(/-|:|\s/,$homloc);
					$group_len += ($end_hom-$beg_hom+1);
				}
				# Go through each homloc in group
				foreach my $homloc (@{$tree->[$i]->{$br}}) {
					# Check if homloc is expressed piC
					if ($homloc =~ /\(.*\)/) {
						# Get piC id
						my($pic) = ($homloc =~ /\((.*)\)/);
						if ($pic =~ /;/) { $pic =~ s/;.*// }
						# Add to hit count
						my($chr_hom,$beg_hom,$end_hom) = split(/-|:|\s/,$homloc);
						my $homloc_len = ($end_hom-$beg_hom+1);
						$hits += $pic_loci->{$pic}->[2]*($homloc_len/$group_len);
					}
				}
			}
		}
		$hits = $hits/$expr if $expr;
		# Link branch to number of homologs and expressed piCs
		push(@{$grouped_brs{$homs}{$expr}{$hits}},$br);
	}
	# Sort branches
	foreach my $homs (sort {$b <=> $a} keys %grouped_brs) {
		foreach my $expr (sort {$b <=> $a} keys %{$grouped_brs{$homs}}) {
			foreach my $hits (sort {$b <=> $a} keys %{$grouped_brs{$homs}{$expr}}) {
				foreach my $br (@{$grouped_brs{$homs}{$expr}{$hits}}) {
					push(@sorted_brs,$br);
					push(@sorted_brs_info,[$homs,$expr,$hits]);
				}
			}
		}
	}
	return \(@sorted_brs,@sorted_brs_info);
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
	my @loc_absent  = (240, 240, 240);
	my @loc_present = (0  , 0  , 0  );
	my @expr_high   = (245, 35 , 0  );
	my @expr_low    = (20 , 51 , 219);
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
			if (scalar(@{$tree->[$i]->{$br}}) == 0) { @rgb = @loc_absent }
			if (scalar(@{$tree->[$i]->{$br}}) >= 1) { @rgb = @loc_present }
			# Add red if piC is expressed
			if (grep {$_ =~ /\(.*\)/} @{$tree->[$i]->{$br}}) {
				# Get total length of homloc group
				my $group_len = 0;
				foreach my $homloc (@{$tree->[$i]->{$br}}) {
					my($chr_hom,$beg_hom,$end_hom) = split(/-|:|\s/,$homloc);
					$group_len += ($end_hom-$beg_hom+1);
				}
				# Get hits
				my $hits = 0;
				foreach my $homloc (@{$tree->[$i]->{$br}}) {
					if ($homloc =~ /\(.*\)/) {
						my($pic) = ($homloc =~ /\((.*)\)/);
						if ($pic =~ /;/) { $pic =~ s/;.*// }
						my($chr_hom,$beg_hom,$end_hom) = split(/-|:|\s/,$homloc);
						my $homloc_len = ($end_hom-$beg_hom+1);
						$hits += $pic_loci->{$pic}->[2]*($homloc_len/$group_len);
					}
				}
				$hits = $hits/scalar(@{$tree->[$i]->{$br}});
				# Make sure that hits values are in range
				if ($hits > $hits_max[$i]) { $hits = $hits_max[$i] }
				if ($hits < $hits_min[$i]) { $hits = $hits_min[$i] }
				# Calculate hits fraction
				@rgb = (80, 150, 64);
			}
			print($out "<td style=\"background-color:rgb($rgb[0],$rgb[1],$rgb[2]);\"></td>");
		}
		print($out "</tr>");
	}
	print($out "</table>");
}

sub draw_global_tree_expression {
	# Take piC tree and its branches
	my($tree, $branches, $outfile) = @_;
	# Option for skipping single piCs
	my $skip_single_clusters = 0;
	# Open output file
	$outfile = FileIO::find_unique_filename($outfile);
	my $out = FileIO::open_outfile($outfile);
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
		# Go through each species (column)
		foreach my $i (0..(scalar(@{$tree})-1)) {
			# piC is expressed
			if (grep {$_ =~ /\(.*\)/} @{$tree->[$i]->{$br}}) {
				# Get total length of homloc group
				my $group_len = 0;
				foreach my $homloc (@{$tree->[$i]->{$br}}) {
					my($chr_hom,$beg_hom,$end_hom) = split(/-|:|\s/,$homloc);
					$group_len += ($end_hom-$beg_hom+1);
				}
				# Get hits
				my $hits = 0;
				foreach my $homloc (@{$tree->[$i]->{$br}}) {
					if ($homloc =~ /\(.*\)/) {
						my($pic) = ($homloc =~ /\((.*)\)/);
						if ($pic =~ /;/) { $pic =~ s/;.*// }
						my($chr_hom,$beg_hom,$end_hom) = split(/-|:|\s/,$homloc);
						my $homloc_len = ($end_hom-$beg_hom+1);
						$hits += $pic_loci->{$pic}->[2]*($homloc_len/$group_len);
					}
				}
				$hits = $hits/scalar(@{$tree->[$i]->{$br}});
				print($out "$hits\t");
			} else {
				# Loc present
				if (scalar(@{$tree->[$i]->{$br}}) >= 1) {
					print($out "0\t");
				} else {
					print($out "\t");
				}
			}
		}
		print($out "\n");
	}
}