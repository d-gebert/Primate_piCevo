#!/usr/bin/perl
use strict;
use warnings;
use lib '/Users/dgebert/Dropbox/Perlmodules/';
use FileIO;

# Constants
$|=1; #Autoflush
# Variables

# Species identifiers
my @species = ('Hsap','Mmul','Mfas','Cjac','Mmur','Ltar');
my %spec_id = ();
foreach my $i (0..$#species) {
	$spec_id{$species[$i]} = $i;
}

# Program name
print("\n--- $0 ---\n");
# Collect command line arguments
my $USAGE = "perl $0 <pic_locs.txt> <repeatmasker_list.txt> <genomes_list.txt>\n";
unless ($ARGV[0] && $ARGV[1]) {
	die("\nUsage: $USAGE\n");
}
# Genome locs
my $pic_locs_file = $ARGV[0];
my @pic_locs_data = FileIO::get_file_data_array($pic_locs_file);
# Repeatmasker files
my $rep_file_list = $ARGV[1];
my @repmask_files = FileIO::get_file_data_array($rep_file_list);
# Genome files
my $gnm_file_list = $ARGV[1];
my @genome_files = FileIO::get_file_data_array($gnm_file_list);

# Get locations
my @pic_locs = ();
foreach my $line (@pic_locs_data) {
	# Skip blank lines
	if ($line =~ /^$/) { next }
	# Get loc for each species
	my @homloc = split(/\t/,$line);
	# Add loc to corresponding species if it exists
	foreach my $i (0..$#species) {
		# Homloc exists for spec i
		if ($homloc[$i]) {
			# Skip if pic not expressed
			unless ($homloc[$i] =~ /\(.*\)/) { next }
			# Get pic id
			my($pic_id) = ($homloc[$i] =~ /\((.*)\)/);
			$homloc[$i] =~ s/\(.*//;
			# Get loc coordinates
			my @loc = split(/:|-|\s/,$homloc[$i]);
			my $id = (keys %{$pic_locs[$i]})+1;
			unshift(@loc,$id);
			push(@loc,$pic_id);
			$pic_locs[$i]{$id} = [@loc];
		}
	}
}

# Go through each species
foreach my $i (0..$#species) {
	# Get repeat loci
	my $rep_data = get_repeatmask_data($repmask_files[$i],$pic_locs[$i]);
	# Get repeat divergence
	my $mean_div = get_mean_repeat_divergence($rep_data);
	# Get genome slices size
	my $pic_size = get_slices_size($pic_locs[$i]);
	# Get repeat share
	my $rep_share = get_total_repeat_share($rep_data,$pic_size);
	# Output
	printf("%s\t%.4f\t%.4f\n",$species[$i],$mean_div,$rep_share);
	# Get repeat class divergence
	my($rep_class_div,$rep_classes) = get_repeat_class_divergence($rep_data,$pic_size);
	# Output results: TE landscape
	my $out1 = FileIO::open_outfile('TE_landscape.'.$species[$i].'_homlocs.txt');
	foreach my $class (sort {$a cmp $b} keys %{$rep_classes}) {
		print($out1 "\t$class");
	}
	print($out1 "\n");
	foreach my $div (sort {$a <=> $b} keys %{$rep_class_div}) {
		print($out1 "$div\t");
		foreach my $class (sort {$a cmp $b} keys %{$rep_classes}) {
			$rep_class_div->{$div}->{$class} = 0 unless $rep_class_div->{$div}->{$class};
			printf($out1 "%.6f\t",$rep_class_div->{$div}->{$class});
		}
		print($out1 "\n");
	}
	# Get repeat loci
	my $rep_data_gnm = get_repeatmask_data($repmask_files[$i]);
	# Get genome size
	my $gnm_size = get_genome_size($genome_files[$i]);
	# Get repeat class divergence
	my($rep_class_div_gnm,$rep_classes_gnm) = get_repeat_class_divergence($rep_data_gnm,$gnm_size);
	# Output results: TE landscape
	my $out2 = FileIO::open_outfile('TE_landscape.'.$species[$i].'_genome.txt');
	foreach my $class (sort {$a cmp $b} keys %{$rep_classes_gnm}) {
		print($out2 "\t$class");
	}
	print($out2 "\n");
	foreach my $div (sort {$a <=> $b} keys %{$rep_class_div_gnm}) {
		print($out2 "$div\t");
		foreach my $class (sort {$a cmp $b} keys %{$rep_classes_gnm}) {
			$rep_class_div_gnm->{$div}->{$class} = 0 unless $rep_class_div_gnm->{$div}->{$class};
			printf($out2 "%.6f\t",$rep_class_div_gnm->{$div}->{$class});
		}
		print($out2 "\n");
	}
}

exit;

################################# subroutines #################################

sub get_repeatmask_data {
	# Take repeatmasker file name
	my($repeatmask_file,$locs) = @_;
	# Storage variable
	my %rep_data = ();
	# If slice locs are given, get positions
	my %chr_poss = ();
	if ($locs) {
		foreach my $loc (sort keys %{$locs}) {
			# Get loc coordinates
			my $chr = $locs->{$loc}->[1];
			my $beg = $locs->{$loc}->[2];
			my $end = $locs->{$loc}->[3];
			# Save positions
			foreach my $pos ($beg..$end) {
				$chr_poss{$chr}{$pos} = 1;
			}
		}
	}
	# Open input file
	my $in = FileIO::open_infile($repeatmask_file);
	# Parse repeatmasker file
	while (my $line = <$in>) {
		$line =~ s/\s+$//; #better chomp
		# Line starts with sw score
		if ($line =~ /^\s*\d+/) {
			# Remove leading whitespace
			$line =~ s/^\s*//;
			# Get line data
			my @d = split(/\s+/,$line);
			my $chr = $d[4];
			my $beg = $d[5];
			my $end = $d[6];
			my $cla = $d[10];
			if ($cla =~ /Simple_repeat/ || $cla =~ /Low_complexity/ || $cla =~ /Satellite/) { next }
			# If slice locs are given, skip repeats outside locs
			my $rep_in = 0;
			if ($chr_poss{$chr}{$beg}) { $rep_in++ }
			if ($chr_poss{$chr}{$end}) { $rep_in++ }
			if ($locs && not $rep_in) { next }
			# If slice locs are given, update repeat coordinates
			if ($locs) {
				# If repeat not completely in pic, update
				if ($rep_in < 2) {
					foreach my $loc (sort keys %{$locs}) {
						# Get loc coordinates
						my $loc_chr = $locs->{$loc}->[1];
						my $loc_beg = $locs->{$loc}->[2];
						my $loc_end = $locs->{$loc}->[3];
						# Same chr/scaff
						if ($loc_chr eq $chr) {
							# Same location
							if ($loc_beg <= $end && $loc_end >= $beg) {
								# Update repeat coordinates if needed
								if ($loc_beg > $beg) { $beg = $loc_beg }
								if ($loc_end < $end) { $end = $loc_end }
							}
						}
					}
				}
				$d[5] = $beg;
				$d[6] = $end;
			}
			# Save repeat data
			push(@{$rep_data{$chr}},\@d);
		}
	}
	return \%rep_data;
}

sub get_mean_repeat_divergence {
	# Take pic and repeat data
	my($rep_data) = @_;
	# Storage variable
	my $div_sum = 0;
	my $rep_count = 0;
	# Go through each chromosome
	foreach my $chr (sort keys %{$rep_data}) {
		# Go through each repeat
		foreach my $rep (@{$rep_data->{$chr}}) {
			# Get repeat coordinates
			my $rep_div = $rep->[1];
			my $rep_beg = $rep->[5];
			my $rep_end = $rep->[6];
			my $rep_name = $rep->[9];
			my $rep_class = $rep->[10];
			# Save divergence
			$div_sum += $rep_div;
			$rep_count++;
		}
	}
	my $mean_div = $div_sum/$rep_count;
	return $mean_div;
}

sub get_slices_size {
	# Take fasta file name
	my($loc_data) = @_;
	# Var init
	my $total_size = 0;
	# Parse fasta file
	foreach my $loc (sort keys %{$loc_data}) {
		# Get loc coordinates
		my $chr = $loc_data->{$loc}->[1];
		my $beg = $loc_data->{$loc}->[2];
		my $end = $loc_data->{$loc}->[3];
		# Increment total size
		$total_size += $end-$beg+1;
	}
	return $total_size;
}

sub get_genome_size {
	# Take fasta file name
	my($fasta_file) = @_;
	# Var init
	my $total_size = 0;
	# Open input file
	my $in = FileIO::open_infile($fasta_file);
	# Parse fasta file
	while (my $line = <$in>) {
		$line =~ s/\s+$//; #better chomp
		# Sequence
		if ($line !~ /^>/) {
			$total_size += length($line);
		}
	}
	return $total_size;
}

sub get_total_repeat_share {
	# Take pic and repeat data
	my($rep_data,$genome_size) = @_;
	# Storage variable
	my $rep_bps = 0;
	# Go through each chromosome
	foreach my $chr (sort keys %{$rep_data}) {
		# Storage variable
		my %rep_pos = ();
		# Go through each repeat
		foreach my $rep (@{$rep_data->{$chr}}) {
			# Get repeat coordinates
			my $rep_div = $rep->[1];
			my $rep_beg = $rep->[5];
			my $rep_end = $rep->[6];
			my $rep_name = $rep->[9];
			my $rep_class = $rep->[10];
			# Count repeat bps
			foreach my $pos ($rep_beg..$rep_end) {
				$rep_pos{$pos} = 1;
			}
		}
		foreach my $pos (keys %rep_pos) {
			$rep_bps += 1;
		}
	}
	my $rep_share = $rep_bps/$genome_size*100;
	return $rep_share;
}

sub get_repeat_class_divergence {
	# Take pic and repeat data
	my($rep_data,$genome_size) = @_;
	# Storage variable
	my %class_div = ();
	my %classes = ();
	# Go through each chromosome
	foreach my $chr (sort keys %{$rep_data}) {
		# Go through each repeat
		foreach my $rep (@{$rep_data->{$chr}}) {
			# Get repeat coordinates
			my $rep_div = $rep->[1];
			my $rep_beg = $rep->[5];
			my $rep_end = $rep->[6];
			my $rep_name = $rep->[9];
			my $rep_class = $rep->[10];
			my $rep_len = $rep_end-$rep_beg+1;
			# Save divergence
			my $div = int($rep_div+0.5);
			$class_div{$div}{$rep_class} += $rep_len;
			$classes{$rep_class} += $rep_len;
		}
	}
	# Calculate relative repeats shares
	foreach my $div (keys %class_div) {
		foreach my $class (keys %{$class_div{$div}}) {
			$class_div{$div}{$class} = $class_div{$div}{$class}/$genome_size*100;
		}
	}
	foreach my $class (keys %classes) {
		$classes{$class} = $classes{$class}/$genome_size*100;
	}
	return \(%class_div,%classes);
}