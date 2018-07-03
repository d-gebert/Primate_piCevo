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
my $USAGE = "perl $0 <mapfile.map>\n";
unless ($ARGV[0]) {
	die("\nUsage: $USAGE\n");
}
# Input files
my $map_file = $ARGV[0];

# Get ping pong signatures
my($len_combos,$pp_nred_rate,$pp_read_rate,$pp_seqs) = get_pp_overlaps($map_file);

## Output results
# Get species id
my($sp_id) = ($map_file =~ /^([^_]+)/);
# Open output file
my $outfile = $map_file.'_pp_partner_lengths.txt';
my $out = FileIO::open_outfile($outfile);
# Title line
foreach my $len_x (sort {$a <=> $b} keys %{$len_combos}) {
	print($out "\t$len_x");
}
print($out "\n");
foreach my $len_x (sort {$a <=> $b} keys %{$len_combos}) {
	print($out "$len_x");
	foreach my $len_y (sort {$a <=> $b} keys %{$len_combos->{$len_x}}) {
		printf($out "\t%.1f",$len_combos->{$len_x}->{$len_y});
	}
	print($out "\n");
}

printf("PP share\nNon-redundant seqs: %.1f%%\nNormalized reads: %.1f%%\n",$pp_nred_rate,$pp_read_rate);

my $outfile2 = $map_file.'_pp_seqs.txt';
my $out2 = FileIO::open_outfile($outfile2);
foreach my $seq (sort keys %{$pp_seqs}) {
	print($out2 ">$pp_seqs->{$seq}\n$seq\n");
}

exit;

################################# subroutines #################################

sub get_pp_overlaps {
	# Take name of map file file
	my($map_file,$norm_hits,$norm_reads) = @_;
	# Options
	$norm_hits  = $norm_hits  ? 1 : 0;
	$norm_reads = $norm_reads ? 1 : 0;
	# Variables 1
	my $max_len = 0;
	my %hits_per_seq  = ();
	my %reads_per_seq = ();
	# Variables 2
	my %overlap = ();
	my @sw = ();
	my $prev_loc = '';
	# Variables 3
	my %len_combos = ();
	my %pp_seqs = ();
	
	## Step 1: Get hit and read counts
	# Open map file
	my $in = FileIO::open_infile($map_file);
	# Parse map file
	while (my $line = <$in>) {
		$line =~ s/\s+$//; #better chomp
		# Hit line
		if ($line !~ /^trans_id/ && $line !~ /^#/) {
			# Get line data
			my @d = split(/\t/,$line);
			my $loc = $d[0];
			my $rds = $d[3];
			my $seq = $d[4];
			my $str = $d[6];
			# Get hits per sequence
			if ($norm_hits) {
				$hits_per_seq{$seq} = 1;
			} else {
				$hits_per_seq{$seq}++;
			}
			# Get reads per sequence
			if ($norm_reads || $rds !~ /^\d+$/) {
				$reads_per_seq{$seq} = 1;
			} else {
				$reads_per_seq{$seq} = $rds;
			}
			# Save max sequence length (-> maximum overlap)
			unless ($max_len) { $max_len = 0 }
			if (length $seq > $max_len) {
				$max_len = length $seq;
			}
		}
	}
	close($in);
	
	## Step 2: Get 5' overlaps
	# Open map file
	$in = FileIO::open_infile($map_file);
	# Read map with sliding window
	while (my $line = <$in>) {
		$line =~ s/\s+$//; #better chomp
		# Hit line
		if ($line !~ /^trans_id/ && $line !~ /^#/) {
			# Get line data
			my @d = split(/\t/,$line);
			my $loc = $d[0];
			my $pos = $d[1];
			my $rds = $d[3];
			my $seq = $d[4];
			my $str = $d[6];
			my $alc = $d[8];
			# Empty sliding window when location changes
			if ($loc ne $prev_loc) { @sw = () }
			# Add hit to sliding window if current hit is on the minus strand
			if ($d[6] =~ /-/) {
				# Sliding window entries comprise coordinate|sequence|allocated reads (if available)
				push(@sw,[$pos,$seq,$alc]);
				# Check if hit(s) must be removed from sliding window
				while (1) {
					if (@sw > 1 && $pos > ($sw[0][0]) + $max_len) {
						shift(@sw);
					} else {
						last;
					}
				}
			}
			# Count overlaps in the sliding window if current hit is on the plus strand
			else {
				foreach my $i (0..@sw-1) {
					# Check whether hits overlap
					if ($sw[$i][0] + (length $sw[$i][1]) > $d[1] && $sw[$i][0] + (length $sw[$i][1]) < $pos + length $seq) {
						# Get overlap length
						my $overlap_len = $sw[$i][0] + (length $sw[$i][1]) - $pos;
						# Post processed map (re-al-locate)
						if ($alc) {
							$overlap{$overlap_len} += ($alc * $sw[$i][2]);
						}
						# Standard map
						else {
							$overlap{$overlap_len} += (($reads_per_seq{$seq} / $hits_per_seq{$seq}) * ($reads_per_seq{$sw[$i][1]} / $hits_per_seq{$sw[$i][1]}));
						}
						# Count length combination if 10 nt overlap
						if ($overlap_len == 10) {
							$len_combos{length($seq)}{length($sw[$i][1])} += (($reads_per_seq{$seq} / $hits_per_seq{$seq}) * ($reads_per_seq{$sw[$i][1]} / $hits_per_seq{$sw[$i][1]})/2);
							$len_combos{length($sw[$i][1])}{length($seq)} += (($reads_per_seq{$seq} / $hits_per_seq{$seq}) * ($reads_per_seq{$sw[$i][1]} / $hits_per_seq{$sw[$i][1]})/2);
							$pp_seqs{$seq}++;
							$pp_seqs{$sw[$i][1]}++;
						}
					}
				}
			}
			$prev_loc = $loc;
		}
	}
	close($in);
	# Get share of pp partner piRNAs
	my $nred_count = 0;
	my $read_count = 0;
	my $pp_nred_count = 0;
	my $pp_read_count = 0;
	foreach my $seq (sort keys %reads_per_seq) {
		if ($pp_seqs{$seq}) {
			$pp_nred_count++;
			#$pp_read_count += ($reads_per_seq{$seq}/$hits_per_seq{$seq})*($pp_seqs{$seq}/2);
			$pp_read_count += $reads_per_seq{$seq};
		}
		$nred_count++;
		$read_count += $reads_per_seq{$seq};
	}
	my $pp_nred_rate = $pp_nred_count/$nred_count*100;
	my $pp_read_rate = $pp_read_count/$read_count*100;
	# Return length combinations and read counts
	return (\%len_combos,$pp_nred_rate,$pp_read_rate,\%pp_seqs);
}