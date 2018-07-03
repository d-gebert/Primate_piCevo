#!/usr/bin/perl
use strict;
use warnings;
use lib '/Users/dgebert/Dropbox/Perlmodules/';
use FileIO;
use FastaIO;

# Constants
my $zscore_th = 2.3264; # Ping-pong signature p=0.01
my $min_rpm = 0;        # Minimum coverage as reads per million
$|=1; #Autoflush
# Variables
my $pp_genes = 0;
my $pi_genes = 0;

# Program name
print("\n--- $0 ---\n");
# Collect command line arguments
my $USAGE = "perl $0 <mapfile.map> <seqfile.fas> <refseqs.fas> [min rpkm]\n";
unless ($ARGV[0] && $ARGV[1]) {
	die("\nUsage: $USAGE\n");
}
# Input files
my $map_file = $ARGV[0];
my $seq_file = $ARGV[1];
my $ref_file = $ARGV[2];
if ($ARGV[3]) { $min_rpm = $ARGV[3] }

# Get total number of reads
my $total_reads = get_total_reads_genic($map_file);

# Get hits per sequence and gene
my($hits_per_seq,$hits_per_gene) = get_hits_per_seq_gene($map_file);

# Get reference sequences
my $reference_seqs = FastaIO::get_fasta_seqs($ref_file);

# Get reads per gene
my $reads_per_gene = get_reads_per_gene($map_file,$hits_per_seq,$hits_per_gene,$total_reads,$reference_seqs);

# Get ping pong signatures
my($overlap,$zscores) = get_pp_overlaps($map_file);

# Get ping pong status
my $pp_status = get_pp_status($reads_per_gene,$overlap,$zscores);

## Output results
# Get species id
my($sp_id) = ($map_file =~ /^([^_]+)/);
# Open output file
my $outfile = FileIO::find_unique_filename($map_file.'_pi_target_genes.txt');
my $out = FileIO::open_outfile($outfile);
# Title line
printf($out "Gene_id\trpkm\tTotal_reads\tSense\tAntisense\tS_rate\tAS_rate\tPP_signature\n");
# Go through each gene
foreach my $gid (sort {$reads_per_gene->{$b}->[0] <=> $reads_per_gene->{$a}->[0]} keys %{$reads_per_gene}) {
	# Get read counts
	my $rpm = $reads_per_gene->{$gid}->[0] ? $reads_per_gene->{$gid}->[0] : 0;
	my $rds = $reads_per_gene->{$gid}->[1] ? $reads_per_gene->{$gid}->[1] : 0;
	my $pls = $reads_per_gene->{$gid}->[2] ? $reads_per_gene->{$gid}->[2] : 0;
	my $mns = $reads_per_gene->{$gid}->[3] ? $reads_per_gene->{$gid}->[3] : 0;
	my $gsy = $reads_per_gene->{$gid}->[4] ? $reads_per_gene->{$gid}->[4] : $gid;
	my $plp = $pls/$rds*100;
	my $mip = $mns/$rds*100;
	# Skip if rpm below threshold
	if ($rpm < $min_rpm) { next }
	# Get ping pong status
	my $pps = $pp_status->{$gid};
	# Count ping pong genes and target genes
	if ($pps == 1) { $pp_genes++ }
	$pi_genes++;
	# Print to output
	printf($out "%s\t%.1f\t%.1f\t%.1f\t%.1f\t%.2f\t%.2f\t%d\t%s\n",$gid,$rpm,$rds,$pls,$mns,$plp,$mip,$pps,$gsy);
}

print("piR genes \[>${min_rpm}rpm\]:\t$pi_genes\nPing-pong genes:\t$pp_genes\n\n");

exit;

################################# subroutines #################################

sub get_total_reads {
	# Take name of map file file
	my($seq_file) = @_;
	# Storage variables
	my $total_reads = 0;
	# Open map file
	my $in = FileIO::open_infile($seq_file);
	# Go through each line
	while (my $line = <$in>) {
		$line =~ s/\s+$//; #better chomp
		# Header line
		if ($line =~ /^>/) {
			my($rds) = ($line =~ />(\d+)/);
			$total_reads += $rds;
		}
	}
	return $total_reads;
}

sub get_total_reads_genic {
	# Take name of map file file
	my($map_file) = @_;
	# Storage variables
	my $total_reads = 0;
	my %reads_per_seq = ();
	# Open map file
	my $in = FileIO::open_infile($map_file);
	# Go through each line
	while (my $line = <$in>) {
		$line =~ s/\s+$//; #better chomp
		# Hit line
		if ($line !~ /^trans_id/ && $line !~ /^\s*$/) {
			# Get line data
			my @d = split(/\t/,$line);
			# Get reads and sequence
			my $rds = $d[3];
			my $seq = $d[4];
			# Save reads per sequence
			$reads_per_seq{$seq} = $rds;
		}
	}
	# Get sum of all reads per seq
	foreach my $seq (keys %reads_per_seq) {
		$total_reads += $reads_per_seq{$seq};
	}
	return $total_reads;
}

sub get_hits_per_seq_gene {
	# Take name of map file file
	my($map_file) = @_;
	# Storage variables
	my %hits_per_seq = ();
	my %hits_per_gene = ();
	# Open map file
	my $in = FileIO::open_infile($map_file);
	# Go through each line
	while (my $line = <$in>) {
		$line =~ s/\s+$//; #better chomp
		# Hit line
		if ($line !~ /^trans_id/ && $line !~ /^\s*$/) {
			# Get line data
			my @d = split(/\t/,$line);
			# Get gene info
			my($gid) = ($d[0] =~ /gene:(\w+)/);
			# Get reads and sequence
			my $seq = $d[4];
			# Save number of hits per sequence / per gene
			$hits_per_seq{$seq}{$gid}++;
			$hits_per_gene{$gid}{$seq}++;
		}
	}
	return \(%hits_per_seq,%hits_per_gene);
}

sub get_reads_per_gene {
	# Take name of map file file
	my($map_file,$hits_per_seq,$hits_per_gene,$total_reads,$reference_seqs) = @_;
	# Storage variable
	my %reads_per_gene = ();
	# Open map file
	my $in = FileIO::open_infile($map_file);
	# Go through each line
	while (my $line = <$in>) {
		$line =~ s/\s+$//; #better chomp
		# Hit line
		if ($line !~ /^trans_id/ && $line !~ /^\s*$/) {
			# Get line data
			my @d = split(/\t/,$line);
			# Get gene info
			my $ref = $d[0];
			my($gid) = ($d[0] =~ /gene:(\w+)/);
			my($gsy) = ($d[0] =~ /gene_symbol:(\w+)/);
			# Get reads and sequence
			my $rds = $d[3];
			my $seq = $d[4];
			# Get strand
			my $str = $d[6];
			# Get ref seq length
			my $len = length($reference_seqs->{$ref});
			# Get normalized reads
			my $nor = ($rds/$hits_per_gene->{$gid}->{$seq})/(keys %{$hits_per_seq->{$seq}});
			my $rpm = $nor/($len/1000 * $total_reads/1000000);
			# Save number of hits per sequence per gene
			$reads_per_gene{$gid}[0] += $rpm;
			$reads_per_gene{$gid}[1] += $nor;
			$reads_per_gene{$gid}[2] += $nor if $str eq '+';
			$reads_per_gene{$gid}[3] += $nor if $str eq '-';
			$reads_per_gene{$gid}[4]  = $gsy;
		}
	}
	return \%reads_per_gene;
}

sub get_pp_overlaps {
	# Take name of map file file
	my($map_file,$norm_hits,$norm_reads) = @_;
	# Options
	$norm_hits  = $norm_hits  ? 1 : 0;
	$norm_reads = $norm_reads ? 1 : 0;
	# Variables 1
	my %max_len = ();
	my %hits_per_seq  = ();
	my %reads_per_seq = ();
	# Variables 2
	my %overlap = ();
	my @sw = ();
	my $prev_loc = '';
	# Variables 3
	my %zscores = ();
	
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
			my($gid) = ($d[0] =~ /gene:(\w+)/);
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
			unless ($max_len{$gid}) { $max_len{$gid} = 0 }
			if (length $seq > $max_len{$gid}) {
				$max_len{$gid} = length $seq;
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
			my($gid) = ($d[0] =~ /gene:(\w+)/);
			my $pos = $d[1];
			my $rds = $d[3];
			my $seq = $d[4];
			my $str = $d[6];
			my $alc = $d[8];
			# Empty sliding window when location changes
			if ($gid ne $prev_loc) { @sw = () }
			# Add hit to sliding window if current hit is on the minus strand
			if ($d[6] =~ /-/) {
				# Sliding window entries comprise coordinate|sequence|allocated reads (if available)
				push(@sw,[$pos,$seq,$alc]);
				# Check if hit(s) must be removed from sliding window
				while (1) {
					if (@sw > 1 && $pos > ($sw[0][0]) + $max_len{$gid}) {
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
						# Post processed map (re-al-locate)
						if ($alc) {
							$overlap{$gid}{$sw[$i][0] + (length $sw[$i][1]) - $pos} += ($alc * $sw[$i][2]);
						}
						# Standard map
						else {
							$overlap{$gid}{$sw[$i][0] + (length $sw[$i][1]) - $pos} += (($reads_per_seq{$seq} / $hits_per_seq{$seq}) * ($reads_per_seq{$sw[$i][1]} / $hits_per_seq{$sw[$i][1]}));
						}
					}
				}
			}
			$prev_loc = $gid;
		}
	}
	close($in);
	
	## Step 3: Calculate ping pong z-scores
	foreach my $gid (sort keys %overlap) {
		# Skip if no ping pong overlaps
		unless ($overlap{$gid}{10}) {
			$zscores{$gid} = 0;
			next;
		}
		# Statistic variables
		my $average_background = 0;
		my $variance = 0;
		my $stddeviation = 0;
		# Get average background (overlaps 1-9+11-20)
		foreach my $len (1..($max_len{$gid}-1)) {
			unless ($overlap{$gid}{$len}) { next }
			if ($len >= 1 && $len <= 9) {
				$average_background += $overlap{$gid}{$len};
			}
			elsif ($len >= 11 && $len <= 20) {
				$average_background += $overlap{$gid}{$len};
			}
		}
		$average_background = $average_background / 19;
		# Get background variance and standard deviation
		foreach my $len (1..($max_len{$gid}-1)) {
			unless ($overlap{$gid}{$len}) { next }
			if ($len >= 1 && $len <= 9) {
				$variance += ($overlap{$gid}{$len} - $average_background)**2;
			}
			elsif ( $len >= 11 && $len <= 20) {
				$variance += ($overlap{$gid}{$len} - $average_background)**2;
			}
		}
		$variance = $variance / 19;
		$stddeviation = sqrt($variance);
		# Get z-score
		if ($stddeviation > 0) {
			$zscores{$gid} = ($overlap{$gid}{10} - $average_background) / $stddeviation;
		} else {
			if ($overlap{$gid}{10} > 0) {
				$zscores{$gid} = 999;
			}
		}
	}
	
	# Return overlaps and z-score
	return \(%overlap,%zscores);
}

sub get_pp_status {
	# Take name of map file file
	my($reads_per_gene,$overlap,$zscores) = @_;
	# Storage variable
	my %pp_status = ();
	# Go through each gene
	foreach my $gid (sort keys %{$reads_per_gene}) {
		# Initialize ping pong status
		$pp_status{$gid} = 0;
		# Check if z-score above threshold
		unless ($zscores->{$gid} && $zscores->{$gid} > $zscore_th) { next }
		# Check if overlaps exist
		unless ($overlap->{$gid}) { next }
		# Check if 10 nt has most overlaps
		my @lens_ranked = ();
		my @lens_values = ();
		foreach my $len (sort {$overlap->{$gid}->{$b} <=> $overlap->{$gid}->{$a}} keys %{$overlap->{$gid}}) {
			push(@lens_ranked,$len);
			push(@lens_values,$overlap->{$gid}->{$len});
		}
		unless ($lens_values[1]) { $lens_values[1] = 0 }
		# Grant ping pong status if 10 nt has unambiguously most overlaps
		if ($lens_ranked[0] == 10 && $lens_values[0] > $lens_values[1]) { 
			$pp_status{$gid} = 1;
		}
	}
	return \%pp_status;
}