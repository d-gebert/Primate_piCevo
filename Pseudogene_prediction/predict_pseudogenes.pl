#!/usr/bin/perl

use strict;
use warnings;
use lib '/Users/dgebert/Dropbox/Perlmodules/';
use FileIO;
use FastaIO;
use BioStat;

# Take time t0
my $t0 = time;

## Constants
# Blast options
my $blast_opts = "-outfmt '7 qseqid sseqid pident evalue qcovs qlen qstart qend slen sstart send sstrand length bitscore' -task dc-megablast";
# Max number of threads
my $n_childs_max = 1;
$|=1; #Autoflush
## Variables

# Collect command line arguments
my $USAGE = "perl $0 <genome.fas.masked> <cds.fas> <genome.gff> [max_cores]\n";
unless ($ARGV[0] && $ARGV[1] && $ARGV[2]) {
	die("\nUsage: $USAGE\n");
}
# Option: Max number of cores
if ($ARGV[3]) { $n_childs_max = $ARGV[3] }

# Input files
my $genome_file = $ARGV[0];
my $cdsseq_file = $ARGV[1];
my $gengff_file = $ARGV[2];

print("Blast on chromosomes...");
# Blast cds on chromosomes and get output file names list
my @blast_files = blast_on_chromosomes($genome_file,$cdsseq_file,$blast_opts,$n_childs_max);
print("....ok\n");

print("Extract cds locations...");
# Get cds locations and positions
my $cds_locs = get_gff_gene_locs($gengff_file);
my $cds_poss = get_cds_positions($cds_locs);
my $cd_genes = get_coding_genes($cds_locs);
print("...ok\n");

print("Extract blast hits...");
# Get blast hits separated by chromosome and strand
my $gene_hits = get_chr_blast_hits(\@blast_files,$cds_poss);
print("......ok\n");

print("Merge overlapping hits...");
# Merge overlapping hits to super-hits
my $super_hits = merge_to_superhits($gene_hits);
print("..ok\n");

print("Combine neighboring hits...");
# Combine neighboring superhits to pseudogene structures
my $unitd_hits = combine_to_pseudo_units($super_hits);
print("ok\n");

print("Choose best parent genes...");
# Get unit properties
my $unit_props = get_unit_properties($unitd_hits);
# Choose best pseudogene parent
remove_redundant_units($unitd_hits,$unit_props);
print("ok\n");

print("Remove short fragments...");
# Eliminate short isolated fragments
remove_short_fragments($unitd_hits,$unit_props);
print("..ok\n");

print("Classify pseudogenes...");
# Add pseudogene type to unit properties
$unit_props = categorize_pseudogenes($unitd_hits,$unit_props,$cd_genes);
print("....ok\n");

print("Prepare output data...");
# Prepare data for output
my($pseudogenes,$pseudoexons) = prepare_data($unitd_hits,$unit_props);
print(".....ok\n");

print("Creating output...");
# Name output file
my($genome_name) = ($genome_file =~ /(\w+)\./);
my $outfile = $genome_name.'.pseudogenes.gff.txt';
# Output pseudogene data
print_gff_output($pseudogenes,$pseudoexons,$outfile);
print(".........ok\n");

# Print elapsed time
print("Elapsed time: ");
print_elapsed_time($t0);

################################# subroutines #################################

sub blast_on_chromosomes {
	# Take genome and cds file names and max thread number
	my($genome_file,$cdsseq_file,$blast_opts,$n_childs_max) = @_;
	# Options
	my $ignore_unreliable_chrs = 1;
	my $min_chr_length = 5_000;
	# Multithreading
	my @childs = ();
	my $n_childs = 0;
	# Get genome sequence and chromosome names
	my $genome = FastaIO::get_fasta_seqs($genome_file);
	my $chroms = FastaIO::get_fasta_headers($genome_file);
	# Get unzipped files/filenames
	if ($cdsseq_file =~ /\.gz$/) { system("gunzip -k -f $cdsseq_file") }
	$genome_file =~ s/\.gz$//;
	$cdsseq_file =~ s/\.gz$//;
	# Get blast file names
	my @blast_files = ();
	foreach my $chr (@{$chroms}) {
		# Ignore unreliable/short chromosomes/scaffolds
		if ($ignore_unreliable_chrs) {
			if ($chr =~ /chrUn_/) { next }
			if ($chr =~ /random/) { next }
			if ($chr =~ /_alt/) { next }
			if (length($genome->{$chr}) < $min_chr_length) { next }
		}
		my $chr_file = $genome_file.'.'.$chr.'.fas';
		my $blast_file_gz = $cdsseq_file.'.'.$chr_file.'.bln.out.gz';
		push(@blast_files,$blast_file_gz);
	}
	# Go through each chromosome
	foreach my $chr (@{$chroms}) {
		# Ignore unreliable/short chromosomes/scaffolds
		if ($ignore_unreliable_chrs) {
			if ($chr =~ /chrUn_/) { next }
			if ($chr =~ /random/) { next }
			if ($chr =~ /_alt/) { next }
			if (length($genome->{$chr}) < $min_chr_length) { next }
		}
		# Fork process
		my $pid;
		if ($n_childs == $n_childs_max) { 
			$pid = wait();
			$n_childs--;
		}
		if (defined($pid = fork())) {
			# Collect child PIDs in parent process
			if ($pid) {
				$n_childs++; 
				push(@childs,$pid);
			}
			# Child process
			else {
				# Open output file for chromosome sequence
				my $chr_file = $genome_file.'.'.$chr.'.fas';
				my $out = FileIO::open_outfile($chr_file);
				# Print (formatted) chromosome sequence to file
				my $fmt_chr = FastaIO::format_sequence($genome->{$chr},80);
				print($out ">$chr\n${$fmt_chr}\n");
				# Name blast output file
				my $blast_file = $cdsseq_file.'.'.$chr_file.'.bln.out';
				my $blast_file_gz = $cdsseq_file.'.'.$chr_file.'.bln.out.gz';
				# Blast and pack if file does not exist
				unless (-e $blast_file || -e $blast_file_gz) {
					system("blastn -query $cdsseq_file -subject $chr_file -out $blast_file $blast_opts");
					system("gzip $blast_file");
				}
				# Delete chromosome sequence file
				unlink($chr_file);
				exit;
			}
		}
		# Unsuccessful fork
		else { warn 'Could not fork' and next }
	}
	# Wait for child processes
	foreach my $pid (@childs) { waitpid($pid,0) }
	# Remove unzipped cds file
	if (-e $cdsseq_file.'.gz') { unlink($cdsseq_file) }
	# Return blast file names array
	return @blast_files;
}

sub get_chr_blast_hits {
	# Take blast output files list
	my($blast_files,$cds_poss) = @_;
	# Storage variable
	my %gene_hits = ();
	# Go through each blast output file
	foreach my $blast_file (@{$blast_files}) {
		# Get chromosome name
		my($chr) = ($blast_file =~ /.([^.]+).fas.bln.out/);
		# Get blast hits for chromosome
		$gene_hits{$chr} = get_blast_hits($blast_file,$cds_poss);
		#last;
	}
	# Return blast hits
	return \%gene_hits;
}

sub get_blast_hits {
	# Take infile name
	my($file,$cds_poss) = @_;
	# E-value threshold
	my $eval_th = 0.0001;
	# Get infile data as array
	my $in = FileIO::open_infile($file);
	# Initialize variables
	my $gene_id = '';
	my $gsymbol = '';
	my %gene_hits = ();
	# Parse blast output
	while (my $line = <$in>) {
		$line =~ s/\s+$//; #better chomp
		# Query name line
		if ($line =~ /^# Query:/) {
			# Get gene id and symbol
			($gene_id) = ($line =~ /GeneID:([^;,\]]+)/);
			($gsymbol) = ($line =~ /\[gene=(\w+)\]/);
			unless ($gsymbol) { ($gsymbol) = ($line =~ /\[protein_id=(\w+)\.*/) }
			unless ($gsymbol) { $gsymbol = 'Unknown' }
		}
		# Results line
		elsif ($line =~ /^[^#]\w+/) {
			# Get info
			my @d = split(/\t/, $line);
			my $chr_id = $d[1];
			my $strand = $d[11];
			my $evalue = $d[3];
			# Skip low e-values
			if ($evalue > $eval_th) { next }
			# Get hit coordinates
			my($beg,$end) = sort {$a <=> $b} ($d[9],$d[10]);
			# Check if hit on known coding genes
			my $region_coding = 0;
			foreach my $pos ($beg..$end) {
				if ($cds_poss->{$chr_id}->{$strand}->{$pos}) {
					$region_coding++;
					last;
				}
			}
			# Save hit if not on coding gene cds
			unless ($region_coding) {
				# Save hits per gene, separated by chromosome and strand
				push(@{$gene_hits{$strand}{$gene_id}},[@d,$gsymbol]);
			}
		}
	}
	return \%gene_hits;
}

sub get_gff_gene_locs {
	# Take infile name and keyword
	my($gff_file) = @_;
	# Get file data
	my $gff = FileIO::open_infile($gff_file);
	# Initialize variables
	my %cds_locs = ();
	my $gene = '';
	my $non_coding = 0;
	# Parse gff file
	while (my $line = <$gff>) {
		# Results line
		if ($line =~ /^[^#]\w+/) {
			# Get info
			my @d = split(/\t/, $line);
			my $chr  = $d[0];
			my $type = $d[2];
			my $str  = $d[6] eq '+' ? 'plus' : 'minus';
			my $info = $d[8];
			# Gene line
			if ($type =~ /gene/i) {
				# Skip non protein coding genes
				if ($info !~ /protein_coding/i) {
					$non_coding = 1;
					next;
				} else {
					$non_coding = 0;
				}
				# Get gene id
				($gene) = ($info =~ /GeneID:([^;,]+)/);
			}
			# CDS line
			elsif ($type =~ /cds/i) {
				# Skip non protein coding genes
				if ($non_coding) { next }
				# Save entries
				push(@{$cds_locs{$chr}{$str}{$gene}},\@d);
			}
		}
	}
	return \%cds_locs;
}

sub get_cds_positions {
	# Take cds locs
	my($cds_locs) = @_;
	# Storage variable
	my %cds_poss = ();
	# Go through each chromosome
	foreach my $chr (sort keys %{$cds_locs}) {
		# Go through each strand of chromosome
		foreach my $str (sort keys %{$cds_locs->{$chr}}) {
			# Go through each query gene
			foreach my $gen (sort keys %{$cds_locs->{$chr}->{$str}}) {
				# Go through each exon
				foreach my $cds (@{$cds_locs->{$chr}->{$str}->{$gen}}) {
					# Get cds coordinates
					my($cds_beg,$cds_end) = ($cds->[3],$cds->[4]);
					# Save each position
					foreach my $pos ($cds_beg..$cds_end) {
						$cds_poss{$chr}{$str}{$pos} = 1;
					}
				}
			}
		}
	}
	return \%cds_poss;
}

sub get_coding_genes {
	# Take cds locs
	my($cds_locs) = @_;
	# Storage variable
	my %cd_genes = ();
	# Go through each chromosome
	foreach my $chr (sort keys %{$cds_locs}) {
		# Go through each strand of chromosome
		foreach my $str (sort keys %{$cds_locs->{$chr}}) {
			# Go through each query gene
			foreach my $gen (sort keys %{$cds_locs->{$chr}->{$str}}) {
				# Save exons
				@{$cd_genes{$gen}} = @{$cds_locs->{$chr}->{$str}->{$gen}};
			}
		}
	}
	return \%cd_genes;
}

sub merge_to_superhits {
	# Take (non-coding) gene hits
	my($hits) = @_;
	# Storage variable
	my %super_hits = ();
	# Go through each chromosome
	foreach my $chr (sort keys %{$hits}) {
		# Go through each strand of chromosome
		foreach my $str (sort keys %{$hits->{$chr}}) {
			# Index of hit start coordinates
			my $beg_i = $str eq 'plus' ? 9 : 10;
			# Go through each query gene
			foreach my $gen (sort keys %{$hits->{$chr}->{$str}}) {
				# Storage for super hit groups
				my %super_group = ();
				# Super hit number
				my $super_n = 0;
				# Previous hit coordinates
				my $prev_beg = 0;
				my $prev_end = 0;
				# Go though each hit
				foreach my $hit (sort {$a->[$beg_i] <=> $b->[$beg_i]} @{$hits->{$chr}->{$str}->{$gen}}) {
					# Get hit coordinates
					my($beg,$end) = sort {$a <=> $b} ($hit->[9],$hit->[10]);
					# Skip redundant hits
					if ($prev_end && $beg <= $prev_end && $end <= $prev_end) { next }
					if ($prev_end && $beg == $prev_beg && $end == $prev_end) { next }
					# Increment super hit number if current and previous hits do not overlap 
					unless ($prev_end && $beg <= $prev_end) { $super_n++ }
					# Link hit to current super hit group
					push(@{$super_group{$super_n}},$hit);
					# Save hit end coordinate for next iteration
					$prev_beg = $beg;
					$prev_end = $end;
				}
				# Merge overlapping hits to super hits
				foreach my $group_n (sort {$a <=> $b} keys %super_group) {
					# Initialize super hit coordinates
					my $beg_s = 1_000_000_000;
					my $end_s = 0;
					# Initialize super hit parameters
					my $q_id_s = $gen;
					my $s_id_s = $chr;
					my $iden_s = 0;
					my $eval_s = 0;
					my $qcov_s = 0;
					my $qlen_s = 0;
					my $qbeg_s = 0;
					my $qend_s = 0;
					my $slen_s = 0;
					my $stra_s = $str;
					my $alen_s = 0;
					my $bits_s = 0;
					my $gsym_s = '';
					my $alen_sum = 0;
					# Go through each overlapping hit
					foreach my $hit (sort {$a->[$beg_i] <=> $b->[$beg_i]} @{$super_group{$group_n}}) {
						# Get hit coordinates
						my($beg,$end) = sort {$a <=> $b} ($hit->[9],$hit->[10]);
						# Get hit parameters
						my $iden = $hit->[2];
						my $eval = $hit->[3];
						my $qcov = $hit->[4];
						my $qlen = $hit->[5];
						my $qbeg = $hit->[6];
						my $qend = $hit->[7];
						my $slen = $hit->[8];
						my $alen = $hit->[12];
						my $bits = $hit->[13];
						my $gsym = $hit->[14];
						# Get super hit coordinates
						if ($beg < $beg_s) {
							$beg_s  = $beg;
							$qbeg_s = $qbeg;
						}
						if ($end > $end_s) {
							$end_s  = $end;
							$qend_s = $qend;
						}
						# Get super hit combined parameters
						$iden_s += $iden*$alen;
						$eval_s += $eval/@{$super_group{$group_n}};
						$qcov_s += $qcov/@{$super_group{$group_n}};
						$qlen_s += $qlen/@{$super_group{$group_n}};
						$slen_s  = $slen;
						$alen_s  = $end_s-$beg_s+1;
						$bits_s += $bits/@{$super_group{$group_n}};
						$gsym_s  = $gsym;
						$alen_sum += $alen;
					}
					# Get finalized unit identity
					$iden_s = $iden_s/$alen_sum;
					# Save super hit
					my $super_hit = [
						$q_id_s,$s_id_s,$iden_s,$eval_s,$qcov_s,$qlen_s,$qbeg_s,$qend_s,$slen_s,$beg_s,$end_s,$stra_s,$alen_s,$bits_s,$gsym_s
					];
					push(@{$super_hits{$chr}{$str}{$gen}},$super_hit);
				}
			}
		}
	}
	return \%super_hits;
}

sub combine_to_pseudo_units {
	# Take super hits
	my($hits) = @_;
	# Storage variable
	my %unitd_hits = ();
	# Maximum intron size, i.e. distance between hits of same unit
	my $max_dist_min = 30_000;
	# Unit number
	my $unit_n = 0;
	# Go through each chromosome
	foreach my $chr (sort keys %{$hits}) {
		# Go through each strand of chromosome
		foreach my $str (sort keys %{$hits->{$chr}}) {
			# Index of hit start coordinates
			my $beg_i = $str eq 'plus' ? 9 : 10;
			# Go through each query gene
			foreach my $gen (sort keys %{$hits->{$chr}->{$str}}) {
				# Get intron sizes
				my @intron_lens = ();
				my $cds_end_pre = 0;
				foreach my $cds (@{$cd_genes->{$gen}}) {
					my $cds_beg = $cds->[3];
					my $cds_end = $cds->[4];
					my $intron_len = $cds_beg-$cds_end_pre+1 if $cds_end_pre;
					push(@intron_lens,$intron_len) if $intron_len;
					$cds_end_pre = $cds_end;
				}
				# Get maximum intron size/distance between (super) hits
				my $max_intron = BioStat::get_maximum(\@intron_lens);
				unless ($max_intron) { $max_intron = 0 }
				my $max_distance = $max_intron*1.5 > $max_dist_min ? $max_intron*1.5 : $max_dist_min;
				# Previous hit end coordinate
				my $prev_end = 0;
				# Go though each hit
				foreach my $hit (sort {$a->[$beg_i] <=> $b->[$beg_i]} @{$hits->{$chr}->{$str}->{$gen}}) {
					# Get hit coordinates
					my($beg,$end) = sort {$a <=> $b} ($hit->[9],$hit->[10]);
					# Increment unit number if distance between hits above threshold
					if ($prev_end && $beg-$prev_end > $max_distance) { $unit_n++ }
					# Link hit to current unit
					push(@{$unitd_hits{$chr}{$str}{$unit_n}},$hit);
					# Save hit end coordinate for next iteration
					$prev_end = $end;
				}
				$unit_n++;
			}
			$unit_n++;
		}
		$unit_n++;
	}
	return \%unitd_hits;
}

sub get_unit_properties {
	# Take hit units
	my($unitd_hits) = @_;
	# Storage variable
	my %unit_properties = ();
	# Go through each chromosome
	foreach my $chr (sort keys %{$unitd_hits}) {
		# Go through each strand of chromosome
		foreach my $str (sort keys %{$unitd_hits->{$chr}}) {
			# Go through each hit unit
			foreach my $unit (sort {$a <=> $b} keys %{$unitd_hits->{$chr}->{$str}}) {
				# Initialize unit property variables
				my $unit_gene = '';
				my $unit_begp = 1_000_000_000;
				my $unit_endp = 0;
				my $unit_iden = 0;
				my $unit_eval = 100;
				my $unit_qcov = 0;
				my $unit_alen = 0;
				my $unit_gsym = '';
				# Go through each hit
				foreach my $hit (@{$unitd_hits->{$chr}->{$str}->{$unit}}) {
					# Get hit properties
					my $iden = $hit->[2];
					my $eval = $hit->[3];
					my $qlen = $hit->[5];
					my $alen = $hit->[12];
					my($begp,$endp) = sort {$a <=> $b} ($hit->[9],$hit->[10]);
					# Get unit properties
					$unit_gene 	= $hit->[0];
					$unit_gsym  = $hit->[14];
					$unit_iden += $iden*$alen;
					$unit_alen += $alen;
					$unit_eval  = $eval < $unit_eval ? $eval : $unit_eval;
					$unit_qcov += $alen/$qlen*100;
					$unit_begp  = $begp < $unit_begp ? $begp : $unit_begp;
					$unit_endp  = $endp > $unit_endp ? $endp : $unit_endp;
				}
				unless ($unit_alen) { next }
				# Get finalized unit identity
				$unit_iden = $unit_iden/$unit_alen;
				$unit_qcov = $unit_qcov <= 100 ? int(($unit_qcov)+0.5) : 100;
				# Save unit properties
				$unit_properties{$unit} = [$unit_gene,$chr,$unit_begp,$unit_endp,$str,$unit_iden,$unit_eval,$unit_qcov,$unit_alen,$unit_gsym];
			}
		}
	}
	return \%unit_properties;
}

sub remove_redundant_units {
	# Take hit units
	my($unitd_hits,$unit_props) = @_;
	# Storage for deleted units
	my %deleted_units = ();
	# Go through each chromosome
	foreach my $chr (sort keys %{$unitd_hits}) {
		# Go through each strand of chromosome
		foreach my $str (sort keys %{$unitd_hits->{$chr}}) {
			# Go through each hit unit
			foreach my $unit_x (sort {$a <=> $b} keys %{$unitd_hits->{$chr}->{$str}}) {
				# Skip deleted units
				if ($deleted_units{$unit_x}) { next }
				# Go through each hit unit
				foreach my $unit_y (sort {$a <=> $b} keys %{$unitd_hits->{$chr}->{$str}}) {
					# Skip if same unit
					if ($unit_x eq $unit_y) { next }
					# Skip deleted units
					if ($deleted_units{$unit_y}) { next }
					# Check if units overlap
					my $units_overlap = 0;
					my $gsy_x = '';
					my $gsy_y = '';
					# Go through each hit
					foreach my $hit_x (@{$unitd_hits->{$chr}->{$str}->{$unit_x}}) {
						# Get hit coordinates
						my($beg_x,$end_x) = sort {$a <=> $b} ($hit_x->[9],$hit_x->[10]);
						$gsy_x = $hit_x->[14];
						# Go through each hit
						foreach my $hit_y (@{$unitd_hits->{$chr}->{$str}->{$unit_y}}) {
							# Get hit coordinates
							my($beg_y,$end_y) = sort {$a <=> $b} ($hit_y->[9],$hit_y->[10]);
							$gsy_y = $hit_y->[14];
							# Check if hits overlap
							if ($end_y >= $beg_x && $beg_y <= $end_x) {
								$units_overlap = 1;
							}
						}
					}
					# Units overlap: 
					if ($units_overlap) {
						# Assign score to compare group
						my $y_score = 0; # Score for identity [5], e-value [6], coverage [7]/length [8]
						$y_score += $unit_props->{$unit_y}->[5] > $unit_props->{$unit_x}->[5] ? 1 : 0;
						$y_score += $unit_props->{$unit_y}->[6] < $unit_props->{$unit_x}->[6] ? 1 : 0;
						$y_score += $unit_props->{$unit_y}->[8] > $unit_props->{$unit_x}->[8] ? 1 : 0;
						# Select best unit
						my $del_unit = $y_score < 2 ? $unit_y : $unit_x;
						# Prefer characterized over uncharacterized
						if ($gsy_x =~ /LOC/ && $gsy_y !~ /LOC/) {
							$del_unit = $unit_x;
						} elsif ($gsy_x !~ /LOC/ && $gsy_y =~ /LOC/) {
							$del_unit = $unit_y;
						}
						if ($gsy_x =~ /Unknown/ && $gsy_y !~ /Unknown/) {
							$del_unit = $unit_x;
						} elsif ($gsy_x !~ /Unknown/ && $gsy_y =~ /Unknown/) {
							$del_unit = $unit_y;
						}
						# Delete worse unit
						delete $unitd_hits->{$chr}->{$str}->{$del_unit};
						delete $unit_props->{$del_unit};
						$deleted_units{$del_unit} = 1;
						# If deleted unit is unit x, break out of loop
						if ($del_unit eq $unit_x) { last }
					}
				}
			}
		}
	}
}

sub remove_short_fragments {
	# Take hit units and unit properties
	my($unitd_hits,$unit_props) = @_;
	# Go through each chromosome
	foreach my $chr (sort keys %{$unitd_hits}) {
		# Go through each strand of chromosome
		foreach my $str (sort keys %{$unitd_hits->{$chr}}) {
			# Go through each hit unit
			foreach my $unit (sort {$a <=> $b} keys %{$unitd_hits->{$chr}->{$str}}) {
				# Get unit properties
				my $unit_eval = $unit_props->{$unit}->[6];
				my $unit_qcov = $unit_props->{$unit}->[7];
				my $unit_alen = $unit_props->{$unit}->[8];
				# Delete unit if below thresholds ($unit_qcov < 10 || $unit_alen < 300 || $unit_eval >= 0.0001)
				if ($unit_qcov < 10 || $unit_alen < 300 || $unit_eval >= 0.0001) {
					delete $unitd_hits->{$chr}->{$str}->{$unit};
					delete $unit_props->{$unit};
				}
			}
		}
	}
}

sub categorize_pseudogenes {
	# Take hit units and unit properties
	my($unitd_hits,$unit_props,$cd_genes) = @_;
	# Go through each chromosome
	foreach my $chr (sort keys %{$unitd_hits}) {
		# Go through each strand of chromosome
		foreach my $str (sort keys %{$unitd_hits->{$chr}}) {
			# Go through each hit unit
			foreach my $unit (sort {$a <=> $b} keys %{$unitd_hits->{$chr}->{$str}}) {
				# Get unit parent gene
				my $unit_gene = $unit_props->{$unit}->[0];
				# Get parent gene exon number
				my %gene_exons = ();
				foreach my $cds (@{$cd_genes->{$unit_gene}}) {
					$gene_exons{$cds->[8]}++;
				}
				my @gene_exons = values %gene_exons ? values %gene_exons : 1;
				my $gene_exons = BioStat::get_median(\@gene_exons);
				# Get (expected) unit exon number
				my $unit_exons = @{$unitd_hits->{$chr}->{$str}->{$unit}};
				my $unit_covrg = $unit_props->{$unit}->[7];
				my $expd_exons = int(($gene_exons*($unit_covrg/100))+0.5);
				# Categorize pseudogene sequence
				my $pgene_type = '';
				if ($unit_exons <= ($expd_exons/2)) {
					$pgene_type = 'processed_pseudogene';
				} else {
					$pgene_type = 'unprocessed_pseudogene';
				}
				#if ($unit_covrg < 34) { $pgene_type = 'pseudogene_fragment' }
				# Add pseudogene type to unit properties
				push(@{$unit_props->{$unit}},$pgene_type);
			}
		}
	}
	return $unit_props;
}

sub prepare_data {
	# Take hit units and unit properties
	my($unitd_hits,$unit_props) = @_;
	# Storage variables
	my %pseudogenes = ();
	my %pseudoexons = ();
	# Go through each chromosome
	foreach my $chr (sort keys %{$unitd_hits}) {
		# Go through each strand of chromosome
		foreach my $str (sort keys %{$unitd_hits->{$chr}}) {
			# Go through each hit unit
			foreach my $unit (sort {$a <=> $b} keys %{$unitd_hits->{$chr}->{$str}}) {
				# Transfer pseudogene data
				@{$pseudogenes{$chr}{$unit}} = @{$unit_props->{$unit}};
				# Transfer pseudo-exon data
				foreach my $hit (@{$unitd_hits->{$chr}->{$str}->{$unit}}) {
					push(@{$pseudoexons{$chr}{$unit}},$hit);
				}
			}
		}
	}
	return \(%pseudogenes,%pseudoexons);
}

sub print_gff_output {
	# Take pseudogenes and pseudoexons
	my($pseudogenes,$pseudoexons,$outfile) = @_;
	# Open output file
	my $out = FileIO::open_outfile($outfile);
	# Go through each chromosome
	foreach my $chr (sort keys %{$pseudogenes}) {
		# Go through each pseudogene (sort by location on chromosome)
		foreach my $pgene (sort {$pseudogenes->{$chr}->{$a}->[2] <=> $pseudogenes->{$chr}->{$b}->[2]} keys %{$pseudogenes->{$chr}}) {
			# Get pseudogene data
			my $pg_beg = $pseudogenes->{$chr}->{$pgene}->[2];
			my $pg_end = $pseudogenes->{$chr}->{$pgene}->[3];
			my $pg_str = $pseudogenes->{$chr}->{$pgene}->[4] eq 'plus' ? '+' : '-';
			my $pg_gid = $pseudogenes->{$chr}->{$pgene}->[0];
			my $pg_gsy = $pseudogenes->{$chr}->{$pgene}->[9];
			my $pg_typ = $pseudogenes->{$chr}->{$pgene}->[10];
			my $pg_ide = $pseudogenes->{$chr}->{$pgene}->[5];
			my $pg_cov = $pseudogenes->{$chr}->{$pgene}->[7];
			my $pg_len = $pseudogenes->{$chr}->{$pgene}->[8];
			my $pg_inf = 'ID='.$pgene;
			$pg_inf .= ';Dbxref=ParentGeneID:'.$pg_gid;
			$pg_inf .= ';Name='.$pg_gsy.';gbkey=Gene;parent_gene='.$pg_gsy;
			$pg_inf .= ';gene_biotype='.$pg_typ.';pseudo=true';
			print($out "$chr\tPseudon\tgene\t$pg_beg\t$pg_end\t.\t$pg_str\t.\t$pg_inf");
			printf($out "\t%.2f\t%d\t%d\n", $pg_ide,$pg_cov,$pg_len);
			# Go through each pseudo-exon
			foreach my $pexon (@{$pseudoexons->{$chr}->{$pgene}}) {
				# Get pseudo-exon data
				my($pe_beg,$pe_end) = sort {$a <=> $b} ($pexon->[9],$pexon->[10]);
				my $pe_ide = $pexon->[2];
				my $pe_inf = 'ID='.$pgene;
				$pe_inf .= ';Dbxref=ParentGeneID:'.$pg_gid;
				$pe_inf .= ';Name='.$pg_gsy.';gbkey=exon;parent_gene='.$pg_gsy;
				$pe_inf .= ';gene_biotype='.$pg_typ.';pseudo=true';
				print($out "$chr\tPseudon\texon\t$pe_beg\t$pe_end\t.\t$pg_str\t.\t$pe_inf");
				printf($out "\t%.2f\n", $pe_ide);
			}
		}
	}
}

sub print_elapsed_time {
	# Take start time t0
	my($t_0) = @_;
	# Get elapsed time
	my $sec = time-$t0;
	my $d = int($sec/(24*60*60));
	my $h = ($sec/(60*60))%24;
	my $m = ($sec/60)%60;
	my $s = $sec%60;
	# Prepare message
	my $message = '';
	if ($d) { $message .= $d == 1 ? "$d day, " : "$d days, " }
	if ($h) { $message .= $h == 1 ? "$h hour, " : "$h hours, " }
	if ($m) { $message .= $m == 1 ? "$m minute, " : "$m minutes, " };
	if ($s) { $message .= $s == 1 ? "$m second" : "$m seconds" };
	$message .= "\n";
	# Return time message
	print $message;
}