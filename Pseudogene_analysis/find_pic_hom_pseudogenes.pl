#!/usr/bin/perl

use strict;
use warnings;
use lib '/Users/dgebert/Dropbox/Perlmodules';
use FileIO;
use FastaIO;
use BioStat;

# Global constants
my $score_th = 0;
my $coverage_th = 50;
my $overlap_th = 0;
my $min_pg_len = 150;
# Global variables
my @homloc_seqs = ();
my @pic_pgenes = ();
my $glob_lin = 0;
my %pgene_lineages = ();
my %pgene_lin_dirs = ();
my %pgene_lin_gsys = ();
# Options
$|=1; #Autoflush

# Species identifiers
my @species = ('Hsap','Mmul','Mfas','Cjac','Mmur','Ltar');
my %spec_id = ();
foreach my $i (0..$#species) { $spec_id{$species[$i]} = $i }
# Genome identifiers
my @genomes = ('GRChg38','rheMac8','macFas5','calJac3','micMur3','otoGar3');
my %gnom_id = ();
foreach my $i (0..$#genomes) { $gnom_id{$genomes[$i]} = $i }

# Program name
print("\n--- $0 ---\n");
# Collect command line arguments
my $USAGE = "perl $0 <homlocs_tree_piCids.txt> <homloc_seq_files_list.txt> <piCpseudo_files_list.txt>\n";
unless ($ARGV[0]&&$ARGV[1]&&$ARGV[2]) {
	die("\nUsage: $USAGE\n");
}
# Input files
my $homloc_tree_file = $ARGV[0];
my $homloc_seq_files = $ARGV[1];
my $pic_pseudo_files = $ARGV[2];

# Get homlocs tree
my($tree, $branches) = get_tree_data($homloc_tree_file);

# Get sequence files
my @seq_files = FileIO::get_file_data_array($homloc_seq_files);
# Get homloc sequences
foreach my $seq_file (@seq_files) {
	# Get genome name and id
	my($gnom_name) = ($seq_file =~ /^([^_]+)_/);
	my $gnom_i = $gnom_id{$gnom_name};
	# Get sequence for the curr. genome
	$homloc_seqs[$gnom_i] = FastaIO::get_fasta_seqs($seq_file);
}

# Get pic pseudogene files
my @pseudo_files = FileIO::get_file_data_array($pic_pseudo_files);
# Get pic pseudogenes
foreach my $pseudo_file (@pseudo_files) {
	# Get species name and id
	my($spec_name) = ($pseudo_file =~ /^(\D+)\d*_/);
	my $spec_i = $spec_id{$spec_name};
	# Get pic pseudogenes for the curr. species
	$pic_pgenes[$spec_i] = get_pic_pseudogenes($pseudo_file);
}

### Main: Search for pseudogenes in homologous piCs
# Go through each homloc lineage
foreach my $br (@{$branches}) {
	# Storage variable for pseudogene positions in current homloc lineage
	my %pseudo_poss = ();
	my %pseudo_dirs = ();
	my %pseudo_gsys = ();
	my $pg_lin_id = 0;
	# Go through each x species
	foreach my $sp_x (0..$#species) {
		# Go through each loc of species x in homloc lineage
		foreach my $loc_x (@{$tree->[$sp_x]->{$br}}) {
			# Skip non-expressed and expressed p5 loci
			if ($loc_x !~ /\(.*\)/) { next }
			if ($loc_x =~ /-p5/) { next }
			# Get pic id(s)
			my($pic_ids) = ($loc_x =~ /\((.*)\)/);
			my @pic_ids = split(/;/,$pic_ids);
			# Go through each pic id
			foreach my $pic (@pic_ids) {
				# Check if pic containes pseudogenes
				unless ($pic_pgenes[$sp_x]->{$pic}) { next }
				# Get pic loc sequence
				my $loc_x_seq = $homloc_seqs[$sp_x]->{$loc_x};
				# Get loc x start coordinate
				my($loc_x_beg) = ($loc_x =~ /:(\d+)-\d+\s/);
				# Go through each pseudogene
				foreach my $pg (sort keys %{$pic_pgenes[$sp_x]->{$pic}}) {
					## Get pseudogene sequence
					my @pge_seqs = ();
					$pg_lin_id++;
					# Go through each pseudo-exon
					foreach my $pg_exon (@{$pic_pgenes[$sp_x]->{$pic}->{$pg}}) {
						# Get pseudo-exon coordinates
						my $pge_beg = $pg_exon->[3];
						my $pge_end = $pg_exon->[4];
						my $pge_len = $pge_end-$pge_beg+1;
						# Get relative coordinates
						my $pge_beg_rel = 0;
						if ($pge_beg > $loc_x_beg) { 
							$pge_beg_rel = $pge_beg-$loc_x_beg;
						} else {
							$pge_len = $pge_len+($pge_beg-$loc_x_beg);
						}
						# Get pseudogene sequence as substring of pic loc
						my $pge_seq = substr($loc_x_seq, $pge_beg_rel, $pge_len);
						push(@pge_seqs,$pge_seq);
						# Save pseudogene positions
						foreach my $pos ($pge_beg..$pge_end) {
							$pseudo_poss{$pg_lin_id}{$sp_x}{$loc_x}{$pos} = 1;
						}
						# Save directionality of pseudo-exon
						my $pge_dir = $pg_exon->[6];
						$pseudo_dirs{$pg_lin_id} = $pge_dir;
						# Save parent gene symbol of pseudo-exon
						my $pge_gsy = $pg_exon->[1];
						$pseudo_gsys{$pg_lin_id} = $pge_gsy;
					}
					## Write pseudogene sequence to file
					my $pg_file = 'Peudogene_seq.fas';
					my $pg_out = FileIO::open_outfile($pg_file);
					foreach my $seq (@pge_seqs) { print($pg_out ">$pg-PG\n$seq\n") }
					close($pg_out);
					## BLAST pseudogene sequence against homologous piC loci
					# Go through each y species
					foreach my $sp_y (0..$#species) {
						# Ignore if same species
						if ($sp_y eq $sp_x) { next }
						# Ignore if no homologs exist
						unless ($tree->[$sp_y]->{$br}) { next }
						# Bool: Pseudogene in homloc (false)
						my $pg_in_homloc_y = 0;
						# Blast pseudogene in piCs of this branch
						foreach my $loc_y (@{$tree->[$sp_y]->{$br}}) {
							# Get homloc sequence
							my $loc_y_seq = $homloc_seqs[$sp_y]->{$loc_y};
							# Write homloc sequence to file
							my $hl_file = 'Homloc_seq.fas';
							my $hl_out = FileIO::open_outfile($hl_file);
							print($hl_out ">$loc_y\n$loc_y_seq\n");
							close($hl_out);
							# Call blast
							my $bl_file = 'Pseudo_vs_homloc.blast.out';
							my $blast_opt = "-outfmt '7 qseqid sseqid pident evalue qcovs qlen qstart qend slen sstart send sstrand length bitscore' -task dc-megablast";
							system("blastn -query $pg_file -subject $hl_file -out $bl_file $blast_opt");
							# Analyze blast output file
							my @blast_data = FileIO::get_file_data_array($bl_file);
							my $hom_pg_y_loc = '';
							foreach my $line (@blast_data) {
								# Results line
								if ($line =~ /^[^#]\w+/) {
									# Get line data
									my @d = split(/\t/, $line);
									my $bitscore = $d[13];
									my $coverage = $d[4];
									# Get subject coordinates
									my @subj_loc = sort {$a <=> $b} ($d[9],$d[10]);
									my $subj_beg = $subj_loc[0];
									my $subj_end = $subj_loc[1];
									# Get homloc coordinates
									my($hom_chr,$hom_beg,$hom_end) = split(/:|-/,$d[1]);
									# Get real hit coordinates
									my $real_beg = $subj_beg+$hom_beg-1;
									my $real_end = $subj_end+$hom_beg-1;
									# Skip hit if score and coverage below thresholds
									if ($bitscore < $score_th || $coverage < $coverage_th) { next }
									# Save positions of subject
									foreach my $pos ($real_beg..$real_end) {
										$pseudo_poss{$pg_lin_id}{$sp_y}{$loc_y}{$pos} = 1;
									}
									#$pg_in_homloc_y = 1;
								}
							}
							# Going to next loc of species y to blast pseudogene sequence against
						}
						# Going to next species y to blast pseudogene sequence against
					}
					# Going to next pseudogene of curr. pic id of curr. loc of curr. species x in curr. homlocs lineage
				}
				# Going to next pic id of curr. loc of curr. species x in curr. homlocs lineage
			}
			# Going to next loc of curr. species x in curr. homlocs lineage
		}
		# Going to next species x in curr. homlocs lineage
	}
	## Check if pseudogene lineages overlap
	# Initialize lineage connections
	my @couples = ();
	foreach my $lin (sort keys %pseudo_poss) { push(@couples,[$lin]) }
	# Go through each lineage
	foreach my $lin_x (sort keys %pseudo_poss) {
		# Compare to each remaining lineage
		foreach my $lin_y (sort keys %pseudo_poss) {
			# Skip if same lineage
			if ($lin_x == $lin_y) { next }
			# Skip if pseudogenes have different orientations
			if ($pseudo_dirs{$lin_x} ne $pseudo_dirs{$lin_y}) { next }
			# Check if lineages overlap
			my $lins_overlap = 0;
			# Go through each species
			foreach my $sp (0..$#species) {
				# Check for identical piC loci in different lineages
				foreach my $loc_x (sort keys %{$pseudo_poss{$lin_x}{$sp}}) {
					foreach my $loc_y (sort keys %{$pseudo_poss{$lin_y}{$sp}}) {
						# Identical locus
						if ($loc_x eq $loc_y) {
							# Check if pseudogene positions overlap
							foreach my $pos (sort keys %{$pseudo_poss{$lin_x}{$sp}{$loc_x}}) {
								# Lineages overlap
								if ($pseudo_poss{$lin_y}{$sp}{$loc_y}{$pos}) {
									$lins_overlap = 1;
								}
							}
						}
					}
				}
			}
			# Connect overlapping lineages
			if ($lins_overlap) {
				push(@couples,[$lin_x,$lin_y]);
			}
		}
	}
	## Combine overlapping pseudogene lineages
	# Go through each lineage couple
	foreach my $c_1 (0..$#couples) {
		# Compare to each other couple
		foreach my $c_2 (0..$#couples) {
			# Skip if same couple index
			if ($c_1 == $c_2) { next }
			# Check if couples overlap
			my $couples_overlap = 0;
			foreach my $lin_c1 (@{$couples[$c_1]}) {
				if (grep {$lin_c1 eq $_} @{$couples[$c_2]}) {
					$couples_overlap = 1;
				}
			}
			# Combine couples if they overlap
			if ($couples_overlap) {
				push(@{$couples[$c_1]},@{$couples[$c_2]});
				@{$couples[$c_2]} = ();
			}
		}
	}
	## Get pseudogene presence/absence
	# Go through each combined lineage group
	foreach my $group (@couples) {
		# Skip empty groups
		unless (@{$group}) { next }
		# Increment global lineage id count
		$glob_lin++;
		# Get grouped lineages
		@{$group} = BioStat::uniq(sort(@{$group}));
		# Initialize pseudogene presence/absence status
		foreach my $sp (0..$#species) { $pgene_lineages{$glob_lin}{$sp} = 0 }
		# Go though each grouped lineage
		foreach my $lin (@{$group}) {
			# Go through each species
			foreach my $sp (0..$#species) {
				# Save presence status of global pseudogene lineages
				if (keys %{$pseudo_poss{$lin}{$sp}}) {
					$pgene_lineages{$glob_lin}{$sp} = 1;
				}
			}
			$pgene_lin_dirs{$glob_lin} = $pseudo_dirs{$lin};
			push(@{$pgene_lin_gsys{$glob_lin}},$pseudo_gsys{$lin});
		}
	}
	# Going to next homlocs lineage
}

## Get sums of different combinations
my %pres_combos = ();
# Go through each lineage
foreach my $lin (sort {$a <=> $b} keys %pgene_lineages) {
	# Get presence/absence pattern
	my $pattern = '';
	foreach my $sp (0..$#species) {
		$pattern .= $pgene_lineages{$lin}{$sp};
	}
	# Get directionality
	my $dir = $pgene_lin_dirs{$lin};
	$pres_combos{$dir}{$pattern}++;
}

# Output
my $outfile = 'pseudo_presence.txt';
my $out = FileIO::open_outfile($outfile);
foreach my $sp (0..$#species) {
	print($out "$species[$sp]\t");
}
print($out "\n");
foreach my $lin (sort {$a <=> $b} keys %pgene_lineages) {
	foreach my $sp (0..$#species) {
		print($out "$pgene_lineages{$lin}{$sp}\t");
	}
	print($out "$pgene_lin_dirs{$lin}\t");
	print($out "@{$pgene_lin_gsys{$lin}}\n");
}

my $outfile2 = 'pseudo_presence_sums.txt';
my $out2 = FileIO::open_outfile($outfile2);
my $outfile2_par = 'pseudo_presence_sums_par.txt';
my $out2_par = FileIO::open_outfile($outfile2_par);
my $outfile2_rev = 'pseudo_presence_sums_rev.txt';
my $out2_rev = FileIO::open_outfile($outfile2_rev);
foreach my $sp (0..$#species) {
	print($out2 "$species[$sp]\t");
	print($out2_par "$species[$sp]\t");
	print($out2_rev "$species[$sp]\t");
}
print($out2 "\n");
print($out2_par "\n");
print($out2_rev "\n");
foreach my $dir (sort keys %pres_combos) {
	foreach my $pat (sort keys %{$pres_combos{$dir}}) {
		my @pattern = split('',$pat);
		foreach my $stat (@pattern) {
			$stat = $stat*$pres_combos{$dir}{$pat};
			print($out2 "$stat\t");
		}
		print($out2 "$dir\n");
	}
}

foreach my $pat (sort keys %{$pres_combos{'par'}}) {
	my @pattern = split('',$pat);
	foreach my $stat (@pattern) {
		$stat = $stat*$pres_combos{'par'}{$pat};
		unless ($stat) { $stat = 'NA' }
		print($out2_par "$stat\t");
	}
	print($out2_par "\n");
}
foreach my $pat (sort keys %{$pres_combos{'rev'}}) {
	my @pattern = split('',$pat);
	foreach my $stat (@pattern) {
		$stat = $stat*$pres_combos{'rev'}{$pat};
		unless ($stat) { $stat = 'NA' }
		print($out2_rev "$stat\t");
	}
	print($out2_rev "\n");
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
	# Return tree and branches
	return \(@tree, @branches);
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
		# Get pic id
		if ($line =~ /^#/) {
			my @d = split(/\t/,$line);
			$pic_id = $d[0];
			$pic_id =~ s/#//;
		} 
		# Get pseudogene
		elsif ($line !~ /^$/) {
			my @d = split(/\t/,$line);
			my $parent_id = $d[0];
			$parent_id =~ s/\,.*//;
			my $parent_gs = $d[1];
			my $pgene_ide = $d[4];
			my $pgene_beg = $d[7];
			my $pgene_end = $d[8];
			my $pgene_str = $d[9];
			my $pgene_dir = $d[11];
			push(@{$pic_pgenes{$pic_id}{$parent_id}},[$parent_id,$parent_gs,$pgene_ide,$pgene_beg,$pgene_end,$pgene_str,$pgene_dir]);
		}
	}
	return \%pic_pgenes;
}

sub blast_pic_gene {
	# Take piC pseudogene info, homologous piCs, and both piC seqs
	my($pic_pseudo_x, $hom_pics_y, $pic_seqs_x, $pic_seqs_y) = @_;
	my @hom_pg_y_pics = ();
	# Thresholds, hit score and coverage
	my $score_th = 500;
	my $coverage_th = 30;
	# Get pseudogene and its piC
	my $pgene_x = $pic_pseudo_x->[2];
	my $pic_x = $pic_pseudo_x->[0];
	# Get pseudogene coordinates
	my $pgene_beg = $pic_pseudo_x->[4];
	my $pgene_end = $pic_pseudo_x->[5];
	# Calculate length of pseudogene
	my $pgene_len = $pgene_end-$pgene_beg+1;
	# Calculate start position of pseudogene substring
	my $pic_x_beg = $pic_seqs_x->{$pic_x}->[1];
	my $str_pos = $pgene_beg-$pic_x_beg+1;
	# Get pseudogene sequence from piC
	my $pgene_seq = substr($pic_seqs_x->{$pic_x}->[0], $str_pos, $pgene_len);
	# Blast pseudogene against homologous piCs
	foreach my $pic_y (@{$hom_pics_y}) {
		# Print pseudogene sequence to query file
		my $qfile = 'qfile';
		my $qfh = FileIO::open_outfile($qfile);
		print $qfh (">$pic_x\_$pgene_x\n$pgene_seq\n");
		close($qfh);
		# Print homologous piC sequence to subject file
		my $sfile = 'sfile';
		my $sfh = FileIO::open_outfile($sfile);
		print $sfh (">$pic_y\n$pic_seqs_y->{$pic_y}->[0]\n");
		close($sfh);
		# Blast output file
		my $ofile = 'ofile';
		# Call blastn
		system("blastn -query $qfile -subject $sfile -out $ofile -outfmt '7 qseqid sseqid pident evalue qcovs qlen qstart qend slen sstart send sstrand length bitscore' -perc_identity 50");
		# Analyze blast output file
		my @blast_data = FileIO::get_file_data_array($ofile);
		foreach my $line (@blast_data) {
			# Results line
			if ($line =~ /^[^#]\w+/) {
				my @d = split(/\t/, $line);
				my $score = $d[13];
				my $coverage = $d[4];
				my $pg_pic_y = $d[1];
				# Consider hit if score and coverage above thresholds
				if ($score > $score_th and $coverage > $coverage_th) {
					# Save homologous piC if not yet saved
					unless (grep { $_ eq $pg_pic_y } @hom_pg_y_pics) {
						push(@hom_pg_y_pics, $pg_pic_y);
					}
				}
			}
		}
		system("rm qfile sfile ofile");
	}
	return @hom_pg_y_pics;
}
