#!/usr/bin/perl
use strict;
use warnings;
use lib '/Users/dgebert/Dropbox/Perlmodules/';
use FileIO;
use FastaIO;

# Constants
my $flank_len = 5_000_000;	# Upstream and downstream bps to look for flank genes
my $open_flank = 5_000_000;	# Expand syntenic region when only one flank found
my $max_genes = 10; 		# Max number of genes per flank
my $syn_score_th = 4; 		# Min number of genes for syntenic region
my $gene_cov_th = 20; 		# Min coverage of flank gene hits
$|=1; #Autoflush
# Variables

# Program name
print("\n--- $0 ---\n");
# Collect command line arguments
my $USAGE = "perl $0 <sp1_piCseqs.fas> <sp1_genome.fas> <sp1_genome.gff> <sp2_genome.fas>\n";
unless ($ARGV[0] && $ARGV[1] && $ARGV[2] && $ARGV[3]) {
	die("\nUsage: $USAGE\n");
}
# Input files
my $sp1_pic_seqs_file = $ARGV[0];
my $sp1_genome_seq_file = $ARGV[1];
my $sp1_genome_gff_file = $ARGV[2];
my $sp2_genome_seq_file = $ARGV[3];

## Phase 1: Extract flanking genes up and downstream of pic
# Get piC locations
my $pic_locs = get_pic_locs($sp1_pic_seqs_file);
# Get gff gene and cds locations
my($gen_locs,$cds_locs) = get_gff_gene_locs($sp1_genome_gff_file);
# Get flanking genes
my $flank_genes = get_flank_genes($pic_locs,$gen_locs,$flank_len,$max_genes);
# Get flanking gene sequences
my $flank_gene_seqs = get_flank_gene_seqs($pic_locs,$cds_locs,$flank_genes,$sp1_genome_seq_file);

## Phase 2: Seek syntenic region with flank genes in second species
# Blast flank genes to genome and get hits
my $gen_hits = get_gene_blast_hits($flank_genes,$flank_gene_seqs,$sp1_pic_seqs_file,$sp2_genome_seq_file,$gene_cov_th);
# Group hits and get locations
my $gen_hit_locs = get_hit_group_locs($gen_hits);
# Get synteny scores for each hit region
my($syn_scores,$identities) = get_synteny_scores($gen_hit_locs,$flank_genes);
# Get syntenic loci and status
my($syn_locs,$syn_stat) = get_syntenic_locs($syn_scores,$identities,$flank_genes,$gen_hit_locs);

## Phase 3: Seek homologous pic sequence in syntenic region
# Blast pic sequences on syntenic loci
my $pic_hits = get_pic_blast_hits($sp1_pic_seqs_file,$sp2_genome_seq_file,$syn_locs);
# Sort and group hits and get best hit groups
my $pic_hits_top = get_best_pic_hits($pic_hits);
# Get pic homolog (best hit groups) coordinates
my $hom_locs = get_hom_pic_locs($pic_hits_top);
# Get pic homolog (best hit groups) properties
my $hom_props = get_hom_pic_props($pic_hits_top);

## Output results
my $outfile = 'homol_pic_locs.txt';
my $out = FileIO::open_outfile($outfile);
print($out "pic_num\t");
print($out "pic_id\tpic_loc\tpic_len\thom_loc\thom_len\thom_rel\taln_cov\taln_ide\taln_bit\thom_stat\tsyn_loc\tsyn_stat\n");
# Go through each pic
foreach my $pic (sort keys %{$pic_hits_top}) {
	# Get pic coordinates and length
	my $pic_chr = $pic_locs->{$pic}->[0];
	my $pic_beg = $pic_locs->{$pic}->[1];
	my $pic_end = $pic_locs->{$pic}->[2];
	my $pic_len = $pic_end-$pic_beg+1;
	# Get syntenic region coordinates and status
	my $syn_chr = $syn_locs->{$pic}->[0];
	my $syn_beg = $syn_locs->{$pic}->[1];
	my $syn_end = $syn_locs->{$pic}->[2];
	my $synteny = $syn_stat->{$pic};
	# Get homologous pic coordinates and length
	my $hom_chr = $hom_locs->{$pic}->[0];
	my $hom_beg = $hom_locs->{$pic}->[1];
	my $hom_end = $hom_locs->{$pic}->[2];
	my $hom_len = $hom_end-$hom_beg+1;
	my $hom_rel = $hom_len/$pic_len*100;
	# Get homolog alignment properties
	my $aln_len = $hom_props->{$pic}->[0];
	my $aln_ide = $hom_props->{$pic}->[1];
	my $aln_sco = $hom_props->{$pic}->[2];
	my $aln_cov = $aln_len/$pic_len*100;
	# Homology status
	my $homstat = 1;
	if ($aln_cov<5 || $hom_len<1500 || $hom_rel<15 || $synteny==0) { $homstat = 0 }
	my($pic_num) = ($pic =~ /^\D+_(\d+)/);
	# Print results
	printf($out "%d\t", $pic_num);
	printf($out "%s\t%s:%d-%d\t%d\t", $pic,$pic_chr,$pic_beg,$pic_end,$pic_len);
	printf($out "%s:%d-%d\t%d\t", $hom_chr,$hom_beg,$hom_end,$hom_len);
	printf($out "%.1f\t%.1f\t%.1f\t%.1f\t%d\t", $hom_rel,$aln_cov,$aln_ide,$aln_sco,$homstat);
	printf($out "%s:%d-%d\t%d\n", $syn_chr,$syn_beg,$syn_end,$synteny);
}

exit;

################################# subroutines #################################

sub get_pic_locs {
	# Take pic seq file name
	my($picseq_file) = @_;
	# Storage variable
	my %pic_locs = ();
	# Get pic infos (from fasta headers)
	my $pic_infos = FastaIO::get_fasta_headers($picseq_file);
	# Extract pic locations
	foreach my $info (@{$pic_infos}) {
		my @d = split(/\s|-/,$info);
		my $pic = $d[0];
		my $chr = $d[1];
		my $beg = $d[2];
		my $end = $d[3];
		# Save pic loc info
		@{$pic_locs{$pic}} = ($chr,$beg,$end);
	}
	return \%pic_locs;
}

sub get_gff_gene_locs {
	# Take infile name and keyword
	my($gff_file) = @_;
	# Get file data
	my $gff = FileIO::open_infile($gff_file);
	# Initialize variables
	my %gff_cds_locs = ();
	my %gff_gen_locs = ();
	my $gene = '';
	my $name = '';
	my $tran = '';
	# Parse gff file
	while (my $line = <$gff>) {
		# Results line
		if ($line =~ /^[^#]\w+/) {
			# Get info
			my @d = split(/\t/, $line);
			my $chr  = $d[0];
			my $type = $d[2];
			my $beg  = $d[3];
			my $end  = $d[4];
			my $info = $d[8];
			# Gene line
			if ($type =~ /gene/i) {
				# Skip non protein coding genes
				if ($info !~ /protein_coding/i) { next }
				# Get gene id (gff & gtf)
				if ($info =~ /GeneID:([^;]+)\;/) {
					($gene) = ($info =~ /GeneID:([^;]+)\;/);
				}
				if ($info =~ /gene_id\s\"([^\"]+)\";/) {
					($gene) = ($info =~ /gene_id\s\"([^\"]+)\";/);
				}
				# Get gene name (gff & gtf)
				if ($info =~ /gene_name\s"([^\"]+)"/) {
					($name) = ($info =~ /gene_name\s"([^\"]+)"/);
				} elsif ($info =~ /Name=(\w+?);/) {
					($name) = ($info =~ /Name=(\w+?);/);
				}
				# Save entries
				push(@d,$name);
				$gff_gen_locs{$chr}{$gene} = \@d;
			}
			# CDS line
			elsif ($type =~ /cds/i) {
				# Get transcript id (gff & gtf)
				if ($info =~ /Genbank:([^;]+)\;/) {
					($tran) = ($info =~ /Genbank:([^;]+)\;/);
				}
				if ($info =~ /transcript_id\s\"([^\"]+)\";/) {
					($tran) = ($info =~ /transcript_id\s\"([^\"]+)\";/);
				}
				# Save entries
				push(@d,$name);
				push(@{$gff_cds_locs{$gene}{$tran}},\@d);
			}
		}
	}
	return \(%gff_gen_locs,%gff_cds_locs);
}

sub get_flank_genes {
	# Take pic locs and gene locs
	my($pic_locs,$gen_locs,$flank_len,$max_genes) = @_;
	# Storage variables
	my %flank_genes = ();
	# Go through each pic and get flanking genes
	foreach my $pic (sort keys %{$pic_locs}) {
		# Get pic coordinates
		my $chr = $pic_locs->{$pic}->[0];
		my $pic_beg = $pic_locs->{$pic}->[1];
		my $pic_end = $pic_locs->{$pic}->[2];
		my $flank_5 = $pic_beg-$flank_len;
		my $flank_3 = $pic_end+$flank_len;
		# Initialize flanking genes arrays
		@{$flank_genes{$pic}[0]} = ();
		@{$flank_genes{$pic}[1]} = ();
		# Go through each gene on same chromosome and get flanking genes
		foreach my $gene (sort { $gen_locs->{$chr}->{$a}->[3] <=> $gen_locs->{$chr}->{$b}->[3]} keys %{$gen_locs->{$chr}}) {
			# Get gene coordinates
			my $gen_beg = $gen_locs->{$chr}->{$gene}->[3];
			my $gen_end = $gen_locs->{$chr}->{$gene}->[4];
			# Save upstream flanking genes
			if ($gen_end > $flank_5 && $gen_end < $pic_beg) {
				push(@{$flank_genes{$pic}[0]},$gene);
			}
			# Save downstream flanking genes
			elsif ($gen_beg > $pic_end && $gen_beg < $flank_3) {
				push(@{$flank_genes{$pic}[1]},$gene);
			}
		}
		# Reduce flanking genes to max
		while (scalar(@{$flank_genes{$pic}[0]}) > $max_genes) { shift(@{$flank_genes{$pic}[0]}) }
		while (scalar(@{$flank_genes{$pic}[1]}) > $max_genes) { pop(@{$flank_genes{$pic}[1]}) }
	}
	return \%flank_genes;
}

sub get_flank_gene_seqs {
	# Take pic locs, flank genes, cds locs, genome
	my($pic_locs,$cds_locs,$flank_genes,$genome_file) = @_;
	# Storage variables
	my %flank_gene_seqs = ();
	# Get reference genome sequence
	my $genome = FastaIO::get_fasta_seqs($genome_file,1);
	# Go through each pic
	foreach my $pic (sort keys %{$pic_locs}) {
		# Get chromosome
		my $chr = $pic_locs->{$pic}->[0];
		# For each flank
		for (my $i=0; $i<2; $i++) {
			# Go through each flanking gene
			foreach my $gene (@{$flank_genes->{$pic}->[$i]}) {
				# Go through each possible transcript and save sequence
				my %trans_seqs = ();
				foreach my $tran (sort keys %{$cds_locs->{$gene}}) {
					# Go through each cds exon and extract sequence
					my $cds_seq = '';
					foreach my $exon (@{$cds_locs->{$gene}->{$tran}}) {
						# Get exon loc
						my $beg = $exon->[3];
						my $end = $exon->[4];
						my $exon_len = $end-$beg+1;
						# Extract cds exon sequence from genome
						$cds_seq .= substr($genome->{$chr},($beg-1),$exon_len);
					}
					$trans_seqs{$tran} = $cds_seq;
				}
				# Choose largest transcript
				my $best_seq = '';
				foreach my $tran (sort {length($trans_seqs{$b}) <=> length($trans_seqs{$a})} keys %trans_seqs) {
					# Save largest peptide and go to next gene
					$best_seq = $trans_seqs{$tran} and last;
				}
				# Save gene cds sequence
				$flank_gene_seqs{$pic}[$i]{$gene} = $best_seq;
			}
		}
	}
	# Clear genome data
	%{$genome} = ();
	# Return flanking gene sequences
	return \%flank_gene_seqs;
}

sub get_gene_blast_hits {
	# Take flank genes, flank gene seqs and second genome file
	my($flank_genes,$flank_gene_seqs,$sp1_pic_seqs_file,$sp2_genome_seq_file,$cov_th) = @_;
	# Blast options
	my $blast_opts = "-outfmt '7 qseqid sseqid pident evalue qcovs qlen qstart qend slen sstart send sstrand length bitscore'";
	# Storage variable
	my %gen_hits = ();
	# Create directory as blast output folder
	my($sp2_genome_name) = ($sp2_genome_seq_file =~ /^([^\.]+)\./);
	my $out_folder = $sp1_pic_seqs_file.'_'.$sp2_genome_name.'_blast';
	unless(-e $out_folder or mkdir $out_folder) {
		die("Error: Unable to create folder \'$out_folder\'\n");
	}
	# Go through each pic
	foreach my $pic (sort keys %{$flank_genes}) {
		# Save flanking gene sequences to fasta file
		my $flank_genes_file = $out_folder.'/'.$pic.'.flanks.fas';
		my $out = FileIO::open_outfile($flank_genes_file);
		# Go through each flank
		for (my $i=0; $i<2; $i++) {
			# Go through each gene
			foreach my $gene (@{$flank_genes->{$pic}->[$i]}) {
				# Get gene sequence and print to file
				my $gene_seq = $flank_gene_seqs->{$pic}->[$i]->{$gene};
				print($out ">$gene\n$gene_seq\n");
			}
		}
		# Blast flank gene sequences to second species genome
		my $blastn_file = $flank_genes_file.'.'.$sp2_genome_name.'.bn.out';
		# If blast file already exists check if file is ok
		my $blastn_file_ok = 0;
		if (-e $blastn_file && -s $blastn_file) {
			my @file_data = FileIO::get_file_data_array($blastn_file);
			if ($file_data[-1] =~ /# BLAST processed \d+ queries/) {
				$blastn_file_ok = 1 if (stat $blastn_file)[7] > 300;
			}
		}
		unless ($blastn_file_ok) {
			system("blastn -query $flank_genes_file -subject $sp2_genome_seq_file -out $blastn_file $blast_opts");
		}
		$gen_hits{$pic} = get_blastn_gene_hits($blastn_file,$cov_th);
	}
	return \%gen_hits;
}

sub get_blastn_gene_hits {
	# Take infile name
	my($file,$cov_th) = @_;
	# Get file data
	my @file_data = FileIO::get_file_data_array($file);
	# Initialize variables
	my %blast_hits = ();
	# Parse tblastn output
	foreach my $line (@file_data) {
		# Results line
		if ($line =~ /^[^#]\w+/) {
			# Get info
			my @d = split(/\t/, $line);
			my $gen_id = $d[0];
			my $chr_id = $d[1];
			my $coverg = $d[4];
			# Save hits
			if ($coverg < $cov_th) { next }
			push(@{$blast_hits{$chr_id}{$gen_id}},\@d);
		}
	}
	return \%blast_hits;
}

sub get_hit_group_locs {
	# Take gene blast hits
	my($gen_hits) = @_;
	# Storage variable
	my %gen_hit_locs = ();
	# Go through each pic
	foreach my $pic (sort keys %{$gen_hits}) {
		# Go through each subject hit chromosome
		foreach my $chr (sort keys %{$gen_hits->{$pic}}) {
			# Go through each gene
			foreach my $gen (sort keys %{$gen_hits->{$pic}->{$chr}}) {
				# Group gene hits
				my $group = 0;
				my $prev_hit_end = 0;
				my $max_distance = 10_000;
				my %grouped_hits = ();
				my %group_score = ();
				my %group_locs = ();
				my %group_ident = ();
				my %group_match = ();
				my %group_mis_m = ();
				my %hit_strands = ();
				my %group_strand = ();
				# Go through each hit
				foreach my $hit (sort {$a->[9] <=> $b->[9]} @{$gen_hits->{$pic}->{$chr}->{$gen}}) {
					# Hit positions
					my $identity = $hit->[2];
					my $hit_beg  = $hit->[9];
					my $hit_end  = $hit->[10];
					my $sstrand  = $hit->[11];
					my $hit_len  = $hit->[12];
					my $bitscore = $hit->[13];
					# Adjust symbols and coordinates according to subject strand
					if ($hit_beg > $hit_end) {
						$hit_beg = $hit->[10];
						$hit_end = $hit->[9];
					}
					# Save strand
					$hit_strands{$sstrand}++;
					# Get distance
					my $distance = $hit_beg-$prev_hit_end;
					# If distance above threshold open new group
					if ($distance > $max_distance) { $group++ }
					# Allocate hit to group
					push(@{$grouped_hits{$group}},$hit);
					# Save group properties
					$group_score{$group} += $bitscore;
					$group_locs{$group}[0] = $hit_beg unless $group_locs{$group}[0];
					$group_locs{$group}[0] = $hit_beg if $hit_beg < $group_locs{$group}[0];
					$group_locs{$group}[1] = $hit_end unless $group_locs{$group}[1];
					$group_locs{$group}[1] = $hit_end if $hit_end > $group_locs{$group}[1];
					# Get matches and mismatches
					my $n_match = int(($hit_len*$identity/100)+0.5);
					my $n_mis_m = $hit_len-$n_match;
					$group_match{$group} += $n_match;
					$group_mis_m{$group} += $n_mis_m;
					# Save end position for next hit
					$prev_hit_end = $hit_end;
				}
				# Get group strand
				foreach my $gr (keys %group_match) {
					foreach my $str (sort {$hit_strands{$b} <=> $hit_strands{$a}} keys %hit_strands) {
						$group_strand{$gr} = $str and last;
					}
				}
				# Calculate identities
				foreach my $gr (keys %group_match) {
					$group_ident{$gr} = $group_match{$gr}/($group_match{$gr}+$group_mis_m{$gr})*100;
				}
				# Choose best group
				my $best_gr = 0;
				foreach my $gr (sort {$group_score{$b} <=> $group_score{$a}} keys %group_score) {
					$best_gr = $gr and last;
				}
				@{$gen_hit_locs{$pic}{$chr}{$gen}} = (@{$group_locs{$best_gr}},$group_ident{$best_gr},$group_strand{$best_gr});
			}
		}
	}
	return \%gen_hit_locs;
}

sub get_synteny_scores {
	# Take flank genes, flank gene seqs and second genome file
	my($gen_hit_locs,$flank_genes) = @_;
	# Score storage variables
	my %syn_scores = ();
	my %identities = ();
	# Get score for each hit region
	foreach my $pic (sort keys %{$gen_hit_locs}) {
		foreach my $chr (sort keys %{$gen_hit_locs->{$pic}}) {
			# Go through each gene
			my @genes_upst = @{$flank_genes->{$pic}->[0]};
			my @genes_down = @{$flank_genes->{$pic}->[1]};
			for my $g (reverse 0..$#genes_upst) {
				if ($gen_hit_locs->{$pic}->{$chr}->{$genes_upst[$g]}) {
					$syn_scores{$pic}{$chr}++;
					my $ident = $gen_hit_locs->{$pic}->{$chr}->{$genes_upst[$g]}->[2];
					$identities{$pic}{$chr} += $ident;
				}
			}
			for my $g (0..$#genes_down) {
				if ($gen_hit_locs->{$pic}->{$chr}->{$genes_down[$g]}) {
					$syn_scores{$pic}{$chr}++;
					my $ident = $gen_hit_locs->{$pic}->{$chr}->{$genes_down[$g]}->[2];
					$identities{$pic}{$chr} += $ident;
				}
			}
			$identities{$pic}{$chr} = $identities{$pic}{$chr}/$syn_scores{$pic}{$chr};
		}
	}
	return \(%syn_scores,%identities);
}

sub get_syntenic_locs {
	# Take synteny scores, identities, flanking genes and gene hit loci
	my($syn_scores,$identities,$flank_genes,$gen_hit_locs) = @_;
	# Storage variable
	my %syn_locs = ();
	my %syn_stat = ();
	# Go through each pic
	foreach my $pic (sort keys %{$syn_scores}) {
		# Get best hit chromosome
		my $syn_chr = '';
		my $top_score = 0;
		my $top_ident = 0;
		foreach my $chr (sort { $syn_scores->{$pic}->{$b} <=> $syn_scores->{$pic}->{$a} } keys %{$syn_scores->{$pic}}) {
			# Get top score
			if ($syn_scores->{$pic}->{$chr} > $top_score) {
				$top_score = $syn_scores->{$pic}->{$chr};
			}
			# Check if chromosome is a topscorer
			if ($syn_scores->{$pic}->{$chr} == $top_score) {
				# Get topscorer chromosome with highest identity
				if ($identities->{$pic}->{$chr} > $top_ident) {
					$top_ident = $identities->{$pic}->{$chr};
					$syn_chr = $chr;
				}
			}
		}
		# Get flank border coordinates
		my @genes_upst = @{$flank_genes->{$pic}->[0]};
		my @genes_down = @{$flank_genes->{$pic}->[1]};
		my($beg_upst,$end_upst,$beg_down,$end_down) = (0,0,0,0);
		# Upstream flanking genes
		for my $g (0..$#genes_upst) {
			my $beg = $gen_hit_locs->{$pic}->{$syn_chr}->{$genes_upst[$g]}->[0];
			my $end = $gen_hit_locs->{$pic}->{$syn_chr}->{$genes_upst[$g]}->[1];
			unless ($beg && $end) { next }
			# Update flank borders
			unless ($beg_upst) { $beg_upst = $beg }
			unless ($end_upst) { $end_upst = $end }
			if ($beg < $beg_upst) { $beg_upst = $beg }
			if ($end > $end_upst) { $end_upst = $end }
		}
		# Downstream flanking genes
		for my $g (0..$#genes_down) {
			my $beg = $gen_hit_locs->{$pic}->{$syn_chr}->{$genes_down[$g]}->[0];
			my $end = $gen_hit_locs->{$pic}->{$syn_chr}->{$genes_down[$g]}->[1];
			unless ($beg && $end) { next }
			# Update flank borders
			unless ($beg_down) { $beg_down = $beg }
			unless ($end_down) { $end_down = $end }
			if ($beg < $beg_down) { $beg_down = $beg }
			if ($end > $end_down) { $end_down = $end }
		}
		# Single flank bool
		my $single_flank = 0;
		# Get coordinates for syntenic region
		my $syn_beg = 0;
		my $syn_end = 0;
		# Both sytenic flanks present: true
		if ($beg_upst && $beg_down) {
			# Syntenic flank regions ovelap: true
			if ($end_upst >= $beg_down && $beg_upst <= $end_down) {
				# Get smallest and largest coordinate
				my @sort_pos = sort {$a <=> $b} ($beg_upst,$end_upst,$beg_down,$end_down);
				$syn_beg = $sort_pos[0];
				$syn_end = $sort_pos[-1];
			}
			# Syntenic flank regions ovelap: false
			else {
				# Plus orientation
				if ($beg_upst < $beg_down) {
					$syn_beg = $end_upst;
					$syn_end = $beg_down;
				}
				# Minus orientation
				elsif ($beg_upst > $beg_down) {
					$syn_beg = $end_down;
					$syn_end = $beg_upst;
				}
			}
		}
		# Both sytenic flanks present: false
		else {
			# Upstream flank present
			if ($beg_upst) {
				$syn_beg = $beg_upst-$open_flank;
				$syn_end = $end_upst+$open_flank;
				# Single flank bool: true
				$single_flank = 1;
			} 
			# Downstream flank present
			elsif ($beg_down) {
				$syn_beg = $beg_down-$open_flank;
				$syn_end = $end_down+$open_flank;
				# Single flank bool: true
				$single_flank = 1;
			}
			if ($syn_beg < 0) { $syn_beg = 1 }
		}
		# Save syntenic region location
		$syn_locs{$pic}[0] = $syn_chr;
		$syn_locs{$pic}[1] = $syn_beg;
		$syn_locs{$pic}[2] = $syn_end;
		
		# Get synteny status
		if ($syn_chr && $syn_scores->{$pic}->{$syn_chr} >= $syn_score_th) {
			if ($single_flank) {
				$syn_stat{$pic} = 1;
			} else {
				$syn_stat{$pic} = 2;
			}
		} else {
			$syn_stat{$pic} = 0;
		}
	}
	return \(%syn_locs,%syn_stat);
}

sub get_pic_blast_hits {
	# Take pic file, second genome file, syntenic loci
	my($sp1_pic_seqs_file,$sp2_genome_seq_file,$syn_locs) = @_;
	# Storage variable
	my %pic_hits = ();
	# Blast options
	my $blast_opts = "-outfmt '7 qseqid sseqid pident evalue qcovs qlen qstart qend slen sstart send sstrand length bitscore' -task dc-megablast";
	# Blast output folder
	my($sp2_genome_name) = ($sp2_genome_seq_file =~ /^([^\.]+)\./);
	my $out_folder = $sp1_pic_seqs_file.'_'.$sp2_genome_name.'_blast';
	# Get pic sequences
	my $pic_seqs = FastaIO::get_fasta_seqs($sp1_pic_seqs_file,1);
	# Get genome sequence
	my $genome = FastaIO::get_fasta_seqs($sp2_genome_seq_file,1);
	# Go through each pic
	foreach my $pic (sort keys %{$pic_seqs}) {
		# Print pic sequence to file
		my $pic_seq_file = $out_folder.'/'.$pic.'.seq.fas';
		my $out_p = FileIO::open_outfile($pic_seq_file);
		print($out_p ">$pic\n$pic_seqs->{$pic}\n");
		# Get sequence of syntenic region
		my $syn_chr = $syn_locs->{$pic}->[0];
		my $syn_beg = $syn_locs->{$pic}->[1];
		my $syn_end = $syn_locs->{$pic}->[2];
		my $syn_len = $syn_end-$syn_beg+1;
		my $syn_seq = substr($genome->{$syn_chr},($syn_beg-1),$syn_len);
		# Print syntenic sequence to file
		my $syn_seq_file = $out_folder.'/'.$pic.'.'.$sp2_genome_name.'.syn.fas';
		my $out_s = FileIO::open_outfile($syn_seq_file);
		print($out_s ">$syn_chr:$syn_beg-$syn_end\n$syn_seq\n");
		# Name blast output file
		my $blastn_file = $pic_seq_file.'.'.$sp2_genome_name.'.syn.dcm.out';
		# If blast file already exists check if file is ok
		my $blastn_file_ok = 0;
		if (-e $blastn_file && -s $blastn_file) {
			my @file_data = FileIO::get_file_data_array($blastn_file);
			if ($file_data[-1] =~ /# BLAST processed \d+ queries/) {
				$blastn_file_ok = 1 if (stat $blastn_file)[7] > 300;
			}
		}
		unless ($blastn_file_ok) {
			system("blastn -query $pic_seq_file -subject $syn_seq_file -out $blastn_file $blast_opts");
		}
		# Get blast hits
		$pic_hits{$pic} = get_blastn_pic_hits($blastn_file);
		# Filter hits repetitive in regions
		$pic_hits{$pic} = filter_repetitive_hits($pic_hits{$pic});
	}
	return \%pic_hits;
}

sub get_blastn_pic_hits {
	# Take infile name
	my($file) = @_;
	# Get file data
	my @file_data = FileIO::get_file_data_array($file);
	# Initialize variables
	my @blast_hits = ();
	# Parse tblastn output
	foreach my $line (@file_data) {
		# Results line
		if ($line =~ /^[^#]\w+/) {
			# Get info
			my @d = split(/\t/,$line);
			my $pic_id = $d[0];
			# Save hits
			push(@blast_hits,\@d);
		}
	}
	return \@blast_hits;
}

sub filter_repetitive_hits {
	# Take blast hits
	my($hits) = @_;
	# Constants
	my $max_reps_hits = 10;
	my $max_rep_share = 0.5;
	# Variables
	my @hits_norep = ();
	# Get pos repetitions
	my $pos_repetitions = get_hits_per_pos($hits);
	# Go through each hit
	foreach my $hit (@{$hits}) {
		# Get pic hit coordinates
		my $pic_beg = $hit->[6];
		my $pic_end = $hit->[7];
		my $hit_len = $hit->[12];
		## Discard gene hits located in a repeat
		my $pos_in_rep = 0; #count
		# Check if region is repetitive (not annotated repeats)
		foreach my $pos ($pic_beg..$pic_end) {
			if ($pos_repetitions->{$pos} >= $max_reps_hits) {
				$pos_in_rep++;
			}
		}
		my $rep_share = $pos_in_rep/$hit_len;
		# Save hit if share in repeat location below threshold
		if ($rep_share < $max_rep_share) {
			push(@hits_norep,$hit);
		}
	}
	return \@hits_norep;
}
## Get hits per position to find repetitive regions
sub get_hits_per_pos {
	# Take blast hits
	my($hits) = @_;
	# Variables
	my %pos_repetitions = ();
	# Go through each hit
	foreach my $hit (@{$hits}) {
		# Get pic hit coordinates
		my $pic_beg = $hit->[6];
		my $pic_end = $hit->[7];
		# Count hits on each position
		foreach my $pos ($pic_beg..$pic_end) {
			$pos_repetitions{$pos}++;
		}
	}
	return \%pos_repetitions;
}

sub get_best_pic_hits {
	# Take pic blast hits
	my($pic_hits) = @_;
	# Storage variable
	my %pic_hits_top = ();
	# Go through each pic
	foreach my $pic (sort keys %{$pic_hits}) {
		# Sort hits
		my %sorted_hits = ();
		foreach my $hit (sort { $a->[9] <=> $b->[9] } @{$pic_hits->{$pic}}) {
			# Hit data
			my @hit_loc = sort {$a <=> $b} ($hit->[9],$hit->[10]);
			my $hit_beg = $hit_loc[0];
			my $hit_end = $hit_loc[1];
			my $hit_str = $hit->[11];
			# Save sorted hit divided by strand
			push(@{$sorted_hits{$hit_str}},$hit);
		}
		# Group sorted hits
		my $group = 0;
		my %grouped_hits = ();
		my $max_distance = 10_000;
		foreach my $strand (sort keys %sorted_hits) {
			# Set previous end position to zero
			my $prev_end = 0;
			# Go through sorted hits
			foreach my $hit (@{$sorted_hits{$strand}}) {
				# Hit data
				my @hit_loc = sort {$a <=> $b} ($hit->[9],$hit->[10]);
				my $hit_beg = $hit_loc[0];
				my $hit_end = $hit_loc[1];
				# Get distance
				my $distance = abs($hit_beg-$prev_end);
				# If distance above threshold open new group
				if ($distance > $max_distance) { $group++ }
				# Allocate hit to group
				push(@{$grouped_hits{$group}},$hit);
				# Save end position for next hit
				$prev_end = $hit_end;
			}
			# New group after strand switch
			$group++;
		}
		# Get group lengths
		my %group_len = ();
		foreach my $gr (sort keys %grouped_hits) {
			# Go through grouped hits
			foreach my $hit (@{$grouped_hits{$gr}}) {
				# Hit length
				$group_len{$gr} += $hit->[12];
			}
		}
		# Get largest group
		my $best_group = 0;
		foreach my $gr (sort {$group_len{$b} <=> $group_len{$a}} keys %group_len) {
			$best_group = $gr and last;
		}
		# Save best group hits
		unless ($grouped_hits{$best_group}) { @{$grouped_hits{$best_group}} = () }
		@{$pic_hits_top{$pic}} = @{$grouped_hits{$best_group}};
	}
	return \%pic_hits_top;
}

sub get_hom_pic_locs {
	# Take best hit groups
	my($pic_hits_top) = @_;
	# Storage variabl
	my %hom_locs = ();
	# Go through each pic
	foreach my $pic (sort keys %{$pic_hits_top}) {
		my $hom_chr = 'N/A';
		my $hom_beg = 0;
		my $hom_end = 0;
		# Go through each hit
		foreach my $hit (@{$pic_hits_top->{$pic}}) {
			# Syntenic loc
			my @syn_loc = split(/:|-/,$hit->[1]);
			my $syn_chr = $syn_loc[0];
			my $syn_beg = $syn_loc[1];
			# Hit loc
			my @hit_loc = sort {$a <=> $b} ($hit->[9],$hit->[10]);
			my $hit_beg = $hit_loc[0]+$syn_beg;
			my $hit_end = $hit_loc[1]+$syn_beg;
			# Get group start and end
			$hom_chr = $syn_chr if $syn_chr;
			unless ($hom_beg) { $hom_beg = $hit_beg }
			unless ($hom_end) { $hom_end = $hit_end }
			if ($hit_beg < $hom_beg) { $hom_beg = $hit_beg }
			if ($hit_end > $hom_end) { $hom_end = $hit_end }
		}
		$hom_locs{$pic}[0] = $hom_chr;
		$hom_locs{$pic}[1] = $hom_beg;
		$hom_locs{$pic}[2] = $hom_end;
	}
	return \%hom_locs;
}

sub get_hom_pic_props {
	# Take best hit groups
	my($pic_hits_top) = @_;
	# Storage variable
	my %hom_props = ();
	# Go through each pic
	foreach my $pic (sort keys %{$pic_hits_top}) {
		# Alignment length
		my $aln_len = get_group_length($pic_hits_top->{$pic});
		# Alignment identity
		my $aln_ident = get_group_identity($pic_hits_top->{$pic});
		# Alignment bit score
		my $aln_score = get_group_score($pic_hits_top->{$pic});
		# Save properties
		@{$hom_props{$pic}} = ($aln_len,$aln_ident,$aln_score);
	}
	return \%hom_props;
}

sub get_group_length {
	# Take blast hit group
	my($hits) = @_;
	# Aligned positions
	my %aln_pos = ();
	# Go through each hit
	foreach my $hit (@{$hits}) {
		# Syntenic loc
		my @syn_loc = split(/:|-/,$hit->[1]);
		my $syn_beg = $syn_loc[1];
		# Hit loc
		my @hit_loc = sort {$a <=> $b} ($hit->[9],$hit->[10]);
		my $hit_beg = $hit_loc[0]+$syn_beg;
		my $hit_end = $hit_loc[1]+$syn_beg;
		# Get alignment positions
		foreach my $pos ($hit_beg..$hit_end) { $aln_pos{$pos}++ }
	}
	# Alignment length
	my $aln_len = keys %aln_pos;
	return $aln_len;
}

sub get_group_identity {
	# Take blast hit group
	my($hits) = @_;
	# Calculate mean identity
	my $ident_len_sum = 0;
	my $a_len_sum = 0;
	my $aln_ident = 0;
	foreach my $hit (@{$hits}) {
		my $ident = $hit->[2];
		my $a_len = $hit->[12];
		$ident_len_sum += ($ident*$a_len);
		$a_len_sum += $a_len;
	}
	$aln_ident = $ident_len_sum/$a_len_sum if $a_len_sum;
	return $aln_ident;
}

sub get_group_score {
	# Take blast hit group
	my($hits) = @_;
	# Get score sum
	my $aln_score = 0;
	foreach my $hit (@{$hits}) {
		$aln_score += $hit->[13];
	}
	return $aln_score;
}