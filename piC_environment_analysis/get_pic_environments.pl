#!/usr/bin/perl
use strict;
use warnings;
use lib '/Users/dgebert/Dropbox/Perlmodules/';
use FileIO;
use FastaIO;
use BioStat;
use GradientRGB;
use GD;

# Constants
my $res = 1_000_000;
my $draw_heatmap = 0;
$|=1; #Autoflush
# Variables
my %rep_content_num = ();
my %rep_content_per = ();
my %rep_content_cla = ();
my %rep_content_div = ();
my %pic_content_num = ();
my %pic_content_num_genic = ();
my %pic_content_num_inter = ();
my %pic_content_num_psgen = ();
my %pic_content_num_gfree = ();
my %gen_content_num = ();
my %psg_content_num = ();
my @rep_classes = ();
my %count = ();
my %count_pp = ();
my %gnm_base_count = ();
my @gfree_pics = ();
my %centromeric_int = ();

# Program name
print("\n--- $0 ---\n");
# Collect command line arguments
my $USAGE = "perl $0 <piCseqs.fas> <genome.fas> <genome.gff> <genome.out> <centromers.txt>\n";
unless ($ARGV[0] && $ARGV[1]) {
	die("\nUsage: $USAGE\n");
}
# Input files
my $pic_seqs_file = $ARGV[0];
my $genome_file = $ARGV[1];
my $gene_gff_file = $ARGV[2];
my $repeat_file = $ARGV[3];
my $centro_file = $ARGV[4];

my $cen_locs = get_centromer_locs($centro_file);
foreach my $chr (sort keys %{$cen_locs}) {
	print("$chr:$cen_locs->{$chr}->[0]-$cen_locs->{$chr}->[1]\n");
}

# Extract file content
print("Extracting data...");
# Get piC locations
my $pic_locs = get_pic_locs($pic_seqs_file);
my $pic_seqs = FastaIO::get_fasta_seqs($pic_seqs_file,1);
# Get genome and chromosome names
my $genome = FastaIO::get_fasta_seqs($genome_file);
my $chroms = FastaIO::get_fasta_headers($genome_file);
# Get repeat locations
my $repeats = get_rm_hits($repeat_file);
# Get gff gene locations
my($gen_locs,$cds_locs) = get_gff_gene_locs($gene_gff_file);
# Get centromer locations
my $cen_locs = get_centromer_locs($centro_file);
print("ok\n");

# Get content
print("Processing data...");
# Go through each chromosome
foreach my $chr (@{$chroms}) {
	# Skip small scaffolds
	my $chr_len = length($genome->{$chr});
	if ($chr_len < $res/5) { next }
	## GC content
	# Slide through chromosome
	for (my $win=0; $win<=$chr_len; $win+=$res) {
		# Get interval
		my $interval = int($win/$res);
		# Get window sequence
		my $win_seq = substr($genome->{$chr}, $win, $res);
		# Get GC/AT share
		while ($win_seq =~ /[AT]/ig){$gnm_base_count{$chr}{$interval}{'AT'}++};
		while ($win_seq =~ /[GC]/ig){$gnm_base_count{$chr}{$interval}{'GC'}++};
	}
	## Repeats
	# Go through each repeat on the chromosome
	foreach my $hit (@{$repeats->{$chr}}) {
		# Get hit coordinates
		my $rep_div = $hit->[1];
		my $rep_beg = $hit->[5];
		my $rep_end = $hit->[6];
		my $rep_fam = $hit->[9];
		my $rep_cla = $hit->[10];
		unless (grep {$rep_cla eq $_} @rep_classes) {
			push(@rep_classes,$rep_cla);
		}
		# Get hit length
		my $rep_len = $rep_end-$rep_beg+1;
		# Get interval
		my $interval = int($rep_beg/$res);
		# Save data for each interval
		$rep_content_div{$chr}{$interval} += $rep_div;
		$rep_content_num{$chr}{$interval}++;
		$rep_content_per{$chr}{$interval} += $rep_len/$res;
		$rep_content_cla{$chr}{$interval}{$rep_cla} += $rep_len/$res;
	}
	## piRNA clusters
	# Go through each pic
	foreach my $pic (sort keys %{$pic_locs}) {
		# Get pic coordinates
		my $pic_chr = $pic_locs->{$pic}->[0];
		my $pic_beg = $pic_locs->{$pic}->[1];
		my $pic_end = $pic_locs->{$pic}->[2];
		my $pic_len = $pic_end-$pic_beg+1;
		unless ($pic_chr eq $chr) { next }
		# Get interval
		my $interval = int($pic_beg/$res);
		# Save data for each interval
		$pic_content_num{$chr}{$interval}++;
		# Get gene content of pic
		my %pic_gen_exon_pos = ();
		my %pic_psg_exon_pos = ();
		# Go through each gene on chromosome
		foreach my $gene (sort { $gen_locs->{$chr}->{$a}->[3] <=> $gen_locs->{$chr}->{$b}->[3]} keys %{$gen_locs->{$chr}}) {
			# Get gene coordinates
			my $gen_beg = $gen_locs->{$chr}->{$gene}->[3];
			my $gen_end = $gen_locs->{$chr}->{$gene}->[4];
			my $gen_inf = $gen_locs->{$chr}->{$gene}->[8];
			# Check if gene lies in pic
			if ($gen_beg <= $pic_end && $gen_end >= $pic_beg) {
				# Go through transcript of the current gene
				foreach my $trans (sort keys %{$cds_locs->{$gene}}) {
					# Go through exon of the current transcript
					foreach my $exon (@{$cds_locs->{$gene}->{$trans}}) {
						# Get exon coordinates
						my $ex_beg = $exon->[3];
						my $ex_end = $exon->[4];
						my $ex_inf = $exon->[8];
						# Check if exon lies in pic
						if ($ex_beg <= $pic_end && $ex_end >= $pic_beg) {
							# Save exon positions on pic
							foreach my $pos ($ex_beg..$ex_end) {
								# Skip if position not in pic
								unless ($pos >= $pic_beg && $pos <= $pic_end) { next }
								$pic_gen_exon_pos{$pos} = 1 if $gen_inf =~ /protein_coding/;
								$pic_psg_exon_pos{$pos} = 1 if $gen_inf =~ /pseudogene/;
							}
						}
					}
				}
			}
		}
		# Get gene exon share
		my $gen_exo_len = 0;
		foreach my $pos (keys %pic_gen_exon_pos) { $gen_exo_len++ }
		my $gen_exon_share = $gen_exo_len/$pic_len*100;
		# Get pseudogene exon share
		my $psg_exo_len = 0;
		foreach my $pos (keys %pic_psg_exon_pos) { $psg_exo_len++ }
		my $psg_exon_share = $psg_exo_len/$pic_len*100;
		# Seperate genic and intergenic pics
		if ($gen_exon_share > 0) {
			$pic_content_num_genic{$chr}{$interval}++;
		} else {
			$pic_content_num_inter{$chr}{$interval}++;
		}
		# Separate pseudogenic and inter(pseudo)genic pics
		if ($gen_exon_share == 0 && $psg_exon_share > 0) {
			$pic_content_num_psgen{$chr}{$interval}++;
		} elsif ($gen_exon_share == 0 && $psg_exon_share == 0) {
			$pic_content_num_gfree{$chr}{$interval}++;
			# Get GC/AT share
			while ($pic_seqs->{$pic} =~ /[AT]/ig){$count{'AT'}++};
			while ($pic_seqs->{$pic} =~ /[GC]/ig){$count{'GC'}++};
			while ($pic_seqs->{$pic} =~ /[AT]/ig){$count_pp{$pic}{'AT'}++};
			while ($pic_seqs->{$pic} =~ /[GC]/ig){$count_pp{$pic}{'GC'}++};
			push(@gfree_pics,$pic);
		}
	}
	## Genes
	# Go through each gene on chromosome
	foreach my $gene (sort { $gen_locs->{$chr}->{$a}->[3] <=> $gen_locs->{$chr}->{$b}->[3]} keys %{$gen_locs->{$chr}}) {
		# Get gene coordinates
		my $gen_beg = $gen_locs->{$chr}->{$gene}->[3];
		my $gen_end = $gen_locs->{$chr}->{$gene}->[4];
		my $gen_inf = $gen_locs->{$chr}->{$gene}->[8];
		# Get interval
		my $interval = int($gen_beg/$res);
		# Save data for each interval
		$gen_content_num{$chr}{$interval}++ if $gen_inf =~ /protein_coding/;
		$psg_content_num{$chr}{$interval}++ if $gen_inf =~ /pseudogene/;
	}
	## Centromer share
	my $cen_beg = $cen_locs->{$chr}->[0];
	my $cen_end = $cen_locs->{$chr}->[1];
	# Slide through chromosome
	for (my $win=0; $win<=$chr_len; $win+=$res) {
		# Get interval
		my $interval = int($win/$res);
		# Get window position
		my $win_beg = $win;
		my $win_end = $win+$res;
		# Check if window in centromer region
		if ($cen_beg < $win_end && $cen_end > $win_beg) {
			# Window fully in centromer
			if ($cen_beg < $win_beg && $cen_end > $win_end) {
				$centromeric_int{$chr}{$interval} = 1;
			} 
			# Window not fully in centromer
			else {
				# Get length of centromeric part of window
				my $cen_part = 0;
				# Window end in centromer
				if ($cen_beg > $win_beg) {
					$cen_part = $win_end-$cen_beg;
				}
				# Window beg in centromer
				elsif ($cen_end < $win_end) {
					$cen_part = $cen_end-$win_beg;
				}
				# Get centromer share of window
				my $cen_share = $cen_part/$res*100;
				# Declare window centromeric if share above threshold
				if ($cen_share < 67) {
					$centromeric_int{$chr}{$interval} = 1;
				} else {
					$centromeric_int{$chr}{$interval} = 0;
				}
			}
		} else {
			$centromeric_int{$chr}{$interval} = 0;
		}
	}
	
}
print("ok\n");

# Extract file content
print("Creating output...");
# Print output
my $pic_out = FileIO::open_outfile('GenesMb_pics.txt');
my $gnm_out = FileIO::open_outfile('GenesMb_gnom.txt');
my $pic_g_out = FileIO::open_outfile('GenesMb_pics_genic.txt');
my $pic_i_out = FileIO::open_outfile('GenesMb_pics_inter.txt');
my $pic_p_out = FileIO::open_outfile('GenesMb_pics_psgen.txt');
my $pic_f_out = FileIO::open_outfile('GenesMb_pics_gfree.txt');
my $tes_out = FileIO::open_outfile('GenesMb_tes.txt');
my $pic_gc_mb_out = FileIO::open_outfile('GenesMb_gc_pics.txt');
my $gnm_gc_mb_out = FileIO::open_outfile('GenesMb_gc_gnom.txt');
print($tes_out "Genes\tTE_div\t");
foreach my $cla (sort @rep_classes) {
	print($tes_out "$cla\t");
}
print($tes_out "\n");
my $pic_r_out = FileIO::open_outfile('GenesMb_tediv_pics.txt');
my $gnm_r_out = FileIO::open_outfile('GenesMb_tediv_gnom.txt');
my $pic_a_out = FileIO::open_outfile('GenesMb_alupc_pics.txt');
my $gnm_a_out = FileIO::open_outfile('GenesMb_alupc_gnom.txt');
my $pic_l_out = FileIO::open_outfile('GenesMb_l1pc_pics.txt');
my $gnm_l_out = FileIO::open_outfile('GenesMb_l1pc_gnom.txt');
my $gnm_p_out = FileIO::open_outfile('GenesMb_pseudo_gnom.txt');
# Go through each chromosome
foreach my $chr (@{$chroms}) {
	foreach my $int (0..int((length($genome->{$chr}))/$res)) {
		unless ($gen_content_num{$chr}{$int}) { $gen_content_num{$chr}{$int} = 0 }
		unless ($psg_content_num{$chr}{$int}) { $psg_content_num{$chr}{$int} = 0 }
		unless ($rep_content_cla{$chr}{$int}{'SINE/Alu'}) { next }
		if ($centromeric_int{$chr}{$int}) { next }
		# Gene density vs piCs
		print($pic_out "$gen_content_num{$chr}{$int}\n" x $pic_content_num{$chr}{$int}) if $pic_content_num{$chr}{$int};
		print($pic_g_out "$gen_content_num{$chr}{$int}\n" x $pic_content_num_genic{$chr}{$int}) if $pic_content_num_genic{$chr}{$int};
		print($pic_i_out "$gen_content_num{$chr}{$int}\n" x $pic_content_num_inter{$chr}{$int}) if $pic_content_num_inter{$chr}{$int};
		print($pic_p_out "$gen_content_num{$chr}{$int}\n" x $pic_content_num_psgen{$chr}{$int}) if $pic_content_num_psgen{$chr}{$int};
		print($pic_f_out "$gen_content_num{$chr}{$int}\n" x $pic_content_num_gfree{$chr}{$int}) if $pic_content_num_gfree{$chr}{$int};
		print($gnm_out "$gen_content_num{$chr}{$int}\n");
		print($gnm_p_out "$psg_content_num{$chr}{$int}\n");
		# Gene density vs repeats
		my $mean_div = $rep_content_div{$chr}{$int}/$rep_content_num{$chr}{$int} if $rep_content_num{$chr}{$int};
		unless ($mean_div) { $mean_div = 0 }
		print($tes_out "$gen_content_num{$chr}{$int}\t$mean_div\t");
		foreach my $cla (sort @rep_classes) {
			unless ($rep_content_cla{$chr}{$int}{$cla}) { $rep_content_cla{$chr}{$int}{$cla} = 0 }
			#push(@{$data_cla{$cla}},$rep_content_cla{$chr}{$int}{$cla});
			print($tes_out "$rep_content_cla{$chr}{$int}{$cla}\t");
		}
		print($tes_out "\n");
		# Repeats vs piCs
		print $pic_r_out "$mean_div\n" x $pic_content_num{$chr}{$int} if $pic_content_num{$chr}{$int};
		print $gnm_r_out "$mean_div\n";
		print $pic_a_out "$rep_content_cla{$chr}{$int}{'SINE/Alu'}\n" x $pic_content_num{$chr}{$int} if $pic_content_num{$chr}{$int};
		print $gnm_a_out "$rep_content_cla{$chr}{$int}{'SINE/Alu'}\n";
		print $pic_l_out "$rep_content_cla{$chr}{$int}{'LINE/L1'}\n" x $pic_content_num{$chr}{$int} if $pic_content_num{$chr}{$int};
		print $gnm_l_out "$rep_content_cla{$chr}{$int}{'LINE/L1'}\n";
		# GC content per mb
		my $gc_share_mb = $gnm_base_count{$chr}{$int}{'GC'}/($gnm_base_count{$chr}{$int}{'AT'}+$gnm_base_count{$chr}{$int}{'GC'})*100;
		print($pic_gc_mb_out "$gc_share_mb\n" x $pic_content_num{$chr}{$int}) if $pic_content_num_gfree{$chr}{$int};
		print($gnm_gc_mb_out "$gc_share_mb\n");
	}
}

print("ok\n");

my $pic_gc_share = $count{'GC'}/($count{'AT'}+$count{'GC'})*100;
my $pic_at_share = $count{'AT'}/($count{'AT'}+$count{'GC'})*100;

my $pic_gc_out = FileIO::open_outfile('GC_contents_pics.txt');
foreach my $pic (sort keys %count_pp) {
	my $gc_share = $count_pp{$pic}{'GC'}/($count_pp{$pic}{'AT'}+$count_pp{$pic}{'GC'})*100;
	print($pic_gc_out "$gc_share\n");
}

print("piC GC share: $pic_gc_share\npiC AT share: $pic_at_share\n");

# Get repeat loci
my $rep_data = get_repeatmask_data($repeat_file,$pic_locs);
# Get repeat divergence
my $mean_div = get_mean_repeat_divergence($rep_data);
# Get genome slices size
my $pic_size = get_slices_size($pic_locs);
# Get repeat share
my $rep_share = get_total_repeat_share($rep_data,$pic_size);

my %pic_locs_gfree = ();
foreach my $pic (@gfree_pics) {
	$pic_locs_gfree{$pic} = $pic_locs->{$pic};
}

# Get repeat loci
my $rep_data_gf = get_repeatmask_data($repeat_file,\%pic_locs_gfree);
# Get repeat divergence
my $mean_div_gf = get_mean_repeat_divergence($rep_data_gf);
# Get genome slices size
my $pic_size_gf = get_slices_size(\%pic_locs_gfree);
# Get repeat share
my $rep_share_gf = get_total_repeat_share($rep_data_gf,$pic_size_gf);

print("piC TE share: $rep_share\npiC TE divergence: $mean_div\n");
print("i.piC TE share: $rep_share_gf\ni.piC TE divergence: $mean_div_gf\n");

exit unless $draw_heatmap;

# Extract file content
print("Creating image...");

## Create images
my $chr = 'chr6';
my @data_gen = ();
my @data_psg = ();
my @data_div = ();
my @data_alu = ();
my @data_sva = ();
my @data_erv = ();
foreach my $int (0..int((length($genome->{$chr}))/$res)) {
	push(@data_gen,$gen_content_num{$chr}{$int});
	push(@data_psg,$psg_content_num{$chr}{$int});
	my $mean_div = $rep_content_div{$chr}{$int}/$rep_content_num{$chr}{$int} if $rep_content_num{$chr}{$int};
	unless ($mean_div) { $mean_div = 0 }
	push(@data_div,$mean_div);
	push(@data_alu,$rep_content_cla{$chr}{$int}{'SINE/Alu'});
	push(@data_sva,$rep_content_cla{$chr}{$int}{'Retroposon/SVA'});
	push(@data_erv,$rep_content_cla{$chr}{$int}{'LTR/ERVK'});
}
# Initialize drawing coordinates
my $width = 1000;
my $hight = $width/4;
my $s = $width/20;
my $block_hight = $width/25;
my $ay1 = $s;
my $ay2 = $s+$block_hight;
my $by1 = $s*6-$s;
my $by2 = $s*4+$s-$block_hight;
my $scale_y1 = $s*6;
my $scale_y2 = $scale_y1+5;
my $text_y1 = $scale_y1+10;
my $a_max = length($genome->{$chr});
# Image object and color allocation
my $image = new GD::Image(($width+4*$s),($hight+2*$s));
my $white = $image->colorAllocate(255, 255, 255);
my $black = $image->colorAllocate(0,   0,   0  );
# Draw region spanning gray blocks
my @rgb_high = (244, 241, 21);
my @rgb_low = (7, 0, 191);
my @gradient = GradientRGB::make_rgb_gradient(11,\@rgb_high,\@rgb_low);
my @colors;
foreach my $i (0..$#gradient) {
	$colors[$i] = $image->colorAllocate($gradient[$i]->[0],$gradient[$i]->[1],$gradient[$i]->[2]);
}
# Genes
my $dist = 0;
foreach my $int (0..int((length($genome->{$chr}))/$res)) {
	my $data_index = $int;
	$int = $int*$res;
	my $ax1 = $int/$a_max*$width;
 	my $ax2 = $int/$a_max*$width;
 	# Calculate hits fraction
 	my $data_min = BioStat::get_minimum(\@data_gen);
 	my $data_max = BioStat::get_maximum(\@data_gen);
 	#print("$data_index\t$data[$data_index]\t$data_min\t$data_max\n");
	my $fract = GradientRGB::calculate_fraction($data_gen[$data_index],$data_min,$data_max);
	my $rgb_i = GradientRGB::fract_to_rgb_gradient_color_index($fract,11,\@rgb_high,\@rgb_low);
	$image->filledRectangle($ax1+$s, $ay1+$dist, $ax2+($s*1.1), $ay2+$dist, $colors[$rgb_i]);
	#if ($data_index == 60) {$image->filledRectangle($ax1+$s, $ay1+$dist, $ax2+($s*1.1), $ay2+$dist, $white)}
}
$image->rectangle($s, $ay1, $width+$s, $ay2, $black);
# Pseudogenes
$dist += 42;
foreach my $int (0..int((length($genome->{$chr}))/$res)) {
	my $data_index = $int;
	$int = $int*$res;
	my $ax1 = $int/$a_max*$width;
 	my $ax2 = $int/$a_max*$width;
 	# Calculate hits fraction
 	my $data_min = BioStat::get_minimum(\@data_psg);
 	my $data_max = BioStat::get_maximum(\@data_psg);
 	#print("$data_index\t$data[$data_index]\t$data_min\t$data_max\n");
	my $fract = GradientRGB::calculate_fraction($data_psg[$data_index],$data_min,$data_max);
	my $rgb_i = GradientRGB::fract_to_rgb_gradient_color_index($fract,11,\@rgb_high,\@rgb_low);
	$image->filledRectangle($ax1+$s, $ay1+$dist, $ax2+($s*1.1), $ay2+$dist, $colors[$rgb_i]);
	#if ($data_index == 60) {$image->filledRectangle($ax1+$s, $ay1+$dist, $ax2+($s*1.1), $ay2+$dist, $white)}
}
$image->rectangle($s, $ay1+$dist, $width+$s, $ay2+$dist, $black);
# TE divergence
$dist += 42;
foreach my $int (0..int((length($genome->{$chr}))/$res)) {
	my $data_index = $int;
	$int = $int*$res;
	my $ax1 = $int/$a_max*$width;
 	my $ax2 = $int/$a_max*$width;
 	my $min = BioStat::get_minimum(\@data_div);
 	my $med = BioStat::get_median(\@data_div);
 	my @data_div_c = @data_div;
 	foreach my $i (0..$#data_div) {
 		$data_div_c[$i] = $med if $data_div[$i] == $min;
 	}
 	# Calculate hits fraction
 	my $data_min = BioStat::get_minimum(\@data_div_c);
 	my $data_max = BioStat::get_maximum(\@data_div_c);
 	#print("$data_index\t$data[$data_index]\t$data_min\t$data_max\n");
	my $fract = GradientRGB::calculate_fraction($data_div[$data_index],$data_min,$data_max);
	my $rgb_i = GradientRGB::fract_to_rgb_gradient_color_index($fract,11,\@rgb_high,\@rgb_low);
	$image->filledRectangle($ax1+$s, $ay1+$dist, $ax2+($s*1.1), $ay2+$dist, $colors[$rgb_i]);
	#if ($data_index == 60) {$image->filledRectangle($ax1+$s, $ay1+$dist, $ax2+($s*1.1), $ay2+$dist, $white)}
}
$image->rectangle($s, $ay1+$dist, $width+$s, $ay2+$dist, $black);
# Alu share
$dist += 42;
foreach my $int (0..int((length($genome->{$chr}))/$res)) {
	my $data_index = $int;
	$int = $int*$res;
	my $ax1 = $int/$a_max*$width;
 	my $ax2 = $int/$a_max*$width;
 	# Calculate hits fraction
 	my $data_min = BioStat::get_minimum(\@data_alu);
 	my $data_max = BioStat::get_maximum(\@data_alu);
 	#print("$data_index\t$data[$data_index]\t$data_min\t$data_max\n");
	my $fract = GradientRGB::calculate_fraction($data_alu[$data_index],$data_min,$data_max);
	my $rgb_i = GradientRGB::fract_to_rgb_gradient_color_index($fract,11,\@rgb_high,\@rgb_low);
	$image->filledRectangle($ax1+$s, $ay1+$dist, $ax2+($s*1.1), $ay2+$dist, $colors[$rgb_i]);
	#if ($data_index == 60) {$image->filledRectangle($ax1+$s, $ay1+$dist, $ax2+($s*1.1), $ay2+$dist, $white)}
}
$image->rectangle($s, $ay1+$dist, $width+$s, $ay2+$dist, $black);
# SVA share
$dist += 42;
foreach my $int (0..int((length($genome->{$chr}))/$res)) {
	my $data_index = $int;
	$int = $int*$res;
	my $ax1 = $int/$a_max*$width;
 	my $ax2 = $int/$a_max*$width;
 	# Calculate hits fraction
 	my $data_min = BioStat::get_minimum(\@data_sva);
 	my $data_max = BioStat::get_maximum(\@data_sva);
 	#print("$data_index\t$data[$data_index]\t$data_min\t$data_max\n");
	my $fract = GradientRGB::calculate_fraction($data_sva[$data_index],$data_min,$data_max);
	my $rgb_i = GradientRGB::fract_to_rgb_gradient_color_index($fract,11,\@rgb_high,\@rgb_low);
	$image->filledRectangle($ax1+$s, $ay1+$dist, $ax2+($s*1.1), $ay2+$dist, $colors[$rgb_i]);
	#if ($data_index == 60) {$image->filledRectangle($ax1+$s, $ay1+$dist, $ax2+($s*1.1), $ay2+$dist, $white)}
}
$image->rectangle($s, $ay1+$dist, $width+$s, $ay2+$dist, $black);
# ERVK share
$dist += 42;
foreach my $int (0..int((length($genome->{$chr}))/$res)) {
	my $data_index = $int;
	$int = $int*$res;
	my $ax1 = $int/$a_max*$width;
 	my $ax2 = $int/$a_max*$width;
 	# Calculate hits fraction
 	my $data_min = BioStat::get_minimum(\@data_erv);
 	my $data_max = BioStat::get_maximum(\@data_erv);
 	#print("$data_index\t$data[$data_index]\t$data_min\t$data_max\n");
	my $fract = GradientRGB::calculate_fraction($data_erv[$data_index],$data_min,$data_max);
	my $rgb_i = GradientRGB::fract_to_rgb_gradient_color_index($fract,11,\@rgb_high,\@rgb_low);
	$image->filledRectangle($ax1+$s, $ay1+$dist, $ax2+($s*1.1), $ay2+$dist, $colors[$rgb_i]);
	#if ($data_index == 60) {$image->filledRectangle($ax1+$s, $ay1+$dist, $ax2+($s*1.1), $ay2+$dist, $white)}
}
$image->rectangle($s, $ay1+$dist, $width+$s, $ay2+$dist, $black);
# piClusters
foreach my $pic (sort keys %{$pic_locs}) {
	# Get pic coordinates
	my $pic_chr = $pic_locs->{$pic}->[0];
	my $pic_beg = $pic_locs->{$pic}->[1];
	my $pic_end = $pic_locs->{$pic}->[2];
	unless ($pic_chr eq $chr) { next }
	# Get mid position of pic
	my $pic_len = $pic_end-$pic_beg+1;
	my $pic_loc = $pic_end-($pic_len/2);
	# Get position on image
	my $pic_pos = $pic_loc/$a_max*$width;
	$image->line($s+$pic_pos, $ay1-40, $s+$pic_pos, $ay2-40, $black);
}
# Scale bar
$image->dashedLine($s, $ay2+$dist+10, $width+$s, $ay2+$dist+10, $black);
$image->line($s, $ay2+$dist+4, $s, $ay2+$dist+16, $black);
$image->line($s+$width, $ay2+$dist+4, $s+$width, $ay2+$dist+16, $black);
for (my $i = 10; $i < int((length($genome->{$chr}))/$res); $i+=10) {
	my $int = $i*$res;
	my $ax1 = $int/$a_max*$width;
 	my $ax2 = $int/$a_max*$width;
	$image->line($s+$ax1, $ay2+$dist+6, $s+$ax1, $ay2+$dist+14, $black);
}
# Color code legend
my $leg_pos = 0;
for (my $fract = 1.0; $fract >= 0; $fract -= 0.1) {
	my $rgb_i = GradientRGB::fract_to_rgb_gradient_color_index($fract,11,\@rgb_high,\@rgb_low);
	$image->filledRectangle($width+2*$s-40, $ay1+$leg_pos, $width+2.5*$s-40, $ay2-28+$leg_pos, $colors[$rgb_i]);
	$leg_pos += 7.5;
}
$image->rectangle($width+2*$s-40, $ay1, $width+2.5*$s-40, $ay2-28+$leg_pos-7.5, $black);

# Create image
my $image_file = $chr.'.png';
open(my $im, '>', $image_file);
binmode($im);
print($im $image->png);
close($im);

print("ok\n");

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

sub get_rm_hits {
	# Take repeatmasker file name
	my($file) = @_;
	# Variables
	my %rmhits = ();
	# Open file
	my $in = FileIO::open_infile($file);
	# Extract rm hits and save with name
	while (my $line = <$in>) {
		$line =~ s/\s+$//; #better chomp
		if ($line =~ /^\s*\d+/) {
			# Get hit data
			$line =~ s/^\s*//;
			my @d = split(/\s+/,$line);
			my $rep_chr = $d[4];
			my $rep_cla = $d[10];
			if ($rep_cla =~ /Simple_repeat/ || $rep_cla =~ /Low_complexity/ || $rep_cla =~ /Satellite/) { next }
			# Save hit
			push(@{$rmhits{$rep_chr}},\@d);
		}
	}
	return \%rmhits;
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
			elsif ($type =~ /exon/i || $type =~ /cds/i) {
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
			my $chr = $locs->{$loc}->[0];
			my $beg = $locs->{$loc}->[1];
			my $end = $locs->{$loc}->[2];
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
						my $loc_chr = $locs->{$loc}->[0];
						my $loc_beg = $locs->{$loc}->[1];
						my $loc_end = $locs->{$loc}->[2];
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

sub get_slices_size {
	# Take fasta file name
	my($loc_data) = @_;
	# Var init
	my $total_size = 0;
	# Parse fasta file
	foreach my $loc (sort keys %{$loc_data}) {
		# Get loc coordinates
		my $chr = $loc_data->{$loc}->[0];
		my $beg = $loc_data->{$loc}->[1];
		my $end = $loc_data->{$loc}->[2];
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

sub get_centromer_locs {
	# Take centromer file name
	my($file) = @_;
	# Get file data
	my @file_data = FileIO::get_file_data_array($file);
	# Storage variable
	my %cen_locs = ();
	# Parse file data
	foreach my $line (@file_data) {
		# Split line
		my @d = split(/\t/,$line);
		my $chr = $d[0];
		my $beg = $d[1];
		my $end = $d[2];
		my $bnd = $d[4];
		# Skip if not centromer
		unless ($bnd eq 'acen') { next }
		# Save coordinates
		unless ($cen_locs{$chr}) {
			$cen_locs{$chr} = [$beg,$end];
		} else {
			if ($beg < $cen_locs{$chr}[0]) {
				$cen_locs{$chr}[0] = $beg;
			}
			if ($end > $cen_locs{$chr}[1]) {
				$cen_locs{$chr}[1] = $end;
			}
		}
	}
	return \%cen_locs;
}