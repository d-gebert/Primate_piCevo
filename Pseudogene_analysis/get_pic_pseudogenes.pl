#!/usr/bin/perl
use strict;
use warnings;
use lib '/Users/dgebert/Dropbox/Perlmodules/';
use FileIO;
use FastaIO;

# Constants
$|=1; #Autoflush
# Variables

# Collect command line arguments
my $USAGE = "perl $0 <piCseqs.fas> <genome.gff>\n";
unless ($ARGV[0] && $ARGV[1]) {
	die("\nUsage: $USAGE\n");
}
# Input files
my $pic_seqs_file = $ARGV[0];
my $gene_gff_file = $ARGV[1];

# Get piC locations
my $pic_locs = get_pic_locs($pic_seqs_file);
# Get pic ids
my $pic_ids = FastaIO::get_fasta_headers($pic_seqs_file,1);
# Get gff gene and cds locations
my($pgene_locs,$pexon_locs) = get_gff_gene_locs($gene_gff_file);

# Output pic pseudogenes
my($genome_name) = ($gene_gff_file =~ /(\w+)\./);
my $outfile = $genome_name.'.piCpseudogenes.txt';
my $out = FileIO::open_outfile($outfile);
# Go through each pic
foreach my $pic (@{$pic_ids}) {
	# Get pic coordinates
	my $chr = $pic_locs->{$pic}->[0];
	my $pic_beg = $pic_locs->{$pic}->[1];
	my $pic_end = $pic_locs->{$pic}->[2];
	my $pic_dir = $pic_locs->{$pic}->[3];
	# Get pic expression direction per position
	my %pic_pos = ();
	if ($pic_dir =~ /mono:/) {
		foreach my $pos ($pic_beg..$pic_end) {
			$pic_pos{$pos} = '+' if $pic_dir =~ /plus/;
			$pic_pos{$pos} = '-' if $pic_dir =~ /minus/;
		}
	} elsif ($pic_dir =~ /bi:/) {
		my($split_a) = ($pic_dir =~ /between (\d+) and \d+/);
		my($split_b) = ($pic_dir =~ /between \d+ and (\d+)/);
		foreach my $pos ($pic_beg..$split_a) {
			$pic_pos{$pos} = '+' if $pic_dir =~ /plus-/;
			$pic_pos{$pos} = '-' if $pic_dir =~ /minus-/;
		}
		foreach my $pos ($split_b..$pic_end) {
			$pic_pos{$pos} = '-' if $pic_dir =~ /-minus/;
			$pic_pos{$pos} = '+' if $pic_dir =~ /-plus/;
		}
	}
	# Bool pic printed: false
	my $pic_printed = 0;
	# Go through each gene on same chromosome
	foreach my $gene (sort { $pgene_locs->{$chr}->{$a}->[3] <=> $pgene_locs->{$chr}->{$b}->[3]} keys %{$pgene_locs->{$chr}}) {
		# Get pseudogene properties
		my $gen_beg = $pgene_locs->{$chr}->{$gene}->[3];
		my $gen_end = $pgene_locs->{$chr}->{$gene}->[4];
		my $gen_str = $pgene_locs->{$chr}->{$gene}->[6] eq '+' ? 'plus' : 'minus';
		my $gen_inf = $pgene_locs->{$chr}->{$gene}->[8];
		my $gen_ide = $pgene_locs->{$chr}->{$gene}->[9];
		my $gen_gsy = $pgene_locs->{$chr}->{$gene}->[-1];
		my $gen_len = $gen_end-$gen_beg+1;
		my($gen_typ) = ($gen_inf =~ /gene_biotype=([^;,]+)/);
		# Pseudogene lies within pic
		if ($gen_beg <= $pic_end && $gen_end >= $pic_beg) {
			my %exon_poss = ();
			# Go through each exon
			foreach my $exon (@{$pexon_locs->{$chr}->{$gene}}) {
				my $ex_beg = $exon->[3];
				my $ex_end = $exon->[4];
				my $ex_ide = $exon->[9];
				my $ex_len = $ex_end-$ex_beg+1;
				# Check if exon region already printed
				my $region_covered = 0;
				foreach my $pos ($ex_beg..$ex_end) {
					if ($exon_poss{$pos}) { $region_covered++ }
					$exon_poss{$pos} = 1;
				}
				if ($region_covered > $ex_len/2) { next }
				# Exon lies within pic
				if ($ex_beg <= $pic_end && $ex_end >= $pic_beg) {
					# Get exon orientation within pic
					my $gen_ori = '';
					my($par_pos,$rev_pos) = (0,0);
					foreach my $pos ($gen_beg..$gen_end) {
						unless ($pic_pos{$pos}) { next }
						my $pos_strand = '';
						if ($gen_str eq 'plus') {$pos_strand = '+'}
						if ($gen_str eq 'minus') {$pos_strand = '-'}
						if ($pic_pos{$pos} eq $pos_strand) { $par_pos++ }
						if ($pic_pos{$pos} ne $pos_strand) { $rev_pos++ }
					}
					if ($par_pos > $rev_pos) { $gen_ori = 'par' }
					if ($rev_pos > $par_pos) { $gen_ori = 'rev' }
					# Print pic info if not already printed
					unless ($pic_printed) {
						print($out "#$pic\t$chr\t$pic_beg-$pic_end\t$pic_dir\n");
						$pic_printed = 1;
					}
					# Print exon info
					printf($out "%s\t%s\t%d\t%d\t%.2f\t%.2f\t%s\t%d\t%d\t%s\t%s\t%s\n", 
						$gene,$gen_gsy,$gen_len,$ex_len,$gen_ide,$ex_ide,$chr,$ex_beg,$ex_end,$gen_str,$gen_typ,$gen_ori);
				}
			}
		}
	}
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
		my @d = split(/\t/,$info);
		my $pic = $d[0];
		my $chr = $d[1];
		my($beg,$end) = split(/-/,$d[2]);
		my $dir = $d[3];
		# Save pic loc info
		@{$pic_locs{$pic}} = ($chr,$beg,$end,$dir);
	}
	return \%pic_locs;
}

sub get_gff_gene_locs {
	# Take infile name and keyword
	my($gff_file) = @_;
	# Get file data
	my $gff = FileIO::open_infile($gff_file);
	# Initialize variables
	my %ggene_locs = ();
	my %pexon_locs = ();
	my $gene = '';
	my $name = '';
	# Parse gff file
	while (my $line = <$gff>) {
		$line =~ s/\s+$//; #better chomp
		# Results line
		if ($line =~ /^[^#]\w+/) {
			# Get info
			my @d = split(/\t/, $line);
			my $chr  = $d[0];
			my $type = $d[2];
			my $beg  = $d[3];
			my $end  = $d[4];
			my $info = $d[8];
			my $iden = $d[9] ? $d[9] : 0;
			$d[9] = $iden;
			# Gene line
			if ($type =~ /gene/i) {
				# Skip non pseudogenic
				if ($info !~ /pseudogene/i) { next }
				# Get gene id
				# GFF
				if ($info =~ /GeneID:([^;,]+)/) {
					($gene) = ($info =~ /GeneID:([^;,]+)/);
				}
				# GTF
				if ($info =~ /gene_id\s\"([^\"]+)\";/) {
					($gene) = ($info =~ /gene_id\s\"([^\"]+)\";/);
				}
				# Get gene name 
				# GFF
				if ($info =~ /Name=([^;,]+)/) {
					($name) = ($info =~ /Name=([^;,]+)/);
				}
				# GTF
				if ($info =~ /gene_name\s"([^\"]+)"/) {
					($name) = ($info =~ /gene_name\s"([^\"]+)"/);
				}
				# If pseudogenes are predicted by Pseudon, add unit id
				if ($d[1] eq 'Pseudon') {
					my($id) = ($info =~ /ID=(\d+);/);
					$gene = $gene.'_'.$id;
				}
				# Save entries
				push(@d,$name);
				$ggene_locs{$chr}{$gene} = \@d;
			}
			# Exon line
			elsif ($type =~ /exon/i) {
				# Save entries
				push(@d,$name);
				push(@{$pexon_locs{$chr}{$gene}},\@d);
			}
		}
	}
	return \(%ggene_locs,%pexon_locs);
}