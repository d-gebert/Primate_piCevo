#!/usr/bin/perl

use strict;
use warnings;
use lib '/Users/dgebert/Dropbox/Perlmodules';
use BIT;

# Program name
print("\n--- $0 ---\n");

# Collect command line arguments
my $USAGE = "perl $0 <spec1a_results.txt> <spec1b_results.txt>\n";
unless (@ARGV>=2) {
	die("\nUsage: $USAGE\n");
}
my $infile_1 = $ARGV[0];
my $infile_2 = $ARGV[1];

# Get piC locations
my $pic_info_1 = get_pic_info($infile_1);
my $pic_info_2 = get_pic_info($infile_2);

# Get identical piCs
my %ident_pics_1 = ();
my %ident_pics_2 = ();
# Compare piCs of both species
foreach my $pic_1 (sort {$a<=>$b} keys %{$pic_info_1}) {
	foreach my $pic_2 (sort {$a<=>$b} keys %{$pic_info_2}) {
		# Get chromosomes
		my $pic_1_chr = $pic_info_1->{$pic_1}->[0];
		my $pic_2_chr = $pic_info_2->{$pic_2}->[0];
		# If piCs on same chromosome
		if ($pic_1_chr eq $pic_2_chr) {
			# Get coordinates
			my $pic_1_beg = $pic_info_1->{$pic_1}->[1];
			my $pic_1_end = $pic_info_1->{$pic_1}->[2];
			my $pic_2_beg = $pic_info_2->{$pic_2}->[1];
			my $pic_2_end = $pic_info_2->{$pic_2}->[2];
			# If coordinates overlap
			if ($pic_1_beg < $pic_2_end && $pic_1_end > $pic_2_beg) {
				# Save identical piCs
				$ident_pics_1{$pic_1} = $pic_2;
				$ident_pics_2{$pic_2} = $pic_1;
			}
		}
	}
}

# Get unique pics
my @uniqe_pics_1;
foreach my $pic_1 (sort {$a<=>$b} keys %{$pic_info_1}) {
	push(@uniqe_pics_1,$pic_1) unless $ident_pics_1{$pic_1};
}
my @uniqe_pics_2;
foreach my $pic_2 (sort {$a<=>$b} keys %{$pic_info_2}) {
	push(@uniqe_pics_2,$pic_2) unless $ident_pics_2{$pic_2};
}

# Output results
my $out = BIT::open_outfile('Intrespec_piCs.txt');
print($out "Hsap1\tHsap2\tsize1\thits1\tdens1\tsize2\thits2\tdens2\n");
foreach my $pic_1 (sort {$a<=>$b} keys %ident_pics_1) {
	my $pic_2 = $ident_pics_1{$pic_1};
	my $size_1 = $pic_info_1->{$pic_1}->[3];
	my $hits_1 = $pic_info_1->{$pic_1}->[4];
	my $dens_1 = $pic_info_1->{$pic_1}->[5];
	my $size_2 = $pic_info_2->{$pic_2}->[3];
	my $hits_2 = $pic_info_2->{$pic_2}->[4];
	my $dens_2 = $pic_info_2->{$pic_2}->[5];
	printf($out "%s\t%s\t%d\t%0.1f\t%0.1f\t%d\t%0.1f\t%0.1f\n", $pic_1,$pic_2,$size_1,$hits_1,$dens_1,$size_2,$hits_2,$dens_2);
}
foreach my $pic_1 (@uniqe_pics_1) {
	my $size = $pic_info_1->{$pic_1}->[3];
	my $hits = $pic_info_1->{$pic_1}->[4];
	my $dens = $pic_info_1->{$pic_1}->[5];
	printf($out "%s\t\t%d\t%0.1f\t%0.1f\n", $pic_1,$size,$hits,$dens);
}
foreach my $pic_2 (@uniqe_pics_2) {
	my $size = $pic_info_2->{$pic_2}->[3];
	my $hits = $pic_info_2->{$pic_2}->[4];
	my $dens = $pic_info_2->{$pic_2}->[5];
	printf($out "\t%s\t\t\t\t%d\t%0.1f\t%0.1f\n", $pic_2,$size,$hits,$dens);
}

exit;

sub get_pic_info {
	# Take infile name
	my($file) = @_;
	# Get infile data as array
	my @in_data = BIT::get_file_data_array($file);
	# Initialize variables
	my %pic_info = ();
	# Parse piC data
	foreach my $line (@in_data) {
		# Result line        
		if ($line =~ /^Cluster/) {
	    	# Get pic info 
	    	my @d = split('\t', $line);
	    	my($pic)  = ($d[0] =~ /Cluster\s(\d+)/);
	    	my($chr)  = ($d[1] =~ /Location:\s(\w+)/);
	    	my($beg)  = ($d[2] =~ /Coordinates:\s(.+)\-/);
	    	my($end)  = ($d[2] =~ /Coordinates:\s.+\-(.+)/);
	    	my($size) = ($d[3] =~ /.+\s(\d+)/);
	    	my($hits) = ($d[5] =~ /.+\s(\d+\.*\d*)/);
			my($dens) = ($d[6] =~ /.+\s(\d+\.*\d*)/);
	    	@{$pic_info{$pic}} = ($chr,$beg,$end,$size,$hits,$dens);
	    }
	}
	return \%pic_info;
}