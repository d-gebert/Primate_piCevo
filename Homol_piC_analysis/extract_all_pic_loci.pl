#!/usr/bin/perl

use strict;
use warnings;
use lib '/Users/dgebert/Dropbox/Perlmodules';
use BIT;
use BioStat;

# Global constants
my $outfile = "piC_loci_list.txt";
# Global variables
my %pic_loc_props = ();
# Options
$|=1; #Autoflush

# Species identifiers
my @species = ('Hsap','Mmul','Mfas','Cjac','Mmur','Ltar');
my %spec_id = ();
foreach my $i (0..$#species) {
	$spec_id{$species[$i]} = $i;
}

# Program name
print("\n--- $0 ---\n");

# Collect command line arguments
my $USAGE = "perl $0 <files_list_p1piCs.txt> <files_list_p5piCs.txt> <files_list_results_p1.txt> <files_list_results_p5.txt>\n";
unless ($ARGV[0]&&$ARGV[1]&&$ARGV[2]&&$ARGV[3]) {
	die("\nUsage: $USAGE\n");
}
my $p1_file_list = $ARGV[0];
my $p5_file_list = $ARGV[1];
my $p1_tabl_list = $ARGV[2];
my $p5_tabl_list = $ARGV[3];

## Get loci of all piCs
my $pic_loci_p1 = collect_pic_loci($p1_file_list);
my $pic_loci_p5 = collect_pic_loci($p5_file_list);

## Get properties of all piCs
my $pic_props_p1 = get_pic_properties_multi($p1_tabl_list);
my $pic_props_p5 = get_pic_properties_multi($p5_tabl_list);

## Link properties to loci (stored in %pic_loc_props)
add_props_to_pic_loci($pic_loci_p1,$pic_props_p1);
add_props_to_pic_loci($pic_loci_p5,$pic_props_p5);

## Output to file
$outfile = BIT::find_unique_filename($outfile);
my $out = BIT::open_outfile($outfile);
# Go through each locus and output pic id and properties
foreach my $pic (sort {$a cmp $b} keys %pic_loc_props) {
	my $loc  = $pic_loc_props{$pic}->[0];
	my $size = $pic_loc_props{$pic}->[1];
	my $hits = $pic_loc_props{$pic}->[2];
	my $dens = $pic_loc_props{$pic}->[3];
	my $dir  = $pic_loc_props{$pic}->[4];
	printf($out "%s\t%s\t%d\t%.2f\t%.2f\t%s\n", $pic, $loc, $size, $hits, $dens, $dir);
}

exit;

################################# subroutines #################################

sub collect_pic_loci {
	# Take file name array ref
	my($file_list) = @_;
	# Initialize subglobal variables
	my %pic_loci = ();
	# Get results files
	my @files = BIT::get_file_data_array($file_list);
	# Go through each file
	foreach my $file (@files) {
		# Get file data
		my @file_data = BIT::get_file_data_array($file);
		# Parse file data
		foreach my $line (@file_data) {
			# Fasta header
			if ($line =~ /^>/) {
				# Get header data
				my @d = split(/\s/,$line);
				my($pic) = ($d[0] =~ />(.*)/);
				my $loc = $d[1].' '.$d[2];
				$pic_loci{$pic} = $loc;
			}
		}
	}
	return \%pic_loci;
}

sub get_pic_properties_multi {
	# Take proTRAC results file list
	my($file_list) = @_;
	# Initialize subglobal variables
	my @properties = ();
	my $sp_i = 0;
	# Get results files
	my @files = BIT::get_file_data_array($file_list);
	# Go though each file
	foreach my $file (@files) {
		# Get piC properties
		my $pic_props = get_pic_properties($file);
		# Include piC properties in multi hash
		foreach my $pic (keys %{$pic_props}) {
			@{$properties[$sp_i]{$pic}} = @{$pic_props->{$pic}};
		}
		# Increment species id number
		$sp_i++;
	}
	return \@properties;
}

sub get_pic_properties {
	# Take proTRAC results file name
	my($file) = @_;
	my %pic_props = ();
	# Get file data
	my @in_data_res = BIT::get_file_data_array($file);
	# Parse proTRAC output
	foreach my $line (@in_data_res) {
		# Cluster info line
		if ($line =~ /^Cluster/) {
			# Get cluster info
			my @d = split(/\t/, $line);
			my($pic)  = ($d[0] =~ /.+\s(\d+)/);
			my($size) = ($d[3] =~ /.+\s(\d+)/);
			my($hits) = ($d[5] =~ /.+\s(\d+\.*\d*)/);
			my($dens) = ($d[6] =~ /.+\s(\d+\.*\d*)/);
			my($dir)  = ($d[11] =~ /.+\s(\w+:.*)/);
			my @bsites = ();
			if ($d[12]) {
				my($bsites) = ($d[12] =~ /sites:\s(.*)/);
				$bsites .= ' ';
				@bsites = split(/\s\(\w+:\s\d+\)\s/, $bsites);
			}
			# Store data in hashe with pic number as key
			@{$pic_props{$pic}} = ($size,$hits,$dens,$dir,\@bsites);
		}
	}
	return \%pic_props;
}

sub add_props_to_pic_loci {
	# Take pic loci and pic properties
	my($pic_loci,$pic_props) = @_;
	# Go through each locus and get (combined) piC properties
	foreach my $pic (keys %{$pic_loci}) {
		# Get species id
		my($spec) = ($pic =~ /^([^-_]+)/);
		my $sp_i = $spec_id{$spec};
		# Get piC ids
		my @pic_ids = split(/\_/,$pic);
		print("$pic\t$spec\t$sp_i\n");
		# Get piC properties for each piC id
		my $hits = 0;
		my $size = 0;
		my $dir = '';
		foreach my $i (1..$#pic_ids) {
			$size += $pic_props->[$sp_i]->{$pic_ids[$i]}->[0];
			$hits += $pic_props->[$sp_i]->{$pic_ids[$i]}->[1];
			$dir = $pic_props->[$sp_i]->{$pic_ids[$i]}->[3];
		}
		my $dens = $hits/$size*1000;
		# Collect properties and link to locus
		@{$pic_loc_props{$pic}} = ($pic_loci->{$pic},$size,$hits,$dens,$dir);
	}
}