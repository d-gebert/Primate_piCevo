#!/usr/bin/perl

use strict;
use warnings;
use lib '/Users/dgebert/Dropbox/Perlmodules/';
use FileIO;
use FastaIO;

# Constants
$|=1; #Autoflush
# Variables
my $univ_predicted = 0;
my %pseudos_x = ();
my %pseudos_y = ();
my %x_unique_pgs = ();
my %y_unique_pgs = ();

# Collect command line arguments
my $USAGE = "perl $0 <piCpseudogenes_1.txt> <piCpseudogenes_2.txt>\n";
unless ($ARGV[0] && $ARGV[1]) {
	die("\nUsage: $USAGE\n");
}
# Input files
my $pic_pseudos_file_x = $ARGV[0];
my $pic_pseudos_file_y = $ARGV[1];

# Get pic pseudogenes
my($pic_info_x,$pic_pgenes_x) = get_pic_pseudogenes($pic_pseudos_file_x);
my($pic_info_y,$pic_pgenes_y) = get_pic_pseudogenes($pic_pseudos_file_y);

# Go through each pseudogene-containing pic
foreach my $pic (sort keys %{$pic_pgenes_x}) {
	# Go through each pseudogene
	foreach my $pg_x (sort keys %{$pic_pgenes_x->{$pic}}) {
		# Look for overlapping pseudogene prediction
		my $pg_overlap = 0;
		# Go through each exon
		foreach my $ex_x (@{$pic_pgenes_x->{$pic}->{$pg_x}}) {
			# Get exon coordinates
			my $ex_x_beg = $ex_x->[7];
			my $ex_x_end = $ex_x->[8];
			# Go through each pseudogene of other prediction
			foreach my $pg_y (sort keys %{$pic_pgenes_y->{$pic}}) {
				# Go through each exon
				foreach my $ex_y (@{$pic_pgenes_y->{$pic}->{$pg_y}}) {
					# Get exon coordinates
					my $ex_y_beg = $ex_y->[7];
					my $ex_y_end = $ex_y->[8];
					# Check if pseudogenes overlap
					if ($ex_y_end >= $ex_x_beg && $ex_y_beg <= $ex_x_end) {
						$pg_overlap++;
					}
				}
			}
		}
		# Count pseudogenes with universal prediction
		if ($pg_overlap) {
			$univ_predicted++;
		} else {
			$x_unique_pgs{$pic}{$pg_x} = 1;
		}
	}
}

# Go through each pseudogene-containing pic
foreach my $pic (sort keys %{$pic_pgenes_y}) {
	# Go through each pseudogene
	foreach my $pg_y (sort keys %{$pic_pgenes_y->{$pic}}) {
		# Look for overlapping pseudogene prediction
		my $pg_overlap = 0;
		# Go through each exon
		foreach my $ex_y (@{$pic_pgenes_y->{$pic}->{$pg_y}}) {
			# Get exon coordinates
			my $ex_y_beg = $ex_y->[7];
			my $ex_y_end = $ex_y->[8];
			# Go through each pseudogene of other prediction
			foreach my $pg_x (sort keys %{$pic_pgenes_x->{$pic}}) {
				# Go through each exon
				foreach my $ex_x (@{$pic_pgenes_x->{$pic}->{$pg_x}}) {
					# Get exon coordinates
					my $ex_x_beg = $ex_x->[7];
					my $ex_x_end = $ex_x->[8];
					# Check if pseudogenes overlap
					if ($ex_x_end >= $ex_y_beg && $ex_x_beg <= $ex_y_end) {
						$pg_overlap++;
					}
				}
			}
		}
		# Count pseudogenes with universal prediction
		if ($pg_overlap) {
			#$univ_predicted++;
		} else {
			$y_unique_pgs{$pic}{$pg_y} = 1;
		}
	}
}

foreach my $pic (sort keys %{$pic_pgenes_x}) {
	# Go through each pseudogene
	foreach my $pg_x (sort keys %{$pic_pgenes_x->{$pic}}) {
		# Count x pseudogenes
		$pseudos_x{$pg_x} = 1;
	}
}

foreach my $pic (sort keys %{$pic_pgenes_y}) {
	# Go through each pseudogene
	foreach my $pg_y (sort keys %{$pic_pgenes_y->{$pic}}) {
		# Count y pseudogenes
		$pseudos_y{$pg_y} = 1;
	}
}

my $x_predicted = keys %pseudos_x;
my $y_predicted = keys %pseudos_y;

print("$x_predicted $univ_predicted $y_predicted\n");

my $outfilex = 'x_unique_pseudos.txt';
my $outx = FileIO::open_outfile($outfilex);
foreach my $pic (sort keys %x_unique_pgs) {
	print($outx "#$pic\n");
	foreach my $pg (sort keys %{$x_unique_pgs{$pic}}) {
		foreach my $ex (@{$pic_pgenes_x->{$pic}->{$pg}}) {
			foreach my $ex_prop (@{$ex}) {
				print($outx "$ex_prop\t");
			}
			print($outx "\n");
		}
	}
}
my $outfiley = 'y_unique_pseudos.txt';
my $outy = FileIO::open_outfile($outfiley);
foreach my $pic (sort keys %y_unique_pgs) {
	print($outy "#$pic\n");
	foreach my $pg (sort keys %{$y_unique_pgs{$pic}}) {
		foreach my $ex (@{$pic_pgenes_y->{$pic}->{$pg}}) {
			foreach my $ex_prop (@{$ex}) {
				print($outy "$ex_prop\t");
			}
			print($outy "\n");
		}
	}
}

exit;

################################# subroutines #################################

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
		# Get pic info
		if ($line =~ /^#/) {
			my @d = split(/\t/,$line);
			$pic_id = $d[0];
			$pic_id =~ s/#//;
			my $chr = $d[1];
			my($beg,$end) = split(/-/,$d[2]);
			my $dir = $d[3];
			$pic_info{$pic_id} = [$chr,$beg,$end,$dir];
		} 
		# Get pseudogene
		elsif ($line !~ /^$/) {
			my @d = split(/\t/,$line);
			my $parent_id = $d[0];
			my $parent_gs = $d[1];
			my $pgene_ide = $d[4];
			my $pgene_beg = $d[7];
			my $pgene_end = $d[8];
			my $pgene_str = $d[9];
			my $pgene_dir = $d[11];
			push(@{$pic_pgenes{$pic_id}{$parent_id}},\@d);
		}
	}
	return \(%pic_info,%pic_pgenes);
}