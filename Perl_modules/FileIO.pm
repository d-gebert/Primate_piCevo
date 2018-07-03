################
#   File I/O   #
################

package FileIO;

use strict;
use warnings;
use Carp;

# Open input file
# Usage: my $in = open_infile($infile);
sub open_infile {
	# Take input file name
    my($file) = @_;
    # Open input file
    my $fh;
    if ($file =~ /.gz$/) {
		open($fh, "gunzip -c $file |") or die("Cannot open file '$file': $!\n");
	} else {
    	open($fh, '<', $file) or die("Cannot open file '$file': $!\n");
    }
    # Return filehandle
    return $fh;
}

# Open output file
# Usage: my $out = open_outfile($outfile);
sub open_outfile {
	# Take output file name
    my($file) = @_;
    # Open output file
    open(my $fh, '>', $file) or die("Cannot open file '$file': $!\n");
    # Return filehandle
    return $fh;
}

# Open output file to append to
# Usage: my $out = open_outfile_append($outfile);
sub open_outfile_append {
	# Take output file name
    my($file) = @_;
    # Open output file
    open(my $fh, '>>', $file) or die("Cannot open file '$file': $!\n");
    # Return filehandle
    return $fh;
}

# Find file name that does not already exist
# Usage: my $filename = find_unique_filename($filename);
sub find_unique_filename {
    # Take file name
    my($filename) = @_;
    # Check if file already exists
    if (-e $filename) {
    	# Filename contains a suffix
    	if ($filename =~ /.+\..+/) {
			my($prefix) = ($filename =~ /(.*)\.[^.]*$/);
			my($suffix) = ($filename =~ /.*(\.[^.]*$)/);
			my $orig_prefix = $prefix;
			# Find a non-existent name
			my $count = 0;
			while (-e $filename) {
				$count++;
				$prefix = $orig_prefix.'_'.$count;
				$filename = $prefix.$suffix;
			}
		}
		# Filename does not contain a suffix
		else {
			my $orig_name = $filename;
			# Find a non-existent name
			my $count = 0;
			while (-e $filename) {
				$count++;
				$filename = $orig_name;
				$filename = $orig_name.'_'.$count;
			}
		}
    }
    # Return unique filename
    return $filename;
}

# Extract file data and save in array
# Usage: my @filedata = get_file_data_array($file);
sub get_file_data_array {
	# Take input file name
    my($file,$ref_opt) = @_;
    my @filedata = ();
    $ref_opt = 0 unless $ref_opt;
	# Open input file
    my $fh = open_infile($file);
	# Extract lines and save in array
    while (my $line = <$fh>) {
    	$line =~ s/\s+$//; #better chomp
    	push(@filedata, $line);
    }
	# Close file
    close($fh) or die("Unable to close: $!\n");
	# Return array containing file data
    if ($ref_opt) {
    	return \@filedata;
    } else {
    	return @filedata;
    }
}

sub get_folder_file_names {
	# Take folder name
    my($folder) = @_;
	# Open folder
    opendir(my $dir, $folder) or die "Cannot open directory \"$folder\"\n\n";
    # Save files and filter out unix files '.', '..' and '.ds_store'
    my @files = grep (!/^\.\.?/, readdir($dir));
    @files = sort {$a cmp $b} @files;
    closedir($dir);
    # Return fiel names
    return @files;
}

1;