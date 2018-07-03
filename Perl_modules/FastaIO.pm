###############
#   FastaIO   #
###############

package FastaIO;

use strict;
use warnings;
use Carp;

# Save fasta data as hash
# Usage: my $sequences = get_fasta_seqs($infile);
sub get_fasta_seqs {
	# Take fasta file name
	my($file,$short) = @_;
	# Variables
	my $name = '';
	my %sequences = ();
	# Open file
	my $in = open_infile($file);
	# Extract sequence and save with name
	while (my $line = <$in>) {
		$line =~ s/\s+$//; #better chomp
		if ($line =~ /^>/) {
			($name) = ($line =~ />(.*)$/);
			if ($short) { $name =~ s/\s.*// }
		} else {
			$sequences{$name} .= $line;
		}
	}
	return \%sequences;
}

sub get_fasta_seqs_red {
	# Take fasta file name
	my($file,$short) = @_;
	# Variables
	my $name = '';
	my %sequences = ();
	my $id = 0;
	# Open file
	my $in = open_infile($file);
	# Extract sequence and save with name
	while (my $line = <$in>) {
		$line =~ s/\s+$//; #better chomp
		if ($line =~ /^>/) {
			$id++;
			($name) = ($line =~ />(.*)$/);
			if ($short) { $name =~ s/\s.*// }
			$sequences{$id}[0] = $name;
		} else {
			$sequences{$id}[1] .= $line;
		}
	}
	return \%sequences;
}

# Save fasta headers as array
# Usage: my $headers = get_fasta_headers($infile);
sub get_fasta_headers {
	# Take fasta file name
	my($file,$short) = @_;
	# Variables
	my $header = '';
	my @headers = ();
	# Open file
	my $in = open_infile($file);
	# Extract header and save in array
	while (my $line = <$in>) {
		$line =~ s/\s+$//; #better chomp
		if ($line =~ /^>/) {
			($header) = ($line =~ />(.*)$/);
			if ($short) { $header =~ s/\s.*// }
			push(@headers,$header);
		}
	}
	return \@headers;
}

# Print sequence with given line length
# Usage: print_sequence($sequence,$length);
sub print_sequence {
	# Take sequence and line length
    my($sequence, $length) = @_;
    # Print sequence in lines of $length
    for (my $pos = 0 ; $pos < length($sequence) ; $pos += $length) {
        print substr($sequence, $pos, $length), "\n";
    }
}

# Format sequence with given line length
# Usage: format_sequence($sequence,$length);
sub format_sequence {
	# Take sequence and line length
    my($sequence, $length) = @_;
    my $fmtd_seq = '';
    # Print sequence in lines of $length
    for (my $pos = 0 ; $pos < length($sequence) ; $pos += $length) {
        $fmtd_seq .= substr($sequence, $pos, $length)."\n";
    }
    # Return formatted sequence reference
    return \$fmtd_seq;
}

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

1;