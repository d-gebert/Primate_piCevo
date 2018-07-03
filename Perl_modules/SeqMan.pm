##############
#   SeqMan   #
##############

package SeqMan;

use strict;
use warnings;
use Carp;

sub rev_com {
	my($seq) = @_;
	
	my $revcom = reverse $seq;
	$revcom =~ tr/ACGTUNacgtun/TGCAANtgcaan/;
	$revcom =~ s/\s//g;
	
	return $revcom;
}

sub rev {
    my($seq) = @_;
    my $revseq = reverse $seq;
    $revseq =~ s/\s//g;
    return $revseq;
}

sub com {
    my($seq) = @_;
    $seq =~ tr/ACGTUNacgtun/TGCAANtgcaan/;
    return $seq;
    $seq =~ s/\s//g;
}

1;