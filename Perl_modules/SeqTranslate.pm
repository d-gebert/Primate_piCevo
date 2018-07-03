##################
#  SeqTranslate  #
##################

package SeqTranslate;

use strict;
use warnings;
use Carp;

sub dna2peptide {

	my($dna,$frame) = @_;

	# Initialize variables
	my $protein = '';
	my $codon;
	if (not $frame) { $frame = 1 }
	if ($frame && $frame>3) { $frame = 1 }
	# Translate each three-base codon into an amino acid, and append to a protein 
	for(my $i=($frame-1); $i < (length($dna) - 2) ; $i += 3) {
   		$codon = substr($dna,$i,3);
    	$protein .= codon2aa($codon);
	}
	return $protein;
}

sub dna2peptide_spaces {

	my($dna,$frame) = @_;

	# Initialize variables
	my $protein = '';
	my $codon;
	if (not $frame) { $frame = 1 }
	if ($frame && $frame>3) { $frame = 1 }
	# Translate each three-base codon into an amino acid, and append to a protein 
	for(my $i=($frame-1); $i < (length($dna) - 2) ; $i += 3) {
   		$codon = substr($dna,$i,3);
    	$protein .= ' '.codon2aa($codon).' ';
	}
	return $protein;
}

sub codon2aa {

    my($codon) = @_;

    $codon = uc $codon;
 
    my(%genetic_code) = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '*',    # Stop
    'TAG' => '*',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '*',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );

    if(exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    }else{
        return '?';
    }
}

sub stop_positions {

    my($regexp, $sequence) = @_;

    # Declare variables
    my @positions = (  );

    # Determine positions of regular expression matches   
    while ( $sequence =~ /$regexp/ig ) {

        push ( @positions, ((pos($sequence) - length($&) + 1)-1));
    }

    return @positions;
}

sub print_aa_dna {

	my($dna, $peptide, $length) = @_;
	
	# print dna sequence with translation
	for ( my $pos = 0 ; $pos < length($dna) ; $pos += $length ) {
		print "P ", $pos/3+1, "	", substr($peptide, $pos, $length), "\n";
		print "N ", $pos+1, "	", substr($dna, $pos, $length), "\n\n";
	}	
}

sub print_aa_dna_label {

	my($dna, $peptide, $length) = @_;
	my $position;
	
	# whitespace until deletion label '^'
	my $label = ' ' x ($position - 1).'^';
	
	# define the remaining whitespace to the end of the sequence
	my $label_rest = ' ' x (length($dna) - length($label));

	# put the parts together
	$label = $label.$label_rest;
	
	# print dna sequence with translation
	for ( my $pos = 0 ; $pos < length($dna) ; $pos += $length ) {
		print "P ", $pos/3+1, "	", substr($peptide, $pos, $length), "\n";
		print "N ", $pos+1, "	", substr($dna, $pos, $length), "\n";
		print "	 ", substr($label, $pos, $length), "\n";
	}	
}

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