#####################
#   BioStatistics   #
#####################

package BioStat;

use strict;
use warnings;

# Calculate percentile value of array values
sub get_percentile_value {
	# Take array list and p-value
	my($array, $p) = @_;
	# Sort array list
	my @ordered_list = sort { $a <=> $b } @{$array};
	# Number of values
	my $N = scalar(@ordered_list);
	# Calculate ordinal rank
	my $n = (1-$p) * $N;
	# Get position in ordered list
	my $i = int($n+0.999999999);
	# Get percentile value
	my $perc_val = 0;
	if ($N % 2) { #odd N
		$perc_val = $ordered_list[$i-1]; 
	} else { #even N
		$perc_val = ($ordered_list[$i-1]+$ordered_list[$i])/2;
	}
	# Return percentile value
	return $perc_val;
}

# Calculate population variance of array values
sub get_population_variance {
	# Take array list
	my($array) = @_;
	# Number of values
	my $N = scalar(@{$array});
	# Get mean value
	my $mean = get_mean($array);
	# Sum squared differences
	my $diff_sum = 0;
	grep { $diff_sum += ($_-$mean)**2 } @{$array};
	# Calculate population variance
	my $variance = $diff_sum/$N;
	# Return population variance
	return $variance;
}

# Calculate sample variance of array values
sub get_sample_variance {
	# Take array list
	my($array) = @_;
	# Number of values
	my $N = scalar(@{$array});
	# Get mean value
	my $mean = get_mean($array);
	# Sum squared differences
	my $diff_sum = 0;
	grep { $diff_sum += ($_-$mean)**2 } @{$array};
	# Calculate sample variance
	my $variance = $diff_sum/($N-1);
	# Return sample variance
	return $variance;
}

# Calculate population standard deviation of array values
sub get_population_standard_deviation {
	# Take array list
	my($array) = @_;
	# Get sample variation
	my $variance = get_population_variance($array);
	# Calculate standard deviation from variance
	my $std_dev = sqrt($variance);
	# Return standard deviation
	return $std_dev;
}

# Calculate sample standard deviation of array values
sub get_sample_standard_deviation {
	# Take array list
	my($array) = @_;
	# Get sample variation
	my $variance = get_sample_variance($array);
	# Calculate standard deviation from variance
	my $std_dev = sqrt($variance);
	# Return standard deviation
	return $std_dev;
}

# Calculate sum of array values
sub get_sum {
	# Take array list
	my($array) = @_;
	# Sum values
	my $sum = 0;
	grep { $sum += $_ } @{$array};
	# Return sum
	return $sum;
}

# Calculate mean of array values
sub get_mean {
	# Take array list
	my($array) = @_;
	# Number of values
	my $N = scalar(@{$array});
	# Sum values
	my $sum = get_sum($array);
	# Calculate mean value
	my $mean = $sum/$N;
	# Return mean
	return $mean;
}

# Calculate meadian of array values
sub get_median {
	# Take array list
	my($array) = @_;
	# Get median, 50th percentile
	my $median = get_percentile_value($array, 0.5);
	$median = get_mean($array) if scalar(@{$array})<3;
	# Return median
	return $median;
}

# Calculate first quartile of array values
sub get_first_quartile {
	# Take array list
	my($array) = @_;
	# Get median, 50th percentile
	my $Q1 = get_percentile_value($array, 0.75);
	# Return median
	return $Q1;
}

# Calculate third quartile of array values
sub get_third_quartile {
	# Take array list
	my($array) = @_;
	# Get median, 50th percentile
	my $Q3 = get_percentile_value($array, 0.25);
	# Return median
	return $Q3;
}

# Get minimum of array values
sub get_minimum {
	# Take array list
	my($array) = @_;
	# Sort array list
	my @ordered_list = sort { $a <=> $b } @{$array};
	# Get minimum value
	my $min = $ordered_list[0];
	# Return minimum value
	return $min;
}

# Get maximum of array values
sub get_maximum {
	# Take array list
	my($array) = @_;
	# Sort array list
	my @ordered_list = sort { $a <=> $b } @{$array};
	# Get maximum value
	my $max = $ordered_list[-1];
	# Return maximum value
	return $max;
}

# Calculate z-score for given value as part of given array
sub get_z_score {
	# Take value and array list
	my($val,$array) = @_;
	# Get mean
	my $mean = get_mean($array);
	# Get standard deviation
	my $std_dev = get_population_standard_deviation($array);
	# Calculate z-score for value
	my $zscore = ($val-$mean)/$std_dev;
	# Return z-score for input value
	return $zscore;
} 

sub pearson_correlation {
	# Take two arrays of same length
	my($array_x,$array_y) = @_;
	# Initialize sum variables
	my $x_sum = 0;
	my $y_sum = 0;
	my $xy_sum = 0;
	my $xx_sum = 0;
	my $yy_sum = 0;
	# Get n
	my $n = scalar(@{$array_x});
	# Get sums
	foreach my $i (0..($n-1)) {
		$x_sum += $array_x->[$i];
		$y_sum += $array_y->[$i];
		$xy_sum += ($array_x->[$i]*$array_y->[$i]);
		$xx_sum += ($array_x->[$i]*$array_x->[$i]);
		$yy_sum += ($array_y->[$i]*$array_y->[$i]);
	}
	unless ($x_sum&&$y_sum) { return 0 }
	# Calculate r
	my $r = ( ($n*$xy_sum) - ($x_sum*$y_sum) ) / sqrt( ( ($n*$xx_sum) - ($x_sum*$x_sum) ) * ( ($n*$yy_sum) - ($y_sum*$y_sum) ) );
	return $r;
}

# Truncate decimal number to desired decimal places
sub truncate_decimal_number {
	# Take value and number of decimal places
	my($val, $n_dec) = @_;
	# Truncate to desired decimal places
	my $trunc_val = int($val*(10**$n_dec)+0.5)/(10**$n_dec);
	# Return truncated value
	return $trunc_val;
}

# Add commas as thounds separators
sub add_thousands_separators {
	# Take number
	my($number) = @_;
	# Add commas as thousands separators to number
	$number = reverse join ",", (reverse $number) =~ /(\d{1,3})/g;
	# Return number
	return $number;
}

# Extract unique array elements
sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

1;