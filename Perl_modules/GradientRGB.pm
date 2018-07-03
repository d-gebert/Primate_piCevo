################
#   Gradient   #
################

package GradientRGB;

use strict;
use warnings;

sub calculate_fraction {
	# Take absolute value, global minimum and maximum
	my($val,$val_min,$val_max) = @_;
	# Calculate fraction from absolute value
	my $fraction = ($val-$val_min) / ($val_max-$val_min);
	# Return fraction value
	return $fraction;
}

sub fract_to_rgb_gradient_color {
	# Take fraction value
	my($val,$steps,$rgb1,$rgb2) = @_;
	# Round value to one decimal place
	my $r_val = int($val*10+0.5)/10;
	# Get rgb gradient
	my @gradient = make_rgb_gradient($steps,$rgb1,$rgb2);
	# Make gradient hash
	my %gradient = ();
	my $i = 0;
	for (my $x = 10; $x >= 0; $x -= 1) {
		my $fract = $x/10;
		$gradient{$fract} = $gradient[$i];
		$i++;
	}
	# Return rgb color according to fraction value
	return @{$gradient{$r_val}};
}

sub fract_to_rgb_gradient_color_index {
	# Take fraction value
	my($val,$steps,$rgb1,$rgb2) = @_;
	# Round value to one decimal place
	my $r_val = int($val*10+0.5)/10;
	# Get rgb gradient
	my @gradient = make_rgb_gradient($steps,$rgb1,$rgb2);
	# Make gradient hash
	my %gradient = ();
	my %index = ();
	my $i = 0;
	for (my $x = 10; $x >= 0; $x -= 1) {
		my $fract = $x/10;
		$gradient{$fract} = $gradient[$i];
		$index{$fract} = $i;
		$i++;
	}
	# Return rgb color according to fraction value
	return $index{$r_val};
}

sub make_rgb_gradient {
	# Take number of steps and max/min rgbs
	my($steps,$rgb1,$rgb2) = @_;
	# Get initial r, g and b
	my($r,$g,$b) = @$rgb1;
	# Initialize gradient with first rgb
	my @grad = ([$r,$g,$b]);
	# Increment rgb for gradient for each step
	foreach my $step (1..$steps) {
		# Calculate new rgb coordinates
		$r += (($rgb2->[0]-$rgb1->[0])/$steps);
		$g += (($rgb2->[1]-$rgb1->[1])/$steps);
		$b += (($rgb2->[2]-$rgb1->[2])/$steps);
		# Add rgb color to gradient array
		push(@grad,[int($r+0.5),int($g+0.5),int($b+0.5)]);
	}
	# Return gradient
	return @grad;
}

sub make_gradient_legend {
	# Use graphics
	use GD;
	# Take number of steps and max/min rgbs
	my($steps,$rgb1,$rgb2) = @_;
	# Get rgb gradient
	my @gradient = make_rgb_gradient($steps,$rgb1,$rgb2);
	# Image object and color allocation
	my $width = 50;
	my $hight = $width * scalar(@gradient);
	my $image = new GD::Image(($width),($hight));
	# Draw color code legend
	my $x1 = 0;
	my $y1 = 0;
	my $x2 = $width;
	my $y2 = $width;
	foreach my $rgb (@gradient) {
		my $color = $image->colorAllocate($rgb->[0],$rgb->[1],$rgb->[2]);
		$image->filledRectangle($x1, $y1, $x2, $y2, $color);
		# Next row coordinates
		$y1 += $width;
		$y2 += $width;
	}
	# Create image
	my $image_file = "gradient.png";
	open(my $im, '>', $image_file);
	binmode($im);
	print($im $image->png);
	close($im);
}

1;