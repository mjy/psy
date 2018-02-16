package Plot;


use strict;
use warnings;
use Data::Dumper;
use Psy::Strings::Composition;

use GD;

 GD::Image->trueColor(1);
=head1 Psy::Graphics::Plot


=cut

=head1 VERSION

Version 0.01

=cut

=head1 SYNOPSIS

Hacks using GD, not  yet integrated  as with other packages


=cut

our $VERSION = '0.01';

=head2 new

Create a new generic plot object.

=cut


sub new {
	my $type = shift;
	my %params = @_;
  	my $self = {};
	bless $self, $type;
	
	$self->_init(%params);
	
    return $self; 
}


sub _init{
	1;
}



=head2 binaryBitmap

Takes a text file matrix composed of 1s and 0s and plots it to a image.  Use with a bare kword output.

	-infile => <foo>  	the input filename
	-outfile => <bar>	 the output filename
	-wordfile => <undef> a list of the actual kwords, if passed pretty pixelation occurs, must have the same number of lines as columns
	
	# note there is now smart path handling yet, see the plot.t file for examples
	
	-max_plot_len => 6000 # the maximum length of the column
	-y_offset => 10  # the number of blank pixel lines between horizontal plots
	
=cut

sub binaryBitmap {
	my $self = shift;

	my %raw_params = @_;
	my %default_params = (
		'-max_plot_len' => 6000,
		'-y_offset' => 10,
		'-wordfile' => '',
	);
	
	# open the file to count the lines
	my $lines = 0;
	my $cols = 0;
	
	open (IN, $raw_params{'-infile'}) || die ("couldn't open the file $raw_params{'-infile'}"); # ERROR HANDLING
		$lines += tr/\n// while sysread(IN, $_, 2 ** 16);
		$lines += 1;
		seek (IN, 0,0); # "rewind" the file
		$cols = length(<IN>) - 1 ; # length of the first line
		print "working with matrix of $cols X $lines\n";
	close (IN);

	my @words;
	# some word properties fu
	if ($raw_params{'-wordfile'}) {
		print "reading wordlist\n";
		open (IN, $raw_params{'-wordfile'}) || die ("couldn't open the file $raw_params{'-wordfile'}"); # ERROR HANDLING
		@words = <IN>;
		close (IN);
	}
	
	my %params = (%default_params, %raw_params);

	my $l = $params{'-max_plot_len'} > $cols ? $cols : $params{'-max_plot_len'} ;  # max length
	my $h = $lines + $params{'-y_offset'}; # the offset for subsequent columns, should equal to # terminals + # pixels whitespace of gap, if smaller than number of 
	

	my $b;
	
	my $ly = 0;
	my $lx = 0;
	my $y = 0;
		
    # create a new image
	my $im = new GD::Image($l, $lines);

  	# make the background transparent and interlaced
	#  $im->transparent($white);
    $im->interlaced('true');

	
    # allocate some colors
  	my $white = $im->colorAllocate(255,255,255);
  	my $black = $im->colorAllocate(0,0,0);       
    my $red = $im->colorAllocate(160,0,0);      
    my $blue = $im->colorAllocate(0,0,160);
    my $green = $im->colorAllocate(255,255,255);
    my $other = $im->colorAllocate(0,160,160);

	my $dark = $im->colorAllocate(128, 128, 128);
	my $light = $im->colorAllocate(255, 255, 255);
	
	my %dot_color = ( 'A' => $dark,
			'C' => $light,
			'G' => $dark,
			'U' => $light);
	my $c = $white;
	


	
	open (IN, $params{'-infile'}) || die ("couldn't open the file $params{'-infile'}"); # ERROR HANDLING
		while (<IN>) {
			chomp;
			my $line = $_;
			$line =~ s/\s//g;
			for (my $x = 0; $x < length($line); $x++) {
				
#				print $x, substr($words[$x],0,1), "\n" if $params{'-wordfile'} ;
#				print  $dot_color{substr($words[$x],0,1)} if $params{'-wordfile'};

				my $color;
				if ( $params{'-wordfile'} ) {
				 #	$c = $dot_color{substr($words[$x],0,1) };
					my $comp = Psy::Strings::Composition->new('-seqs' => [($words[$x])]);

					
					# print "comp:", $comp->sumPct('C', 'G'), " ", $words[$x];z
					
					$color = $dot_color{substr($words[$x],0,1)} ;
					
					# $color = int(255 - (192 * ( $comp->sumPct('A', 'U'))));
					# print " color: $color\n\n";
					#	print $color,"\n";
					#	print  $im->colorAllocate($color, $color, $color);
					#print  "\n";
					$comp = undef;
				}
			
				$color ||= $white;

				substr($line, $x, 1) == 0 ? $im->setPixel($lx, ($ly * $h) + $y, $black) : $im->setPixel($lx, ($ly * $h) + $y ,  $im->colorAllocate($color, $color, $color)); 
	
				$lx++;	
	
				if ($lx > $l) {
					$lx = 0;
					$ly++;
				}
			}
			
			# new line
			$y++;
			$lx = 0;
			$ly = 0;
		}
			
	close (IN);

	
	open (OUT, ">$params{'-outfile'}") || die ('blah'); # ERROR HANDLING
	# make sure we are writing to a binary stream
    binmode OUT;
	local (*STDOUT) = *OUT;
	 

    # Convert the image to PNG and print it on standard output
    print OUT $im->png;
	close (OUT);
	
}

1;
	
__DATA__
				if ($y <13) {
					$c = $green;
				}
				elsif ($y > 12 and $y < 52) {
					$c = $blue;
				}
				elsif ($y > 51 and $y < 77) {
					$c = $red;
				}

				elsif ($y > 76 and $y < 79) {
				 $c = $other;
				}
				else {
				 $c = $white;
				}


