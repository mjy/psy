package blkmapped;

use strict;
use warnings;
use Data::Dumper;

use vars qw(@ISA);
use Psy::Dna::Alphabet;
use Psy::Helpers::Fasthmmer;
use Psy::Strings::Strings;

@ISA = qw(output);

our $OUT_DIR = "analyses/blkmapped";
our $LEGAL_CHARS =  Psy::Dna::Alphabet->new(); # defaults to the RNA alphabet

=head1 NAME

Psy::Output::Blkmapped

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.02';

=head1 SYNOPSIS

*** under development ***

Maps a block pattern from a matrix to another matrix, splitting the new matrix into whitespace delimited blocks.

Used with hmmer or stockholm derived matrixes that have exquivalent masks.

Output subclass for a Psy object.
Referenced externally only.
Inherited via AUTOLOAD through Psy::Output::Output, see therein for general usage.
	
=cut

sub process {
	my $self = shift;
	my %raw_params = @_;
	my %default_params = (
		'-legal_alphabet' => $LEGAL_CHARS,
		'-u2t' => 1 		# translate U to T?
		# '-clean' => 1, 		# translate non legal_alphabet chars to N 
		# '-replace' => 'N'  	## should go with the alphabet object ultimately
	);
	# requires -mx AND a fastmatrix object
		
	$default_params{'-slice'} = $raw_params{'-mx'}->origSlice if not defined $default_params{'-slice'};	
	my %params = (%default_params, %raw_params);

	my $path = $OUT_DIR;
    if (defined $params{'-sub_path'}) {$path = $params{'-sub_path'}};
	$params{'-file_name'} = 'out_mapped_matrix.txt' if not defined $params{'-file_name'};

 	print "OUTPUT PATH: $path/$params{'-file_name'}\n";
	
	my $mx = $params{'-mx'};
	 my $h = $params{'-fastmatrix'};	
	
	# get the hmmer file
	$h->trimByRf;
	$h->fixGaps;
	my $lt = $h->longestTer;
	# $h->deleteFixedGapPositions;

	my @lengths = $mx->blkLengths('-mode' => 'all');

	# gather needed data
	my $data;					   
	# $data = &output::mxData(%params); # get a basic data object

	# handle the hmmer or incoming data
	
	# split the incoming
	$data->{terminals} = [$h->loopBySeqLength] ;
	$data->{rowlabel} = sub { my $t = shift; &strings::padright($t, " ", $lt + 1) };
	$data->{row} = sub {my $t = shift; join " ", &strings::string2Blocks('-str' => $h->seq($t), '-lengths' => \@lengths)};
	
	# format the mask	
	$data->{'mask'} = &strings::padright("", " ", $lt + 1).$mx->{structures}->{'original'}->sliceMask(%params, '-gap' => ' ');  
   		
	# print Dumper($data);

	# attach the model data 
	
	# get each taxon split one at a time from $h and map to $data
	# foreach my $t ($params{'-slice'}->loop("Taxa")) {
	
	my $tt = Template->new($params{'-config'}); 
	$tt->process(\*DATA, $data, "$path/$params{'-file_name'}") || die 'died trying to process blkmapped template'; 
	1;
}

__DATA__
[% mask %]
[% FOREACH ter = terminals %][% rowlabel(ter) %][% row(ter) %]
[% END %]

