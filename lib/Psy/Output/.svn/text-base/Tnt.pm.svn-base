package Tnt;

use strict;
use warnings;
use Data::Dumper;

use vars qw(@ISA);

use Psy::Dna::Alphabet;

# use column; ## ??

@ISA = qw(output);

=head1 NAME

Psy::Output::Tnt

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.02';

=head1 SYNOPSIS

Output subclass for a Psy object.
Referenced externally only.
Inherited via AUTOLOAD through Psy::Output::Output, see therein for general usage.

This is hacked together relatively quickly, the block choice could be simplified likely.

Requires:
	-mx

Optional: 
	-file_name (output file name)
	-u2t
	-legal_alphabet
	-mode => < bracketed | [unbracketed] | all >
	-replace
	
	
=cut

our $OUT_DIR = "analyses/tnt";
our $LEGAL_CHARS =  Psy::Dna::Alphabet->new; # defaults to the RNA alphabet
$LEGAL_CHARS->remove('-'); # we don't want to include - by defaults

sub process {
	my $self = shift;
	my %raw_params = @_;
	my %default_params = (
		'-u2t' => 1, # translate U to T when == 1
		'-legal_alphabet' => $LEGAL_CHARS,	
		'-slicecharmode' => 'all', # used in mx->totalChars # required (jan31/07), at the point of using this param we always want all chars, i.e. data are removed before this
		'-mode' => 'unbracketed',
		'-replace' => '?' 		   # used with alphabet, replace al non LEGAL_CHARS with ? for tnt
		);

	# requires -mx

	$default_params{'-slice'} ||= $raw_params{'-mx'}->origSlice;	
	my %params = (%default_params, %raw_params);
	
	my $path = $OUT_DIR;

	$params{'-file_name'} ||= 'psy';
	
 	print "PATH: $path/$params{'-file_name'}\n";

	my $mx = $params{'-mx'};
	

	# remove blocks that aren't requested, via -mode, defaults to unbracketed only
	if ($params{'-mode'} eq 'unbracketed') {
		$params{'-slice'}->prune('-mode' => 'bracketed', '-mx' => $mx)
	}
	elsif ($params{'-mode'} eq 'bracketed') {
		$params{'-slice'}->prune('-mode' => 'unbracketed', '-mx' => $mx); # remove unbracketed blocks (they are by default not aligned) for tnt;
	}
	
	my $data;					   
	$data = &output::mxData(%params); # get a basic data object

	#print Dumper($data);
	
	my $tt = Template->new($self->ttConfig) || die $Template::ERROR, "\n";  	
	$tt->process('tnt.tt', $data, "$path/$params{'-file_name'}.tnt") || die 'died trying to process Tnt template', $tt->error(); 
	1;
}


## need to do safenames here

# ccode
#  /1-[0       /1-[1       /1-[2       /1-[3      /1-[4       /1-[5   
#  /1-[6       /1-[7       /1-[8       /1-[9      /1-[10      /1-[11   
# ;

1;

__DATA__
nstates dna;
xread
'some version bit here'
[% total_chars %] [% total_ters %]
[% FOREACH ter = terminals %][% justifiedrowlabel(ter) %][% row(ter) %]
[% END %];
