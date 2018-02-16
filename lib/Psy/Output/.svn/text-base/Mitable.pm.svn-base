package Psy::Output::Mitable;

use strict;
use warnings;
use Data::Dumper;

use vars qw(@ISA);
use Psy::Dna::Alphabet;

@ISA = qw(output);

our $OUT_DIR = "analyses/tables/mi";
our $LEGAL_CHARS =  Psy::Dna::Alphabet->new(); # defaults to the RNA alphabet

=head1 NAME

Psy::Output::Mitable

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

NOT DONE

Output subclass for a Psy object.
Referenced externally only.
Inherited via AUTOLOAD through Psy::Output::Output, see therein for general usage.
	
=cut

=head2 process

=cut

sub process {
	my $self = shift;
	my %raw_params = @_;
	my %default_params = (
		'-max_seq_length' => 100,
		'-include_structure_line' => 0, # doesn't make sense in all combinations
		'-strip_gaps' => 1,
		'-legal_alphabet' => $LEGAL_CHARS,
		'-u2t' => 1, 		# translate U to T?
		'-clean' => 1, 		# translate non legal_alphabet chars to N 
		'-replace' => 'N'  	## should go with the alphabet object ultimately
	);
	# requires -mx
		
	$default_params{'-slice'} = $raw_params{'-mx'}->origSlice if not defined $default_params{'-slice'};	
	my %params = (%default_params, %raw_params);

	my $path = $OUT_DIR;
    if (defined $params{'-sub_path'}) {$path = $params{'-sub_path'}};
	$params{'-file_name'} = 'mi.txt' if not defined $params{'-file_name'};

 	print "PATH: $path/$params{'-file_name'}\n";
	
	# gather needed data
	my $data;					   
	$data = &output::mxData(%params); # get a basic data object

	## FIX
	$data->{basepairs} = sub {my $blk = shift}; # just return an array indexing (not position) of each, so we can loop each

	my $tt = Template->new($params{'-config'}); 
	$tt->process(\*DATA, $data, "$path/$params{'-file_name'}") || die 'died trying to process mitable template'; 
	1;
}

__DATA__

block
[% FOREACH blk = fiveprimeblocks %] [% blk %] foo

[% END %]



