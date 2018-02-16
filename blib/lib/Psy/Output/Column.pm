package Psy::Output::Column;

use strict;
use warnings;
use Data::Dumper;

use Psy::Dna::Alphabet;

use vars qw(@ISA);
@ISA = qw(output);

our $OUT_DIR = "analyses/stats/column";
our $LEGAL_CHARS = Psy::Dna::Alphabet->new(); # defaults to the RNA alphabet

=head1 NAME

Psy::Output::Column

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Output subclass for a Psy object.
Referenced externally only.
Inherited via AUTOLOAD through Psy::Output::Output, see therein for general usage.
	

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

	$default_params{'-slice'} ||= $raw_params{'-mx'}->origSlice;	
	my %params = (%default_params, %raw_params);

	my $path = $OUT_DIR;
    if (defined $params{'-sub_path'}) {$path = $params{'-sub_path'}};
	$params{'-file_name'} ||= 'out.txt';

 	print "PATH: $path/$params{'-file_name'}\n";

	# doesn't use template toolkit (yet)
	
	## fix
	my $data;
	
	my $tt = Template->new($params{'-config'}); 
	$tt->process(\*DATA, $data, "$path/$params{'-file_name'}") || die 'died trying to process FASTA template'; 
	1;
}

1;

__DATA__



