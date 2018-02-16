package Stockholm;

use strict;
use warnings;
use Data::Dumper;

use vars qw(@ISA);
use Psy::Dna::Alphabet;

@ISA = qw(output);

=head1 NAME

Psy::Output::Stockholmn

=head1 VERSION

Version 0.02

=cut

our $VERSION = '0.02';

=head1 SYNOPSIS

Outputs to Stockholm format, see http://www.cgb.ki.se/cgb/groups/sonnhammer/Stockholm.html for details.

Referenced externally only.
Inherited via AUTOLOAD through Psy::Output::Output, see therein for general usage.
	
=cut


our $OUT_DIR = "analyses/stockholm";
our $LEGAL_CHARS =  Psy::Dna::Alphabet->new(); # defaults to the RNA alphabet

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
	$params{'-file_name'} ||= 'out.stk';

	# gather needed data
	my $data;					   
	$data = &output::mxData(%params); # get a basic data object

	$data->{mask} = $params{'-mx'}->{structures}->{original}->sliceMask(%params, '-bracket_chars' => {'left' => '<', 'right' => '>'});
 	$data->{alignedMask} =  $params{'-mx'}->alignedMask(%params);

	my $tt = Template->new($params{'-config'}); 
	$tt->process(\*DATA, $data, "$path/$params{'-file_name'}") || die 'died trying to process STOCKHOLM template'; 
	1;

}

# hmmbuild foo.hmm out.fas


__DATA__
# STOCKHOLM 1.0
#=GF ID    test\
#=GF SOURCE psy
[% FOREACH ter = terminals %][% justifiedrowlabel(ter) %][% row(ter) %]
[% END %]
#=GC SS_cons [% mask %]
#=GC RF [% alignedMask %]
//
