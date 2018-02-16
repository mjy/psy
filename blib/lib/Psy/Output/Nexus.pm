package Nexus;

use Data::Dumper;

use vars qw(@ISA);
@ISA = qw(output); # required

=head1 NAME

Psy::Output::Nexus

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS



Output subclass for a Psy object.
Referenced externally only.
Inherited via AUTOLOAD through Psy::Output::Output, see therein for general usage.


=cut


our $OUT_DIR = "analyses/nexus";
our $LEGAL_CHARS =  Psy::Dna::Alphabet->new(); # defaults to the RNA alphabet


sub process {
	my $self = shift;
	my %raw_params = @_;
	my %default_params = (
		'-max_seq_length' => 100,
		'-include_structure_line' => 0, # doesn't make sense in all combinations
		'-strip_gaps' => 1,
		'-legal_alphabet' => $LEGAL_CHARS,
		'-u2t' => 1, 			# translate U to T?
		'-clean' => 1, 			# translate non legal_alphabet chars to N 
		'-replace' => 'N',  	# should go with the alphabet object ultimately
		'-mode' => 'unbracketed'
	);
	# requires -mx

	$default_params{'-slice'} ||= $raw_params{'-mx'}->origSlice;	
	my %params = (%default_params, %raw_params);

	my $mx = $params{'-mx'}; # easier reference

    # should do this prior to calling perhaps...
    # remove blocks that aren't requested, via -mode, defaults to unbracketed only
	if ($params{'-mode'} eq 'unbracketed') {
		$params{'-slice'}->prune('-mode' => 'bracketed', '-mx' => $mx)
	}
	elsif ($params{'-mode'} eq 'bracketed') {
		$params{'-slice'}->prune('-mode' => 'unbracketed', '-mx' => $mx); 
	}

    $params{'-path'} ||= $OUT_DIR; 
	$params{'-file_name'} = 'out.nex' if not defined $params{'-file_name'};

 	print "PATH: $params{'-path'}/$params{'-file_name'}\n";
	
	# gather needed data
	my $data;					   
	$data = &output::mxData(%params); # get a basic data object

	my $tt = Template->new($self->ttConfig) || die $Template::ERROR, "\n"; 
	$tt->process('nexus.tt', $data, "$params{'-path'}/$params{'-file_name'}") || die $tt->error(), " died trying to process Nexus template\n"; 
	1;	
}

1;

