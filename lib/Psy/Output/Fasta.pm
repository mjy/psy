=head1 Psy::Output::Fasta

fasta

=head1 VERSION

Version 0.01

=cut

$VERSION = '0.01';

=head1 SYNOPSIS

Basic FASTA output.  By default returns all the data in the slice (including bracketed regions), unaligned, in a multi-line format. 

=cut

package Fasta;

use strict;
use warnings;
use Data::Dumper;

use vars qw(@ISA);
use Psy::Dna::Alphabet;

@ISA = qw(output);

our $OUT_DIR = "analyses/fasta";
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
	$params{'-file_name'} = 'out.fas' if not defined $params{'-file_name'};

 	print "PATH: $path/$params{'-file_name'}\n";
	
	# gather needed data
	my $data;					   
	$data = &output::mxData(%params); # get a basic data object
	
	# add some custom formatted data for FASTA format
	if ($params{'-include_structure_line'} == 1) {
		$data->{mask} = ">mask";
	   	$data->{mask} .= $params{'-mx'}->{structures}->{'original'}->sliceMask(%params); ## make structure variable not fixed?
	}
	else {$data->{mask} = "";}
		
	my $tt = Template->new($self->ttConfig); 
	$tt->process('fasta.tt', $data, "$path/$params{'-file_name'}") || die 'died trying to process FASTA template'; 
	1;
}

__DATA__
[% mask %]
[% FOREACH ter = terminals %]>[% rowlabel(ter) %]
[% multi_row_seq.$ter %]

[% END %]

