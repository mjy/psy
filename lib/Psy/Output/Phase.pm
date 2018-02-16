package Phase;

use strict;
use warnings;
use Data::Dumper;

use vars qw(@ISA);
use Psy::Dna::Alphabet;
use Psy::Io::Io;

@ISA = qw(output);

=head1 NAME

Psy::Output::Phase

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Output subclass for a Psy object.
Referenced externally only.
Inherited via AUTOLOAD through Psy::Output::Output, see therein for general usage.

Outputs a PHASE2 legal file and a stub controlfile.
	
=cut


our $OUT_DIR = "analyses/phase2";
our $LEGAL_CHARS =  Psy::Dna::Alphabet->new(); # defaults to the RNA alphabet

sub process {
	my $self = shift;
	my %raw_params = @_;
	my %default_params = (
		'-max_seq_length' => 100,
		'-include_mask' => 1, # doesn't make sense in all combinations
		'-strip_gaps' => 0,
		'-legal_alphabet' => $LEGAL_CHARS,
		'-u2t' => 1,					# translate U to T?
		'-clean' => 1,					# translate non legal_alphabet chars to N 
		'-replace' => 'N',				## should go with the alphabet object ultimately
		'-mask_type' => '2', 			# 1 - untouched characters, 2- 
		'-phase_pseudonot_mode' => 1,	## allows for pseudonots in PHASE 1 by reordering helices	
		'-outgroup_index' => 0,
		'-include_class' => 1
	);

	# requires -mx	
	if (not defined $raw_params{'-slice'}) {
			print "no slice passed for Phase, using default\n";
	 		$raw_params{'-slice'} = $raw_params{'-mx'}->origSlice};
		
	my %params = (%default_params, %raw_params);

	# strip unbracketed blocks from slice
	$params{'-slice'}->remove('Blocks', $params{'-mx'}->loopBlocks(%params, '-mode' => 'bracketed'));	

	my $path = $OUT_DIR;
	if (defined $params{'-sub_path'}) {$path = $params{'-sub_path'}};

	$params{'-file_name'} = "data.rna" if not defined $params{'-file_name'};

	#&io::confirmDir("analyses/phase2/results/"); 	# need to setup a results directory
	#chdir('../../..');

	print "PATH: $path/$params{'-file_name'}\n";

	# gather needed data
	## two modes - PHASE 1.1 / PHASE 2.0 (phancy mask)?
	
	my $data;					   
	$data = &output::mxData(%params); # get a basic data object

	$data->{datafile} = $params{'-file_name'};
	$data->{randseed} = int(rand(999));
	$data->{rnamodel} = 'RNA7A';

	$data->{outgroup} = $params{'-mx'}->{'terminals'}->{$params{'-outgroup_index'}}->label;

	# $data->{results_path} = '';
 	if ($params{'-include_mask'} == 1) { $data->{mask} .= $params{'-mx'}->{structures}->{'original'}->phase2Mask(%params)}; ## make structure variable not fixed?	
 	if ($params{'-include_class'} == 1) { $data->{class} .= $params{'-mx'}->{structures}->{'original'}->phase2Class(%params)}  ## make structure variable not fixed?
	
	my $tt = Template->new(); # don't use config, as it trims needed whitespace from end
	$tt->process(\*DATA, $data, "$path/$params{'-file_name'}") || die 'died trying to process PHASE template'; 

	print "output path: ", $self->basePath , "\n";
	
	# generate a MCMC control file
	my $tc = Template->new($self->ttConfig, 'TRIM' => 0); # don't use config, as it trims needed whitespace from end
	$tc->process("mcmcphase.control.tt", $data, "$path/mcmc.control" ) || die $tt->error(), "\n", 'died trying to process mcmcphase_control template';

	1;
}

__DATA__
[% total_ters %] [% total_chars %] STRUCT

[% mask %]

[% FOREACH ter = terminals %][% rowlabel(ter) %]
[% multi_row_seq.$ter %]

[% END %]

[% class %]


