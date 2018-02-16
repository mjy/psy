package MrBayes;

use strict;
use warnings;
use Data::Dumper;

use vars qw(@ISA);
 # use Psy::Dna::Alphabet;
 # use Psy::Strings::Strings;

@ISA = qw(output);

=head1 Psy::Output::MrBayes


=cut


=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Translation to MrBayes format.

Required parameters
	-mx

Optional parameters:

	
=cut

our $OUT_DIR = "analyses/mrbayes";
our $LEGAL_CHARS =  Psy::Dna::Alphabet->new(); 


sub process {
	my $self = shift;
	my %raw_params = @_;
	my %default_params = (
		# add parameters here;
	);

	# requires -mx
	
	my %params = (%default_params, %raw_params);

	$params{'-legal_alphabet'} ||= $LEGAL_CHARS;

	$params{'-path'} ||= $OUT_DIR;
	
	$params{'-file_name'} = 'mrbayes.nex' if not defined $params{'-file_name'};
	$params{'-slice'} ||= $params{'-mx'}->origSlice;
	
 	print "PATH: $params{'-path'}/$params{'-file_name'}\n";
	
	# simplify 
 	my $mx = $params{'-mx'};
	
	# gather needed data
	my $data;					   
	$data = &output::mxData(%params); # get a basic data object
		
	# get the data

	print Dumper($data->{'pairs_in_block'}(9));
	

	
#	my %formats = ( 
#		'tnt' => 'kword_tnt.tt',
#		'mrbayes' => 'kword_mrbayes.tt',
#	   	'mesquite' => 'kword_mesquite.tt',
#		'bare' => 'kword_bare.tt',
#	);

	
	my $tt = Template->new($self->ttConfig) || die  $Template::ERROR , "\n"; 
	$tt->process('mrbayes.tt', $data, "$params{'-path'}/$params{'-file_name'}") || die $tt->error(), " died trying to process MB template\n"; 
1; # must return true
}



__DATA__
#Nexus

BEGIN TAXA;
	TITLE foo;
	DIMENSIONS NTAX=[% total_ters %];
	TAXLABELS
	[% FOREACH ter = terminals %] [% rowlabel(ter) %] [% END %]
;

END;

BEGIN CHARACTERS;
	TITLE  Blockmeta_[% c %];
	DIMENSIONS NCHAR=[% total_chars %];
	FORMAT DATATYPE = STANDARD GAP = - MISSING = ? SYMBOLS = "  0 1";
CHARSTATELABELS [% SET i = 1 %][% FOREACH blk = trans_blk_loop %][% FOREACH c = blk_chrs(blk) %] 
	[% i %] blk[% blk %]_[% c %][% SET i = (i+1) %],[% END %][% END %]
	;
	
	MATRIX
[% FOREACH ter = terminals %]
[% justifiedrowlabel(ter) %][% FOREACH blk = trans_blk_loop %][% FOREACH c = blk_chrs(blk) %][% chars(c, ter, blk) %][% END %]  [% END %][% END %]
;

END;



