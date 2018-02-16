package Pairedparsimony;

use strict;
use warnings;
use Data::Dumper;
use Carp;
#use diagnostics;

use vars qw(@ISA);
use Psy::Dna::Alphabet;

@ISA = qw(output);

our $OUT_DIR = "analyses/pairedparsimony/";
our $LEGAL_CHARS =  Psy::Dna::Alphabet->new('-type' => 'custom', '-alphabet' => 'ACGU-'); 

=head1 NAME

Psy::Output::Pairedparsimony

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.02';

=head1 SYNOPSIS

Translates the dataset following the 20 state model of Smith et. al 2003. 

Output subclass for a Psy object.  
Referenced externally only.
Inherited via AUTOLOAD through Psy::Output::Output, see therein for general usage.

Required parameters
	-mx

Optional parameters
	-out_format => [nexus | tnt]
	-mode => [paired_only | all_aligned ]
	-mi_char_weights (does nothing at present)
	
=cut

sub process {
	my $self = shift;
	my %raw_params = @_;
	my %default_params = (
		'-out_format' => 'nexus',
		'-mode' => 'paired_only', # 'paired_only, all_aligned
		'-mi_char_weights' => 0 # include a mi character weighting
	);

	$default_params{'-slice'} ||= $raw_params{'-mx'}->origSlice; 	
	my %params = (%default_params, %raw_params);
	 
	$params{'-path'} ||= $OUT_DIR;

	$params{'-file_name'} ||= 'pairedp.nex';
 	
	print "output path: $params{'-path'}\n";

	# gather needed data
	my $data;					   
	$data = &output::mxData(%params); # get a basic data object

	my ($totalchars);


	my %translate = qw(
			A A
			C R
			G N
			U D
			AA C
			AC Q
			AG E
			AU G
			CA H
			CC I
			CG L
			CU K
			GA M
			GC F
			GG P
			GU S
			UA T
			UC W
			UG Y
			UU V
			X- -
			XX -
			-X -
			-- -
			A- A
			-A A
			C- R
			-C R
			G- N
			-G N
			U- D
			-U D
			XA -
			AX -
			UX -
			XU - 
			XG -
			GX -
			XC -
			CX -
		);

	# TNT can only handle characters A-V?
	my	%translate_tnt = qw(
			A A
			C B
			G C
			U D
			AA E
			AC F
			AG G
			AU H
			CA I
			CC J
			CG K
			CU L
			GA M
			GC N
			GG O
			GU P
			UA Q
			UC R
			UG S
			UU T
			X- -
			XX -
			-X -
			-- -
			A- A
			-A A
			C- B
			-C B
			G- C
			-G C
			U- D
			-U D
			XA -
			AX -
			UX -
			XU -
			XG -
			GX -
			XC -
			CX -
		);
 
		## ask joe about translation of XA etc.
		
		my $mx;
	   	$mx	= $params{'-mx'}; # shorten for ease of use
			
		foreach my $t ($params{'-slice'}->loop("Taxa")) {
			foreach my $blk ($params{'-mx'}->{'structures'}->{'original'}->loopHelices(%params) ) { 
				my @bps;	
				my @clean;
				@bps = $mx->{structures}->{'original'}->loopBasePairs(%params, '-blk' => $blk, '-ter' => $t);
				foreach (@bps) {push @clean, $LEGAL_CHARS->clean('-replace' => 'X', '-str' => uc($_))};	
	
				foreach my $bp (@clean) {
					#print $bp, " ", $translate{$bp}, "\n";
					if ($params{'-out_format'} eq 'tnt') {
						$data->{rowdata}->{$t} .= $translate_tnt{$bp};	
					}
					elsif ($params{'-out_format'} eq 'nexus') {
						$data->{rowdata}->{$t} .= $translate{$bp};	
					}
					else {
					 die 'not a legal out_format'
					}
				}
			}
			$data->{rowdata}->{$t} =~ s/\-/\?/gi;
		}
	
	$data->{total_chars} = length $data->{rowdata}->{0}; # override because we have fewer characters ## a kludge
	
	# print Dumper($self->ttConfig); #;params{'-config'});

	my %formats = (  # maps params to templates
		'tnt' => 'paired_tnt.tt',
	   	'nexus' => 'paired_nexus.tt',
	);

	my $tt = Template->new($self->ttConfig)  || die  $Template::ERROR , "\n";  # || die $tt->error(), "\n";  #params{'-config'}); 
	$tt->process($formats{$params{'-out_format'}}, $data, "$params{'-path'}/$params{'-file_name'}") || croak $tt->error(), " died trying to process Pairedparasimony template"; 
	1;
}

# The template below is for testing only, the true templates are in /templates

__DATA__
#Nexus

Begin data;
Dimensions
	ntax = [% total_ters %] 
	nchar = [% total_chars %]
;
Format
	datatype = standard
	missing = ?
	gap = -
	symbols = "A R N D C Q E G H I L K M F P S T W Y V X"
;


Matrix
[% FOREACH ter = terminals %][% justifiedrowlabel(ter) %]  [% rowdata.$ter %]
[% END %]
;
end;

