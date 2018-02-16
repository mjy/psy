package Pairmeta;

use strict;
use warnings;

use Data::Dumper;

use vars qw(@ISA);
use Psy::Dna::Alphabet;
use Psy::Strings::Composition;
use Psy::Strings::Strings;

@ISA = qw(output);

=head1 NAME

Psy::Output::Pairmeta

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Output subclass for a Psy object.  
Referenced externally only.
Inherited via AUTOLOAD through Psy::Output::Output, see therein for general usage.

Required parameters
	-mx

Optional parameters
	'-out_format' => < {table} | ccode | nexus_weights >
	'-chars' => < [	mi cramersV chi2 chi2Sig ] >
	
=cut

our $OUT_DIR = "analyses/pairmeta";
our $LEGAL_CHARS =  Psy::Dna::Alphabet->new(); # defaults to the RNA alphabet

sub process {
	my $self = shift;
	my %raw_params = @_;
	my %default_params = (
		'-strip_gaps' => 1,
		'-legal_alphabet' => $LEGAL_CHARS,
		'-u2t' => 1, 		# translate U to T?
		'-clean' => 1, 		# translate non legal_alphabet chars to N 
		'-replace' => 'N',  	## should go with the alphabet object ultimately
		'-out_format' => 'table'
	);
		   	
	my %params = (%default_params, %raw_params);

	$params{'-slice'} ||= $params{'-mx'}->origSlice;
	
	my $path = $OUT_DIR;
	$params{'-file_name'} = "pairmeta_$params{'-out_format'}.txt" if not defined $params{'-file_name'};

 	print "PATH: $path/$params{'-file_name'}\n";

	# simplify 
	my $mx = $params{'-mx'};
	
	# a ordered list (could use an ordered hash module ultimately
	my @chars = qw/
	   	mi
		cramersV
		chi2
		chi2Sig	
	/;

	$params{'-chars'} ||= \@chars;
		
	# gather needed data (used sparingly in Blkmeta, but used)
	my $data;					   
	$data = &output::mxData(%params); # get a basic data object

#	$data->{cell} = sub { my ($y, $x) = @_; return $data->{mx}->{$y}->{$x} }; 

	my $cols = sub {
		my ($p, $blk) = @_;
		my @cols = $mx->origStructure->pairPos('-pos' => $p, '-blk' => $blk, '-mx' => $mx);
		my $l = $mx->column('-col' => $cols[0]);
	   	my $r = $mx->column('-col' => $cols[1]);
		
		return ($l, $r);
	};

	# you can loop by character or by block
	$data->{char_loop} = $params{'-chars'};

	# everything in place, just need to map characters to the function that returns their value
	# maps characters to subs,  $_[0] is a composition object loaded with the pertinent seq data, &f2d is just sprintf to 2 decimals
	
	# kludge but fast
	my %char_subs;
	if ( $params{'-out_format'} eq 'ccode' ) { # note that you can't have a weight of zero, you need to deactivate
		 %char_subs = (
			'mi' => sub { my $s = $_[0]->mi; $s == 0 ? '] ' : sprintf("%.0f", 100 * $s) }, 	
			'chi2' => sub {   $_[0]->chi2 },
			'chi2Sig' => sub { my $s = $_[0]->chi2Sig; $s == 0 ? '] ' : $s }, # need to deactivate zero characters here
			'cramersV' => sub { my $s = $_[0]->cramersV; $s == 0 ? '] ' : sprintf("%.0f", 100 * $s) }
		);
	}
	else {
		 %char_subs = (
			'mi' => sub { &f2d($_[0]->mi) }, 	
			'chi2' => sub { &f2d($_[0]->chi2) },
			'chi2Sig' => sub { $_[0]->chi2Sig },
			'cramersV' => sub { &f2d($_[0]->cramersV) }
		);
	}

	my $alphabet = Psy::Dna::Alphabet->new('-type' => 'rna');
	$alphabet->excluded('-');

	# add character references to our template data
	$data->{'chars'} = sub {
		my ($chr, $pair, $blk) = @_;
		my $sp = Psy::Analysis::Seqpair->new( &{$cols}($pair, $blk)); 
		return $char_subs{$chr}($sp) # return the composition method specific to the character in question
	};
	
	$data->{blockmeta_ver} = $VERSION;

#	if (($params{'-out_format'} eq 'mesquite_chr_by_blk') || ($params{'-out_format'} eq 'tnt') ) {
#		$data->{nchar_meta} = ($#{$data->{char_loop}} + 1) * ($#t_blks +1); 
##	}	
#	elsif ( $params{'-out_format'} eq 'mesquite_blk_by_chr' ) {
#		$data->{nchar_meta} = $#t_blks +1; 
#	}
	#
#	$data->{nchar_nuc} = 0;
#	map { $data->{nchar_nuc} += length $data->{mx}->{0}->{$_} } @ut_blks;
#	
#	($params{'-out_format'} eq 'tnt') && ($data->{nchar_total} = $data->{nchar_nuc} + $data->{nchar_meta});

	my %formats = ( 
		'table' => 'pairmeta_table.tt',
		'ccode' => 'pairmeta_ccode.tt',
		'nexus_weights' => 'pairmeta_nexus_weights.tt' 
	);

#$formats{$params{'-out_format'}},	
	my $tt = Template->new($self->ttConfig); 
	$tt->process($formats{$params{'-out_format'}}, $data, "$path/$params{'-file_name'}") || die  "died trying to process a blockmetanex template.", die $tt->error(); 
	1;
}

# the actual output is not DATA but found in the /template directory, DATA is used as a testbed

__DATA__
Pair meta data via Psy.
columns: [% FOREACH c = char_loop %] [% c %] [% END %]

[% FOREACH blk = helices %]
	block [% blk %]
	[% FOREACH pair = pairs_in_block(blk) %][% pair %] [% FOREACH c = char_loop %][% chars(c, pair, blk) %] [% END %]
	[% END %]
[% END %] 
	


