package Blockmeta;

use strict;
use warnings;

use Data::Dumper;

use vars qw(@ISA);
use Psy::Dna::Alphabet;
use Psy::Strings::Composition;
use Psy::Strings::Strings;

@ISA = qw(output);

=head1 NAME

Psy::Output::Blkmetanex

=head1 VERSION

Version 0.02

=cut

our $VERSION = '0.02';

=head1 SYNOPSIS

Output subclass for a Psy object.  
Referenced externally only.
Inherited via AUTOLOAD through Psy::Output::Output, see therein for general usage.
	
This output was originally concieved as a way to plot descriptive states for a given block on a tree (some nice phenetic stuff).
After some though/development on the system I found that the Lutzoni lab ARC software does more or less what I wanted to do.

Miadlikowska, Lutzoni, Goward, Zoller, Posada (2003) "New approach to an old problem: Gap-rich regions from ITS and rDNA large-subunit are incorporated into phylogenetic analyses to resolve the Peltigera canina species complex.", Mycologia 95(6):1181-1203.

Therefor, {Miadlikowska et al., 2003 should be given credit for much of the ideas coded herin, though with very few exceptions the implementation is mine (Matt).  I have retained the character order they present for ease of comparison.  

"F. Kauff, J. Miadlikowska & F. Lutzoni (2003), "ARC - a program for Ambiguous Regions Coding",
distributed by the authors (http://www.lutzonilab.net/pages/download.shtml), Dept. of Biology,
Duke University, USA" hmm- did he code this?

Their ARC characters are:

(1) sequence length,
(2)	A frequency	(no. of 'A' divided by sequence length)
(3)	C frequency	(no. of 'C' divided by sequence length)
(4)	G frequency	(no. of 'G' divided by sequence length)
(5)	T frequency	(no. of 'T' divided by sequence length)
(6)	AA frequency	(no. of 'AA' divided by sequence length-1)
(7)	AC frequency	(no. of 'AC' divided by sequence length-1)
(8)	AG frequency	(no. of 'AG' divided by sequence length-1)
(9)	AT frequency	(no. of 'AT' divided by sequence length-1)
(10)	CC frequency	(no. of 'CC' divided by sequence length-1)
(11)	CG frequency	(no. of 'CG' divided by sequence length-1)
(12)	CT frequency	(no. of 'CT' divided by sequence length-1)
(13)	GG frequency	(no. of 'GG' divided by sequence length-1)
(14)	GT frequency	(no. of 'GT' divided by sequence length-1)
(15)	TT frequency	(no. of 'TT' divided by sequence length-1)
(16)	A pairing	(no. of 'AA' divided by total no. of 'A')
(17)	C pairing	(no. of 'CC' divided by total no. of 'C')
(18)	G pairing	(no. of 'GG' divided by total no. of 'G')
(19)	T pairing	(no. of 'TT' divided by total no. of 'T')
(20)	A distribution	(no. of spaces between 'A' divided by total no. of 'A')
(21)	C distribution	(no. of spaces between 'C' divided by total no. of 'C')
(22)	G distribution	(no. of spaces between 'G' divided by total no. of 'G')
(23)	T distribution	(no. of spaces between 'T' divided by total no. of 'T')
(24+x)	motif presence [1] or absence [0] ### NOTE THAT THIS IS DEFINED a-priori ###


Noteable differences from their methods are as follows:
1) Psy has no amino acid functionality.
2) Psy can not read ARC formated input (yet)
3) The motif presence/absence (character 24+x) has been extended by allowing for automatic generation of motifs based on sampleing the data (kwords).
4) Several addition characters have been added (length, counts rather than pct).
5) Outputting is much more flexible, as it is integrated into overall slice/matrix functinos available to Psy.  Available output formats are now 2 types of nexus files formatted for use in Mesquite and at TNT readable format.
 
n) Pair counts are sliding in Psy, I need ot check if they are so in ARC.  e.g. 'AAAA' has 3 'AA' pairs, not 2.


Required parameters
	-mx
	

Optional parameters
	'-out_format' => < mesquite_chr_by_blk | mesquite_blk_by_chr | [tnt] | original_arc >

If a -slice is passed it will be used for taxa partition.  In absence of -plan it will also be used in place of origSlice (see below).

If no -plan is given each block? in the origSlice will be individually translated to metachars.

A -plan can be made as follows: 

	my $p;
	$p->{0}->{blks} = ([0 1 2 3 4]); # block 0 will be a fusion of blocks 0..4.  Required.  
	$p->{0}->{name} = 'foo';  # block 0 label will be 'foo'.  Not necessary, defaults to the index if not set.
	$p->{0}->{type} = 'trans'; # block 0 will be translated to meta characters or left as nucleotide < 'trans' | 'orig' >.  Not required, defaults to 'trans' if not set.
	$p->{1} ... 

	'-plan' => $p
	
	And used like: 
	
	my $o = output->new();
	$o->Blockmeta('-mx' => $mx, '-plan' => $p);

If $mx is a matrix, a plan can also be generated like:
	my $p = $mx->plan;


=cut

our $OUT_DIR = "analyses/blkmeta";
our $LEGAL_CHARS =  Psy::Dna::Alphabet->new(); # defaults to the RNA alphabet

sub process {
	my $self = shift;
	my %raw_params = @_;
	my %default_params = (
		'-gapchar' => '-',
		'-legal_alphabet' => $LEGAL_CHARS,
		'-u2t' => 1, 	# translate U to T? (MUST DO THIS)
		'-out_format' => 'tnt',
		'-translation_mode' => 'basic', # used with plans
	);
		   	
	my %params = (%default_params, %raw_params);

	$params{'-slice'} ||= $params{'-mx'}->origSlice;
	
	my $path = $OUT_DIR;
	$params{'-file_name'} = 'blkmeta.nex' if not defined $params{'-file_name'};

 	print "PATH: $path/$params{'-file_name'}\n";

	# simplify 
	my $mx = $params{'-mx'};
	
	# a ordered list (could use an ordered hash module ultimately
	my @chars = qw/
		pct_a
	   	pct_t 
	   	pct_g 
	   	pct_c
	    pct_aa 
	    pct_ac 
	    pct_ag 
	    pct_at 
	    pct_cc 
	    pct_cg 
	    pct_ct 
	    pct_gg 
	    pct_gt 
	    pct_tt 
		pct_a_pair 
		pct_c_pair 
		pct_g_pair
		pct_t_pair 
		pct_a_dist 
		pct_c_dist
		pct_g_dist 
		pct_t_dist 
		blk_length 
		presence/;

	$params{'-chars'} ||= \@chars;


	# check for a plan passed, build the default if not
	my $plan;
	if (not $params{'-plan'}) {
		$plan = $params{'-mx'}->plan('-mode' => $params{'-translation_mode'}, %params);
	}
	else {
		$plan = $params{'-plan'};
	}
			
	# gather needed data (used sparingly in Blkmeta, but used)
	my $data;					   
	$data = &output::mxData(%params); # get a basic data object
			
	# fuse the data and override the Output->cell method
	foreach my $t ( $params{'-slice'}->loop('Taxa')) {
		foreach my $blk ( keys %{$plan} ) {
			my $s = slice->new;
			$s->blocks( @{$plan->{$blk}->{blks}} );
			$s->taxa( $params{'-slice'}->loop('Taxa') );
			# for 'trans' get the data (strip gaps and translate u2t by default, do NOT pass %params)
			if 	($plan->{$blk}->{type} eq 'trans') {;
				$data->{mx}->{$t}->{$blk} = $mx->rawRowData('-slice' => $s, '-ter' => $t, '-u2t' => $params{'-u2t'}, '-gapchar' => $params{'-gapchar'} );
			}	
			else {
				$data->{mx}->{$t}->{$blk} = $mx->rawRowData('-slice' => $s, '-ter' => $t, '-u2t' => $params{'-u2t'});
			}
		}
		
	}
	$data->{cell} = sub { my ($y, $x) = @_; return $data->{mx}->{$y}->{$x} }; 
	
	# return the number of meta chrs you are using.
	$data->{'total_meta_chrs'} = $#{$params{'-chars'}}; 
	
	# you can loop by character or by block
	$data->{char_loop} = $params{'-chars'};

	# make a list of translated and untranslated blocks
	my @t_blks;
	my @ut_blks;
	map { ( $plan->{$_}->{type} eq 'trans' ) && push @t_blks, $_; } sort {$a <=> $b} ( keys %{$plan} );
	map { ( $plan->{$_}->{type} eq 'orig'  ) && push @ut_blks, $_; } sort {$a <=> $b} ( keys %{$plan} );	
	$data->{trans_blk_loop} = ([@t_blks]);
	$data->{orig_blk_loop} = ([@ut_blks]);

	# everything in place, just need to map characters to the function that returns their value
	# maps characters to subs,  $_[0] is a composition object loaded with the pertinent seq data, &f2d is just sprintf to 2 decimals
	my %char_subs = (
		pct_a => sub { &f2d( $_[0]->pct('A')) }, 
	   	pct_t => sub { &f2d( $_[0]->pct('T')) }, 
	   	pct_g => sub { &f2d( $_[0]->pct('G')) },
	   	pct_c => sub { &f2d( $_[0]->pct('C')) },
	    pct_aa => sub { &f2d( $_[0]->pctPair('-pair' => 'AA')) },
	    pct_ac => sub { &f2d( $_[0]->pctPair('-pair' => 'AC')) },
	    pct_ag => sub { &f2d( $_[0]->pctPair('-pair' => 'AG')) },
	    pct_at => sub { &f2d( $_[0]->pctPair('-pair' => 'AT')) },
	    pct_cc => sub { &f2d( $_[0]->pctPair('-pair' => 'CC')) },
	    pct_cg => sub { &f2d( $_[0]->pctPair('-pair' => 'CG')) },
	    pct_ct => sub { &f2d( $_[0]->pctPair('-pair' => 'CT')) },
	    pct_gg => sub { &f2d( $_[0]->pctPair('-pair' => 'GG')) },
	    pct_gt => sub { &f2d( $_[0]->pctPair('-pair' => 'GT')) },
	    pct_tt => sub { &f2d( $_[0]->pctPair('-pair' => 'TT')) },
		pct_a_pair => sub { &f2d( $_[0]->pairDivChar('-pair' => 'AA')) },
		pct_t_pair => sub { &f2d( $_[0]->pairDivChar('-pair' => 'TT')) },
		pct_g_pair => sub { &f2d( $_[0]->pairDivChar('-pair' => 'GG')) },
		pct_c_pair => sub { &f2d( $_[0]->pairDivChar('-pair' => 'CC')) },
		pct_a_dist => sub { &f2d( $_[0]->ltrDist('-ltr' => 'A')) },
		pct_t_dist => sub { &f2d( $_[0]->ltrDist('-ltr' => 'T')) },
		pct_g_dist => sub { &f2d( $_[0]->ltrDist('-ltr' => 'G')) },
		pct_c_dist => sub { &f2d( $_[0]->ltrDist('-ltr' => 'C')) },
		blk_length => sub {  $_[0]->len },  
		presence => sub { $_[0]->len > 0 ? '1' : '0'},
	);

	my $alphabet = Psy::Dna::Alphabet->new('-type' => 'dna');
	$alphabet->excluded('-'); # should be redundant as they are stripped earlier

	# add character references to our template data
	$data->{'chars'} = sub {
		my ($chr, $t, $blk) = @_;
	   	my $c = Psy::Strings::Composition->new( 
			'-seqs' => ([$data->{'mx'}->{$t}->{$blk}]),
			'-alphabet' => $alphabet,
			'-u2t' => $params{'-u2t'}
		);
	    return $char_subs{$chr}($c) # return the composition method specific to the character in question
	};
	

	# use fixed states and mimic the original ARC if requested
	
	if ($params{'-out_format'} eq 'original_arc') {
		# 0123456789abcdefghijklmnopqrstuvwxyz
		my $state_labels = '0123456789ABCDEFGHIJKLMNOPQRSTUV'; # use 32 states so we are tnt legal
		my $trans_chars = {};
			foreach my $blk ( @t_blks) { # loop the translated blocks and build up a list of all the states
				foreach my $c (@{$params{'-chars'}}) { # loop each character
					my $s = 0; # count of the number of states present for this block

					$trans_chars->{$blk}->{$c}->{'0.00'} = '?'; # zeroed values are not shared states
					$trans_chars->{$blk}->{$c}->{'0'} = '?'; # neither are length 0 or other zeroed states
					
					foreach my $t ($params{'-slice'}->loop('Taxa')) { # loop each taxon
						# add the state for a given taxon to a hash containing all states and assign a character
						# print join "\t", ($c,$t,$blk,"\n");
						# print $data->{'chars'}($c, $t, $blk);
										
						my $v = $data->{'chars'}($c, $t, $blk);
					
						if ($v != 0.00) {
							$state_labels =~ /.{$s}(.)/;
			
							unless ($trans_chars->{$blk}->{$c}->{$v}) { # assign a state for a value of $v
								$trans_chars->{$blk}->{$c}->{$v} = $1;
								$s += 1 
							}
						}
					}

					# if $s > length($state_labels) then set all characters in this column to ? by nuking all keys, could speed things slightly by escaping the loop
					if ($s > length($state_labels)) {
						push @{$data->{'too_many_states'}}, "($blk\:$s\:$c)"; 
						$trans_chars->{$blk}->{$c} = undef;
					}
				}
			}
			#	print Dumper($trans_chars);	

		# link up to the new method 
			
			$data->{'chars_orig'} = sub {
			my ($chr, $t, $blk) = @_;
			my $c = Psy::Strings::Composition->new( 
				'-seqs' => ([$data->{'mx'}->{$t}->{$blk}]),
				'-alphabet' => $alphabet,
				'-u2t' => $params{'-u2t'}
			);
			return $trans_chars->{$blk}->{$chr}->{$char_subs{$chr}($c)} # return the composition method specific to the character in question
			};
		}

	$data->{blockmeta_ver} = $VERSION;

	if (($params{'-out_format'} eq 'mesquite_chr_by_blk') or ($params{'-out_format'} eq 'tnt')  or ($params{'-out_format'} eq 'original_arc') ) {
		$data->{nchar_meta} = ($#{$data->{char_loop}} + 1) * ($#t_blks +1); 
	}	
	elsif ( $params{'-out_format'} eq 'mesquite_blk_by_chr' ) {
		$data->{nchar_meta} = $#t_blks +1; 
	}
	
	# count the number of characters that are standard nucleotide characters
	$data->{nchar_nuc} = 0;
	my $ind = @{[$params{'-slice'}->loop('Taxa')]}[0]; # grab a legal taxon index
	map {$data->{nchar_nuc} += ($data->{mx}->{$ind}->{$_} ? length $data->{mx}->{$ind}->{$_} : 0)} @ut_blks;
	
	(($params{'-out_format'} eq 'tnt') or ($params{'-out_format'} eq 'original_arc')) && ($data->{nchar_total} = $data->{nchar_nuc} + $data->{nchar_meta});

	my %formats = ( 
		'tnt' => 'blkmeta_tnt.tt',
	   	'mesquite_chr_by_blk' => 'blkmeta_mesq_chrXblk.tt',
		'mesquite_blk_by_chr' => 'blkmeta_mesq_blkXchr.tt',
		'original_arc' => 'blkmeta_tnt_orig_arc.tt'
	);
	
	my $tt = Template->new($self->ttConfig) || die $Template::ERROR, "\n"; 
	$tt->process($formats{$params{'-out_format'}}, $data, "$path/$params{'-file_name'}") || die "died trying to process a blockmeta template.", $tt->error(); 
	1;
}

# the actual output is not DATA but found in the Template directory, DATA is used as a testbed

__DATA__
nstates cont ; 
xread
'
TNT formatted data
'
[% nchar_total %] [% total_ters %]

[% FOREACH blk = trans_blk_loop %]
&[cont]
[% FOREACH ter = terminals %][% FOREACH c = char_loop %][% justifiedrowlabel(ter) %][% chars(c, ter, blk) %][% END %][% END %]
[% END %]

&[dna]
[% FOREACH ter = terminals %]	[% justifiedrowlabel(ter) %] [% FOREACH blk = orig_blk_loop %] [% cell(ter, blk) %] [% END %]
[% END %]
'[block order [% FOREACH blk = orig_blk_loop %] [% blk %] [% END %] ]'
;
cc-.; 
proc/; 




