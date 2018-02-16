package Kword;

use strict;
use warnings;
use Data::Dumper;

use vars qw(@ISA);
 # use Psy::Dna::Alphabet;
use Psy::Strings::Strings;

@ISA = qw(output);

=head1 Psy::Output::Kword

=cut


=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Kword parsimony analysis.  Very likely been done before, but code is all original (and flexible!)

Required parameters
	-mx

Optional parameters:
	-data => %hash  # where keys are pointers to slices, blocks of each slice are concatonated together, if not passed then bracketed blocks only are converted
	
	-translation_mode => see Matrix::plan < [basic] | all | all_fused >
	
	-weight => < both | wrd_len | wrd_freq > ### NOT IMPLMENTED!! ###
	
	-out_format => < [tnt] |  mesquite | bare >
	-tax_slice => $slice # a slice which defines the taxa to include
	-kword_size \@ - the kwords sizes you want to explicitly include, all other sizes are discarded, if omitted all possible lengths (block length -1) are generated
	-no_uninformatives => < [1] | 0 > strips parsimony uninformative characters or not
	-write_word_list => 'filename' | undef  # if a filename is passed an additional file listing all words, one word per line is created, practical mostly for single partition dumps of bare kwords

  -zero_length_missing => < [0] | 1 > if 1 then sequences with no data are coded as ? rather than 0

=cut

our $OUT_DIR = "analyses/kword";
our $LEGAL_CHARS =  Psy::Dna::Alphabet->new('-type' => 'custom', '-alphabet' => 'ACGTU', '-replace' => '?'); # REQUIRED FOR kword!!


sub process {
	my $self = shift;
	my %raw_params = @_;
	my %default_params = (
    '-zero_length_missing' => 0,
		'-weight' => 'word_len',
		'-out_format' => 'tnt',
		'-tax_slice' => $raw_params{'-mx'}->origSlice,
		'-gapchar' => '-',
 		'-no_uninformatives' => 1,
		'-out_format' => 'tnt',
		'-translation_mode' => 'basic',
		'-write_word_list' => undef,
	);

	# requires -mx
	
	my %params = (%default_params, %raw_params);

	# hmm- revisit this, used in fusing slices	
	if (not $params{'-data'}) {
		my %all;
	 	map {my $s = slice->new; $s->blocks($_); $all{$_} = $s} $params{'-mx'}->loopBlocks('-mode' => 'bracketed');	
		$params{'-data'} = \%all;
	}

	$params{'-legal_alphabet'} ||= $LEGAL_CHARS;

	$params{'-path'} ||= $OUT_DIR;
	
	$params{'-file_name'} = 'kword' if not defined $params{'-file_name'};
	$params{'-slice'} ||= $params{'-mx'}->origSlice;
	
 	print "PATH: $params{'-path'}/$params{'-file_name'}\n";
	
	# simplify 
 	my $mx = $params{'-mx'};
	
	# gather needed data
	my $data;					   
	$data = &output::mxData(%params); # get a basic data object

	# check for a plan passed, build the default if not
	my $plan;
	if (not $params{'-plan'}) {
		$plan = $params{'-mx'}->plan('-mode' => $params{'-translation_mode'}, %params);
	}
	else {
		$plan = $params{'-plan'};
	}
			
	# get the data

	# fuse the data and OVERRIDE the Output->cell method
	foreach my $t ( $params{'-slice'}->loop('Taxa')) {
		foreach my $blk ( keys %{$plan} ) {
			my $s = slice->new;
			$s->blocks( @{$plan->{$blk}->{blks}});
			$s->taxa( $params{'-slice'}->loop('Taxa') );
			$data->{mx}->{$t}->{$blk} = $mx->rawRowData('-slice' => $s, '-ter' => $t);
			$data->{mx}->{$t}->{$blk} =~ s/-//g; # differs from Blockmeta here
		}
	}
	
	# loop through the kword blocks, doing foo on each cell
	# find the min-max length for each block (speedups vs. search)
	my $kw;	
	foreach my $blk (keys %{$plan}) {
		if ($plan->{$blk}->{kword} == 1) { 
			# my $min = 1; assume min is 1
			my $max = 0;
			foreach my $t ( $params{'-slice'}->loop('Taxa'))  {
				my $foo = length ($data->{mx}->{$t}->{$blk});
				$max = $foo if $foo > $max;	
				# find the min-max length for each block
			}		
			$kw->{$blk}->{max} = $max;	
		}
	}
	
	# print Dumper($kw);
	# set the kwords to find for each block (array)- if not passed {do all}
	
	my %kw_lengths; # keeps track of the actual lengths SEARCHED FOR, not necessarily found, this is however a set of possible ranges
	foreach my $blk (keys %{$kw} ) {
		if (not $params{'-kword_size'}) {
			$kw->{$blk}->{kw_size} = [ 1..$kw->{$blk}->{max}];
			map {	$kw_lengths{$_} = undef; } @{$kw->{$blk}->{kw_size}} 
		}
		else { # ... not proofed
			# only allow kwords that are possible for the block size
			foreach my $v (@{$params{'-kword_size'}}) {				
				if ($v <= $kw->{$blk}->{max}) {
					push @{$kw->{$blk}->{kw_size}}, $v;
					$kw_lengths{$v} = undef;
				}
			}
		}  
   	}
	
	
	
	# now we know what kwords to find- go for it!
	# add character references to our template data
	$data->{'chars'} = sub {
		my ($chr, $t, $blk) = @_;
 
    # if the blk is all missing data, and the appropriate flag is set, score as missing
    if ((not $data->{mx}->{$t}->{$blk} =~ /[^\?]/) and ($params{'-zero_length_missing'} == 1)) { # the first regex checks for the existence of at least one other non "?" (all - are ? at this point)
      return "?"; 
    }

  	$data->{mx}->{$t}->{$blk} =~ /($chr)/; # checks for the presence of the word
    $1 ? "1" : "0"; 
	};
	
	# create the kword character matrix- map it to the blocks, this is the master character variable now
	foreach my $blk (keys %{$kw} ) { # loop the blocks
		foreach my $l (@{$kw->{$blk}->{kw_size}}) { # loop the word size
			foreach my $t ( $params{'-slice'}->loop('Taxa'))  { # loop the taxa
				my %words = &kword( $params{'-legal_alphabet'}->clean('-str' => $data->{mx}->{$t}->{$blk}, %params), $l );  # get the words
				$kw->{$blk}->{chars}->{$_} += $words{$_} for keys %words; # add the words to the hash for that block
			}
		}

		# nuke words containing '?' ** note that '-' are previously removed!
		foreach my $k (keys %{$kw->{$blk}->{chars}} ) {
			delete $kw->{$blk}->{chars}->{$k} if $k =~ /\?/g
		}	
		
		# trim down the characters if we don't want uniformative characters
		if ( $params{'-no_uninformatives'} == 1) { # nuke all the uninformatives (only method working at present) 
			foreach my $k (keys %{$kw->{$blk}->{chars}} ) { # single instances can be eliminated right away 
				delete $kw->{$blk}->{chars}->{$k} if $kw->{$blk}->{chars}->{$k} == 1;
			}		
			
			# multiple hits in same block need to be eliminated as autapomorphies as well
			foreach my $k (keys %{$kw->{$blk}->{chars}} ) { # loop through words
				my $tot = 0;
				foreach my $t ( $params{'-slice'}->loop('Taxa'))  { # loop through taxa
					$tot++ if $data->{'chars'}($k, $t, $blk) eq "1";	
				}
				
				delete $kw->{$blk}->{chars}->{$k} if $tot == $params{'-slice'}->total('Taxa'); # everyone has it
				delete $kw->{$blk}->{chars}->{$k} if $tot == 1; #  1 is an autapomorphy
				delete $kw->{$blk}->{chars}->{$k} if $tot == ($params{'-slice'}->total('Taxa') - 1); # 0 is an autapomorphy	
			}
		}
	}
	
	# calculate the total characters
	$data->{total_chars} = 0;
	foreach my $k ( keys %{$kw} ) {
		$data->{total_chars} += $#{ [keys %{$kw->{$k}->{chars}}] } + 1 > 0 ? ($#{ [keys %{$kw->{$k}->{chars}}] } + 1) : 0; # ($#{ [keys %{$kw->{$k}->{chars}}] } + 1) if ($#{ [keys %{$kw->{$k}->{chars}}] } > 0);
	}

	$data->{trans_blk_loop} = [sort {$a eq $b} keys %{$kw}] ;
	
	# alphabetize the k-words
	$data->{blk_chrs} = sub {my $blk = shift; return [ sort {length $a <=> length $b || $a cmp $b } keys %{$kw->{$blk}->{chars}}  ] };
		
	# add a refence to the k-word sizes
	$data->{kword_size} = [sort {$a <=> $b} keys %kw_lengths] ; 
	
	my %formats = ( 
		'tnt' => 'kword_tnt.tt',
		'mrbayes' => 'kword_mrbayes.tt',
	  'mesquite' => 'kword_mesquite.tt',
		'bare' => 'kword_bare.tt',
	);

	$params{'-file_name'} = "$params{'-file_name'}.tnt" if $params{'-out_format'} eq 'tnt';
	$params{'-file_name'} = "$params{'-file_name'}.nex" if $params{'-out_format'} eq 'mesquite';
	$params{'-file_name'} = "$params{'-file_name'}.nex" if $params{'-out_format'} eq 'mrbayes';
	
	my $tt = Template->new($self->ttConfig) || die  $Template::ERROR , "\n"; 
	$tt->process($formats{$params{'-out_format'}}, $data, "$params{'-path'}/$params{'-file_name'}") || die $tt->error(), " died trying to process Kword template\n"; 

	if ($params{'-write_word_list'}) {
		print "generating the kword list\n";
		open (WORDS, ">$OUT_DIR/$params{'-write_word_list'}") || die ("couldn't open the wordlist file $params{'-write_word_list'}"); # ERROR HANDLING
			foreach my $blk (@{$data->{trans_blk_loop}}) {
				foreach my $w (@{$data->{blk_chrs}($blk)}) {
					print WORDS "$w\n";
				}
			}
		close (WORDS);
		
	}
	1;
}

# weighting schemes
# 1- inversely proportional to length
# 1- inversely proportional to frequency (remember 2 is the least common)
# 1- (1- (max_length -  length / max_length)) 

# DATA used in testing only, see templates

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


