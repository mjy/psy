package KwordProfile;

# a hack to generate profiles, same code is used in Kword.pm

use strict;
use warnings;
use Data::Dumper;

use vars qw(@ISA);
 # use Psy::Dna::Alphabet;
use Psy::Strings::Strings;

@ISA = qw(output);

=head1 Psy::Output::KwordProfile

fasta

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Required parameters
	-mx

Optional parameters:
	-data => %hash  # where keys are pointers to slices, blocks of each slice are concatonated together
	-tax_slice => $slice # a slice which defines the taxa to include
	-kword_size => \@ - the kwords sizes you want to explicitly include, all other sizes are discarded
	-plan
	-mode
	
=cut

our $OUT_DIR = "analyses/kwordprofile";
our $LEGAL_CHARS =  Psy::Dna::Alphabet->new('-type' => 'custom', '-alphabet' => 'ACGTU', '-replace' => '?'); # REQUIRED FOR kword!!


sub process {
	my $self = shift;
	my %raw_params = @_;
	my %default_params = (
		'-weight' => 'word_len',
		'-out_format' => 'tnt',
		'-tax_slice' => $raw_params{'-mx'}->origSlice,
		'-gapchar' => '-',
 		'-no_uninformatives' => 1
	);
	# requires -mx

	my %params = (%default_params, %raw_params);
	
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
	
	my $data;					   
	$data = &output::mxData(%params); # get a basic data object

	# check for a plan passed, build the default if not
	my $plan;
	if (not $params{'-plan'}) {
		$plan = $params{'-mx'}->plan('-mode' => 'basic', %params);
	}
	else {
		$plan = $params{'-plan'};
	}
			
	# get the data

	# fuse the data and override the Output->cell method
	foreach my $t ( $params{'-slice'}->loop('Taxa')) {
		foreach my $blk ( keys %{$plan} ) {
			my $s = slice->new;
			$s->blocks( @{$plan->{$blk}->{blks}} );
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
	
   	foreach my $blk (keys %{$kw} ) {
		if (not $params{'-kword_size'}) {
			$kw->{$blk}->{kw_size} = [ 1..$kw->{$blk}->{max}];
		}
		else { # ... not proofed
			# only allow kwords that are possible for the block size
			foreach my $v (@{$params{'-kword_size'}}) {				
				(push @{$kw->{$blk}->{kw_size}}, $v) if ($v <= $kw->{$blk}->{max});
			}
		}  
   	}
	
	# now we know what kwords to find- go for it!
	# add character references to our template data
	$data->{'chars'} = sub {
		my ($chr, $t, $blk) = @_;
		$data->{mx}->{$t}->{$blk} =~ /($chr)/;
		# print "$t $blk ($chr) " , $data->{mx}->{$t}->{$blk}, " [$1]\n";
		$1 ? 1 : 0;
	};
	
	# create the kword character matrix- map it to the blocks, this is the master character variable now
	foreach my $blk (keys %{$kw} ) { # loop the blocks
		foreach my $l (@{$kw->{$blk}->{kw_size}}) { # loop the word size
			foreach my $t ( $params{'-slice'}->loop('Taxa'))  { # loop the taxa
				my %words = &kword( $params{'-legal_alphabet'}->clean('-str' => $data->{mx}->{$t}->{$blk}, %params), $l      );  # get the words
				$kw->{$blk}->{chars}->{$_} += $words{$_} for keys %words; # add the words to the hash for that block
			}
		}
	
		foreach my $k (keys %{$kw->{$blk}->{chars}} ) {	# nuke words containing '?' ** note that '-' are previously removed!
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
					$tot++ if $data->{'chars'}($k, $t, $blk) == 1;	
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

	# calculate the total characters / word size
	
	my $ws;

	foreach my $blk (keys %{$kw} ) { # loop the blocks
		foreach my $k (keys %{$kw->{$blk}->{chars}} ) { # loop the words
			$ws->{$blk}->{length($k)} += 1;
		}
	}
	

	# foreach block print the number of words of a given size, prints to STDOUT
	
	foreach my $blk (sort {$a <=> $b} keys %{$ws} ) {
		print "$blk\n";
		foreach my $k (sort {$a <=> $b} keys %{$ws->{$blk}} ) {
			print $k, "\t", $ws->{$blk}->{$k}, "\n"
		}
	}
	

	#	print Dumper( $kw );
	#	print Dumper($ws);
	
	# $data->{trans_blk_loop} = [sort {$a <=> $b} keys %{$kw}] ;
	# $data->{blk_chrs} = sub {my $blk = shift; return [ sort {length $a <=> length $b || $a cmp $b } keys %{$kw->{$blk}->{chars}}  ] };
	1;
}

