package Psy::Strings::Strings;

use warnings;
use strict;
use Carp;
use Data::Dumper;

use vars qw(@EXPORT_OK @ISA);

=head1 NAME

Psy::Strings::Strings 

=head1 VERSION

Version 0.01

=cut


our $VERSION = '0.01';

=head1 SYNOPSIS

Quick summary of what the module does.


Functions for use on strings.

    use Psy::Strings::Strings;


=head1 EXPORT

=cut 


require Exporter;
@ISA = qw(Exporter);

our @EXPORT = qw(f2d kword lengthLongestString regexSplit reverseComp gappedRe delFixedPositions string2Blocks lineWordMeta uniqueChars wordsFromLineMeta padleft padright alphabetize arrayAsRange);

			
=head1 FUNCTIONS

=cut 



=head2 maxBound

returns the maximum length of a string (first match)

=cut


=head2 regexSplit

Splits a string on -regexes, inserting '-insert' between them. 

Starts at begining of string, each successive regex must be found after the first.
	
	&regexSplit('-str' => $s, '-regexes' => \@r, '-insert' => '  ')
	
I'm sure this could be done much more efficiently with look-ahead/behind.
	
=cut

sub regexSplit {
	my %raw_params = @_;
	my %default_params = (
			'-insert' => '  '
		);
	my %params = (%default_params, %raw_params);
	$params{'-str'} || return;
	
	my $s =	$params{'-str'};
	my $n=0;
	
	for my $r ( @{$params{'-regexes'}} ) {
		print "$r - ";
		if ( $s =~ /(.{$n})(.*)($r)(.*)/i ) {
			print " [$1] [$2] [$3] [$4] $n";

			$s = $1.$2.$params{'-insert'}.$3.$params{'-insert'}.$4;
			$n = length ($1) + length ($3) + 2;
		}
		print "\n";
	}
	
	# cheat - collapse >2 spaces to nothing	
	$s =~ s/\s{3,}/  /g;
	
	return $s;
}



=head2 reverseComp

A simple reverse compliment a string (DNA).  Only handles 'ACGTU' and '[]' and '\?' for regexes.  Mayhaps move to DNA?

=cut


sub reverseComp { 
	my $str = shift;	
	
	$str =~ s/[Aa]/1/ig;
	$str =~ s/[uU]/2/ig;
	$str =~ s/[gG]/3/ig;
	$str =~ s/[cC]/4/ig;
	$str =~ s/[tT]/5/ig;
	
	$str =~ s/1/T/ig;
	$str =~ s/2/A/ig;
	$str =~ s/3/C/ig;
	$str =~ s/4/G/ig;
	$str =~ s/5/A/ig;

	# regex fixing
	$str =~ s/\[/6/ig;
	$str =~ s/\]/7/ig;
	$str =~ s/6/\]/ig;
	$str =~ s/7/\[/ig;

	$str =~ s/\?\\/\\\?/ig; 
	
	return scalar(reverse($str));
}


=head2 gappedRe

Returns a single regex that should search against an set of  input patterns

=cut


sub gappedRe { # 
	my %raw_params = @_;
	my %default_params = ();
	my %params = (%default_params, %raw_params);	
	my $r;
	
#	print Dumper(%params);
	
	foreach my $blk (@{$params{'-strings'}}) {
#		print "$blk - ";
		$r .= "(-*?";
		map {$blk =~ /.{$_}(.)/; $1 ? ($r .= "$1\-*") : $r.=''} (0..length($blk));
		$r .= "-*?)";
	}
	return $r;
}


=head2 delFixedPositions

Given an array -strings and a character -char, deletes all fixed positions within that string for that character

=cut


sub delFixedPositions { 
	my %params = @_;
	my @del_pos = (0..(length(@{$params{'-strings'}}[0])-1));
	foreach my $s ( @{$params{'-strings'}}) {
		my @del_pos1;
		foreach my $p (@del_pos) {
			$s =~ /.{$p}(.)/;
			($1 eq $params{'-char'}) and push @del_pos1, $p;
		}
		@del_pos = @del_pos1; # we need only check positions that are still fixed
	}
	
	foreach my $p (@del_pos) {
		map { substr($_, $p, 1, '') } @{$params{'-strings'}};
		map {$_--} @del_pos; # we have removed a position!
	}
	return @{$params{'-strings'}}
}


=head2 string2Blocks

Given an array of integers return a block of each length

=cut


sub string2Blocks { 
	# requires -lengths, -str
	my %params = @_;
	my @ar;

	my $pos = 0;
	foreach my $i ( @{$params{'-lengths'}}) {
		# print "$i, $pos\n";
		push @ar, substr($params{'-str'}, $pos, $i);
   		$pos += $i;
		#	$pos > length $params{'-str'} and return 0;
	}
	return @ar;
}


=head2 pieces

duplicates string2Blocks in a cleaner way...doh

=cut


sub pieces { # 
	my $str = shift;
	my @sizes = @_;
		my @r;
	foreach my $c (@sizes) {
		push @r, substr($str,0,$c,'');
	}
	return @r;
}

=head2 lineWordMeta

returns a hash of 'word index' => { 'word_end', 'word_start', 'char_sum_start', 'char_sum_end', 'word_length' }

parameters are '-line' , '-offset' (added to char_sum_start and char_sum_end), '-regex' (the word/pattern to search for)	
with help from http://perlmonks.org/index.pl?node_id=299194

=cut


sub lineWordMeta { 
	# ?? *? +? {}?
	
	my %raw_params = @_;
	my %default_params = (
			'-line' => '',
   			'-offset' => 0,
		   	'-regex' => '\S+' 
		); # '\b(\S)+'	
	my %p = (%default_params, %raw_params);

	my $hash;
	
	my ($i, $j) = (0) x 2;
	while ($p{'-line'}  =~ /($p{'-regex'})/g) { # this was the suggested original form (?=($p{'-regex'}))/g
		$hash->{$i}->{word_start} = $-[1];
	   	$hash->{$i}->{word_end} = $+[1];
		$hash->{$i}->{char_sum_start} = $j + $p{'-offset'};
		$hash->{$i}->{char_sum_end} = $j + $+[1] - $-[1] - 1 + $p{'-offset'};
		$hash->{$i}->{word_length} =  $+[1] - $-[1];
		
		$j += $+[1] - $-[1];
		
	   	$hash->{$i}->{fc} = substr($1, 0, 1);
	    $i++;
	};
	# print Dumper($hash);
	return $hash;	
}


=head2 wordsFromLineMeta

see lineWordMeta
returns a list of words (ordered) as demarked by a hash produced by &lineWordMeta

=cut


sub wordsFromLineMeta {
	my $line = shift;
	my $hash = shift;
	my @words;
	
    # "0" by itself is a special case apparently
	foreach my $i (sort {$a <=> $b} keys %{$hash}) {
		#print "$i  $hash->{$i}->{char_sum_start}  $hash->{$i}->{char_sum_end} ";
        die "!! ERROR: '$line' \n doesn't work with substring bounds (start, length): ", $hash->{$i}->{word_start}, ", ", $hash->{$i}->{word_length}, "\n" if not length($line) >= $hash->{$i}->{word_length} and not $hash->{$i}->{word_start} > 0;
		my $w = substr($line, $hash->{$i}->{word_start}, $hash->{$i}->{word_length});
		#	(not $w) && return;
		push @words, $w; # $hash->{$i}->{word_end} - $hash->{$i}->{word_start});
	}
	#print join " ", @words;
	return @words;
}


=head2 arrayAsRange



=cut


sub arrayAsRange {  # need to change reference to this
	my @in = @_;
	my @array = ( sort {$a <=> $b}  (@in) );
	my $j = shift @array;
	my ($periods, $i);
	my $string = "($j";

	foreach $i (@array) {
		if ($i == $j + 1) {
			$periods = 1;
		} 
		else {
			if ($periods) {
				$string .= "..$j, $i";
			}
			else {
				$string .= ", $i";		
			}
			$periods = 0;
		}
		$j = $i;
	}
	if ($periods) { 
		$string .= "..$j)"
		} 
	else { 
		$string .= ")";
	}
	return $string;
}


=head2 alphabetize

Returns a string containing the alphabetized (read 'sort') letters of the input string.

	alphabetize('cba'); # returns abc


=cut

sub alphabetize {
	my $in = shift;
	my @str = split //, $in;
	return join "", sort @str;
}

=head2 posOffsets

Returns an array containing the offset positions for each instance of $chars in $string 
zero = 1

=cut

sub posOffsets () {
	my (
		$chars,
		$string,
	) = @_;
	
	my @offsets;
	while ( $string =~ /$chars/g ) {
	    	push @offsets, $-[0];
	}	
	return @offsets; 	
}

=head2 totalUniqueChars

Pilfered from "http://www.newts.org/~troc/perl/uniqchar.perl"
returns the number of unique characters in a string	
=cut


sub totalUniqueChars () { 
	my ($str) = shift;
	use integer;
	my %c;
	for (split //, $str) {$c{$_}=undef} scalar keys %c 
}


=head2 uniqueChars

returns an array of the unique characters seen in a string

=cut


sub uniqueChars { 
	my $str = shift;
	my @seen;		
	while ((length $str) != 0) {
    	$str =~ /(.)/g;
		my $c = $1;
#		print "found '$c'\n";
		push @seen, $c; # $seen{$c} = undef;
		$str =~ s/[$c]//g; # need for chars like ?
	}
#    print Dumper(@seen);
	return sort(@seen);
}


=head2 uniqueCharsC

returns ord() instead of chars

=cut


sub uniqueCharsC {
	my %seen = ();		
	my $str = shift;
	while ($str =~ /(.)/g) {
    	$seen{ord($1)}++;
	}
	return sort(keys %seen);
}


=head2 padleft

needs a list (text to pad, character to pad with (assuming single character!!), total length

=cut


sub padleft {
	my  ($intxt, $padchar, $maxchars) = @_;
 
	my $length = length ($intxt) || 0;
	for (my $i=0; $i<$maxchars-$length; $i++) {
		$intxt = $padchar.$intxt;
	}
	return $intxt;
}


=head2 padright

Returns a right padded string of a certain length
	&padright('foo', ' ', '24');

=cut


sub padright {	my  (
		$intxt,		# text to pad
		$padchar,	# pad char
		$maxchars	# total string length
	) = @_;

	my $length = length ($intxt) || 0;
	$intxt .= $padchar x ($maxchars - $length);
	return $intxt;
}


=head2 lengthLongestString

length of the longest string in an array

=cut


sub lengthLongestString () {
	my @array = @_;
	my $longest = 0;
	($#array == 0) && return -1;
	foreach my $i (@array) {
		$i || next; # handles length zeros
		$longest= length($i) if length($i) > $longest;
	}
	return $longest;
}


=head2 numMax

the max number in an array

=cut


sub numMax { 
	my $cur = shift;
    foreach my $num (@_) {
    	$cur = $num if $num > $cur;
     }
	return $cur;
}


=head2 kword

Returns a hash count of all strings of length $len occuring in $str, stepping by +1

This is equivalent to overlapping kwords? 

=cut


sub kword {
	my ($str, $len) = @_;
	(not $len) && return undef;
	my %r;
	for (my $i = 0;  $i < length($str) - $len +1; $i++) {
		$r{substr($str,$i,$len)}++
	}
	return %r;
}



=head2 f2d

Short form.  Return a sprintf to 2 decimals.

=cut

sub f2d {
	return  sprintf("%.2f", shift) ;
}



=head1 AUTHOR

'Matt, C<< <m{j}yoder@{tee}[aye](em)(you).domain4unis> >>

=head1 BUGS


Please report any bugs or feature requests to
C<bug-psy-strings-strings@rt.cpan.org>, or through the web interface at

<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Psy>.

I will be notified, and then you'll automatically be notified of progress on
your bug as I make changes.

=head1 ACKNOWLEDGEMENTS

Some here.

=head1 COPYRIGHT & LICENSE


Copyright 2005 'Matt, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1; # End of Psy::Strings::Strings
