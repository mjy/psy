
package Psy::Strings::Composition;

use warnings;
use strict;

use Data::Dumper;

use vars qw(@ISA);
use Psy::Dna::Alphabet;
use Psy::Strings::Strings; #  qw(delFixedPositions);

=head1 NAME

Psy::Strings::Composition 

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Composition based functions for sequences of letters (single string(.

Integrated with Psy::Dna::Alphabet object.

Concatenates all strings in passed to it and treats them as one.

    use Psy::Strings::Composition;

    my $foo = Psy::Strings::Composition->new(
		'-seqs' => \@strings
	);
    ...

	Letters in both excluded and alphabet are treated as excluded.
	
=head1 object

=cut

=head2 new

Takes reference to array of seqs.


=cut

sub new {
	my $type = shift;
	my %params = @_; 
	my $self = {};
	
	bless $self, $type;

	$self->_init(%params);
	
	return $self; 
}

sub _init {
	my $self = shift;
	my %raw_params = @_;
	my %default_params = (
							'-u2t' => 0,      # u2t two states - 0, leave/convert to U/ 1 leave/convert to t
							'-no_other' => 1  # delete all characters that would be other as well 
						);	 	 
							
	my %params = (%default_params, %raw_params);
	$params{'-alphabet'} ||= Psy::Dna::Alphabet->new('-type' => 'rna');

	my $alph;
	$alph = $params{'-alphabet'};

	my $seq = '';
	map {$seq .= $_ } @{$params{'-seqs'}};
	$self->inLen(length $seq);	## this would be original length only
	($params{'-u2t'} == 1) ? $seq =~ s/[uU]/T/gi : $seq =~ s/[tT]/U/gi;
	
	# print "$seq - ";
	
	my $ex = join '', $alph->excluded;
	$ex && ($seq =~ s/[$ex]//g); # if the alphabet has excluded letters, we can get rid of those 
	
	
	$self->seq($seq); # store the original sequence
	length($seq) == 0 and return; # no further work needed
	
	foreach my $l (split //, $self->{seq}) {
		$self->{ltr}->{$l}->{count}++;
	}

	(not $params{'-no_other'}) && ($self->{ltr}->{other}->{count} = 0);
	
	## BELOW IS WORKING BUT REDUNDANT, as EXCLUDED CHARS are excluded before you hit this.
	# simplify the hash now
	# there are three categories: in | excluded | other,
	# !in && !ex = other
	# iin && !ex = next
	# !in && ex = delete and next 
	foreach my $l ($self->ltrs) {
		$alph->in($l) && next; # quotemeta done on alaphabet side
		if ($alph->inEx($l)) {delete $self->{ltr}->{$l}; next};

		(not $params{'-no_other'}) && ($self->{ltr}->{other}->{count} += $self->{ltr}->{$l}->{count});
		delete $self->{ltr}->{$l}
	}

	foreach my $l ($self->ltrs) {
		$self->{len} += $self->{ltr}->{$l}->{count};
	}
	
	foreach my $l ($self->ltrs) {
		$self->pct($l, $self->count($l) / $self->len);
	}
	
	$self->ranked(sort {$self->count($a) <=> $self->count($b) } ($self->ltrs));

	# original seq vs counted seq
}

=head2 seq

Returns a string containing the concatanated input

=cut


sub seq {
	my $self = shift;
	if (@_) { $self->{seq} = shift}
	return $self->{seq}
}


=head2 ltrs

all the letters *counted*

=cut


sub ltrs { 
	my $self = shift;
	return sort (keys %{$self->{ltr}});
}


=head2 ranked

Returns array of ltr(s) ordered from least to most.

=cut


sub ranked { 
	my $self = shift;
#	print Dumper($self);
	if (@_) {$self->{ranked} = \@_}
	return () if not $self->{ranked};
	return @{$self->{ranked}} ;
}


=head2 inLen

length of incoming letter# s

=cut

sub inLen { 
	my $self=shift;
	if (@_) { $self->{inlen} = shift} ;
	return $self->{inlen}
}


=head2 len

Returns integer length of final processed letters (minus excluded).

=cut

sub len { 
	my $self=shift;
	if (@_) { $self->{len} = shift} ;
	return $self->{len} || 0
}	

=head2 count

Returns the count of letter $l

=cut

sub count { 
	my $self = shift;
	my $l = shift;
	if (@_) { $self->{ltr}->{$l}->{count} = shift };
	return  $self->{ltr}->{$l}->{count} || 0;
	
}

=head2 pct

Return ltr by percent for letter $l

=cut

sub pct { # return ltr by percent
	my $self = shift;
	my $l = shift;
	if (@_) { $self->{ltr}->{$l}->{pct} = shift };
	return $self->{ltr}->{$l}->{pct} || 0;
}


=head2 sumPct

Return sum of pct ltr for letters in @l

=cut

sub sumPct { # return ltr by percent
	my $self = shift;
	my @letters = @_;
	my $sum = 0;
	for my $l (@letters) {
		$sum += $self->{ltr}->{$l}->{pct} || 0 # could use some error checking here
	}
	return $sum;
}



=head2 least

Return an array of the least common ltr(s) (in case of ties > 1). 

=cut

sub least { 
	my $self = shift;
	my @most;
	my @ranked = $self->ranked;
	@most = shift (@ranked);
	foreach my $l (@ranked) {
		if ($self->count($l) ==  $self->count($most[0]) ) { push @most, $l }
		else { return @most; }
	}
}

=head2 most

Return the most common ltr(s)

=cut

sub most {  
	my $self = shift;
	my @most;
	my @ranked = $self->ranked;
	@most = pop (@ranked);
	foreach my $l (@ranked) {
		if ($self->count($l) == $self->count($most[0])) { push @most, $l }
		else { return @most; }
	}
}


=head2 pctPair

Sensu: "F. Kauff, J. Miadlikowska & F. Lutzoni (2003), "ARC - a program for Ambiguous Regions Coding",
distributed by the authors (http://www.lutzonilab.net/pages/download.shtml), Dept. of Biology,
Duke University, USA"
 

=cut

sub pctPair {
	my $self=shift;
	my %params = @_;
	$params{'-pair'} || die 'no -pair to pctPair';
    my %h = &kword($self->seq, 2);
	$h{$params{'-pair'}} || return 0;
 	return 0 if ($self->len - 1) == 0; ## mayhaps -1 or undef
	return sprintf("%.2f", ($h{$params{'-pair'}} / ($self->len - 1))   )
}


=head2 pairDivChar

Sensu: "F. Kauff, J. Miadlikowska & F. Lutzoni (2003), "ARC - a program for Ambiguous Regions Coding",
distributed by the authors (http://www.lutzonilab.net/pages/download.shtml), Dept. of Biology,
Duke University, USA"
 
From their example: (no. of 'AA' divided by total no. of 'A')

=cut


sub pairDivChar {
	my $self=shift;
	my %params = @_;
	$params{'-pair'} || die 'no -pair to pctPair';
    my %h = &kword($self->seq, 2);
	$h{$params{'-pair'}} || return 0;
    return 0 if	 $self->count(substr($params{'-pair'},1,1)) == 0;
 	return  sprintf("%.2f", ( $h{$params{'-pair'}} / $self->count(substr($params{'-pair'},1,1))) );
}


=head2 ltrDist

Sensu: "F. Kauff, J. Miadlikowska & F. Lutzoni (2003), "ARC - a program for Ambiguous Regions Coding",
distributed by the authors (http://www.lutzonilab.net/pages/download.shtml), Dept. of Biology,
Duke University, USA"

From above: 'A distribution	(no. of spaces [= chars] between 'A' divided by total no. of 'A')'

=cut 


sub ltrDist {
	my $self=shift;
	my %params = @_;
	$params{'-ltr'} || die 'no -ltr to pctPair';
	
	# count by regex
	my $r = "$params{'-ltr'}";
    $r .= "(.*?)$params{'-ltr'}" x ($self->count($params{'-ltr'}) -1 );
	# print "regex $r\n";
	my $s = $self->seq;
 	my $space = 0;
	map {$space += length ($_)} $s =~ /$r/;
 	return 0 if $self->count($params{'-ltr'}) == 0 ;
	
	# print "space: $space\n";
 	return sprintf("%.2f",  ($space / $self->count($params{'-ltr'})) );
}


=head1 AUTHOR

'Matt, C<< <m{j}yoder@{tee}[aye](em)(you).domain4unis> >>

=head1 BUGS

Please report any bugs or feature requests to
C<bug-psy-strings-composition@rt.cpan.org>, or through the web interface at
<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Psy>.
I will be notified, and then you'll automatically be notified of progress on
your bug as I make changes.

=head1 ACKNOWLEDGEMENTS

None yet.

=head1 COPYRIGHT & LICENSE

Copyright 2005 'Matt, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1; # End of Psy::Strings::Composition
