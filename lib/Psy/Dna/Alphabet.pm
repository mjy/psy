package Psy::Dna::Alphabet;

use warnings;
use strict;
use Data::Dumper;

# our @ISA = qw(Psy);

=head1 NAME

Psy::Dna::Alphabet - The great new Psy::Dna::Alphabet!

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Quick summary of what the module does.

Perhaps a little code snippet.

    use Psy::Dna::Alphabet;

    my $foo = Psy::Dna::Alphabet->new(
										'-type' => 'custom',
									   	'-alphabet' => 'ACGT');

	Valid types:	
	-type  (valid: dna, rna, iupac, custom)
										
	Valid params for type custom:
	-alphabet (string)
	-excluded (string)
	
=cut


=head2 new

return/reset the alphabet (included)

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
	my %default_params = ('-type' => 'rna');
		
	my %params = (%default_params, %raw_params);
	my @empty = \();
	
	if ($params{'-type'} eq 'custom') {
		$self->type('custom');
		$self->alphabet($params{'-alphabet'});
		$self->excluded($params{'-excluded'});	
	}
	
	elsif ($params{'-type'} eq 'dna') {
		$self->alphabet('ATGC-?N');
		$self->type('dna');
	}
	
	elsif ($params{'-type'} eq 'iupac') {
		$self->alphabet('ACTGURYSWKMBDHVN.-?') ; # period is identity?
		$self->type('iupac');

	}
	else { # default to rna 
		$self->alphabet('AUCG-?N');
		$self->type('rna');
	}
	
	return 1
}

sub _verify {
	my $self=shift;
	my %params = @_;
	foreach my $p (keys %params) { # check matrix parameters here
		print $p;
	}
	return %params;
}

=head2 type

Accessor for alphabet type.

=cut

sub type { 
	my $self = shift;
	if (@_)  {$self->{type} = shift }
	return $self->{type};
}

=head2 alphabet

return/reset the alphabet (included)

=cut

sub alphabet { 
	my $self = shift;
	if (@_)  {$self->{alphabet} = shift; $self->{len} = $self->_lengthAlphabet }
	return split //,($self->{'alphabet'});
}

=head2 strOk

Returns true if the string passed in is composed only of 'in' characters

=cut

sub strOk { 
	my $self = shift;
	my $str = shift;
	my $ok = quotemeta($self->{alphabet});
	return 1 if ($str =~ /(?)[^$ok]/);
	0;

}

=head2 in

Returns true unless $l is not found in alphabet (included).

=cut

sub in {
	my ($self, $l)  = @_;
	$l = quotemeta($l);
	($self->{alphabet} =~ /$l/) && (return 1);
	return 0;	
}

=head2 inEx

Returns true unless $l is not found in excluded alphabet.

=cut

sub inEx {
	my ($self, $l)  = @_;
	$l = quotemeta($l);
	($self->{alphabet} =~ /$l/) && (return 1);
	return 0;	
}

=head2 alphabetLength

Returns length of alphabet (included)

=cut

sub alphabetLength {
	my $self = shift;
	return $self->{len};		
}

=head2 excluded

Return/reset explicitly excluded characters.  Returns a reference to an array.

=cut

sub excluded { 
	my $self = shift;
	if (@_)  {$self->{excluded} = shift}
	if (defined $self->{excluded}) {return split //,$self->{excluded}}
	else {return ''};
}

=head2 add

Add a char(class) to the alphabet

=cut

sub add { 
	my $self = shift;
	my $chars = shift;
	$self->remove($chars); # so they don't double up
	$self->{alphabet} .= $chars ; 
	$self->{len} = $self->_lengthAlphabet ;
	return 1;
}

=head2 remove

Remove a char(class) from the alphabet

=cut

sub remove { # 
	my $self = shift;
	my $chars = shift;
	$self->{alphabet} =~ s/[$chars]//gi;
	$self->{len} = $self->_lengthAlphabet;
	return 1;
}

=head2 addExclude

Add a char(class) to the explictily excluded characters 

=cut

sub addExclude { #
	my $self = shift;
	if (@_) {$self->{excluded} .= shift};
	return 1;
}

=head2 removeExclude

Add a char(class) to the explicitly excluded characters

=cut

sub removeExclude {
	my $self = shift;
	my $chars = shift;
	$self->{alphabet} =~ s/[$chars]//gi;
	return 1;
}

sub _lengthAlphabet { ## hmm- redundant?
	my $self = shift;
	return length $self->{alphabet};
}

=head2 possibleBPs

Returns all legal BPs, with -o flag included all BPs including 'o' combinations

=cut

sub possibleBPs {
	my $self = shift;
	my %bps;
	foreach my $chr1 ($self->alphabet) {
		foreach my $chr2 ($self->alphabet) {
			my $bp = join "", sort ($chr1,$chr2);
			$bps{$bp} = undef;
		}
	}
	return sort (keys %bps);
}

=head2 clean

Cleans the string of all characters (case invariant) not in the alphabet, replacing the characters of -str with -replace
-tu = 1 additionally swaps all Ts to Us if $self->type = 'rna' ## REDUNDANT??

=cut

sub clean {
	my $self = shift;
	my %raw_params = @_;
	my %default_params = ('-replace' => '',  '-tu' => 1);
	
	return undef if not defined $raw_params{'-str'};

	my %params = (%default_params, %raw_params);
	my $chars = quotemeta (join "", $self->alphabet );
	
	(($self->type eq 'rna') && ($params{'-tu'} == 1)) && ($params{'-str'} =~ s/[tT]/U/gi );

	$params{'-str'} =~ s/[^$chars]/$params{'-replace'}/gi;
	$params{'-str'} = uc($params{'-str'});

	return $params{'-str'};
}

	
=head1 AUTHOR

'Matt, C<< <m{j}yoder@{tee}[aye](em)(you).domain4unis> >>

=head1 BUGS

Please report any bugs or feature requests to
C<bug-psy-dna-alphabet@rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Psy>.
I will be notified, and then you'll automatically be notified of progress on
your bug as I make changes.

=head1 ACKNOWLEDGEMENTS

=head1 COPYRIGHT & LICENSE

Copyright 2005 'Matt, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1; # End of Psy::Dna::Alphabet
