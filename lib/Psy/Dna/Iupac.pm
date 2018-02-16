package Psy::Dna::Iupac;

use warnings;
use strict;
use Data::Dumper;

use vars qw(@ISA @EXPORTER @EXPORT_OK $VERSION);

@ISA = qw(Psy);


=head1 NAME

Psy::Dna::Iupac

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

A static meta data construct that facilitates too/from iupac translations

Case is not checked! lowercase input will not return true


    use Psy::Dna::Iupac;

    my $foo = Psy::Dna::Iupac->new();
    ...

=head1 EXPORT

=cut

=head3 notes

http://www.cellbiol.com/cgi-bin/oligo/html/iupac.txt

ambiguous_dna_values = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "M": "AC",
    "R": "AG",
    "W": "AT",
    "S": "CG",
    "Y": "CT",
    "K": "GT",
    "V": "ACG",
    "H": "ACT",
    "D": "AGT",
    "B": "CGT",
    "X": "GATC",
    "N": "GATC",
    }
ambiguous_rna_values = {
    "A": "A",
    "C": "C",
    "G": "G",
    "U": "U",
    "M": "AC",
    "R": "AG",
    "W": "AU",
    "S": "CG",
    "Y": "CU",
    "K": "GU",
    "V": "ACG",
    "H": "ACU",
    "D": "AGU",
    "B": "CGU",
    "X": "GAUC",
    "N": "GAUC",
    }

ambiguous_dna_complement = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "M": "K",
    "R": "Y",
    "W": "W",
    "S": "S",
    "Y": "R",
    "K": "M",
    "V": "B",
    "H": "D",
    "D": "H",
    "B": "V",
    "X": "X",
    "N": "N",
    }


=cut

=head2 new

Return a new object.

=cut

sub new  {
	my $type = shift;
	my $self = {};

	my %iupac = (
		'A'    => 'A',
		'C'    => 'C',
		'T'    => 'T',
		'G'    => 'G',
		'AG'   => 'R',
		'CT'   => 'Y',
		'GT'   => 'K',
		'AC'   => 'M',
		'CG'   => 'S',
		'AT'   => 'W',
		'CGT'  => 'B',
		'AGT'  => 'D',
		'ACT'  => 'H',
		'ACG'  => 'V',	
		'ACGT' => 'N',
	);
	
	foreach my $key (keys %iupac) {
		$self->{$key} = $iupac{$key};
	}

	bless $self, $type;
	return $self; 		
}
	
sub _verify {
	my $self=shift;
	my %params = @_;
	foreach my $p (keys %params) { # check matrix parameters here
		print $p;
	}
	return %params;
}

=head2 nuc

Returns the nuc(s) corresponding to the Iupac code, should be integrated into the dataobject for speed

=cut

sub nuc { 
	my $self = shift;
	my $in = shift;

	return 0 if not defined $in;
	return 0 if not grep $_ eq $in, values %{$self};
	my %t = reverse %{$self};
	return $t{$in};
}


=head2 iupac

Returns the corresponding IUPAC code, input must be alphabeticall!!!

=cut


sub iupac { #
	my $self = shift;
	my $in = shift;
	
	$in =~ s/U/T/g;

	return '' if not defined $in;	
	return $self->{$in} if $self->legalnuc($in);
	return 'X';
}


=head2 tstv

A table/report test

=cut

sub tstv {
	my $self = shift;
	my $in = shift;
	return "--" if not $self->legalnuc($in);
	
	$in =~ s/U/T/g;

	if (grep $_ eq $in, qw/A C G T AG CT/) {
		return 'TS';
	}
	return 'TV';
}


=head2 legalnuc 

Returns 0/1 depending on  whether $_ contains a legal nucleotide character.

=cut

sub legalnuc {
	my $self = shift;
	my $in = shift;
	
	$in =~ s/U/T/g;

	#$in = &str_alphabetize($in);	# don't assume its input string is in the right order	
	return 0 if not defined $in;
	if (grep $_ eq $in, keys %{$self} ) {
		return 1;		
	}
	return 0;
}


=head2 legaliupac

Checks that $in is a legal iupac code

=cut

sub legaliupac { 
	my $self = shift;
	my $in = shift;
	return 'NULL IN' if not defined $in;
	if (grep $_ eq $in, values %{$self} ) {
		return 1;		
	}
	return 0;	
}


=head2 loopIupac

Returns an array of the legal IUPAC characters

=cut

sub loopIupac {
	my $self= shift;
	return sort values %{$self};
}

=head1 AUTHOR

'Matt, C<< <m{j}yoder@{tee}[aye](em)(you).domain4unis> >>

=head1 BUGS

Please report any bugs or feature requests to
C<bug-psy-dna-iupac@rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Psy>.
I will be notified, and then you'll automatically be notified of progress on
your bug as I make changes.

=head1 ACKNOWLEDGEMENTS

=head1 COPYRIGHT & LICENSE

Copyright 2005 'Matt, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1; # End of Psy::Dna::Iupac
