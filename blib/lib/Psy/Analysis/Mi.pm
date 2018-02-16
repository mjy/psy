package Psy::Analysis::Mi;

use warnings;
use strict;
use Data::Dumper;
#use Inline::Files -backup;

use vars qw(@ISA);
use Psy::Dna::Alphabet;
use Psy::Io::Io;

@ISA = qw(output);

=head1 NAME

Psy::Analysis::Mi - The great new Psy::Analysis::Mi!

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

    use Psy::Analysis::Mi;

    my $foo = Psy::Analysis::Mi->new();
    ...


=cut
	
our $OUT_DIR = "analyses/mi";
our $LEGAL_CHARS = Psy::Dna::Alphabet->new(); # defaults to the RNA alphabet

=head2 new

=cut


	sub new {
		my $type = shift;
		my @seqs = @_; # reference to array of seqs, column to return (starts at 1)

		# requires two strings, if only one is passed the other is assumed to be '-' throughout
				
		my $self = {};
		
		bless $self, $type;

		$self->_init(@seqs);
		
		return $self; 
	}


	#
	# M(X;Y) = sum (over i,j) p(xi, yj)log  ( P(xi,yj) / (P(xi)P(yj) )
	
__DATA__


=head1 AUTHOR

'Matt, C<< <m{j}yoder@{tee}[aye](em)(you).domain4unis> >>

=head1 BUGS

Please report any bugs or feature requests to
C<bug-psy-analysis-mi@rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Psy>.
I will be notified, and then you'll automatically be notified of progress on
your bug as I make changes.

=head1 ACKNOWLEDGEMENTS

=head1 COPYRIGHT & LICENSE

Copyright 2005 'Matt, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1; # End of Psy::Analysis::Mi
