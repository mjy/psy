package Psy::Misc::Colour;

use warnings;
use strict;

=head1 NAME

Psy::Misc::Colour

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Quick summary of what the module does.

Perhaps a little code snippet.

    use Psy::Misc::Colour;

    my $foo = Psy::Misc::Colour->new();
    ...

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 FUNCTIONS

=head2 hexFromRange

=cut

sub hexFromRange { # returns the  hex of value scaled to range from 0 to 255 (roughly computed)
	my ($i, $n) = @_;
		
	$i > $n && die "out of range in hexFromRange";
	return unpack ("H2", pack "c", int(255 * ($i/$n)) )
}

=head2 textFromIndex

=cut


sub textFromIndex { 
	my $i = shift;

	my @colors = qw(
					Red Green Blue Cyan Magenta Yellow Gray Black Orange 
					LIghtRed LightGreen LightBlue LightCyan LightMagenta LightYellow LightGray White Purple 
					DarkRed DarkGreen DarkBlue DarkCyan DarkMagenta Brown DarkGray Violet SeaGreen SlateBlue DarkYellow
			); 
	return $colors[$i % $#colors] || "no color";
}


=head1 AUTHOR

'Matt, C<< <m{j}yoder@{tee}[aye](em)(you).domain4unis> >>

=head1 BUGS

Please report any bugs or feature requests to
C<bug-psy-misc-colour@rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Psy>.
I will be notified, and then you'll automatically be notified of progress on
your bug as I make changes.

=head1 ACKNOWLEDGEMENTS

=head1 COPYRIGHT & LICENSE

Copyright 2005 'Matt, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1; # End of Psy::Misc::Colour
