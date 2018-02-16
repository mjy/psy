package Psy::Io::Io;

use warnings;
use strict;
use Carp;

use Data::Dumper;
use Psy::Psy qw($PSY_BASE_DIR);

=head1 NAME

Psy::Io::Io 

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Quick summary of what the module does.

Perhaps a little code snippet.

    use Psy::Io::Io;

    my $foo = Psy::Io::Io->new();
    ...

=head1 EXPORT

=cut


use vars qw(@EXPORT_OK);
use Exporter;
@EXPORT_OK = qw(confirmDir);

=head1 FUNCTIONS

=head2 confirmDir

Confirms that the path exists by cd/building successive folders.

Creates the path string (split on '/') and leaves dir in the last folder created.
Ignores preceeding '/', i.e. paths built from current directory.

=cut


sub confirmDir { 
	my $path = shift;	
	return 'no path passed to confirmdir' if not defined $path; 
#	$path = substr($path, 1, length($path)) if  (substr($path,0,1) eq '/'); # strip preceeding / if passed
	
	# chdir($Psy::ROOT_DIR);
	
	my @dir = split (/\//, $path);

	foreach (@dir) {
		# if (length ($_) == 0) { chdir ($Psy::PSY_BASE_DIR) ; next};
		# print "[ $_ ]";
		chdir $_ or mkdir $_;
		chdir $_;
	}

	return $path;
	# thanks to Petruchio, jc of perlmonks for the hints here
	# rd $dir /qy 	#rd $dir /sy
	# By the way, for a non-recursive io_confirmdir solution: perl -le 'sub io_confirmdir ($) { chdir $_[0] or mkdir $_[0] } io_confirmdir "/home/ray/zee"'
	#http://secu.zzu.edu.cn/book/Perl/Perl%20Bookshelf%20%5B3rd%20Ed%5D/cookbook/ch09_09.htm
	
	# Actually, it seems that File::Find now contains an rmtree sub.
	# grinder - your recursive delete problem is best cured with File::Path qw/rmtree mkpath/
}


=head1 AUTHOR

'Matt, C<< <m{j}yoder@{tee}[aye](em)(you).domain4unis> >>

=head1 BUGS

Please report any bugs or feature requests to
C<bug-psy-io-io@rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Psy>.
I will be notified, and then you'll automatically be notified of progress on
your bug as I make changes.

=head1 ACKNOWLEDGEMENTS

=head1 COPYRIGHT & LICENSE

Copyright 2005 'Matt, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1; # End of Psy::Io::Io
