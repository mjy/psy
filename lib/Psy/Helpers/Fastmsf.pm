package Psy::Helpers::Fastmsf;

use warnings;
use strict;
use Data::Dumper;

use Psy::Helpers::Fastmatrix;
our @ISA = ('Psy::Helpers::Fastmatrix');

=head1 NAME

Psy::Helpers::Fastmsf

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

A fastmatrix subclass for msf formatted files as output by clustalw 1.81

Translates periods to dashes (no periods currently allowed in tax labels).

=cut

=head2 new

=cut


sub _init {
	my $self = shift;
	my %params = @_;
	
	my $current_ter;
	
	# read in, doesn't assume anything about taxa, just appends all seq data to the preceeding '>taxon'
	open ( IN, $params{'-infile'} ) || die print "cannot open $params{'-infile'} \n ";
		while (<IN>) {	
			if ($_ =~ '//') {
				my $i = 0;
				while (<IN>) {		
					my $row = $_;
					chomp $row;
					next if $row eq '';
					$row =~ s/\./\-/g; # translate . to /

					$row =~ /(\w*)\s*(.*)/i;
				
					$self->{mx}->{$1} .= $2;
				}
			}
			else {
				next
			}
		}	
		map {$self->{mx}->{$_} =~ s/\s*//g} keys %{$self->{mx}}; # strip the msf whitespace
	close IN;
}

1;


=head1 AUTHOR

'Matt, C<< <m{j}yoder@{tee}[aye](em)(you).domain4unis> >>

=head1 BUGS

=head1 ACKNOWLEDGEMENTS

=head1 COPYRIGHT & LICENSE

Copyright 2005 'Matt, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1; # End of Psy::Helpers::Fastfasta
