package Psy::Helpers::Fastoneline;

use warnings;
use strict;
use Data::Dumper;

use Psy::Strings::Strings;
use Psy::Helpers::Fastmatrix;
# use Psy::Helpers::Gbname; ## might want to remove this

our @ISA = ('Psy::Helpers::Fastmatrix'); # required

=head1 NAME

Psy::Helpers::Fastoneline!

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

=cut


=head2 new

=cut

sub new {
	my $type = shift;
	my %params = @_;
  	my $self = {};
# 	$self->{Interleaves} = {};	
	bless $self, $type;
	
	$self->_init(%params);

    return $self; 
}

sub _init {
	my $self = shift;
	my %params = @_;
	
	my $current_ter;
	$params{'-infile'} || die 'no file to open';
	# read in, doesn't assume anything about taxa, just appends all seq data to the preceeding '>taxon'
	open ( IN, $params{'-infile'} ) || die print "cannot open [ $params{'-infile'} ] \n";
		my $i = 0;
		while (<IN>) {	
			my $row = $_;
			chomp $row;
			$row =~ s/^(\s+)//; # strip space from begining of the row, Justin Case
			($row =~ /^\[/)  && next;
			(length $row == 0) && next;
			$i++;
			
			$row =~ /(\w*)\s*(.*)/;
			$current_ter = $1;
			my $data = $2;
			$data ||= '';

			$self->seq($current_ter, $data);
			
		}	
		print "read $i terminals\n";
	close IN;
	return $i;
}

1;


=head1 AUTHOR

'Matt, C<< <m{j}yoder@{tee}[aye](em)(you).domain4unis> >>

=head1 BUGS

Please report any bugs or feature requests to
C<bug-psy-helpers-fastfasta@rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Psy>.
I will be notified, and then you'll automatically be notified of progress on
your bug as I make changes.

=head1 ACKNOWLEDGEMENTS

=head1 COPYRIGHT & LICENSE

Copyright 2005 'Matt, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1; # End of Psy::Helpers::Fastfasta
