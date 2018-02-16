package Psy::Helpers::Fasthmmer;

use warnings;
use strict;

use Data::Dumper;
use	Psy::Helpers::Fastmatrix;

our @ISA = ('Psy::Helpers::Fastmatrix');

=head1 NAME

Psy::Helpers::Fasthmmer - The great new Psy::Helpers::Fasthmmer!

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

# &_test;

Quick summary of what the module does.

Perhaps a little code snippet.

    use Psy::Helpers::Fasthmmer;

    my $foo = Psy::Helpers::Fasthmmer->new();
    ...

=cut


=head2 new

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
	my $self= shift;
	my %params = @_;

	my $current_ter;
	
	# read in, doesn't assume anything about taxa, just appends all seq data to the preceeding '>taxon'
	open ( IN, $params{'-infile'} ) || die print "cannot open $params{'-infile'} \n ";
		my $i = 0;
		while (<IN>) {	
			my $row = $_;
			chomp $row;
			
			$row =~ s/^(\s+)//; # strip space from begining of the row, Justin Case
				if ($row =~ /# STOCKHOLM/) {
   				while (<IN>) {
					my $row = $_;
					chomp $row;
					if ($row =~ /^(#=GC\sRF)(\s*)(\S+)/) { $self->RF($3) }
					if ($row =~ /^(\w+)(\s*)(\S+)/) { 
						$self->{mx}->{$1} .= $3;
						$i++;
					}
					next;
				}
			}
		}	
		
	print "read $i terminals\n";
	close IN;

	return 1
}

=head2 fixGaps

Replaces the xX. with -

=cut

sub fixGaps { 
	my $self = shift;
	map {$self->{mx}->{$_} =~ s/[xX\.]/-/gi;} ($self->loopTer);
	return 1
}

=head2 RF

# returns the RF line

=cut

sub RF { 
	my $self = shift;
	if (@_) {$self->{RF} = shift; }
	return $self->{RF}
}


=head2 trimByRf

Hmmer built alignments add a lot of gapped space on the front and back ends of the aligment, who knows why, strip this space

=cut

sub trimByRf {
	
	my $self = shift;
	my $rf = $self->RF;
	$rf =~ /^(\.*)x(.*)x(\.*)/;
    # get the lengths of the space on the front and back end- 
	# this should sum to the differences in observed length from the aligment the model was built on to the aligment that is returned with hmmalign
	my $f = length($1);
	my $b = length($3);
	my $l = $self->numCols;
	map {substr($self->{mx}->{$_},0, $f,'')} ($self->loopTer); # strip that space
	map {substr($self->{mx}->{$_}, $l-$f-$b, $b,'')} ($self->loopTer);
}


=head2 RFBoundries 

Return the boundaries of an RF line based on an array containing the number of fixed positions for each block
Not used and not complete (see blkmapped.pm)

=cut

sub RFBoundries { 
	my $self = shift;
	my %raw_params = @_;
	my %default_params = ('-char' => 'x');
	my %params = (%default_params, %raw_params);
	my $s = $self->RF;
	my @r = ();
	
	my $t = $s;
	
	print "\n\ntr: " , ($t =~ tr/x//);
	# print "$s\n";
	
	my $d;
	for my $i (@{$params{'-counts'}}) {
		$d += $i;
	}
	print "count: ", $d, "\n";
	
	$s =~ m/^([^$params{'-char'}]+)/gi; # match the longest string of characters that are not the char of interest starting at the begining
	push @r, length $1;
	print "\n\nFirst: ", length $1, " ", $1, "\n";
	substr($s,0,$r[0],''); # cut the first string off
	
	my @n = @{$params{'-counts'}};
	my $c = $params{'-char'};

	# find first
	# start loop
	#  grab first
	#  if next is . take positions
	
	foreach my $i (@n) {
		print $i;	
		# thanks to the perlmonks for help on this (see my nodes)
		my $regex = '\A' . ("[^$c]*?$c" x $i);
		print "$regex\n"; 
	
		if ($s =~ m{ ($regex) }xms ) {    # PBP orthodoxy :)
			push @n, length $1;
			substr($s, 0, length $1,"");
			print " $1\n\n";

			if ($s =~ m/^((\.)+)/gi) {
				push @n, length $1;
				substr($s, 0, length $1,"");
			}
		}

		else {
			die "danger! unxpected shortness of RF!!" 
        }
	}

	push @n, length $s; # capture the last piece
   	return @r; 	
}

sub _test {
	my $foo =fasthmmer->new('-infile' => 'result.foo');
	#print $foo->RF;
	my @bar = ([1,4,7]);
	$foo->nTransSplit('-counts' => @bar);
	#print Dumper($foo)
}


=head1 AUTHOR

'Matt, C<< <m{j}yoder@{tee}[aye](em)(you).domain4unis> >>

=head1 BUGS

Please report any bugs or feature requests to
C<bug-psy-helpers-fasthmmer@rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Psy>.
I will be notified, and then you'll automatically be notified of progress on
your bug as I make changes.

=head1 ACKNOWLEDGEMENTS

=head1 COPYRIGHT & LICENSE

Copyright 2005 'Matt, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1; # End of Psy::Helpers::Fasthmmer
