package Psy::Helpers::Fastgenbank;

use warnings;
use strict;
use Data::Dumper;
use Psy::Helpers::Fastmatrix;
use Carp;
use Psy::Helpers::Gbname; 

my @ISA = qw(Psy::Helpers::Fastmatrix); ## necessary?

=head1 NAME

Psy::Helpers::Fastgenbank 

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Quick summary of what the module does.

Perhaps a little code snippet.

    use Psy::Helpers::Fastgenbank;
    my $foo = Psy::Helpers::Fastgenbank->new();

	Translate genbank to Fasta: 
	
    $foo->out_unBoundMx('-outfile' => 't/out/out.fas');

=cut

=head2 new

=cut


sub new {
	my $type = shift;
	my %params = @_;
  	my $self = {};

	bless $self, $type;
	
	return $self if $self->_init(%params);

	croak('failed to initialize fastgenbank object');
}


sub _init {
	my $self = shift;
	my %params = @_;
	
	my ($curseq, $accession, $common_name, $locus, $seq);
	my @classification;
	$accession = '';
	
	open ( IN, $params{'-infile'} ) || die print "cannot open $params{'-infile'} \n ";
		my $i = 0;
		
		# you can only use the ACCESSION index once a full record is parsed, as the LOCUS line does 
		# not always share the ACCESSION number (though in many cases it does)
		# if you want to access by LOCUS line a few small tweaks in the write section can be made
		
		while (<IN>) {	
			my $row = $_;
			chomp $row;	
			$row =~ s/^(\s+)//; # strip space from begining of the row, Justin Case
						
			if ($row =~ /^LOCUS/) { 
				# write the (previous) sequence
				unless ($accession eq "") { # don't write the first time round!
					
					#	$self->accession($curseq, $seq); # load the seq;
					$self->seq($accession, $seq); # load the seq;
					$self->locus($accession, $locus);
					$self->commonName($accession, $common_name);
					$self->classification($accession, @classification);

					# wipe previous holders
					($accession, $locus, $common_name, $seq) = ('') x 4;
					@classification = ();
				}
				
				$i++;
				my @locus = split (/\s+/, $row);
				$locus = $locus[1];
			}

			elsif ($row =~ /(^ACCESSION)(.*)/) { # the sequence data
				my @accession = split (/\s+/, $row);
				$accession = $accession[1];	
			}
			
			elsif ($row =~ /(ORGANISM)(.*)/) { # the common name and classification
			$common_name = $2;
				$common_name =~ s/^\s*//;
				$common_name =~ s/\s+/_/g;
				$common_name =~ s/\.//g;
				print "$common_name\n"; ## debuging 
				
				my $class = '';
				while (<IN>) {
					my $row = $_;
					chomp $row;
					if ($row =~ /\./) {
						$class .= $row;
						$class =~ s/[\s*|\.]//g;
						@classification = split /;/, $class;
						# print join ":", @class;
						# print "\n";
						last;	
					}
					$class .= $row;
				}
			}
			elsif ($row =~ /(^ORIGIN)(.*)/) { # the sequence data	
				while (<IN>) {
					my $row = $_;
					chomp $row;
					if ($row =~ /\/\//) {
						$seq =~ s/\s*//g;	
						# print "$seq\n";
						last
					};
					$row =~ /(^\s*\d*\s)(.*)/;
					$seq .= $2;			
				}
			}
		}
		print "read $i records\n";
	close IN;
	1;	#print Dumper($self);
}


=head1 accessors

=cut 


=head2 seq
=cut

sub seq { 
	my $self = shift;
	my $t = shift;
	if (@_) { $self->{$t}->{'seq'} = shift };
	return $self->{$t}->{'seq'};
}


=head2 classification

=cut


sub classification { 
	my $self = shift;
	my $t = shift;
	if (@_) { $self->{$t}->{classification} = [@_] };
	return $self->{$t}->{classification};
}


=head2 commonName

=cut

sub commonName { 
	my $self = shift;
	my $t = shift;
	if (@_) { $self->{$t}->{commonname} = shift };
	return $self->{$t}->{commonname};
}


=head2 locus 

=cut

sub locus { 
	my $self = shift;
	my $t = shift;
	if (@_) { $self->{$t}->{locus} = shift };
	return $self->{$t}->{locus};
}


=head1 functions

=cut 



=head2 found 

Determines (by ACCESSION) wether a sequence is present in the fastgenbank object

=cut

sub found { 
	my $self = shift;
	my $t = shift; # ACCESSION number
	$self->{$t} && return 1;
	return 0
}


=head2 seqObj 

Returns a pseudo - object for use in gbname (mimics the necessary bioperl seqobj structure)

=cut

sub seqObj { 
	my $self = shift;
	my $t = shift; ## a genbank display accessor
	my $obj = {};
	$obj->{species}->{_classification} = $self->classification($t) || die "$t doesn't have a classification";
	$obj->{primary_seq}->{display_id} = $t;
	$obj->{species}->{_common_name} = $self->commonName($t);

	bless $obj, "gbname"; # it just matters that its an object, not what kind
	return $obj;
}


=head2 loopTer 

Returns an array of terminals

=cut

sub loopTer {
	my $self = shift;
	return sort keys %{$self};
}


=head2 out_unBoundMx 

OVERRIDES fastmx

=cut

sub out_unBoundMx { 
	my $self = shift;
	my %params = @_;

	# 	print $params{'-outfile'};
	my $of = $params{'-outfile'} || "out.fas";


	my $gbname = Psy::Helpers::Gbname->new;

	open (OUT, ">$of" ) or die "can't open outputfile\n";
	
	my $i=0;
	
	# print join " ", $self->loopTer;
	
	foreach my $t ($self->loopTer) {	
		print OUT ">", $gbname->name('-in' => $self->seqObj($t)), "\n";
		print OUT $self->seq($t);
		print OUT "\n\n";
		print ".";
		$i++;	
	}
	close OUT;
	print "wrote: $i seqs\n";
	
	$gbname->summaryHigher;
	$gbname->summaryNoHigher;
	
1;
}


=head2 out_namesViaAccFile  

Reads a single column textfile containing genbank accessors and returns a formated

=cut

sub out_namesViaAccFile {
	my $self = shift;
	my %params = @_;
	my $of = $params{'-outfile'} || "formatted_names.txt";
	my $if = $params{'-infile'} || "names.txt";
	
	my $gbname = gbname->new;

	open (OUT, ">$of" ) or die "can't open outputfile\n";
	
	open (IN, "$if") || die "can't open infile";
		while (<IN>) {
			chomp;
			# print Dumper $self->commonName($_);
			if ($self->found($_)) {print OUT $gbname->name('-in' => $self->seqObj($_))}
			else {print OUT "$_ not found"};
			print OUT "\n";
		}
	close IN;

	close OUT;

	print "\n names dumped to: $of\n";
}


=head1 AUTHOR

'Matt, C<< <m{j}yoder@{tee}[aye](em)(you).domain4unis> >>

=head1 BUGS

Please report any bugs or feature requests to
C<bug-psy-helpers-fastgenbank@rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Psy>.
I will be notified, and then you'll automatically be notified of progress on
your bug as I make changes.

=head1 ACKNOWLEDGEMENTS

=head1 COPYRIGHT & LICENSE

Copyright 2005 'Matt, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1; # End of Psy::Helpers::Fastgenbank
