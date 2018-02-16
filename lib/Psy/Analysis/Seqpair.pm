package Psy::Analysis::Seqpair;

use warnings;
use strict;
use Statistics::Distributions;
use Data::Dumper;

use Psy::Dna::Alphabet;
use Psy::Strings::Strings;

use Psy::Psy;
our @ISA = qw(Psy);

use vars qw(@ISA @EXPORTER @EXPORT_OK $VERSION);


=head1 NAME

Psy::Analysis::Seqpair -

functions for 2 sequences (strings really) of same length

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

An object and and methods for use on two strings of equal length (generally RNA sequences).

This object does not reference terminals.

Alphabet is passed throughout as \@

    use Psy::Analysis::Seqpair;
	@seqs = ('ACGT', 'ACGT');
    my $sp = Psy::Analysis::Seqpair->new(\@seqs);
   	$sp->mi; 
	...

=head1 METHODS

=cut


=head2 new

=cut

sub new {
	my $type = shift;
	my @seqs = @_; # reference to array of seqs, column to return (starts at 1)
	my $self = {};
	
	bless $self, $type;

	$self->_init(@seqs);
	
	return $self; 
}


sub _init {
	my $self = shift;
	my @seqs = @_;

	$self->{legal_chars} = Psy::Dna::Alphabet->new;
	
	# translate non alphabet chars to ?
	my $i = 0;
	foreach my $seq (@seqs) {	
		$self->{seq}->{$i} =  $self->{legal_chars}->clean('-str' => $seq, '-replace' => '?');
		$i++;
	}
	
	(length $seqs[0] != length $seqs[1]) && die 'not equal length';
	
	$self->seqLen(length ($seqs[0]) );
	
	# store all basepairs for faster reference
	# foreach my $l ($self->{legal_chars}->alphabet) {
	# }
	
	for (my $i =0; $i < length $seqs[0] ; $i++ ) { 
			$self->{seq}->{0} =~ /.{$i}(.)/;
			$self->{pairs}->{$i} .= $1;
			$self->{nucs}->{1}->{$1}->{t}++;
			
			$self->{alpha}->{$1} = undef;
			# $self->{tnucs}->{$1}++;
			
			$self->{seq}->{1}  =~ /.{$i}(.)/;
			$self->{pairs}->{$i} .= $1;
			$self->{nucs}->{2}->{$1}->{t}++;

			$self->{alpha}->{$1} = undef;

			# $self->{tnucs}->{$1}++;
				
			$self->{uniquepairs}->{	$self->{pairs}->{$i} }->{t}++;
	}

	# convert to percentages too
	foreach my $s ('1', '2') {
		foreach my $k (keys %{$self->{nucs}->{$s}} ) {
			#	print " $k \n";
			$self->np($s, $k, ($self->nt($s, $k) / $self->seqLen));	
			# print $self->np($s, $k); 
		}
	}

	foreach my $pr ($self->up) {
		$self->pp($pr, $self->pt($pr) / $self->seqLen);
	}
}


=head2 seqLen 

length of sequences

=cut


sub seqLen { 
	my $self = shift;
	if (@_) { $self->{seqlen} = shift}
	return $self->{seqlen};
}

#	sub pi { # length of sequences
#		my $self = shift;
#		if (@_) { $self->{seqlen} = shift}
#		return $self->{seqlen};
#	}


=head2 obsAlpha 

Returns the array of observed letters in the two strings. ## T is converted to U??

=cut


sub obsAlpha {
	my $self = shift;
	return (sort keys %{$self->{alpha}});
}


=head2 bp 

Return the bp from a given position (starts at zero)

=cut


sub bp {
	my $self = shift;
	my $i = shift;
	return $self->{pairs}->{$i};
}


=head2 up 

Returns array of unique pairs observed in the two strings

=cut


sub up { 
	my $self = shift;
	return sort keys %{$self->{uniquepairs}};	
}



=head2 np 

?? returns the probability of nucleotide n for sequence (0/1)

	$foo->np(1,'A');


=cut


sub np { 
	my $self = shift;
	my $s = shift;
	my $n = shift;
	
	if (@_) { $self->{nucs}->{$s}->{$n}->{p} =  shift}
	return $self->{nucs}->{$s}->{$n}->{p};
}


=head2 nt 

Returns the total for a seq (0/1) for a given letter.

	$foo->nt(1,'A');


=cut

sub nt {
	my $self = shift;
	my $s = shift;
	my $n = shift;
	
	if (@_) { $self->{nucs}->{$s}->{$n}->{t} =  shift}
	return $self->{nucs}->{$s}->{$n}->{t};
}


=head2 pp

Returns probability of pair $p.

=cut


sub pp {
	my $self = shift;
	my $p = shift;
	if (@_) { $self->{uniquepairs}->{$p}->{p} = shift}
	return $self->{uniquepairs}->{$p}->{p} 
}



=head2 pt

returns count (total) of pair $p

=cut

sub pt {
	my $self = shift;
	my $p = shift;
	if (@_) { $self->{uniquepairs}->{$p}->{t} = shift}
	return $self->{uniquepairs}->{$p}->{t} || 0
}



sub _verify {
	my $self=shift;
	my %params = @_;
	foreach my $p (keys %params) { # check matrix parameters here	
	}
	return %params;
}


=head2 mi

Returns the mutal information (ln) for the object. !!! The mi should range from 0 (perfectly independent) to (2).

=cut


sub mi {
	my $self = shift;
	my $mi = 0;
	foreach my $pr ($self->up) {
		my $c1 = substr($pr, 0,1);
		my $c2 = substr($pr, 1,1);
					
		$mi +=  ( $self->pp($pr) / log(2) * log ( $self->pp($pr) / ( $self->np(1, $c1) * $self->np(2, $c2) )));
	}
	return $mi;
}


=head2 chi2VecTot

Given an array of nuc pairs, returns the total pairs (usable for both rows and cols - see chi2)

=cut

sub chi2VecTot { 
	my $self = shift;
	my @np = @_;
	my $t = 0;
	foreach my $np (@np) {
		$t += $self->pt($np)
	}
	return $t;
}


=head2 chi2

calculates the chi2 for basepairs

=cut


sub chi2 { 
	my $self = shift;
	my @alpha = $self->obsAlpha;

	my $chi2 = 0;
	my $t;	
	my $i = 1;
	foreach my $bp1 (@alpha) {
		my $j = 1;
		foreach my $bp2 (@alpha) {
			push @{$t->{row}->{$i}}, "$bp1$bp2";
			push @{$t->{col}->{$j}}, "$bp1$bp2";
			$j++;
		}
		$i++;
	}	
	
	$i = 1;
	foreach my $bp1 (@alpha) {
		my $j = 1;
		foreach my $bp2 (@alpha) {
			my $obs = $self->pt("$bp1$bp2");
			my $exp = ($self->chi2VecTot(@{$t->{col}->{$j}}) * $self->chi2VecTot(@{$t->{row}->{$i}}) / $self->seqLen); # row tot * col tot / tot-tot

			$exp == 0 && next;	
			#print "{ $j $bp1 $bp2 }\n";

			$chi2 += ( ($obs - $exp)**2) / $exp;
			$j++;
		}
		$i++;
	}	
	return $chi2;
}


=head2 cramersV 

Returns Cramers V, see planetmath.org for comments/code basis.
V ranges from 0-1, the closer to 0 the smaller the association of categorical variables, if X==Y then V==1.

=cut


sub cramersV { 
	my $self = shift;
	my $chi2 = $self->chi2;
	my @alpha = $self->obsAlpha;
	
	return ($chi2 / ($self->seqLen * $#alpha))**.5
}


=head2 contigencyTable

Print to stdout a contingency table and some other stats.

=cut


sub contigencyTable {
	my $self = shift;
	my @alpha = $self->obsAlpha;
	
	# $#alpha == 0 && die "both vectors identical";
			
	print "   ";
	print join "\t", @alpha;
	print "\n";

	my $t;
		
	my $i = 1;
	foreach my $bp1 (@alpha) {
		print "$bp1 |";
		my $j = 1;
		foreach my $bp2 (@alpha) {
			push @{$t->{row}->{$i}}, "$bp1$bp2";
			push @{$t->{col}->{$j}}, "$bp1$bp2";

			print $self->pt("$bp1$bp2");
			print "\t";
			$j++;
		}	
		print $self->chi2VecTot( @{$t->{row}->{$i}});
		$i++;
		print "\n";
	}	
	
	# last row of contingency table, note the total "observations" is == the sequence length!!

	print "   ";
	foreach my $j (sort keys %{$t->{col}}) {
		print $self->chi2VecTot( @{$t->{col}->{$j}}), "\t";
	}
	print $self->seqLen;
	print "\n\n";
	
	my $c2 = $self->chi2;

	if ($#alpha == 0) {print "no freedom!\n"; return 0} 
	my $chis = Statistics::Distributions::chisqrdistr( ($#alpha)**2, 0.01) ; # alpha indexes from zero, so its already -1, and its always symetrical, so square it

	print "\ndof: ", ($#alpha)**2, "; chi sqr: $c2; ";
	print "crit val: $chis\n";
	
	$c2 < $chis ? print "accept null (independant)\n" : print "reject null (dependant)\n";

	print "Cramer's V: ", $self->cramersV, "\n";
}

=head2 chi2Sig

Returns 0 if the chi2 value is not significant (independance), 1 if there is dependance.

=cut

sub chi2Sig {
	my $self=shift;
	my $c2 = $self->chi2;

	my @alpha = $self->obsAlpha;
	
	if ($#alpha == 0) { return 0} # no degrees of freedom 
	my $chis = Statistics::Distributions::chisqrdistr( ($#alpha)**2, 0.01) ; # alpha indexes from zero, so its already -1, and its always symetrical, so square it	
	$c2 < $chis ? 0 : 1;
}



1;

__END__

=head1 AUTHOR

'Matt, C<< <m{j}yoder@{tee}[aye](em)(you).domain4unis> >>

=head1 BUGS

Please report any bugs or feature requests to
C<bug-psy-analysis-seqpair@rt.cpan.org>, or through the web interface at
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

1; # End of Psy::Analysis::Seqpair
