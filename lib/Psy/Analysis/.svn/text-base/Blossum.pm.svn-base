package Psy::Analysis::Blossum;

use warnings;
use strict;
use Data::Dumper;

=head1 NAME

Psy::Analysis::Blossum 

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

uses the package column below, switch to columns object

Perhaps a little code snippet.

    use Psy::Analysis::Blossum;

    my $foo = Psy::Analysis::Blossum->new();
    ...

=cut

=head1 EXPORT

=cut

#	use Exporter;
#	use vars qw(@ISA @EXPORTER @EXPORT_OK $VERSION);
#	@EXPORT_OK = qw(version rootpath projectname starttime);
	
=head1 OBJECT
	
=cut

=head2 new

=cut

	sub new {
		my $type = shift;
		my $rows = shift; # reference to an array of rows

		my $len = length (@{$rows}[0]);
		foreach my $row (@{$rows}) {
			die "input data not a matrix in package blossum\n" if length ( $row )  != $len;
		}
	
		my $self = {};
		
		my $k = 0;
		my ($l1, $l2);
		my $total_rows =  $#{$rows} + 1;
		
		# get data
		for (my $r = 0; $r < $len; $r++) {
			my $cur_col = column->new($r, $rows);
		
			for (my $i=0; $i< $total_rows; $i++) {

				$self->{position}->{$r}->{alphabet}->{ ${$cur_col->{rows}}[$i] }->{total} += 1;

				$self->{alphabet}->{${$cur_col->{rows}}[$i]} = defined; # record the letter to the alphabet master set

				for (my $j=$i+1; $j< $total_rows; $j++) {
					($l1, $l2) = sort (${$cur_col->{rows}}[$i], ${$cur_col->{rows}}[$j]);
					
					$self->{nucs}->{$l1}->{total} += .5; 
					$self->{nucs}->{$l2}->{total} += .5;

					$self->{pairs}->{$l1.$l2}->{total} += 1;
					$k++;
				}
			}
		}

		# calculate (need to force to count as reals, not text - somethings not quite right (apparently with the numbers on bioneq
		# based on http://apps.bioneq.qc.ca/twiki/bin/view/Knowledgebase/BLOSUM
		
		foreach my $r (keys %{$self->{position}} ) {
			foreach my $alpha (keys %{$self->{alphabet}} ) {
				$self->{position}->{$r}->{alphabet}->{$alpha}->{q_ij} = $self->{position}->{$r}->{alphabet}->{ $alpha }->{total} / $total_rows; # q alphabet by position from data
			}
		}

		foreach my $key (keys %{$self->{nucs}} ) {
			$self->{nucs}->{$key}->{percent} = $self->{nucs}->{$key}->{total} / $k; # p_ij  
			$self->{alignment_probability} = $self->{alignment_probability} * $self->{nucs}->{$key}->{percent};
		}
		
		foreach my  $key (keys %{$self->{pairs}} ) {
			$self->{pairs}->{$key}->{percent} = $self->{pairs}->{$key}->{total} / $k; # q_ij ( q paired matrices calculated from observed frequencies)
			
			($l1, $l2) = split //, $key;

			if ($l1 eq $l2) {
				$self->{pairs}->{$key}->{pr_occurance} = ( $self->{nucs}->{$l1}->{percent} * $self->{nucs}->{$l2}->{percent} ); 	# e_ij
			}
			else {
				$self->{pairs}->{$key}->{pr_occurance} = 2 * $self->{nucs}->{$l1}->{percent} * $self->{nucs}->{$l2}->{percent};  # e_ij
			}

			$self->{pairs}->{$key}->{odds} = ( $self->{pairs}->{$key}->{percent} / $self->{pairs}->{$key}->{pr_occurance} );		# q_ij / e_ij
			$self->{pairs}->{$key}->{log_odds} = log ($self->{pairs}->{$key}->{odds}) / log(2); 
			$self->{pairs}->{$key}->{blossum_odds} = int ($self->{pairs}->{$key}->{log_odds} * 2);	
		}	
		bless $self, $type;
		return $self; 
	}


=head2 display_matrix

=cut


	sub display_matrix { # not working
		my $self = shift;
		my $type = shift;
	
		foreach my $l1 (keys %{$self->{nucs}} ) {
			foreach my $l2 (keys %{$self->{nucs}} ) {

				($l1, $l2) = sort ($l1, $l2);
				
				print "[ "; # $l1 $l2
				print $self->{pairs}->{$l1.$l2}->{$type} if defined $self->{pairs}->{$l1.$l2}->{$type};
				print "]";
				print "\t";
			}

			print "\n";
		}
	}

1;



package column; # depreciate to columns ultimately

		sub new {
			my $type = shift;
			my $self = {};
			my ($pos, $seqs) = @_; # reference to array of seqs, column to return (starts at 1)

			#my $self=>{rows} = [];
			
			#my $i=0;
			foreach my $seq (@{$seqs}) {
				$seq =~ /.{$pos}(.)/;
				#print $1;
				push @{$self->{rows}} , $1;
				#$i++;
			}
			
			bless $self, $type;
			return $self; 
		}		
1;


__END__


=head1 AUTHOR

'Matt, C<< <m{j}yoder@{tee}[aye](em)(you).domain4unis> >>

=head1 BUGS

Please report any bugs or feature requests to
C<bug-psy-analysis-blossum@rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Psy>.
I will be notified, and then you'll automatically be notified of progress on
your bug as I make changes.

=head1 ACKNOWLEDGEMENTS

=head1 COPYRIGHT & LICENSE

Copyright 2005 'Matt, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1; # End of Psy::Analysis::Blossum
