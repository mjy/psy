# package Psy::Helpers::Bayes2gelr;

=head1 Bayes2gelr - Bayes to Gelman's R

(Hopefully) as implemented in Miller et al (2004)
# Bayes to Gelman's R, as implemented in Miller et al. 2004 - code by Matt Yoder 2004
# 
# "when the within-sequence variance is much less than the between-sequence variance, the separate
# sequences are not sampling from the same distribution and have not converged"
#
# root R hat - is a measure of how much better (scale) we could do if we sampled to inifinity.  1 is no better (what's less than 1???), numbers > 1 means we could do much better with > n
#


=head1 VERSION

# ver 0.02
# 	- all_stats() function reports all stats for all parameters
# 	- unsure as to why R est values can go below 1??
#
# ver 0.01 
# 	- basic objects completed, trial runs finished and apparently working
# 	- this whole thing is simple stats though, probably much more efficiently implemented
# 	- MrBayes file object might be useful for further manipulations/merging etc.

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

use Psy::Helpers::Bayes2gelr;
  
my $foo = converge->new("d_151.p", "d_152.p");
$foo->all_stats();

print $foo->val_B('TL{all}') , "\n";
print $foo->val_W('TL{all}') , "\n";
print $foo->varhat('TL{all}') , "\n";
print $foo->rootRhat('pi(G){1}'), "\n";
print $foo->psi('Gen');
print Dumper($foo->{data});

=cut
package converge;

use Data::Dumper;

use strict;
use warnings;

sub new {
	my $type = shift;
	my @infiles = @_;
	my $self = {};

	die "too few input files" if $#infiles < 1;
	print "loading ", $#infiles +1, " files\n";
	my $i=0;
	my $generations; # test to see all files have the same number of generations
	
	foreach my $file (@infiles) {
		 $self->{data}->{$i} = bayesPfile->new($file);
		 
		 if ($i == 0) {
			$generations = $self->{data}->{$i}->stats('Gen','count');
		 }
		 
		 if ($i>1) { # after the first file is loaded need to check that each file has the same number of samples
			 if ($generations != $self->{data}->{$i}->stats('Gen','count') ) {
				die " number of generations not equal in input files! \n";
			 }
		 }
		 $i++;
	}

	print "total runs: " , $i, "\n";
	$self->{total_runs} = $i;
	print "number generations: ", $generations, "\n";
	$self->{number_generations} = $generations;
	
	$self->{run_time} = localtime($^T);
	
	bless $self, $type;
   	return $self;	
}

sub psi { 
	my $self = shift;
	my $param = shift;
	my $psi;
	foreach my $d (values %{ $self->{data}} ) {	# loop through each bayesPfile object
			$psi += $d->psiDotJ($param);	
	}

	$psi = ( (1 / $self->{total_runs}) * $psi);
	return $psi;
}

sub val_B { # between sample variation (simple statistics really) - zero for identical samples?
	my $self = shift;
	my $param = shift;
	
	my $val_B;

	my $psi = $self->psi($param); # don't keep re-calculating this

	foreach my $d (values %{ $self->{data}} ) {	# loop through each bayesPfile object
		$val_B += ( ($d->psiDotJ($param) - $psi) ) **2	
	}

	$val_B = ( $self->{number_generations} / ($self->{total_runs} - 1 ) ) * $val_B;
	
	return $val_B;
}

sub val_W { # between sample variation (simple statistics really) - zero for identical samples?
		my $self = shift;
		my $param = shift;

		my ($tmp, $psidotJ);
		my $val_W = 0;

		
		my $psi = $self->psi($param); # don't keep re-calculating this
		
		foreach my $d (values %{ $self->{data}} ) {	# loop through each bayesPfile object
			#print "$val_W ";
			$tmp = 0;
			$psidotJ = $d->psiDotJ($param);
			
			foreach my $val ($d->valsByParam($param)) {
				$tmp += ($val - $psidotJ)**2 ; 
			}
			
			$val_W += $tmp / ($self->{number_generations} - 1) ;
		}

		return ($val_W / $self->{total_runs});
}

sub varhat { # between sample variation (simple statistics really) - zero for identical samples?
		my $self = shift;
		my $param = shift;
		
		return ( ((($self->{number_generations} - 1)/$self->{number_generations})*$self->val_W($param)) + ($self->val_B($param) / $self->{number_generations})   );
}

sub rootRhat {
	my $self = shift;
	my $param = shift;
	my $vh =  $self->varhat($param);
	my $W =  $self->val_W($param);
	return 0 if $W==0;
	return sqrt( $vh / $W )
}

sub all_stats {
	my $self = shift;

	foreach my $var (values %{$self->{data}->{0}->{vars} }) {
			print "$var\n";
			
			printf ("\tB\t %6.8f\n",$self->val_B($var));
			printf ("\tW\t %6.8f\n", $self->val_W($var));
			printf ("\tvar est\t %6.8f\n",$self->varhat($var));
			printf ("\tR est\t %6.8f\n",$self->rootRhat($var));

			print  "\n"  ;
	}
}

1;

=head2 bayesPFile

A second package

=cut

package bayesPfile;

use Statistics::Descriptive;

use strict;
use warnings;

sub new {
	my $type = shift;
	my $self = {};
	$self->{infile} = shift;

	my @vars;
	my @vals;
	my ($i, $t) = (0) x 2;

	#$t=1;
	
	if (not defined $self->{infile}) { die " must provide a filename to load from ";}

	open (IN, $self->{infile} ) || die " can't open columns file ", $self->{infile},"\n";
	
		while (<IN>) {
			/^\[ID:\s(.+)]/ && do { 	# matches first line
				$self->{id} = $1;
				next
			};

			/^Gen/ && do {				# store the parameters 
				#$self->{vars}->{0} = "iter"; # index the generation 
				@vars = split;
				my $i=0;
				foreach my $var (@vars) {
					$self->{vars}->{$i} = $var; 
					$i++;
				}
				next
			};

			# if its not the first too lines, its a data line
			# assumes each parameter is sampled only once per generation (safe in MrBayes/PHASE case)	
			@vals = split;
			my $i=0;
			foreach my $val (@vals) {
				push @{$self->{vals}->{$self->{vars}->{$i}}} , $val; 
				$i++;
			}
			$t++;
		}

		die "No lines in input file!\n" if $t == 0;
		$self->{total_lines} = $t;
		print "  total lines: ", $self->{total_lines} , "\n";

		
	close (IN);
	
	bless $self, $type;
   	return $self;
}

sub psiDotJ { # returns Miller's 2004 psi.j
	my $self = shift;
	my $param = shift;
	
	return ( (1/$self->{total_lines}) * $self->stats($param, 'sum') );
}


sub valsByParam { # returns an array of vals
	my $self = shift;
	my $param = shift;
	
	return  @{ $self->{vals}->{$param} };
}

sub valByGenIndex {  # returns a single value for a given parameter for a single generation as indexed from 0
	my $self = shift;
	my ($gen, $param) = @_;
	
	die "generation too great" if $gen > $#{$self->{vals}->{$self->{vars}->{0} } } ; 
	return ${ @{$self->{vals}->{$param}} } [$gen];
}

sub valByGen {  # returns a single value for a given parameter for a single generation as indexed from column 0 in .p files (SLOW, shouldn't be used)
	my $self = shift;
	my ($gen, $param) = @_;
	my $k;
	print  join " ", @{$self->{vals}->{$param}};
	
	my $i= 0;
	foreach my $val (@{$self->{vals}->{Gen}}) {
		last if $val == $gen;
		$i++;
	}
	return ${ @{$self->{vals}->{$param}} } [$i];
}

sub stats() { # return a hash with some statistics for a given parameter name (not column)
	my $self = shift;
	my ($param, $stat) = @_;

 	die "parameter undefined\n" if not defined $param;
	die "not a valid parameter?\n" if not defined $self->{vals}->{$param};
	
	my $sta=Statistics::Descriptive::Full->new();	
	$sta->add_data( @{$self->{vals}->{$param}} );

	if (not defined $stat) { # if you don't ask for one, you get them all
		return {"sum" => $sta->sum(), "min" => $sta->min() , "nax" => $sta->max(),  "mean" => $sta->mean(), "standard_deviation" => $sta->standard_deviation, "variance" => $sta->variance };
	}
	elsif ($stat eq "sum") {
		return $sta->sum();
	}
	
	elsif ($stat eq "min") {
		return $sta->min();
	}
	
	elsif ($stat eq "max") {
		return $sta->max();
	}
	
	elsif ($stat eq "mean") {
		return $sta->mean();
	}
	elsif ($stat eq "count") { # might not work
		return $sta->count();
	}
	elsif ($stat eq "standard_deviation") { # might not work
		return $sta->standard_deviation();
	}
	else {
		print "->stat() returned no value!\n";
		return undef;
	}	
}

=head1 AUTHOR

'Matt, C<< <m{j}yoder@{tee}[aye](em)(you).domain4unis> >>

=head1 BUGS

Please report any bugs or feature requests to
C<bug-psy-helpers-bayes2gelr@rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Psy>.
I will be notified, and then you'll automatically be notified of progress on
your bug as I make changes.

=head1 ACKNOWLEDGEMENTS

=head1 COPYRIGHT & LICENSE

Copyright 2005 'Matt, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1; # End of Psy::Helpers::Bayes2gelr
