
# phase22tracer.pl

# summarizes phase2 mcmcphase output to a file readable by tracer 1.2

# version 0.01  
# free software Copyright Matt Yoder 2005

use strict;
use warnings;

my $usage = "\nusage: 'perl phase_to_tracer.pl <phase filename root> <output file name> <step size>\n"; #; [burnin]

my $root_filename = shift @ARGV or die $usage, $!;
my $tracer_output = shift @ARGV or die $usage, $!;
# my $bunin = shift @ARGV or die $usage, $!;

my $step_size;
if (@ARGV) {
		$step_size = shift
	}
	else {
		$step_size = 100; print " WARNING- arbitrarily setting step size to 100\n";
	}

open (MP, "$root_filename.mp") or die $!;
open (FO, ">$tracer_output") or die $!;

open (PLOT, "$root_filename.plt") or die $!;
	my @plot = <PLOT>;
close PLOT;

open (BL, "$root_filename.bl") or die $!;
my @bl = <BL>;
close BL;

my $counter = 0;

while (<MP>) {
	my @row = split (/\s+/, $_);
	#print @row;
	shift @row unless $#row == 0;
	
	if ($counter == 0) {
		print FO "Gen\tlnLk\ttree_prior\tmdl_prior\tTL"; # case matters, keep 'Gen' 
		my $i = 0;
		foreach my $col (@row) {
			print FO "\tparam$i";
			$i++
		}
		print FO "\n";
	}
		
	print FO $counter*$step_size, "\t";

	my $totBL = 0;
	my @blvals = split (/\s+/, $bl[$counter]);
	shift @blvals; # funny formating of PHASE .bl
	map ($totBL += $_ , @blvals);
	
	my @bar =	split (/\s+/, $plot[$counter]);
	$counter == 0 and shift @bar; ## temporary fix of misformated? .plt files
	print FO "$bar[0]\t$bar[1]\t$bar[2]\t"; # log likelihoods are negative in Phase2 
	
	print FO "$totBL\t";
	
	print FO join "\t", @row;
	print FO "\n";

	$counter++;
}

close MP;
close FO;
