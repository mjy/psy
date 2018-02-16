#!/usr/local/bin/perl -w

# version 1.03 / original script by Chris Elsik /  Modified July 26/04 by M. Yoder

#usage: './phase_to_nexus.pl <phase .samples file (treefile)> <phase .output file> <name of nexus file to be made (don't include .nex)> <burnin>'
my $usage = "\nusage: perl phase_to_nexus.pl <phase .samples file (treefile)> <phase .output file> <name of nexus file to be made (do not include .nex)> <burnin>\n";

my $phase_tree_file = shift @ARGV or die $usage;
my $phase_output = shift @ARGV or die $usage;
my $nexus_output = shift @ARGV or die $usage;
my $burnin = shift @ARGV or die $usage;

open (OUT, ">$nexus_output.nex") or die $usage;

print OUT "#Nexus;\n";
print OUT "Begin trees;\n";

open (TREE, "$phase_tree_file") or die $usage;

my $tree_counter = 0;
while (<TREE>) {
    $tree_counter++;
    my $tree = $_;
    print OUT "tree PHASE_" . $tree_counter . " = [&U] " . $tree;
} 
close TREE;

print OUT "End;\n";


print OUT "contree $burnin-$tree_counter/ strict=no indi=yes majrule=yes le50=yes treefile\=phase_consensus.tre;\n";


print OUT "\n\n" . "[!\n";

open (DATA, "$phase_output") or die $!;

LINE: while (<DATA>) {
    if (m/^species/) {
	while (<DATA>) {
	    if (m/^invariants/) {
		print OUT $_;
		next LINE;
	    }
	}
    }
    else {
	print OUT $_;
    }
}

close DATA;

print OUT "\n]\n";

close OUT;

print "WARNING!  Your supplied burnin value of $burnin is greater than the number of trees in the treefile ($tree_counter).  The nexus file will not execute properly\n" if $burnin > $tree_counter;
print "Done.\n";

exit;

