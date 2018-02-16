#!/usr/local/bin/perl -w

# version 0.1
#
# doesn't really convert tree formats, just appends necessary info for PAUP 

#usage: './phase_to_nexus.pl <phase .samples file (treefile)> <phase .output file> <name of nexus file to be made (don't include .nex)> <burnin>'
my $usage = "\nusage: perl phylip2nexus.pl <phylip trees (tree lines only)> <name of nexus file to be made (do not include .nex)>\n";

my $input = shift @ARGV or die $usage;
my $output = shift @ARGV or die $usage;

open (OUT, ">$output.nex") or die $usage;

print OUT "#Nexus;\n";
print OUT "Begin trees;\n";

open (TREE, "$input") or die $usage;

my $tree_counter = 0;
while (<TREE>) {
    $tree_counter++;
    my $tree = $_;
    print OUT "tree converted_" . $tree_counter . " = [&U] " . $tree;
} 
close TREE;

print OUT "End;\n";


close OUT;

print "Done!\n";

exit;

