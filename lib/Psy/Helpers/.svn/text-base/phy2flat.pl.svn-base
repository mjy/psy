
# phylip 2 flatfile
# version 0.01 - Matt Yoder
# converts interleaved phylip files to one big interleave

use strict;
use warnings;

my $usage = " died\nusage: perl phy2flat.pl <infile> <outfile>\n";
my $infile = shift @ARGV or die "$usage";
my $outfile = shift @ARGV or die "$usage";

my @row;
my @mx;
my $row;
my $curtaxa = 1;

my ($numtaxa, $numchars);

	open ( IN, $infile ) || die print "cannot open $infile \n $usage ";
		my $i = 0;
		my $havetaxa = 0;
		
		while (<IN>) {
			$row = $_;
			chomp $row;
			$row =~ s/^(\s+)//; # strip space from begining of the row
			
			if ($i == 0) {
				@row = split /\s/, $row;
				$numtaxa = $row[0]; $numchars = $row[1];
				print "numtaxa $numtaxa\nnumchars $numchars\n";
			}
			
			elsif ($_ !~ /(^(\[)|^(#=)|^(\s*#=)|^(\s+\[)|^([\n\r])|^(\s+[\n\r]))/) {
				if ($havetaxa == 0) {
					@row = split /(\s+)/, $row;
					$mx[$i][0] = shift @row;
					$mx[$i][1] .= join "", @row;
				}
				
				else {
					$row =~ s/(\s+)//g;
					$mx[$curtaxa][1] .= $row;
					$curtaxa++;
				}
			}
			
			if ($curtaxa > $numtaxa) {$curtaxa = 1};
			if ($i > $numtaxa) {$havetaxa = 1} ;
			
			$i++;
			#print "$i\n";
		}
		print "total lines: $i\n";
	close (IN);

	open (OUT, ">$outfile") or die "can't open file to output to $usage";
		for (my $j=1; $j<$numtaxa+1;$j++) {
			print OUT "$mx[$j][0]\t$mx[$j][1]\n";
		}

	close (OUT);

	1;

exit;
	
