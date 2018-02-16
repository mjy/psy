use Test::Simple 'no_plan';

use strict;
use warnings;
use lib ('lib');

use Psy::Helpers::Fastoneline;

ok (my $t = Psy::Helpers::Fastoneline->new('-infile' => 't/in/t/oneline_matrix')) ;

# my @r = qw/aaactccaTCTAAGGCT AaATATGACCa CGAGACCGATA GCGAACAAGTAC CGTGAGGGAAA gttgaaaag [GA]ACTTTGAAGAGAG AGTTCAAG AGTACGTGA/;

# paired : H15 GTGAAACCGTTCAGGG
# unpaired (before D2) : GUAAA

my @d2 = qw/CTTTGAAGAGAG
			C[CT]T[GC]AGA 
			GGGA[GT]A 
			GGGG
			TGCACT
			TCTCCCCT
			CGTCGCGAC
			CCGT[TG]GGG
			GGGCCG
			GCGACGCT 
			/;

my @d3d5 = qw/GAAGGTC
				GCGTACACG
				CTGGCACTC
				TCTCATCTG
				TGGGTGAGA
				GAGTGCCAAG
			 /;

my @rna18S = qw/AAATCcTTTaaCGAGG/; 

# ok ($t->out_splitRegex( '-outfile' => 't/out/splitRegex.out', '-regexes' => \@r) );

ok ($t->out_splitRegexAligned('-outfile' => 't/out/splitRegexAligned.out', '-regexes' => \@d3d5) );


