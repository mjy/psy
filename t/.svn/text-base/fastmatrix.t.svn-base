use Test::More tests => 11;

use strict;
use warnings;
use Data::Dumper;

use lib ('lib');

use Psy::Helpers::Fastmatrix;

my $t = Psy::Helpers::Fastmatrix->new;

# some test data
$t->{'mx'}->{'seq_1'} = 'ACGTACGTACGT';
$t->{'mx'}->{'seq_2'} = '--------ACGT';
$t->{'mx'}->{'seq_3'} = 'ACxXACGT----';
$t->{'mx'}->{'U39957'} = 'CTGATCAAAGATTAAAAATTAATTAAATAGGTGCAATAAATTAATTGTACAAAGGTAGGATAATAATTTGGTTTTTTATTAAAATCTTGAATGAATGATTAAATGAAATAAAAACTGTCTCAATTTAATTAAATTTAAAAATTTTTTTTTGGGGTAAAAAAAAAAATTTTTTTTAAAAGGCGAGAAGACCCTATAGGATTTTATATAATAAAAATAATTATATTATTAAAATTAAATTATAATAAATTATATTTAATTGGGGTAATTTAAAAAATAAAAAAATTTTTTTTTAATATTTACATAAATAAATGAAATAAATTAATAAATAAATAATATTTAAATAAAATAAATTACCTTAGGGATAACAGCGTTATTTTTTTTAGATAGATCATATTGATACAAAAGATTGCGACCTCGATGTTGAATTA';

is( ($t->loopTer)[0] , 'U39957', '->loopTer');
is( $t->seq('seq_1') , 'ACGTACGTACGT', '->seq ');
ok( $t->seqLen('seq_1') == 12, , '->seqLen');
ok( $t->seqLen('seq_2') == 4, '->seqLen');
ok( $t->seqLen('seq_3') == 6, '->seqLen');

ok( ($t->loopBySeqLength)[0] eq 'U39957', '->loopBySeqLength longest'); 
ok( ($t->loopBySeqLength)[2] eq 'seq_3', '->loopBySeqLength shortest'); 
ok( ($t->loopByGBnameThenSeqLen)[0] eq 'seq_1', '->loopByGbnameThenSeqLen longest'); 
ok( ($t->loopByGBnameThenSeqLen)[3] eq 'U39957', '->loopByGbnameThenSeqLen shortest'); 

print Dumper($t);
ok ($t->clustalAlign, '->clustalAlign');
print Dumper($t);
ok ($t->muscleAlign, '->muscleAlign');
print Dumper($t);




