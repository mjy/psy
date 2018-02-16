use Test::Simple tests => 5;

# this test requires network/ncbi access

use strict;
use warnings;
use Data::Dumper;	

use lib('lib');

use Psy::Helpers::Gbname;

my $s = 'U39957';

ok(my $foo = new Psy::Helpers::Gbname, 'made a new object');

ok(my $seqObj = $foo->seqObjLookup($s), 'retrieved a genbank seq object');

ok($foo->higher('-seq_obj' => $seqObj) eq 'Diapriidae', 'got default higher name');
ok($foo->name('-in' => $seqObj, '-species_label' => 'none') eq 'SSU39957_DIAPRIIDAE', 'got default higher name');

my %bar = $foo->seqFamilyGroupNames($seqObj);
ok($bar{'family'} eq 'Diapriidae' , 'got family name');

# U39957 CTGATCAAAGATTAAAAATTAATTAAATAGGTGCAATAAATTAATTGTACAAAGGTAGGATAATAATTTGGTTTTTTATTAAAATCTTGAATGAATGATTAAATGAAATAAAAACTGTCTCAATTTAATTAAATTTAAAAATTTTTTTTTGGGGTAAAAAAAAAAATTTTTTTTAAAAGGCGAGAAGACCCTATAGGATTTTATATAATAAAAATAATTATATTATTAAAATTAAATTATAATAAATTATATTTAATTGGGGTAATTTAAAAAATAAAAAAATTTTTTTTTAATATTTACATAAATAAATGAAATAAATTAATAAATAAATAATATTTAAATAAAATAAATTACCTTAGGGATAACAGCGTTATTTTTTTTAGATAGATCATATTGATACAAAAGATTGCGACCTCGATGTTGAATTA
