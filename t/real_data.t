use Test::More ('tests' => 1);

use strict;
use warnings;
use Data::Dumper;	
# use diagnostics;

use lib qw(lib/);

use Psy::Psy;
use Psy::Output::Output;
use Psy::Output::Tnt;
use Psy::Output::Blockmeta;
use Psy::Output::Html; # gets around Psy::Output::Autoload
use Psy::Output::Kword;

#my %params = ( 
#	'-matrix_file' => 'Death_Star_3.11.Nex',
#	'-matrix_label' => 'my_label',
#	'-helix_index' => 'Ichs2_Index.3.11.txt',
#	'-path' => 't/in/ichs/',
#	'-number_terminals' => 289,  
#	'-project_name' => 'ichs2'
#);


my %params2 = ( 
	'-matrix_file' => 'Diapriidae_28S_D2_0_00.Nex',
	'-matrix_label' => 'my_label',
	'-helix_index' => 'Diapriidae_HI_0_00.txt',
	'-path' => 't/in/diap/',
	'-number_terminals' => 152,  
	'-project_name' => 'diap'
);

my %params3 = ( 
	'-matrix_file' => 'Chalcidoidea.28S.1.02.NEX',
	'-matrix_label' => 'my_label',
	'-helix_index' => 'Chalcidoidea.Block.Index.1.01.txt',
	'-path' => 't/in/chal/',
	'-number_terminals' => 526,  
	'-project_name' => 'chal'
);

my $foo = Psy->new(%params2);
my $mx = $foo->mx('-matrix_label' => 'my_label');


# ok( $o->Html('-mx' => $mx ) ); 

# PATH ISSUES HERE!
# ok($mx->alignSliceBlocks, 'Clustal test (need clustalW in path, see docs)');
my $s = $mx->origSlice;

my @outgroups = (8, 35, 54, 57, 112, 120, 130, 131, 139..152); # 57 - proctotrupidae
my @good_outgroups = (120, 112, 130, 131);
my @weirdness = (80);
$s->remove('Taxa', @outgroups, @weirdness);
$s->taxa(@good_outgroups);


$s->pruneUninformativeTaxa('-mx' => $mx);
my $o = output->new();

my $p;
for my $blk ($mx->origSlice->loop('Blocks')) {
#	if ($mx->blk($blk)->bracketed == 1) {
		$p->{$blk}->{blks} = ([$blk]);
		$p->{$blk}->{kword} = 1;
		$p->{$blk}->{type} = 'trans'
#	}
}

delete $p->{85};


# nuke 85

 $o->Blockmeta('-mx' => $mx, '-slice' => $s, '-out_format' => 'tnt',  '-plan' => $p); # ,
# $o->Kword('-mx' => $mx, '-plan' => $p, '-kword_size' => [3..300], '-slice' => $s,  '-plan' => $p );
