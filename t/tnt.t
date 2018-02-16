use Test::More 'no_plan';

use strict;
use warnings;

use Data::Dumper;	

use lib qw(lib/);

# and use some test data

use Psy::Psy;
use Psy::Output::Tnt;

my %params = ( 
	'-matrix_file' => 'matrix',
	'-matrix_label' => 'my_label',
	'-helix_index' => 'stem_index',
	'-path' => 't/in/t/',
	'-number_terminals' => 8,  
	'-project_name' => 'test'
);

ok(my $foo = Psy->new(%params), 'made a new object');
ok(my $mx = $foo->mx('-matrix_label' => 'my_label'));


my $o = output->new();
ok( $o->Tnt('-mx' => $mx)); 
ok( $o->Tnt('-mx' => $mx, '-mode' => 'unbracketed', '-file_name' => 'unbracketed'));
ok( $o->Tnt('-mx' => $mx, '-mode' => 'bracketed', '-file_name' => 'bracketed')); 
ok( $o->Tnt('-mx' => $mx, '-mode' => 'all', '-file_name' => 'all')); 



 $o->Tnt('-mx' => $mx,  '-mode' => 'all', '-file_name' => 'unaligned');
 $mx->alignSliceBlocks('-align_method' => 'muscle');
 $o->Tnt('-mx' => $mx, '-mode' => 'all', '-file_name' => 'aligned');
