use Test::More 'no_plan';

use strict;
use warnings;

use Data::Dumper;	

use lib qw(lib/);

# and use some test data
use Psy::Psy;
use Psy::Output::Output;
use Psy::Output::Nexus;

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
ok( $o->Nexus('-mx' => $mx, '-file_name' => 'vanilla.nex')); 
ok( $o->Nexus('-mx' => $mx, '-file_name' => 'spaced.nex', '-blk_spacer' => " ")); 
ok( $o->Nexus('-mx' => $mx, '-file_name' => "bracketed.nex", '-bracket_blocks' => 1, '-number_bracketed' => 0, '-blk_spacer' => " ", '-mode' => 'all', '-original_indexing' => 0)); 

