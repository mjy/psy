use Test::More 'no_plan';

use strict;
use warnings;

use Data::Dumper;	

use lib qw(lib/);

# and use some test data

use Psy::Psy;

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
ok($mx->blk(0)->bracketed == 0);

ok($mx->blk(2)->bracketed == 1);
is($mx->colInBlk('-col' => 31), 7);
is_deeply($mx->column('-col' => 0),  'CAGGGGGGG');
is_deeply($mx->column('-col' => 27), 'CAGGGGGGG');
is_deeply($mx->column('-col' => 38), 'GGGGAAGAA');
is($mx->ter(0)->label, 'aus',  '->ter');

ok($mx->blk2Fastmatrix, '->blk2Fastmatrix called');

ok($mx->alignBlock('-blk' => 10), 'you need "clustalw" in your path for this test to succeed ');
ok($mx->alignBlock('-blk' => 10, '-align_method' => 'muscle'), 'you need "muscle" in your path for this test to succeed ');


# not setup/working
my %params2 = ( 
	'-matrix_file' => 'phase_simulated.rna',
	'-matrix_label' => 'my_label',
	'-path' => 't/in/t/',
	'-fileformat' => 'PHASE_simulate',
	'-project_name' => 'test',
	'-struct_size' => 100,
);


my %params3 = ( 
	'-matrix_file' => 'interleave_single_block',
	'-matrix_label' => 'my_label3',
	'-path' => 't/in/t/',
	'-number_terminals' => 13,  
	'-project_name' => 'test2',
);

ok(my $foo2 = Psy->new(%params3), 'made a new object');
ok(my $mx2 = $foo2->mx('-matrix_label' => 'my_label3'));

is($mx2->ter(0)->label, 'foo',  '->ter');
is($mx2->totalChars, 256 ,  '->totalChars');
is($mx2->totalBlocks, 4 ,  '->totalBlocks');

# ok(my $bar = Psy->new(%params2), 'made a new object');
# ok(my $mx1 = $bar->mx('-matrix_label' => 'my_label'));

#print $mx->column('-col' => 38);
# print Dumper($mx1->origSlice);
# my $s = $mx->origSlice;
# $s->remove('Blocks', 0..5);
# print Dumper($mx->origSlice);

