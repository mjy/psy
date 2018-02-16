use Test::More ('tests' => 23);

use strict;
use warnings;

use Data::Dumper;	

use lib qw(lib/);

use Psy::Matrix::Slice;

ok(my $s = slice->new, 'new slice');
ok( $s->blocks(0..9), 'added some blocks');
ok( $s->taxa(0..9), 'added some taxa');

ok( $s->total('Blocks') == 10, 'have 10 blocks');

is( $s->containsBlk(1), 1, 'contains 1');
is( $s->containsBlk(50), 0, 'does not contain 50');
is( $s->containsBlk(-1), 0, 'does not contain -1');



ok( $s->remove('Blocks', (0..4)), 'removed 5');
ok( $s->total('Blocks') == 5, 'have 5 blocks');
is( $s->containsBlk(1), 0, 'does not contain 1');

ok( $s->final('Blocks') == 9, 'final');
ok( $s->first('Blocks') == 5, 'first');
ok( $s->lengthLastBlk == 1, 'lengthLastBlk');


ok( my $s1 = $s->clone, 'cloned a block');




# print Dumper($s1);

# and some tests on a matrix

use Psy::Psy;

my %params = ( 
	'-matrix_file' => 'matrix',
	'-matrix_label' => 'my_label',
	'-helix_index' => 'stem_index',
	'-path' => 't/in/t/',
	'-number_terminals' => 8,  
	'-project_name' => 'test'
);

my $foo = Psy->new(%params);
my $mx = $foo->mx('-matrix_label' => 'my_label');
my $s2 = $mx->origSlice;

ok(my $c = $s2->collapse('-mx' => $foo->mx('-matrix_label' => 'my_label')) );
my @r = keys %{$c};
ok($#r == 5);
ok( $s2->describeCollapsed($c) ); # not really an object method
ok( $s2->decribeCoded );

my $s3 = slice->new;
$s3->blocks(0, 2, 4, 6, 8, 10);
$s3->taxa(0..9);

is( $s3->excludedBlkStart('-mx' => $mx, '-blk_index' => 7), undef, 'not present'); # unbracketed only

is( $s3->excludedBlkStart('-mx' => $mx, '-blk_index' => 8), 11, 'excludedBlkStart'); # unbracketed only
is( $s3->excludedBlkEnd('-mx' => $mx, '-blk_index' => 8), 18, 'excludedBlkEnd');

is( $s3->unexcludedBlkStart('-mx' => $mx, '-blk_index' => 8), 14,  'unxcludedBlkStart' ); # includes bracketed
is( $s3->unexcludedBlkEnd('-mx' => $mx, '-blk_index' => 8), 21,  'unexcludedBlkEnd' );

