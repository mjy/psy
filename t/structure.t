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

my $foo = Psy->new(%params);

my $mx = $foo->mx('-matrix_label' => 'my_label');
ok(my $s = $mx->origStructure);
is($s->label('-blk' => 6, '-label_name' => 'label1'),'REC2');
is($s->phase2Mask('-slice' => $mx->origSlice), '((.(.((..(....<<<<.<<<<.>>>>>>>>..((.()))..))))))............');
is_deeply([$s->loopStems('-slice' => $mx->origSlice)], [qw/0 3 5 7 8 9/]);
is_deeply([$s->loopHelices('-slice' => $mx->origSlice)], [qw/0 3 7/]);
is_deeply([$s->loopBasePairs('-mx' => $mx, '-blk' => 0, '-ter' => 0)], [qw/CG CG GC UU UU GC/]);
is_deeply([$s->pairPos('-mx' => $mx, '-blk' => 0, '-pos' => 1)], [qw/0 39/]);

is_deeply( [$s->pseudoKnotPairs(qw/0 100 1 99 10 88 11 99/)], [qw/0 1 10/], 'pseudoKnotPairs');

# print Dumper($s->loopPairs('-mx' => $mx, '-blk' => 0));
# print $s->loopStems('-slice' => $mx->origSlice);
is_deeply( $s->nucPairs('-mx' => $mx)->{'pairs'}, [qw/1:40 2:39 4:38 6:35 7:34 10:33 11:28 12:27 13:26 14:25 16:24 17:23 18:22 19:21 29:43 30:42 32:41/], 'nucPairs, pairs'); 
is_deeply( $s->nucPairs('-mx' => $mx, '-zero_beers' => 1)->{'pairs'}, [qw/0:39 1:38 3:37 5:34 6:33 9:32 10:27 11:26 12:25 13:24 15:23 16:22 17:21 18:20 28:42 29:41 31:40/] , 'nucPairs, pairs, zero'); 
is_deeply( $s->nucPairs('-mx' => $mx, '-zero_beers' => 1)->{'nucs'}, [qw/0 39 1 38 3 37 5 34 6 33 9 32 10 27 11 26 12 25 13 24 15 23 16 22 17 21 18 20 28 42 29 41 31 40/], 'nucPairs, nucs, zero'); 

is_deeply( $s->nucPairs('-mx' => $mx, '-zero_beers' => 1)->{'nucs_ns'}, [qw/2 4 7 8 14 19 30 35 36/]);
is_deeply( $s->nucPairs('-mx' => $mx)->{'nucs_ns'}, [qw/3 5 8 9 15 20 31 36 37/], 'nupairs nucs_ns');

# test with a short slice

my $slice = slice->new;

$slice->blocks(0..7);
$slice->taxa(0..8);

is_deeply( $s->nucPairs('-mx' => $mx, '-slice' => $slice, '-zero_beers' => 1)->{'nucs_ns'}, [qw/0 1 2 3 4 5 6 7 8 9 14 19 28 29 30 31/]);
is_deeply( $s->nucPairs('-mx' => $mx, '-slice' => $slice, '-zero_beers' => 1)->{'nucs'}, [qw/10 27 11 26 12 25 13 24 15 23 16 22 17 21 18 20/]);
is_deeply( $s->nucPairs('-mx' => $mx, '-slice' => $slice, '-zero_beers' => 1)->{'pairs'}, [qw/10:27 11:26 12:25 13:24 15:23 16:22 17:21 18:20/], 'short slice pairs');

# non contiguous slices

my $slice2 = slice->new;
$slice2->blocks(0,1,2,3,4,6,7,8);
$slice2->taxa(0..8);
is_deeply( $s->nucPairs('-mx' => $mx, '-slice' => $slice2, '-zero_beers' => 1)->{'nucs_ns'}, [qw/2 4 7 8  10 11 12 13 14 15 16 17 18   19   20 21 22 23    27 28/], 'non contig nucs_ns');
is_deeply( $s->nucPairs('-mx' => $mx, '-slice' => $slice2, '-zero_beers' => 1)->{'nucs'}, [qw/0 31 1 30 3 29 5 26 6 25 9 24/]);
is_deeply( $s->nucPairs('-mx' => $mx, '-slice' => $slice2, '-zero_beers' => 1)->{'pairs'}, [qw/0:31 1:30 3:29 5:26 6:25 9:24/], 'short slice pairs');

ok($s->nucPairs('-mx' => $mx, '-slice' => $slice2));



# print Dumper($s->nucPairs('-mx' => $mx, '-slice' => $slice2));

#print Dumper($s->nucPairs('-mx' => $mx));
# $mx->_debugBlockRanges;
# print Dumper($s);


