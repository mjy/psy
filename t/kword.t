use Test::More 'no_plan';

use strict;
use warnings;

use Data::Dumper;	

use lib qw(lib/);

# and use some test data

use Psy::Psy;
use Psy::Output::Kword;

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

# REMEMBER: By default kword only translates bracketed blocks.  To translate everything use '-translation_mode' => 'all'

ok( $o->Kword('-mx' => $mx, '-format' => 'tnt' )); # '-kword_size' => [1..4] 
ok( $o->Kword('-mx' => $mx, '-out_format' => 'tnt'));  #  '-kword_size' => [1..4] 
ok( $o->Kword('-mx' => $mx, '-out_format' => 'mesquite', '-file_name' => 'test'));  #  '-kword_size' => [1..4] 
ok( $o->Kword('-mx' => $mx, '-out_format' => 'tnt', '-file_name' => 'test2', 'translation_mode' => 'all_fused'   ));  #  '-kword_size' => [1..4] 
ok( $o->Kword('-mx' => $mx, '-out_format' => 'tnt',  '-file_name' => 'test3', '-translation_mode' => 'all' ));  #  '-kword_size' => [1..4] 
ok( $o->Kword('-mx' => $mx, '-out_format' => 'mrbayes',  '-file_name' => 'test4_mb', '-translation_mode' => 'all' ));  #  '-kword_size' => [1..4] 
ok( $o->Kword('-mx' => $mx, '-out_format' => 'tnt',  '-file_name' => 'test5', '-translation_mode' => 'all', '-zero_length_missing' => 1 ));  #  '-kword_size' => [1..4] 


