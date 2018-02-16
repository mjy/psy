use Test::More 'no_plan';

use strict;
use warnings;

use Data::Dumper;	

use lib qw(lib/);

# and use some test data

use Psy::Psy;
use Psy::Output::Kword;

my %params = ( 
	'-matrix_file' => 'kword_test',
	'-matrix_label' => 'my_label',
	'-helix_index' => 'stem_index2',
	'-path' => 't/in/t/',
	'-number_terminals' => 3,  
	'-project_name' => 'test'
);


ok(my $foo = Psy->new(%params), 'made a new object');
ok(my $mx = $foo->mx('-matrix_label' => 'my_label'));

my $o = output->new();

ok( $o->Kword('-mx' => $mx, '-out_format' => 'mesquite', '-file_name' => 'kword_test_eg' ));  #



