use Test::More 'no_plan';

use strict;
use warnings;

use Data::Dumper;	

use lib qw(lib/);

# and use some test data

use Psy::Psy;
use Psy::Output::KwordProfile;

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

ok( $o->KwordProfile('-mx' => $mx,  '-mode' => 'all_fused', '-kword_size' => [1..4]    ));    

ok( $o->KwordProfile('-mx' => $mx,  '-mode' => 'all_fused'    ));    



