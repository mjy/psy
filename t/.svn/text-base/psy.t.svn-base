use Test::More ('tests' => 2);

use strict;
use warnings;
use Data::Dumper;	

use lib qw(lib/);
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

