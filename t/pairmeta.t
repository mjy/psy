use Test::More ('tests' => 4);

use strict;
use warnings;
use Data::Dumper;	
# use diagnostics;

use lib qw(lib/);

use Psy::Psy;
use Psy::Output::Output;

use Psy::Output::Pairmeta; # gets around Psy::Output::Autoload

my %params = ( 
	'-matrix_file' => 'matrix',
	'-matrix_label' => 'my_label',
	'-helix_index' => 'stem_index',
	'-path' => 't/in/t/',
	'-number_terminals' => 8,  
	'-project_name' => 'test'
);

#my %params2 = ( 
#	'-matrix_file' => 'phase_simulated.rna',
#	'-matrix_label' => 'my_label',
#	'-path' => 't/in/t/',
#	'-fileformat' => 'PHASE_simulate',
#	'-project_name' => 'test',
#	'-struct_size' => 2,
#);

ok(my $foo = Psy->new(%params), 'made a new Psy object');

my $mx = $foo->mx('-matrix_label' => 'my_label');

my $o = output->new();

ok( $o->Pairmeta('-mx' => $mx )); # defaults to '-out_format' => 'table'

ok( $o->Pairmeta('-mx' => $mx, '-out_format' => 'nexus_weights' ));
ok( $o->Pairmeta('-mx' => $mx, '-out_format' => 'ccode' ));


