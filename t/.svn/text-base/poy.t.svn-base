use Test::More ('tests' => 10);

use strict;
use warnings;
use Data::Dumper;	
# use diagnostics;

use lib qw(lib/);

use Psy::Psy;
use Psy::Output::Output;

use Psy::Output::Poy; # gets around Psy::Output::Autoload

my %params = ( 
	'-matrix_file' => 'matrix',
	'-matrix_label' => 'my_label',
	'-helix_index' => 'stem_index',
	'-path' => 't/in/t/',
	'-number_terminals' => 8,  
	'-project_name' => 'test'
);


my $foo = Psy->new(%params), 'made a new Psy object';
my $mx = $foo->mx('-matrix_label' => 'my_label');

my $o = output->new();
 ok( $o->Poy('-mx' => $mx ) == 9);

# ok( $o->Poy('-mx' => $mx , '-mode' => 'do_all') == 9); # third block excluded for length 1
# ok( $o->Poy('-mx' => $mx , '-mode' => 'fso_all') == 11);
#ok( $o->Poy('-mx' => $mx , '-mode' => 'do_unbrack_fso_brack') == 10 );
#ok( $o->Poy('-mx' => $mx , '-mode' => 'do_brack_fso_unbrack') == 10);

#ok( $o->Poy('-mx' => $mx , '-mode' => 'fso_unbrack') == 7); # one block is excluded
#ok( $o->Poy('-mx' => $mx , '-mode' => 'do_unbrack') == 6 );
#ok( $o->Poy('-mx' => $mx , '-mode' => 'do_brack') == 3);
#ok( $o->Poy('-mx' => $mx , '-mode' => 'fso_brack') == 4 );

#ok( $o->Poy('-mx' => $mx , '-collapse' => 1) );
