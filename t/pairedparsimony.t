use Test::More ('tests' => 3);

use strict;
use warnings;
use Data::Dumper;	
#  use diagnostics;

use lib qw(lib/);

use Psy::Psy;
use Psy::Output::Output;
use Psy::Output::Pairedparsimony; # gets around Psy::Output::Autoload

my %params = ( 
	'-matrix_file' => 'matrix',
	'-matrix_label' => 'my_label',
	'-helix_index' => 'stem_index',
	'-path' => 't/in/t/',
	'-number_terminals' => 8,  
	'-project_name' => 'test'
);

ok(my $foo = Psy->new(%params), 'made a new Psy object');
my $mx = $foo->mx('-matrix_label' => 'my_label');

my $o = output->new();
ok($o->Pairedparsimony('-mx' => $mx));
ok($o->Pairedparsimony('-mx' => $mx, '-out_format' => 'tnt', '-file_name' => 'tnt.tnt'));
