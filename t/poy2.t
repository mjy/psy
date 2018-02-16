use Test::More ('tests' => 1);

use strict;
use warnings;
use Data::Dumper;	
# use diagnostics;

use lib qw(lib/);

use Psy::Psy;
use Psy::Output::Output;

use Psy::Output::Poy; # gets around Psy::Output::Autoload

my %params = ( 
	'-matrix_file' => 'Encyrtidae_28S.T.08.Nex',
	'-matrix_label' => 'my_label',
	'-helix_index' => 'Encyrtidae_Index.T.05.txt',
	'-path' => 't/in/ency2/',
	'-number_terminals' => 134,  
	'-project_name' => 'ency2'
);

my $foo = Psy->new(%params);
my $mx = $foo->mx('-matrix_label' => 'my_label');

my $o = output->new();
ok($o->Poy('-mx' => $mx, '-collapse' => 1));



