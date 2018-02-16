use Test::More 'no_plan';

use strict;
use warnings;
use Data::Dumper;	
# use diagnostics;

use lib qw(lib/);

use Psy::Psy;
use Psy::Output::Output;

use Psy::Output::Fasta; # gets around Psy::Output::Autoload

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
ok( $o->Fasta('-mx' => $mx), 'exported data in Fasta format');


ok( $o->Fasta('-mx' => $mx,'-file_name' => 'test2'), 'exported data in Fasta format');
ok( $o->Fasta('-mx' => $mx, '-file_name' => 'test3'), 'exported data in Fasta format');


