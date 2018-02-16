use Test::More 'no_plan';

use strict;
use warnings;

use Data::Dumper;	

use lib qw(lib/);

# and use some test data

use Psy::Psy;
use Psy::Output::MrBayes;

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
ok($o->MrBayes('-mx'=> $mx));

my $s = slice->new; # has no helices
$s->blocks(0,1,2,3,4);
$s->taxa(0..8);
ok($o->MrBayes('-mx'=> $mx, '-slice' => $s, '-file_name' => 'mrbayes_short'));


my $s2 = slice->new; # has one helix
$s2->blocks(0..6);
$s2->taxa(0..8);
ok($o->MrBayes('-mx'=> $mx, '-slice' => $s2, '-file_name' => 'mrbayes_less_short'));

