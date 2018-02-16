use Test::More ('tests' => 7);

use strict;
use warnings;
use Data::Dumper;	
# use diagnostics;

use lib qw(lib/);

use Psy::Psy;
use Psy::Matrix::Column;

use Psy::Output::Poy; # gets around Psy::Output::Autoload

my %params = ( 
	'-matrix_file' => 'matrix',
	'-matrix_label' => 'my_label',
	'-helix_index' => 'stem_index',
	'-path' => 't/in/t/',
	'-number_terminals' => 8,  
	'-project_name' => 'test'
);

my $foo = Psy->new(%params);
my $mx = $foo->mx('-matrix_label' => 'my_label');

ok (my $c = column->new($mx->blockData('-blk' => 0)) );
ok ($c->totalSeqs == 9, 'totalSeqs');
ok ( %{$c->returnColHash(1)}->{'C'} == 4, 'returnColHash');

ok ($c->returnColVec(1) eq 'CCCCUUUUU', 'returnColVec');
ok ($c->pluralConcensus('-column' => 1) eq 'CU', 'pluralConcensus');
is ($c->primerRegex, '[ACG]C[AG][AC]GG', 'primerRegex');
is($c->primerRegex('-reverse_comp' => 1) ,'CC[GT][CT]G[CGT]' , 'primerRegex reverse complimented');
