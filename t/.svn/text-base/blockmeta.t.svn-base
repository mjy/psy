use Test::More ('tests' => 2);

use strict;
use warnings;
use Data::Dumper;	
# use diagnostics;

use lib qw(lib/);

use Psy::Psy;
use Psy::Output::Output;

use Psy::Output::Blockmeta; # gets around Psy::Output::Autoload

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

#ok( $o->Blockmeta('-mx' => $mx ));

ok( $o->Blockmeta('-mx' => $mx, '-file_name' => 'original_arc', '-out_format' => 'original_arc' ));
#ok( $o->Blockmeta('-mx' => $mx, '-file_name' => 'mesq_chr_by_blk', '-out_format' => 'mesquite_chr_by_blk' ));
#ok( $o->Blockmeta('-mx' => $mx, '-file_name' => 'mesq_blk_by_chr', '-out_format' => 'mesquite_blk_by_chr' ));
