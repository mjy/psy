use Test::More 'no_plan';

use strict;
use warnings;

use Data::Dumper;	

use lib qw(lib/);

# and use some test data

use Psy::Psy;
use Psy::Output::Kword;
use Psy::Graphics::Plot;

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
my $g = Plot->new();

ok( $o->Kword('-mx' => $mx, '-out_format' => 'bare', '-file_name' => 'plot', '-translation_mode' => 'all_fused', '-write_word_list' => 'word_list' ) );  #  '-kword_size' => [1..4] 



#use Cwd;
#my $dir = cwd;
#print $dir;

ok( $g->binaryBitmap('-infile' => 'analyses/kword/plot', '-outfile' => 'foo.png'));
ok( $g->binaryBitmap('-infile' => 'analyses/kword/plot', '-outfile' => 'foo2.png', '-max_plot_len' => 10));
ok( $g->binaryBitmap('-infile' => 'analyses/kword/plot', '-outfile' => 'foo3.png', '-max_plot_len' => 10, '-y_offset' => 2));
ok( $g->binaryBitmap('-infile' => 'analyses/kword/plot', '-outfile' => 'foo4.png', '-wordfile' => 'analyses/kword/word_list'));
