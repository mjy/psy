use Test::More 'no_plan';

use strict;
use warnings;

use Data::Dumper;	

use lib qw(lib/);

# and use some test data

use Psy::Helpers::Crawler;
use Psy::Output::Kword;

my %params = ( 
	'-root_folder' => 't/in/t/crawler',
	'-output_folder' => 't/out/crawler',
);

ok(my $c = Crawler->new(%params), 'made a new crawler');

print Dumper($c);
is_deeply([$c->datasets], [qw/001 002 003/], 'detect folders');
is_deeply([$c->datafiles('-folder' => '001')], [qw/data_bar.nex data_foo.nex/], 'return files');
is_deeply($c->treefiles('-folder' => '001'), [qw/tree.tre/], 'return trees');

$c->run('-convert', '-analyze');

is($c->return_analysis(1), '0_data_bar.nex_1.tnt');
# is_deeply([$c->analyses], [qw/0_data_bar.nex.tnt 1_data_foo.nex.tnt 2_data_bar.nex.tnt 3_data_foo.nex.tnt/], 'return analyses');

# print $c->analysis_loop;
# print "\n\n\n";
print Dumper( $c->analysis_params(5));

# print Dumper($c->{'analysis_loop'});
# print Dumper($c->{'result'});

# $c->run('-analyze');

