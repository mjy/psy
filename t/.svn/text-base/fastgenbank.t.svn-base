use Test::Simple tests => 5;

use strict;
use warnings;

use lib('lib');

use Psy::Helpers::Fastgenbank;

ok(my $o = Psy::Helpers::Fastgenbank->new('-infile' => 't/in/t/sequences.gb' ), 'made new, loaded file');

my @foo = $o->loopTer;
ok( $#foo > 0, ' read some records');

ok($o->out_unBoundMx('-outfile' => 't/out/out.fas'), 'call out_unBoundMx');
ok($o->found('DQ022280') == 1, 'found a sequence');
ok($o->found('FOOBALISHOUS') == 0, 'did not find a sequence');
