use Test::More 'no_plan';

use strict;
use warnings;

use lib('lib');
use Data::Dumper;

use Psy::Strings::Composition; #  qw(delFixedPositions);
use Psy::Strings::Strings;

my @foo = qw(-----A ---A-C --C--G G----G);

ok(my $c = Psy::Strings::Composition->new('-seqs' => \@foo));
is($c->seq, '-----A---A-C--C--GG----G');
is_deeply( [$c->ltrs], [qw/- A C G/] );
is_deeply( [$c->ranked], [qw/A C G -/] );

is_deeply( [$c->ranked], [qw/A C G -/] );
is($c->inLen, 24 );
is($c->count('A'), 2);

is($c->pct('A'), 0.0833333333333333);
is_deeply( [$c->least], [qw/A C/]);
is_deeply( [$c->most], [qw/-/]);

my $s = 'AA---CT-AA--A-GCGT';
my $c1 = Psy::Strings::Composition->new('-seqs' => ([$s]));	

is($c1->pctPair('-pair' => 'AA'), 0.12);
is($c1->pairDivChar('-pair' => 'AA'), '0.40');
is($c1->ltrDist('-ltr' => 'A'), '1.60');

is($c1->sumPct('A'), 0.277777777777778);
is($c1->sumPct('A', 'U'), 0.388888888888889);
is($c1->sumPct('G', 'C'), 0.222222222222222);
is($c1->sumPct('G', 'C', 'A'), 0.500000000000000);

