use Test::More qw(no_plan);

use strict;
use warnings;
use Data::Dumper;
# use diagnostics;

use lib('lib');

use Psy::Analysis::Seqpair; #  qw(delFixedPositions);

# some test data
my @foo = qw(ACGT TGCA);

ok(my $o = Psy::Analysis::Seqpair->new(@foo));
ok($o->seqLen == 4);
is( $o->bp(1), 'CG', 'bp');
is_deeply([$o->obsAlpha], [qw/A C G U/], '->obsAlpha');
is_deeply([$o->up], [qw/AU CG GC UA/], '->up');
is($o->np(1, 'A'), .25, '->np');
is($o->nt(1, 'C'), 1, '->nt');
is($o->pp('AC'), undef, '->pp');
is($o->pp('AU'), .25, '->pp');


