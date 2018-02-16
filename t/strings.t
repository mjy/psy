use Test::More ('no_plan');

use strict;
use warnings;

use lib('lib');

use Psy::Strings::Strings; #  qw(delFixedPositions);

# some test data
my @foo = qw(-----A ---A-C --C--G G----G);

my @bar  = delFixedPositions('-strings' => \@foo, '-char' => '-');

my @blorf = (1,2,94,5,12..34);

ok($bar[0] eq '---A');
ok($bar[1] eq '--AC');
ok($bar[2] eq '-C-G');
ok($bar[3] eq 'G--G');

# different test data
my @blk = qw(ACCA AC - AC AC-C);
my $expanded = 'A--CCA-A---C---AC-AC-C';
my $r= gappedRe('-strings' => \@blk);
$expanded =~ m/$r/i;
ok($1 eq 'A--CCA-');
ok($3 eq '-');
ok($5 eq 'AC-C');

my $s = 'ABCDEFGHIJKLMNOPQRS123';
my @r = qw/ABC HIJ PQRS 1 23/;

is(&regexSplit('-str' => $s, '-regexes' => \@r) , '  ABC  DEFG  HIJ  KLMNO  PQRS  1  23  ');

@r = qw/ABC HIJ PQRS SomethingNotFoundHere 1 23/;
is(&regexSplit('-str' => $s, '-regexes' => \@r) , '  ABC  DEFG  HIJ  KLMNO  PQRS  1  23  ');


is(&arrayAsRange(@blorf), '(1..2, 5, 12..34, 94)');

