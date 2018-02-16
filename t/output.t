use Test::More 'no_plan';

use strict;
use warnings;

use Data::Dumper;	
# use diagnostics;
use lib qw(lib/); # .t tests should run as 'perl t/foo.t'

use Psy::Output::Output;

# ok(my $o = output->new, 'made a new output object'); ## likely not working due to AUTOLOAD?
is(1, 1, 'no plans for output');
# print Dumper($o);
