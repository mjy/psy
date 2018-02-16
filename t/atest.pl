

use Data::Dumper;
use strict;
use warnings;

my $foo = 'ACTATGA_CA_TGGACTTG_1232132_a123123';
$foo =~ /(.+)_(.*)\z/;
print $2;

