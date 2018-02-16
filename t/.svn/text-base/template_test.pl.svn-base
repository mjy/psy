use strict;
use warnings; 
use Template;

my ($foo, $bar);
my $tt = Template->new();
$foo->{'blorf'}->{florf} = [0..4];
$bar= "BLKBJLFG";
$tt->process('t.txt', $foo) || die $tt->ERROR, "\n";
