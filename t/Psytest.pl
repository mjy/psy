
=head1 Psytest.pl

# uncomment the pointer to the pertinent lib

=cut


use warnings;
#use String::Random;

use lib ('../lib');
use lib ('../lib/Psy');

use Psy::Psy;
# use Psy::Output::Poy;
use Template;

use Psy::Output::Output;
use Data::Dumper;

use Psy::Strings::Strings;
use Psy::Dna::Alphabet;

# use diagnostics;  

# $Psy::debug = 1;

my $Psy = Psy->new(	'-matrix_file' => 'Hymenoptera_18S_0.01.Nex',
	'-matrix_label' => 'my_label',
	'-helix_index' => 'Hyme_18S_index_0.00.txt',	#	'-helix_index' => 'Evaniodea.index.0.02.txt',
	'-path' => 'in/hym18s/', 
	# '-path' => 'D:/lib/JRNA/data/hym18S/',			 	# '-path' => '../../data/evan/',
	'-number_terminals' => 285,  # n-1
	'-project_name' => 'hym18S'
);

my $mx =  $Psy->mx('-matrix_label' => 'my_label');

my $slice = slice->new;  
$slice->blocks(0..291);
$slice->taxa(0..285);

$slice->describeCollapsed($slice->collapse('-mx' => $mx));


#$slice->taxa(0..15, 17..36, 40..82, 86..121, 123..139, 141..151, 153..168, 170..181, 183..185, 187..235, 237..239, 242..246);	# taxa with 28S data

# print join " ", $mx->sumAlignedPosByBlk;

# $mx->debug_blockLengths;

# $slice->pruneUninformativeTaxa('-mx' => $mx, '-cutoff' => .5);

# my %foo = $mx->fragmentHash('-blk' => 4);
# print Dumper %foo;

my $out = output->new('-project_name' => 'psy_hym18S');

# use helpers::fasthmmer;

# my $h = fasthmmer->new('-infile' => '/lib/Psy/test/in/hym18S/out_all.sto');

# $out->blkmapped('-mx' => $mx, '-slice' => $slice, '-fastmatrix' => $h );

$out->Poy('-mx' => $mx, '-slice' => $slice);

# $out->tnt('-mx' => $mx, '-slice' => $slice);
# $out->stockholm('-mx' => $mx);
# $out->blockmetanex('-mx' => $mx, '-slice' => $slice);

# $out->fasta('-mx' => $mx, '-slice' => $slice, '-strip_gaps' => 1);

# print $mx->{structures}->{original}->phase2Mask('-slice' => $slice);

#print "\n";
# print Dumper($foo);

#my $var1 = '..................................................((((((((((.(((((((((((((....))
#)))..)))))))).......))))))))))..<<<<<<<AAAAA...AAAA...............aa.....aaaaaaa
#..AA.......aaAAAAAAAAAAAAAA......aaaaaaaaaaaaaa>>>>>>>...<<<............AA..AA..
#AAA..aaAAAA.........aaaa...aa.........aaa..>>>((((((...((((((((((())))))........
#.)))))..)))))).............................';

#print "\n$var1"; 

## TRY SHORTER LINE LENGTHS IN PHASE 1.1

# you fill up %dispatch and just say $dispatch{$foo}->()

# build up hash of default arg values
#    my %default_args = ( foo => 'bar', baz => 'quux' );

# now join the two, giving precedence to %raw_args
#    my %args = ( %default_args, %raw_args );

# output ultimately $Psy->output->fasta (default all the way through)

#print $Psy->mx('-matrix_label' => 'my_label')->lengthLongestTer;

# $mx->removeBlock('-blk' => 211);
# print $mx->cell('-blk' => 211, '-ter' => 0);

#print $mx->lastBlockIndex;
#$mx->_debugBlockRanges;
#$mx->deleteBlock('-blk' => 186);
#$mx->_calcBlockPositions; 
#$mx->_debugBlockRanges;

#print $mx->terLabel2Index('-label' => 'DERV005');

#print Dumper($mx->{blocks}->{0});

#print $mx->{_block_accession};
#print keys %{$mx->{blocks}};
#print Dumper($mx->{blocks}->{211});
#print "FOO ", $mx->{blocks}->{211}->literalEnd;
#print $mx->column('-col' => 0);
#print Dumper($mx->{structures});
#print "column in block: ", $mx->colInBlk('-col' => 0, '-mode' => 'includebracketed'), "\n";
# $mx->{alphabet}->add('-?NBCDG');
#print join " ", $mx->{alphabet}->possibleBPs;

#print join " " ,$mx->loopBlocks;

#print Dumper( $mx->{structures}->{'original'}->nucPairs('-mx' => $mx) );

#$out->nexus('-mx' => "FOO");

#print "foo: $foo\n";;
#print strings::uniqueChars(join "", $mx->sliceData);

#my $bar = "Template";

#my $tt = $bar->new();

#$tt->process('template_test', &output::mxData('-mx' => $mx), \*STDOUT) || die $tt->error();; 

# print join " ", $mx->{structures}->{'original'}->loopHelices('-slice'=>$mx->orig_slice);

# $mx->generateSafenames;
# print Dumper ($mx->{structures}->{original}->alignedHelices('-blk'=> 42, '-mx' => $mx, '-ter' => 55) );

#print $mx->rawRowData(-blk=>0, -ter=>0, -gapchar=>'AUAA|-');
#print join " ", $mx->loopTaxa;

 # >mx('-matrix_label' => 'my_label')->{structures}
#print Dumper( $Psy->mx('-matrix_label' => 'my_label')->total_chars('-slicecharmode'=> 'all') );

#&output::matrix( 
#			'-mx' => $Psy->mx('-matrix_label' => 'my_label')
#		);
			
# helix_slice;
#print Dumper($Psy->{structures});

# =cut
 
__END__
