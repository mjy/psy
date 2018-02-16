
# please read the __README included with this distribution prior to using this code

# jrna - scripts copyright Matt Yoder, 2004, gnu public license applies

package Jrna;
use vars qw($VERSION);
$VERSION = '0.777';

use strict;
use warnings;


use FindBin qw($Bin); 		# sets $Bin to the root directory

srand ();

####################################
# Jrna package variables (globals) #
####################################

my $TIME = localtime($^T);

my $rootdir = $Bin;

#ultimately in stem descriptor obj?
my @bonds;	# array of compensatory bond changes for stem regions bonds[non-prime block number][basepairindex][statistic], not to be confused with @bond!!, used in displaying stems on the webpage output
my %sindex;	# should be a member of a stems descriptor object
my %stems;	# key is a block index that points to the prime end
			
# mtrx obj?
my @mx;			# $mx[terminal][block][0] -> the actuall data in the matrix
my %ters;		# terminal hash #->name

# mtrx meta data obj?
my @deschead;			# description headers $deschead[block][0] - first row; [1] - second row
my @nbd = ();
my @nbd1 = (); 			# [0] - block contents;	[1] - column start; [2] - column end (needs to be eliminated at some point)
my $defaultcolstart = 20;	# the first column length (taxa+whitespace), used various places, should be nuked

# i/o
my @fileroot = ("blkd", "mx", "sl", "stm", "bpf"); 	# path and filename roots for the various interleave output in out_mx
my %out_modes = (			# not implemented - future global "formating" reference
			"raw", 0,			
			"nexus", 1, 
			"phylip", 2, 
			"fasta", 3,	# not implemented
			"tnt", 4,
			"mrbayes", 5,
			"phase", 6,
			"poy", 7,
			"inaase", 8,
			"css", 9,
			"xml", 10	#not implemented
		);

# alphabet obj?
	     #0   1    2    3    4    5     6    7
my @rna  = qw(A U G C N - ? o);  	#(index order is important-> should be a hash likely) - double check this on reporting - is N out of place?

my @iupac = qw(A C T G U R Y S W K M B D H V N \. \-);  # used as a regex to search against, thus \. \- (merge to IUPAC object) 

my %basepairs = ( #  hash of all possible base-pair combinations for aucg?-no, where o is other
	"aa", 0,  "au", 1,  "ag", 2,  "ac", 3,  "an", 4,  "a-", 5,  "a?", 6,  "ao", 7,
	"ua", 8,  "uu", 9,  "ug", 10, "uc", 11, "un", 12, "u-", 13, "u?", 14, "uo", 15,
	"ga", 16, "gu", 17, "gg", 18, "gc", 19, "gn", 20, "g-", 21, "g?", 22, "go", 23,
	"ca", 24, "cu", 25, "cg", 26, "cc", 27, "cn", 28, "c-", 29, "c?", 30, "co", 31,
	"na", 32, "nu", 33, "ng", 34, "nc", 35, "nn", 36, "n-", 37, "n?", 38, "no", 39,
	"-a", 40, "-u", 41, "-g", 42, "-c", 43, "-n", 44, "--", 45, "-?", 46, "-o", 47,
	"?a", 48, "?u", 49, "?g", 50, "?c", 51, "?n", 52, "?-", 53, "??", 54, "?o", 55,
	"oa", 56, "ou", 57, "og", 58, "oc", 59, "on", 60, "o-", 61, "o?", 62, "oo", 63
);

# colum obj ultimately
my @ca;

# misc globals to clean up
my ($cur_block, $stemtot, $numpairsinstems, $totalstemlength ) = (0) x 4; 

my $orig_interleave = interleave->new;	# an interleave object for the original interleave -> sink to mx object ultimately

# ---------------------------------------------------------------------------------------------

sub calc_colstats () { 	# column stats algorithms
	my $slice = shift;

	$slice->describe;
	
	my (@ccol); # reset the descriptive array
	my ($a, $c, $g, $n, $u, $dash, $question, $other) = (0) x 8;
	
	# build the columns (remove to seperate fn)
	foreach my $x ($slice->loop("Blocks")) {
		if ($nbd[$x][3] == 0) {							# only examine bracketed blocks
			for (my $p=0; $p<$nbd[$x][0]; $p++) {  				# loop through all positions in the string
				foreach my $y ($slice->loop("Taxa")) {			# loop through the taxa - remember smaller subsets may have all dashes, while larger don't
			 		$mx[$y][$x][0] =~ /.{$p}(.)/; 			# get the character a position $p			
					$ccol[$nbd[$x][6]+$p] .= $1;			# concat
				}	
			}
		}
	}

	foreach my $blk ($slice->loop("Blocks")) {					# loop through the blocks matrix
		if ($nbd[$blk][3]==0) {									# only examine bracketed blocks
			for (my $i = $nbd[$blk][6]; $i < $nbd[$blk][7]+1; $i++) { 		# loop through all positions in the string
				$a = ($ccol[$i] =~ tr/A//);
				$c = ($ccol[$i] =~ tr/C//);
				$g = ($ccol[$i] =~ tr/G//);
				$u = ($ccol[$i] =~ tr/T|U//);
				$n = ($ccol[$i] =~ tr/N//);
				$dash = ($ccol[$i] =~ tr/-//);
				$question = ($ccol[$i] =~ tr/\?//);
				$other = $slice->total("Taxa") - ($a+$c+$g+$n+$u+$dash+$question); # $t
		
				#print " $other ";
				#$sum = $other + $a+$c+$g+$t+$n+$u+$dash+$question; 	
		
				#store values in a matrix for later use;
				$ca[$i][0] = $a;
				$ca[$i][1] = $u;
				$ca[$i][2] = $g;
				$ca[$i][3] = $c;
				$ca[$i][4] = $n;
				$ca[$i][5] = $dash;
				$ca[$i][6] = $question;
				$ca[$i][7] = $other;

				for (my $z=0; $z<8; $z++) { # zero out undefs for stats purposes, there is probably a better way to do this
					$ca[$i][$z] = 0 unless defined $ca[$i][$z];
				}
				
				#DEBUG
				# print "$blk : $i : $ca[$i][0] $ca[$i][1] $ca[$i][2] $ca[$i][3] $ca[$i][4] $ca[$i][5] $ca[$i][6] $ca[$i][7]  \n "; 
			}
		}
	}
}

sub calc_stem_bp_stats() {
	my (
		$slice,
		$summary_mode	# 0-include all taxa; 1-exclude "??" from calculations
	) = @_;

	#dependancies:
	# prior call to stems_initialize

	my ($stem1, $stem2, $mask, $curstem, $bpindex);
	my $a = 0;	
	my @wrkstems;
	my @tmp;

	@bonds = (); # scope ok, reset to at each call

	my $total_taxa = $slice->total("Taxa");
	#bonds[non-prime block number][columposition][basepairindex][statistic], {statsitics: 0-total; 1- ?} note BP index is NOT the same as columnum!!!! has to be rebuilt because of stem expansion

	foreach my $blk ($slice->loop("Blocks")) {
		if (exists $stems{$blk}) {
			push @wrkstems, $blk;
		}
	}

	print "working stems: || @wrkstems |\n";
	#could probably build $bonds as %bonds, but I'm not smart enought to do hashofhashofhash
	
	foreach $curstem (@wrkstems) { # loop through block indicies

	$bpindex = 0; # recalculated index

		#loop through taxa
		foreach my $t ($slice->loop("Taxa")) {
		
			#get the raw data
			@tmp = &str_aligned_stems($mx[$t][$curstem][0], $mx[$t][$stems{$curstem}][0], $nbd[$curstem][8], $nbd[$stems{$curstem}][8]);
			$stem1 = $tmp[0];
			$stem2 = $tmp[1];
			$mask =  $tmp[2];
			
			#print "| $stem1, $stem2, $mask |\n";
			
			$bpindex = 0; #recalculated index
			for (my $p = 0; $p < length $mask; $p ++ ) { # loop through positions
		 					
				#figure out the BP match
				$mask =~ /.{$p}(.)/;
			
				if ($1 eq '.') { # mark only bonded areas 
				}
				
				else {
					$stem1 =~ /.{$p}(.)/;
					$a = $1;	
					$stem2 =~ /.{$p}(.)/;
					$a .= $1;	
		
					$a = lc ($a); # perhaps not needed
		
					# exchange non aucg-?n characters with 'o' ########
					$a =~ s/(?)([^aucg\-\?n])/o/gi;	#[check me]	
			
					# increment totals based on index
					# print $a;
					# print $basepairs{$a};

					$bonds[$curstem][$bpindex][$basepairs{$a}][0]++; # note that not all values are initialized!
					$bpindex++;
				}	
			}	
			#print "bpindex: $bpindex\n";
		}
	}

	# build percentage and other stats here; note: (#bonds[non-prime block number][columposition][basepairindex][statistic])
	# calculate relative percentages of each base pair combination (all data mode)

	foreach $curstem (@wrkstems) { #loop through block indicies	
		for (my $p=0; $p < $#{$bonds[$curstem]} + 1; $p++) { #loop through positions ------- $bonds[$blk][ END OF THIS INDEX] --------------
			foreach $a (keys %basepairs) {
				#DEBUG
				my $qq_tot = 0;
				if (defined $bonds[$curstem][$p][$basepairs{$a}][0]) { # only calculate % for those with > 0
					if ($summary_mode == 1) {
						$qq_tot = $bonds[$curstem][$p][$basepairs{'??'}][0] unless not defined $bonds[$curstem][$p][$basepairs{'??'}][0];
						$bonds[$curstem][$p][$basepairs{$a}][1] = $bonds[$curstem][$p][$basepairs{$a}][0] / ($total_taxa - $qq_tot);  
					}
					else {
						$bonds[$curstem][$p][$basepairs{$a}][1] = $bonds[$curstem][$p][$basepairs{$a}][0] / ($total_taxa);  
					}
				}
				else {
					$bonds[$curstem][$p][$basepairs{$a}][1] = 0;
				}
			}
		}
	}

	# calculate covariation for each position/basepair
	my $cvpercent = .03;
	my $bpcomp1;
	my $bpcomp2; #left and right strings of existing bps
	
	my ($p1, $p2, $tmpcomp1, $tmpcomp2);
	
	# for each stem 
	foreach $curstem (@wrkstems) { #loop through block indicies	
	#build the composition hash
		for (my $p=0; $p < $#{$bonds[$curstem]} +1; $p++) { #plus 1 right here

			$bpcomp1 = "";
			$bpcomp2 = "";

			#build "composition" strings
			foreach $a (keys %basepairs) {
				if ($bonds[$curstem][$p][$basepairs{$a}][1] >= $cvpercent) {
					$a =~ /(.)(.)/;
					$bpcomp1 .= $1;
					$bpcomp2 .= $2;
				}
			}
		
			# check each basepair combination (existing only?) versus bphash where presence > $cvpercent
			foreach $a (keys %basepairs) {
				$tmpcomp1 = $bpcomp1;
				$tmpcomp2 = $bpcomp2;
				if ($bonds[$curstem][$p][$basepairs{$a}][1] > $cvpercent) { # check this row
					
					#print "$a\n";
					$a =~ /(.)(.)/;

					my $re3 = quotemeta($1);
					my $re4 = quotemeta($2);
				
					#print "A: (|$tmpcomp1|\t$c1\t$re3)\t\t(|$tmpcomp2|\t$c2\t$re4)\n";
								
					$tmpcomp1 =~ s/$re3//gi;
				       	$tmpcomp2 =~ s/$re4//gi;

					$tmpcomp1 =~ s/\?|N|-|o//gi;
					$tmpcomp2 =~ s/\?|N|-|o//gi;
					
					#print "B: (|$tmpcomp1|\t$c1\t$re3)\t\t(|$tmpcomp2|\t$c2\t$re4)\n\n";
					
					if ((length $tmpcomp1 > 0) and (length $tmpcomp2 > 0) ) {
						#print "$curstem $a $p covarying!\n";
						$bonds[$curstem][$p][$basepairs{$a}][2] = 1;
					}
				}
			}
		}
	}
	
	# debug
	#	foreach $curstem (@wrkstems) {
	#		print "| $curstem: $#{$bonds[$curstem]} \n";	
	#			for (my $i=0; $i<9; $i++) {	#print the first ith positions - will need to be changed to last index of $i, which = length
	#				for (my $j=0; $j<24; $j++) {
	#					print "$curstem";
	#				print " $j: ";
	#				if (defined $bonds[$curstem][$i][$j][0]) { print "$bonds[$curstem][$i][$j][0] "} 
	#			        else { print  "undef "};
	#			      
	#		}
	#		print "\n";
	#		}
	#	}	
	
}

sub slice_calc_mx_totals () { # slice method ultimately?
	# need to divide slice/non slice methods
	my $slice = shift;

	# returns this hash
	my %tot = ( # available totals 
		"numbrakchars", 0,			# number of bracketed characters
		"numbrakblks", 0,			# number of bracketed blocks
		"inaaseblks", 0,			# number of Inaase blocks -? what's this for?
		"longestbrakblk", 0,		# length longest bracketed block
		"totchars", 0,				# total columns
		"numunbrakchars", 0			# number non-bracketed characters 
	);
	
	foreach my $x ($slice->loop("Blocks")) {
		if ($nbd[$x][3] == 1) { # note this doesn't include the and $nbd[$b][0] < $imaxblklen  
			$tot{"numbrakchars"} += $nbd[$x][0];
			$tot{"numbrakblks"}++;
			$tot{"inaaseblks"} = $tot{"inaaseblks"}."$nbd[$x][1]-$nbd[$x][2] "; # the complete string of inaase data
			if ($tot{"longestbrakblk"} < $nbd[$x][0]) {$tot{"longestbrakblk"} = $nbd[$x][0] }
		}

		else {
			$tot{"numunbrakchars"} += $nbd[$x][0];
		}	
	}
	$tot{"totchars"} = $tot{"numunbrakchars"} + $tot{"numbrakchars"};
	return %tot;


sub slice_calc_stem_totals () {
	my $slice = shift;
	my %tot = (
		"bdas", 0,			#bdas = blocks defined as stems
		"length_bdas", 0,
		"num_non_prime_bdas", 0,
		"num_prime_bdas", 0,
		"num_complete_pairs", 0,
		"num_complete_helices", 0
	);

	foreach my $blk ($slice->loop("Blocks")) {
		if (grep ($_ == $blk, %stems)) {
			$tot{"bdas"}++;
			$tot{"length_bdas"} += $nbd[$blk][0];
		}

		if (grep ($_ == $blk, keys %stems)) {
			$tot{"num_prime_bdas"}++;
		
			if (grep ($_ == $stems{$blk}, values %stems)) {
				$tot{"num_complete_pairs"} += ($nbd[$blk][8] =~ tr/\(//);	
				$tot{"num_complete_helices"}++;
			}	
		}

		if (grep ($_ == $blk, values %stems)) {
			$tot{"num_non_prime_bdas"}++;				
		}
	}	
	return %tot;
}

sub mx_colcompgraph () { # returns an array corresponding to percentage of a u g c n - ? other and a graph/logo of said composition
	# dependencies: a prior call to &calc_colstats with containing block
	my (
		$slice,
		$col,		# column number
		$length		# graph length (in pixels/characters)
	) = @_;
	
	$length = $slice->total("Taxa") if $slice->total("Taxa") < $length;
	
	my @out;
	my $textgraph;
	my @tmparray;
	my @mar;
	
	my $sta = Statistics::Descriptive::Sparse->new();

	for (my $i=0; $i<8; $i++) {  # why nine previously? should be 8??????
		if (defined $ca[$col][$i]) {  # this might need to be "if ( $ca[$col][$i] ) "
			push @tmparray, $ca[$col][$i];
		}	
		else {
			push @tmparray, 0;
		}
	}
	
	$sta -> add_data(@tmparray); #@{$ca[$col]}[0..8]);
	
	my $sum = $sta->sum();
	if ($sum==0) {die "sum is zero! - died in mx_colcompgraph on $col"}
	
	my ($mores, $ones, $subtracts) = (0) x 3;
	for ( my $i = 0; $i < 8; $i++) {
		if (int(($ca[$col][$i] / $sum) * $length) == 0 and $ca[$col][$i] > 0) { # if there is a tiny percentage plot a least one
			$out[$i] = 1;
			$ones++;
		}
		elsif ($ca[$col][$i] == 0) {
			$out[$i] = 0;
		}
		else {
			$out[$i] = int(($ca[$col][$i]/$sum) * $length);
			$mores += int(($ca[$col][$i]/$sum) * $length);
		}
	}

	my $gh = Statistics::Descriptive::Full->new();
	if ($mores+$ones-$length > 0) { #must fix length by pseudo subtracting
		for (my $r = 0; $r < $mores+$ones-$length; $r++) {
			$gh->add_data(@out);
			$out[$gh->maxdex]--;
			$gh=Statistics::Descriptive::Full->new();
		}
	}
	
	elsif ($mores+$ones-$length < 0) { #must fix length by psuedo adding
		for (my $r = 0; $r < abs($mores+$ones-$length); $r++) {	
			$gh->add_data(@out);
			$out[$gh->maxdex]++;
			$gh=Statistics::Descriptive::Full->new();
		}
	}

	for (my $i=0; $i < 8; $i++) {
		$textgraph .= ($rna[$i] x $out[$i]);
	}
		
	push @out, $textgraph; # just for fun, toss on another value (index is 8)
	return @out;
}


sub tax_stats_bp_comp () {
	# returns a two dimensional array of stats describing block composition for a given taxon/fragment
	#  [0 all data; 1 - bracketed only; 2- stems only; 3- non stems, non brackets]
	#  [0-a; 1-u; 2-g; 3-c; 4-other; 5-length]
		
	my (
		$tax_index,
		$slice,	
	) = @_;

	my @tmpstats; 
	
	my $fragment;
	foreach my $blk ($slice->loop("Blocks")) { 	# all data
		$fragment .= $mx[$tax_index][$blk][0];
	}
	@{$tmpstats[0]}[0..8] = &str_rnafragmentcomp($fragment);
	#my %tmp = &str_unique_chars($fragment);
	@{$tmpstats[0]}[9] = &str_unique_chars($fragment); 
	
	$fragment = "";
	foreach my $blk ($slice->loop("Blocks")) {	#bracketed only
		if ($nbd[$blk][3] == 1) { 
			$fragment .= $mx[$tax_index][$blk][0];
		}
	}
	@{$tmpstats[1]}[0..9] = &str_rnafragmentcomp($fragment);
	#@{$tmpstats[3]}[9] = &str_unique_chars($fragment);
	
	$fragment = "";
	foreach my $blk ($slice->loop("Blocks")) {	# helicies only
		if (grep $_ eq $blk, %stems) {
			$fragment .= $mx[$tax_index][$blk][0];
		}
	}
	@{$tmpstats[2]}[0..9] = &str_rnafragmentcomp($fragment);
	#@{$tmpstats[3]}[9] = &str_unique_chars($fragment);
	
	$fragment = "";
	foreach my $blk ($slice->loop("Blocks")) { 	#non-stems, non brackets
		$fragment .= $mx[$tax_index][$blk][0] unless (grep $_ eq $blk, %stems); 
	}
	@{$tmpstats[3]}[0..9] = &str_rnafragmentcomp($fragment);
	#@{$tmpstats[3]}[9] = &str_unique_chars($fragment);
	
	return @tmpstats;
}

sub out_column_stats () {
		my $input = shift;
		my $slice = shift;
		my %params = @_;
		
		my @tmpblks;
		my $col_obj;
		
		my $head_check = 1;
		my @empty = ();
			
		chdir($rootdir);
		&io_confirmdir("analyses/stat/column/$input->{Modelroot}");

		$params{'-exclude'} = \@empty if not defined $params{'-exclude'};
	
	open (OUT, ">$input->{Modelroot}.txt") || die "couldn't open file output in out_column_stats - \n";
		select (OUT);
		
		print &str_file_header();
		
		# make sure a global alphabet object is predefined (column object alphabet is local to object)
		if ( not defined $params{'-alphabet'} ) {
			my $mega;
			foreach my $blk ($slice->loop("Blocks")) {	
				foreach my $tax ($slice->loop("Taxa")) {
					$mega .= $mx[$tax][$blk][0];
				}
			}
			my @alpha = &str_unique_chars($mega) ;
			$params{'-alphabet'} = \@alpha;
			print "[no alphabet was passed - using all chars present in slice: (", @{$params{'-alphabet'}}, ")]\n";
		};

		foreach my $blk ($slice->loop("Blocks")) {
			next if $nbd[$blk][3] == 1;	
			foreach my $tax ($slice->loop("Taxa")) {
				push @tmpblks, $mx[$tax][$blk][0];
			}
			my $col_obj = columns->new(@tmpblks);

			if ($head_check == 1)  {
				$head_check = 0; 
				$col_obj->table(
								%params,
								'-offset' => $nbd[$blk][6],
								'-header' => 'true'
							);
			}
			else{
				$col_obj->table(
								%params,					
								'-offset' => $nbd[$blk][6],
								'-header' => 'false'			# overrides %params!
							);
			}
			@tmpblks=();
			print "\n" if $params{'-blockgaps'} eq 'true';
		}
		select (*STDOUT);
		close OUT;
		print "\nstat results saved at analyses/stat/column/",$input->{Modelroot},"/",$input->{Modelroot},".txt \n";
}


sub out_block_stats () { 
	my (
		$input,
		$slice,
		$char	# a regex expression, the count/stats of which will be returned
	
	) = @_;

	die "must provide a regex for out_block_stats\n" if $char eq "";
	
	my %count;
	my $c;

	chdir($rootdir);
	&io_confirmdir("analyses/stat/block/$input->{Modelroot}");

	open (STAT, ">$input->{Modelroot}.txt") || die "couldn't open nexus file output too in out_stat - \n";

	print STAT "stat regex: $char\n\n";
	
	foreach my $tax ($slice->loop("Taxa")) {
		print STAT "$ters{$tax}\t";
		foreach my $blk ($slice->loop("Blocks")) {	
			$c = ( $mx[$tax][$blk][0] =~ s/$char//gi);
			$c = 0 if not $c;
 			print STAT $c, "\t";
			push @{$count{$blk}} , $c;
		}	
		print STAT "\n";	
	}

	my $sta_out;

	foreach my $blk ($slice->loop("Blocks")) {
		my $sta = Statistics::Descriptive::Sparse->new();
		$sta -> add_data(@{$count{$blk}});
			$sta_out->{"min"}->{$blk} = $sta->min();
			$sta_out->{"max"}->{$blk} = $sta->max();
			$sta_out->{"sum"}->{$blk} = $sta->sum();
			$sta_out->{"mean"}->{$blk} = $sta->mean();
			$sta_out->{"std_dev"}->{$blk} = $sta->standard_deviation();
			$sta_out->{"var"}->{$blk} = $sta->variance();
	}

	my @report_order = qw(min max sum mean std_dev var); # so the stats reports in a standard order

	print STAT "\n";	
	foreach my $key (@report_order) {
			print STAT "$key\t";
		foreach my $blk ($slice->loop("Blocks")) {
			print STAT $sta_out->{$key}->{$blk};
			print STAT "\t";	
		}
		print STAT "\n";
	}
	
	close STAT;
	print "\nout_stat results saved at analyses/stat/block/",$input->{Modelroot},".txt \n";
}

sub out_slice_stats () { # acts the same as blockstats, but assumes the slice represents one block
	my (
		$input,
		$slice,
		$char,	# a regex expression, the count/stats of which will be returned
		$mode	# 0 stats only; 1- include merged slice column 
	) = @_;

	die "must provide a regex for out_block_stats\n" if $char eq "";
	
	my @count;
	my $c;

	my %merged;

	foreach my $tax ($slice->loop("Taxa")) {
		foreach my $blk ($slice->loop("Blocks")) {
			$merged{$tax} .= $mx[$tax][$blk][0];
		}
	}
		
	chdir($rootdir);
	&io_confirmdir("analyses/stat/slice/$input->{Modelroot}");

	open (STAT, ">$input->{Modelroot}.txt") || die "couldn't open nexus file output too in out_stat - \n";

	print STAT "stat regex: $char\n\n";
	
	foreach my $tax ($slice->loop("Taxa")) {
		print STAT "$ters{$tax}\t";

			$c = ( $merged{$tax} =~ s/$char//gi);
			$c = 0 if not $c;
 			print STAT $c, "\t";
			print STAT "\t$merged{$tax}" if $mode == 1;
			push @count , $c;
	
		print STAT "\n";	
	}

	my $sta_out;
	my $sta = Statistics::Descriptive::Sparse->new();
	$sta -> add_data(@count);
	$sta_out->{"min"} = $sta->min();
	$sta_out->{"max"} = $sta->max();
	$sta_out->{"sum"} = $sta->sum();
	$sta_out->{"mean"} = $sta->mean();
	$sta_out->{"std_dev"} = $sta->standard_deviation();
	$sta_out->{"var"} = $sta->variance();
	
	my @report_order = qw(min max sum mean std_dev var); # so the stats reports in a standard order

	print STAT "\n";	
	foreach my $key (@report_order) {
		print STAT "$key\t$sta_out->{$key}\n";
	}
	
	close STAT;
	print "\nout_stat results saved at analyses/stat/slice/",$input->{Modelroot},".txt \n";
}


sub out_web () {
	my $input = shift;	# an input object
	my (
		$slice,			# restrict output to this slice
		$mode, 			# see intlv_working
		$size,			# see intlv_build, not required for mode == 1
		$interleave,	# model interleave, if not included is based on $orig_intlv ** will need to modify ultimately to explicitly pass
	) = @_;
		
	$slice->describe;

	my $working_interleave;
	
	if (defined $interleave) {
		$working_interleave = $interleave->Jrna::intlv_working($slice, $mode, $size);
	}
	else {
		$working_interleave = $orig_interleave->Jrna::intlv_working($slice, $mode, $size);
	}
	
	$working_interleave->describe;
	
	my $colcounttmp = 1; 
		
	#directory functions - confirmat that they exist

	chdir($rootdir);
	&io_confirmdir("models/$input->{Modelroot}");
	
	foreach my $path ( qw(blkd bpf mx sl stm) ) {
		chdir($rootdir);
		&io_confirmdir("models/$input->{Modelroot}/$path");
	}
	chdir($rootdir);

	&calc_colstats($slice);
	print "\ncolumn stats computed\n";

	print "stem base-pair frequencies computed\n";
	
	print "\ngenerating html\n";
	#interleave overview
	open (ILVOOUT, ">models/$input->{Modelroot}/interleave_index.html") || die " - can't open interleave overview file - \n";			
		&html_head (*ILVOOUT, "model - interleave overview", "../../style/model_index.css", "../../style/nav.css" );
			print &descr_mx(*ILVOOUT, $input, $working_interleave,$slice, 1);
		&html_foot (*ILVOOUT);
	close ILVOOUT;
	
	my $count =1;
	foreach my $i ($working_interleave->interleaves) {
		print "writing web_out - interleave $i, $count of ", $working_interleave->total_interleaves, "\n";
		$count++;
		
		# generate a slice object to pass to each descr_
		my $curslice = slice->new;
		$curslice->blocks($working_interleave->blocks($i));
		$curslice->taxa(keys %{$slice->taxa});
		
		## bpf
		open (BPFOUT, ">models/$input->{Modelroot}/bpf/bpf$i.html") || die " - can't open stm file - \n";
			&descr_blkheader(*BPFOUT, $i);
			&html_head (*BPFOUT, "Stems - interleave $i, $deschead[$i][1]", "../../../style/bpf_style.css", "../../../style/nav.css");
			&out_web_nav (*BPFOUT, $working_interleave, $i, 4);
			&descr_stembpfreq (*BPFOUT, $curslice, 1, 1);
			&out_web_nav (*BPFOUT, $working_interleave, $i, 4);
			&html_foot (*BPFOUT);
		close BPFOUT;
		
		## stems
		open (STEMOUT, ">models/$input->{Modelroot}/stm/stm$i.html") || die " - can't open stm file - \n";
			&html_head (*STEMOUT, "Stems - interleave $i, $deschead[$i][1]", "../../../style/stm_style.css", "../../../style/base_style.css", "../../../style/nav.css");
			&descr_blkheader(*STEMOUT, $i);
			&out_web_nav (*STEMOUT, $working_interleave, $i, 3);
			&descr_stems(*STEMOUT, 1, $curslice, $i, $working_interleave); 
			&out_web_nav (*STEMOUT, $working_interleave, $i, 3);
			&html_foot (*STEMOUT);
		close STEMOUT;
		
		## matrix
		open (MXBLKOUT, ">models/$input->{Modelroot}/mx/mx$i.html") || die " - can't open mx file - \n";
			&html_head (*MXBLKOUT, "Matrix data - interleave $i, $deschead[$i][1]", "../../../style/jrnastyle.css","../../../style/nav.css");
			&descr_blkheader(*MXBLKOUT, $i);
			&out_web_nav (*MXBLKOUT, $working_interleave, $i, 1);
			&out_mx (*MXBLKOUT, $curslice, 3, 0, 1, 1, 1, " ", 0, 2, 2, $defaultcolstart, 2, 0, 0);
			&out_web_nav (*MXBLKOUT,  $working_interleave, $i, 1);
			&html_foot (*MXBLKOUT);
		close MXBLKOUT;

		## sequence lines
		open (SLOUT, ">models/$input->{Modelroot}/sl/sl$i.html") || die " - can't open sl file - \n";
			&html_head (*SLOUT, "Column composition - interleave $i, $deschead[$i][1]", "../../../style/sl_style.css", "../../../style/nav.css");
			&descr_blkheader(*SLOUT, $i);
			&out_web_nav (*SLOUT, $working_interleave, $i, 2);
			$colcounttmp = &descr_seqline (*SLOUT, $curslice, 1, 40, 10, 1, $colcounttmp);    # doesn't number on zero mode correctly $intlvdesc[$i][0], 
			#print "colcount: $colcounttmp\n";
			&out_web_nav (*SLOUT,  $working_interleave, $i, 2);
			&html_foot (*SLOUT);
		close SLOUT;

		# block descriptions
		open (MXBLKDESCOUT, ">models/$input->{Modelroot}/blkd/blkd$i.html") || die " - can't open blkd file - \n";
			#multiple calls!
			&html_head (*MXBLKDESCOUT, "Block description - interleave $i, $deschead[$i][1]", "../../../style/blkd_style.css","../../../style/nav.css");
			&descr_blkheader(*MXBLKDESCOUT, $i);
			&out_web_nav (*MXBLKDESCOUT, $working_interleave, $i, 0);

			foreach my $blk ($curslice->loop("Blocks")) {
				&descr_blkstats (*MXBLKDESCOUT, $blk, $curslice, $working_interleave, 1);		#may need to doublecheck
			}
	
			&out_web_nav (*MXBLKDESCOUT, $working_interleave, $i, 0);
			
			&html_foot (*MXBLKDESCOUT);
		close MXBLKDESCOUT;

		#print "done\n";
	}		

	sub out_web_nav {
		local (*FH) = shift;
		my ($interleave, $intlv, $pos) = @_;

		print FH '<div class="nav">'; 
		print FH &html_blknav($interleave, $intlv, $pos, 0); 
		print FH &html_blknavlink( $intlv, $pos, 1, "", 0);
		print FH '</div>'; 
		print FH "\n";
	}
}

sub str_nuc_at_stem_pos () {
	# return the bp's at the nth stem positio in given the two helical pairs
	my (
		$pos,
		$helix,
		$mask
	) = @_;

	my @stempositions;
	my $nuc;
	
	if ($mask =~ /\)/g) {
		$helix = scalar reverse $helix;
		$mask = scalar reverse $mask;
		#print $mask,$helix;
		@stempositions = &str_pos_offsets('\)', $mask)
	}
	else {
		@stempositions = &str_pos_offsets('\(', $mask)
	}

	return (substr($helix, $stempositions[$pos], 1), $stempositions[$pos])	
}
# interleave package methods to be ported to interleave package

sub intlv_numseq () { # method for an interleave object; returns a string that column-numbers a given range of blocks or columns for a given interleave object
	# requires:
	my $self = shift;
	my (	
		$interleave, 	# specific interleave to number
		$mode,			# numbering mode: 0 - number individual columns; 1- number blocks	
		$format,		# string output format: 0 - normal; 1- ??; 3 - css; 4 - colored columns;
		$consecutive,	# 0 - renumber consecutively; 1 - number according to original format	- NOT IMPLEMENTED
		$gapped,		# gapped (insert gaps between description blocks)
		$gapchar,		# gapchar (used to be gapsize)
		$bracketed,		# bracketed - inlcude or execlude bracketed blocks
		$startcol	
	) = @_;
		
	my ($numrows, $longest) = (0) x 2;
	my ($t, $outstring);
	
	if ($mode == 0) { 
		foreach my $blk ($self->blocks($interleave)) {
			$longest = $nbd[$blk][2] if $nbd[$blk][3]==0;
		}
	}
	elsif ($mode == 1) {
		$longest = $self->last_block_in_interleave($interleave);
	}
	else { die "illegal mode [$mode] in intlv_numseq"}
	
	$numrows = length("$longest");
		
	if ($consecutive == 0) {

	}
	
	elsif ($consecutive == 1) {	
		for (my $i=0; $i < $numrows; $i++) {  # print this many rows of data (for sequence numbering or color composition blocks
					
		$outstring .= &css_tagswitch($format, 0, 0, "numrow");
		$outstring .= &css_tagswitch($format, 1, 0, "rightholder");
		
		#print &str_padright("["," ", $startcol);  #put $i. to check number of rows
		$outstring = $outstring.&str_padright("["," ", $startcol);
		
		$outstring .= &css_tagswitch($format, 1, 1, "rightholder");
		$outstring .= &css_tagswitch($format, 1, 0, "colnums");
		
		foreach my $blk ($self->blocks($interleave)) {			
			if ($mode == 0) {
				if ($nbd[$blk][3]==0) { 	#number or color code this block this block				
					for (my $k = $nbd[$blk][6]; $k < $nbd[$blk][7]+1 ; $k++) {  #using unbracketed blocks only now, number individual positions	
						#if mode is standard (0,3)
						$t = &str_padleft($k,"0",$numrows);
						$t =~ /.{$i}(.)/; 			# there's likely a better way to choose either the ith character or print a zero.
						#print $1;
						$outstring = $outstring.$1;

						# if mode is color columns (4) 
						# print a spanned space with color according to r
						
						# if mode is links (5)	
					}
			
					if ($gapped == 1) {
						$outstring = $outstring.$gapchar
					}
				}

				elsif ($nbd[$blk][3]==1 and $bracketed == 1) { #gap this block
					$outstring = $outstring.&str_padright("[", "-", $nbd[$blk][0]+1); #+1 is correct!!!

					#print ']';
					$outstring = $outstring.']';
					
					if ($gapped == 1) {
						#print &str_padright(" "," ",$gapsize-1);
						$outstring = $outstring.$gapchar #&str_padright(" "," ",$gapsize-1);
					}
				}	
	
				else {
				} #do nothing
			}

			elsif ($mode ==1) { # number blocks by their original position 

				$t = &str_padleft($blk, "0" , $numrows);
				$t =~ /.{$i}(.)/; 	
				my $t1 = &str_padleft("$1"," ", $nbd[$blk][0]);
				
				if ($nbd[$blk][3]==0) { 							
					$outstring = $outstring.$t1;
					
					if ($gapped == 1) {
						$outstring = $outstring.$gapchar;
					}
				}

				elsif ($nbd[$blk][3]==1 and $bracketed == 1) { # also bracket this block to indicate original bracketing
					$outstring = $outstring."[$t1]";

					if ($gapped == 1) {				
						$outstring = $outstring.$gapchar
					}
				}	
	
				else {
				} #do nothing
			}
		}

				
		#print "]\n";
		$outstring = $outstring."]";

		$outstring .= &css_tagswitch($format, 1, 1, "colnums");

		$outstring .= &css_tagswitch($format, 0, 1, "numrow");

		$outstring .= "\n";
	}

}
	return $outstring;
}	


sub slice_stems {	# a slice method that returns a list of blocks to exclude given a certain criterion
	my $self = shift;
	my $mode = shift;

	my @excludedblks;
	
	#pre-parse	
	if ($mode == 1) { # include only helices
		foreach my $blk (keys %{$self->blocks}) {
			 push (@excludedblks, $blk) unless grep ($_ == $blk, %stems); 
		}
	}
	elsif ($mode == 2) { # exclude all non-helices
		foreach my $blk (keys %{$self->blocks}) {
			 push (@excludedblks, $blk) unless not grep ($_ == $blk, %stems);
		}
	}
	elsif ($mode == 3) { # exclude everyting but 5' strands
		foreach my $blk (keys %{$self->blocks})  {
			 push (@excludedblks, $blk) unless grep ($_ == $blk, keys %stems);
		}
	}
	
	elsif ($mode == 4) { # exclude everyting but 3' strands
		foreach my $blk (keys %{$self->blocks}) {
			 push (@excludedblks, $blk) unless grep($_ == $blk, values %stems);
		}
	}
	elsif ($mode == 5) { # exclude everything but a random 1/2 of the stems <--- NOT WORKING	
	}

	return @excludedblks;	
}

sub descr_taxa  () {
	my (
		$slice
	) = @_;

	foreach my $t ($slice->loop("Taxa")) {
		my @stats = &tax_stats_bp_comp($t, $slice);
		print "$ters{$t}\t", join " ", @{$stats[0]}[0..9] , "\n";
	}	
}


sub blkstats () { # returns an array with a number of results from Statistics::Descriptive
	my (@blk) = @_;
	my @out;

	my $sta = Statistics::Descriptive::Sparse->new();

	$sta -> add_data(@blk);
	
	push @out, (
		 $sta->min(),
		 $sta->max(),
		 $sta->mean(),
		 $sta->standard_deviation(),
		 $sta->variance()
		);
	return @out;
}


sub descr_average_length () { # descr
	my $slice = shift;
	
	#!! variance is not reported correctly? must pass 0, not undefs to the stats object?
	
	my (@all, @bracketed, @unbracketed, @helices);
	my (@s_all, @s_bracketed, @s_unbracketed, @s_helices);
	
	foreach my $taxon ($slice->loop("Taxa")) {	
		push @all, length (&str_rawseq($slice, $taxon, 0));
		push @unbracketed, length (&str_rawseq($slice, $taxon, 1));
		push @bracketed, length (&str_rawseq($slice, $taxon, 2));
		push @helices, length (&str_rawseq($slice, $taxon, 3));	
	}

	@s_all = &blkstats(@all);
	@s_bracketed = &blkstats(@bracketed);
	@s_unbracketed = &blkstats(@unbracketed);
	@s_helices = &blkstats(@helices);

	print "\n";
	print "blocks: "; 
	print "region\tmin\tmax\tmean\tstd dev\tvariance\n";
	print  "all\t"; printf("%3.1f\t", $_) for @s_all;
	print  "\nbrak\t";  printf("%3.1f\t", $_) for @s_bracketed;
	print  "\nunbrak\t";  printf("%3.1f\t", $_) for @s_unbracketed;
	print  "\nhelix\t";  printf("%3.1f\t", $_) for @s_helices;
	print "\n";
	
}

sub descr_mx() { # describes slice/interleave combination for given input 
	# requires
	local (*FH) = shift;
	my (	
		$input,
		$interleave,
		$slice,
		$mode,	# 0-txt; 1-css
	) = @_;

	my %totmx = &slice_calc_mx_totals($slice);
	my %totslice = &slice_calc_stem_totals($slice);
	
	$interleave = $orig_interleave if not $interleave; # a cheat - if the interleave isn't defined use the original one
	
	local *STDOUT = *FH;
	my $thekey;
	my $i;

	my $totals_str = join " ", (
		"total taxa:",	$slice->total("Taxa"), "\n",
		"total interleaves:", $interleave->total_interleaves, "\n",
		"total blocks:", $slice->total("Blocks"), "\n",
		"total bracketed blocks:", $totmx{"numbrakblks"},"\n",
		"total unbracketed blocks:", $slice->total("Blocks")-$totmx{"numbrakblks"}, "\n",
		"total bracketed chars:", $totmx{"numbrakchars"},"\n",
		"total unbracketed chars: ", $totmx{"numunbrakchars"},"\n",
		"total chars:", $totmx{"totchars"}, "\n",
		"longest bracketed block:", $totmx{"longestbrakblk"}, "\n",
		"total blocks defined as helicies:", $totslice{"bdas"}, "\n",
		"total length of defined non-helicies:", $totslice{"length_bdas"}, "\n",
		"total five-prime strands:", $totslice{"num_non_prime_bdas"}, "\n",
		"total three-prime strands:", $totslice{"num_prime_bdas"}, "\n",
		"total complete helices:", $totslice{"num_complete_helices"},"\n",
		"total complete pairs in blocks defined as helices:", $totslice{"num_complete_pairs"},"\n"
	);
	$totals_str =~ s/\n/\<br\>\n/g if $mode == 1;
	
	#summary statistics
	if ($mode == 0 ) {
		# other things to print here:
		# - interleave description
		# - ordered terminal list/index
		# - block starting, ending and length

		print $totals_str;
	}
	
	else { # css mode, returns to *FH
		print '<div class="page_header">';
			print '<a href="index.php">return home</a>'; #may want to change this to html
		print '</div>';

		print '<div class="section">';
			print 'model: <span class="hmod">',$input->{Modelroot}, '</span><span class="hfiles">datafile: ', "$input->{Infile}<br>stem index file: $input->{Stemindex}",'</span>';
			print "generated with jrna version: $VERSION<br>"; 
		print '</div>', "\n";

		print '<div class="section">', "\n";			
				print '<span class="header">Totals</span>', "\n";	
				print $totals_str;
		print '</div>', "\n";

		print '<div class="section">';			
			print '<span class="header">Taxa</span>', "\n";	
			map  {print "[$_ - $ters{$_}] "} ($slice->loop("Taxa"));
		print '</div>', "\n";

		print '<div class="section">', "\n";
		print '<span class="header">Interleaves</span>', "\n";			
			foreach my $i ($interleave->interleaves) {
				
				print '<div class="iblk">', "\n";
						print '<span class="col_half">[', $i, '] ', "$deschead[$i][0]: $deschead[$i][1] </span>";
						print '<span class="col_half"><span style="font-size: smaller;">', &html_blknavlink($i, 999, 999, "root", 1), '</span></span>';

						print '<span class="col_whole">Blocks: ', join " ", $interleave->blocks($i),'</span>'; 
				
					# print '<span class="colnum">';if (not defined $intlvdesc[$i][2]){print '0'} else {print $intlvdesc[$i][2]} ; print '</span>'; # was 1 in second >=1
					# print '<span class="colnum">';if (not defined $intlvdesc[$i][3]){print '0'} else {print $intlvdesc[$i][3]} ; print '</span>'; # was 1 in second >=1 ; UNINITIALIZED!!!
					
					# print '<span class="shead">'; print "$i $deschead[$i][0] $deschead[$i][1]  unbrakchars: $intlvdesc[$i][2] brakchars: $intlvdesc[$i][3]";
					# print '</span>';
				
				print '</div>', "\n";
			}
		print '&nbsp;</div>'; print "\n";
	}	
	
	return 1;
}

sub descr_stembpfreq () {
	#requires:
	local (*FH) = shift; 	# filehandle
	my (
		$slice,		# slice 
		$mode,		# 0- text only; 1- css;
		$summary_mode   # 0- include all taxa; 1- exclude "??" from calculations
	) = @_;

	my $num_quest;
	
	&calc_stem_bp_stats($slice, $summary_mode);

	local *STDOUT = *FH;

	my @displayedbp = ("gc", "cg", "ua", "au", "gu", "ug", "aa", "ac", "ag", "ca", "cc", "cu", "ga", "gg", "uc", "uu", "--"); # re-include "??" for web
	my $intlvhead = 0;
	my ($blk, $bp);
	my @wrkstems;

	my $total_taxa = $slice->total("Taxa");

	
	foreach $blk ($slice->loop("Blocks")) {
		if (exists $stems{$blk}) {
			push @wrkstems, $blk;
		}
	}
	
	if ($mode == 0) { #text mode
		if ($#wrkstems == -1) {
			print "no stems begin in this interleave\n";
			return;
		}
		
		# print initial header
		print "\t\t\#_tax\t";
		foreach $bp (@displayedbp) { 
			print "$bp\t";
		}
		
		print "\n";
	
		foreach $blk ($slice->loop("Blocks"))  {
			for (my $p=0; $p < $#{$bonds[$blk]} + 1  ; $p++) { # loop each character  
			# print intlv title?		
				# print stem/header title
				if ($p == 0) {
					
					if ($intlvhead >= 20) { #reprint the header line
						print "\t\t\#_tax";
						foreach  $bp (@displayedbp) {
							print "\t$bp";
						}
			
						print "\n";
						$intlvhead = 0;
					}
					print "$sindex{$blk}";
				}
				
				#print position
				print "\t",$p+1;

				#print taxa used
				if (not defined $bonds[$blk][$p][$basepairs{'??'}][0]) {$num_quest = 0} else { $num_quest = $bonds[$blk][$p][$basepairs{'??'}][0]} 
				print "\t", $slice->total("Taxa") + 1 - $num_quest ;
				
				#print stats
				foreach $bp (@displayedbp) { #loop through basepairs to display
					print "\t";
					printf ("%3.1f", $bonds[$blk][$p][$basepairs{$bp}][1]*100);
					unless (not defined ($bonds[$blk][$p][$basepairs{$bp}][2])) {
						print '*' if $bonds[$blk][$p][$basepairs{$bp}][2] == 1;				
					}
				}
				print "\n";
				$intlvhead++;
			}
		}
	}
		
	else { #css mode
	
		sub desc_stembpfreq_bp_head_row() {
			my @displayedbp = @_;
			print '<div class="brow">';
				print '<span class="col1">&nbsp;</span>';	
				print '<span class="col2">&nbsp;</span>';
				print '<span class="col2">&nbsp;</span>';
			
				foreach my $bp (@displayedbp) {
					print '<span class="col">';
					print $bp;
					print '</span>';
				}
			print '</div>';		
			print "\n";
		}
	
	if ($#wrkstems == -1) {
		print '<div class="hrow">';
		print "<p>no stems begin in this interleave</p>\n";
		print '</div>';
		return;
	}
		
	# print initial header
		print '<div class="SOMECLASS">';
		print "\n";	

		#header rows
		print '<div class="hrow">';
			print '<span class="col1">Helix</span>';	
			print '<span class="col2">Base</span>';
			print '<span style="float: left; text-align: center; width: 65%;">base pair composition %</span>';
		print '</div>';
		print "\n";
		
		print '<div class="hrow">';
			print '<span class="col1">&nbsp;</span>';
			print '<span class="col2">pair</span>';
			print '<span class="col2"># seqs</span>';
			
			print '<span style="float: left; width: 24%; padding-right: 1%; border-right: 1px solid; text-align: center;">complementary</span>'; #width should be .col % x 6
			print  '<span style="float: left; width: 52%; text-align: center;">non-complementary</span>';
		print '</div>';
		print "\n";
			
		# bp header row
		&desc_stembpfreq_bp_head_row(@displayedbp);
		
		foreach $blk (@wrkstems) {
			for (my $p=0; $p < $#{$bonds[$blk]} +1 ; $p++) { #loop through positions  ----------- $bonds[$blk][ END OF THIS INDEX] --------------

				# print stem/header title
				if ($p == 0) {
					if ($intlvhead >= 20) { #reprint the bp header line
						# bp header row
						&desc_stembpfreq_bp_head_row(@displayedbp);
						$intlvhead = 0;
					}

					#print first data row columns here
					print '<div class="drowb">';
					print '<span class="col1">';
						print $sindex{$blk};
					print '</span>';
				}
				
				else {
					print '<div class="drow">';
						print '<span class="col1">&nbsp;</span>';
				}

					print '<span class="col2">';
						print $p+1;
					print '</span>';
					
					print '<span class="col2">';
						if (not defined $bonds[$blk][$p][$basepairs{'??'}][0]) {$num_quest = 0} else { $num_quest = $bonds[$blk][$p][$basepairs{'??'}][0]} 
						 print $slice->total("Taxa") + 1 - $num_quest ;
					print '</span>';
					
			#print stats
				foreach $bp (@displayedbp) { #loop through basepairs to display
					print '<span class="d';
				
					#if (($bonds[$blk][$p][$basepairs{$bp}][1] == 0) or (not defined $bonds[$blk][$p][$basepairs{$bp}][2])) {
						#	print '0">';
						#}

					#else {
					
						if (not defined $bonds[$blk][$p][$basepairs{$bp}][2]) {
							print '1" style="background: rgb('; # not covarying - blue
							print 240 - int($bonds[$blk][$p][$basepairs{$bp}][1]*240);
							print ",";
							print 240 - int($bonds[$blk][$p][$basepairs{$bp}][1]*240);
							print ",";
							print 240; 
							print ');">';
						}
						
						else { # ($bonds[$blk][$p][$basepairs{$bp}][2] == 1) { #bp's covarying - red
							print '1" style="background: rgb(';
							print 240;
							print ",";
							print 240 - int($bonds[$blk][$p][$basepairs{$bp}][1]*240);
							print ",";
							print 240 - int($bonds[$blk][$p][$basepairs{$bp}][1]*240);
							print ');">';
						}

					unless ($bonds[$blk][$p][$basepairs{$bp}][1] == 0){
						printf ("%3.1f", $bonds[$blk][$p][$basepairs{$bp}][1] *100);
						
					}
					print '&nbsp;';
					print '</span>';				
				}

				print '</div>';
				print "\n";
				$intlvhead++;
			} #loop rows
		} #loop blocks
		print '</div>'; #end completely block
	}
}

sub descr_seqline () {   # output sequence line - one of the three major block pages; returns colcount (equivalent to endposition)
	#requires:
	local (*FH) = shift; 	# filehandle
	my (
		$slice,
		$mode,			# mode (0) - screen; (1) - css
		$bgl,			# length in letters of bar graph 
		$interval,		# left seq interval	
		$blknummode,	# block mode (0) - absolute (1) - paup
		$startpos		# startposition - if in block mode (1) and the interleve starts with a bracketed block we need an absolute starting point
	) = @_;

	local *STDOUT = *FH;
	
	my ($mid, $bnm1, $bnm2, $navdirs, $seqnumspace, $colcount) = (0) x 6;
	my @barcomp;
	my ($nbm1, $nbm2, $midhit, $gc, $tmp);
	my $inblkpos = 1;	

	#print "\nstartpos: $startpos\n";
	
	# mode of numbering sequences 	- 0 - include bracketed columns in numbering
	if ($blknummode == 0) {  
		$nbm1 = 1;
		$nbm2 = 2;
	}
	else {			#	- 1 - don't include bracketed columns in numbering
		$nbm1 = 6; 
		$nbm2 = 7;
	}

	foreach my $blk ($slice->loop("Blocks")) {
		$mid =int ($nbd[$blk][0] / 2);	#find mid point of this block
		if ($mid==0) {$mid=1};

		$midhit = 1;

		if ($blknummode == 1 and $nbd[$blk][3] == 1) {
			$colcount = $startpos;   	
			$seqnumspace = $colcount % $interval; 
		}
		
		else {
			$colcount = $nbd[$blk][$nbm1];   		
			$seqnumspace = $nbd[$blk][$nbm1] % $interval; 
		}

		$inblkpos = 1;
		
		for  (my $p = $nbd[$blk][1]; $p < $nbd[$blk][2]+1; $p++) {  #this as nbm of paup doesn't return zero (temp)
		
			if ($mode == 0) { # !!!!!!!!!! doesn't work right now replace with logic from mode 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
			#different lines for blocked or not - text mode only
			if ($nbd[$blk][3]==0) {
				$gc = "|";
			}
			else {
				$gc = "*";
			}

				
				#block 1
				if ($colcount==$interval) {
					print "$p";
					$seqnumspace=0;
				}
	
				#block 2
				print "\t";
				print "$gc";
			
				#block 3
				if ($nbd[$blk][$nbm2]==$p) {
					print "_";
				}
				else {
					print " ";
				}
			
				if ($nbd[$blk][$nbm1]==$nbd[$blk][$nbm2] or $p==$nbd[$blk][$nbm1] or $p==$nbd[$blk][$nbm2]) {
					print "$p";
				} 
			
				#block 4
				print "\t";
				if ($nbd[$blk][3]==0) {
					@barcomp = &mx_colcompgraph($p,$bgl);
					print "$barcomp[8]";
				}

				else {
					print " " x $bgl;
				}
			
					
				#block 5
				if ($midhit==$mid) {
					print "\t";
					print " - $blk";					
				}
				print "\n";
			}

			elsif ($mode==1) { # linked css output
									
			#row
				print "\<div class=\"slr\">";

			#block 1
				print "\<span class=\"sl1\">";
					if ($seqnumspace == $interval) { #chance to write out
						if ($blknummode == 0) {
							print $colcount;
							$seqnumspace = 0;
						}
						else {
							if ($nbd[$blk][3] == 0) {
								print $colcount;
								$seqnumspace = 0;
							}
							else { # bracketed block that's equal but in 1 mode
								print '&nbsp;';
							}
						}
					}
					else {
						print '&nbsp;';
					}				
				print "\</span>";
				
			#block 2 (all css borders - the line and dividing lines)
				print "\<span class=\"sl2";		
					if ($nbd[$blk][3]==1) {
						if ($inblkpos == $nbd[$blk][0]) {			#crummy test (length vs counter), but it works
							print "bd";
						}
						else {
							print "b";
						}
					}
				
					else {
						if ($inblkpos == $nbd[$blk][0]) {
							print "nbd";
						}
						else {
							print "nb";
						}
					}
					print "\">";
					print '&nbsp;';	
				print "\</span>";
				

			#block 3
				print "\<span\ class =\"sl3\">";	
					if ($inblkpos == $nbd[$blk][0] || $inblkpos == 1) { #chance to print start/end
						if ($blknummode == 0) {
							print " |$p|";
						}
						else {
							if ($nbd[$blk][3] == 0) {
								#print " |";
								# print $nbd[$b][6] + $inblkpos -1; same as below
								print " $colcount";
								#print "|";
							}
							else {
								print '&nbsp;';	
							}
						}
					}
					else {
						print '&nbsp;';	
					}
				
				print "\</span>";
				
			#block 4
				print "\<span class=\"sl4\">";
					if ($nbd[$blk][3] == 0) {
						if ($blknummode == 0) { #
							@barcomp = &mx_colcompgraph($slice, $nbd[$b][6] + $inblkpos -1 ,$bgl);
						}
						else {
							@barcomp =  &mx_colcompgraph($slice, $colcount, $bgl);
						}

					$tmp = &css_colorRNAstr ($barcomp[8]);
					print "$tmp";
					}

					else {
						print '&nbsp;';
					}											
				print "\</span>";

				
			#block 5
				print "\<span class=\"sl5\">";
					if ($midhit==$mid) {
						print "\t";
						print " - $blk - $sindex{$blk}";
				
					}
					else {
						print '&nbsp;';
					}	
				print "\</span>";

				
			print "\<\/div\>";
			print "\n";

			} # end css output

			if ($blknummode == 0) { #easy don't do anything
				$seqnumspace++;
				$colcount++;
			}
			else {
				if ($nbd[$blk][3]==0) {
					$seqnumspace++;
					$colcount++;
				}
			}
			$inblkpos++;
			$midhit++;	
		}
	}		

	return $colcount; # return reference to the start of the next descr_seqline
	#*FH = *STDOUT;
}

sub descr_stems () {
	#requires: 
	local (*FH) = shift @_;
	my (	$mode,		# mode: 0- text; 1- css  - both not implemented
		$slice,		# a slice object
		$intlv,		# specific interleave 
		$interleave	# interleave object
	) = @_;
	
	#dependancies:
	# prior call to &stems_initialize

	my @bond;		 # $bond[$block][position]
	my @wrkstems;
	my @cws;
	my @collen;
	my @outblks;

	my $l;
	
	my ($tmp, $a, $reversed, $mask, $headtxt);	
	my ($i, $t, $curstate, $colwidth, $calc_on, $calc_headers, $headtxt_done) = (0) x 7;
	my $header = 15;

	local (*STDOUT) = *FH;
	
	# build a shorter array that holds only those stems to print (exmaine interleave blocks only)
	foreach $i ($interleave->blocks($intlv)) {	
		if (exists $stems{$i}) {
			push @wrkstems, $i;
		}
	}

	if ($#wrkstems == -1) {
		print '<div class="hrow">';
		print "<p>no stems begin in this interleave</p>\n";
		print '</div>';
		return;
	}
	
	#print "@wkrstems";

	# need a default - no stems in this block!! code to escape here	
	foreach $t ($slice->loop("Taxa")) { # taxon loop

		my (@row1, @row2, @bond);
		
		# pre parse the data, building the stem pieces and rows, and bond rows
		foreach my $blk (@wrkstems) {
			#build the arrays for the rows, calulate bonds here 
			my @rows = &str_aligned_stems($mx[$t][$blk][0], $mx[$t][$stems{$blk}][0], $nbd[$blk][8], $nbd[$stems{$blk}][8]);
			push @row1, $rows[0];
			push @row2, $rows[1];
			$mask = $rows[2];	#grab the mask	
			
			#build a header line
			if ($calc_headers  == 0) {
					#calculate an array of column lengths for each column
					push @collen, &num_max (length $rows[1], length "$sindex{$blk}-$sindex{$stems{$blk}}");	#? col+1	
			}
			
			#calculate bond, need to recalc for each taxon
			for (my $p = 0; $p < length $mask; $p++) { 		
				$mask =~ /.{$p}(.)/;

				if ($1 eq '.') { #no bond here
					$bond[$blk][$p] = 0;	# not compensated
				}
						
				else { # bonded
					$rows[0] =~ /.{$p}(.)/;
					$a = $1;	
					$rows[1] =~ /.{$p}(.)/;
			
					$a .= $1;	
										
					$a = lc ($a); #make sure its lowercase for matching purposes
															
					if ($a =~ /(au)|(ua)|(gc)|(cg)/ ) {
						$bond[$blk][$p] = 2;	# a fully compensated bond
					}
						
					elsif ($a =~ /(gu)|(ug)/) {
						$bond[$blk][$p] = 1;	# partially compensated
					}
					else {
						$bond[$blk][$p] = 0;	# not compensated
					}
				}
			}
		}
		
		if ($headtxt_done == 0) { #build the header row (s?)
			$headtxt .= '<div class="hrow">';  # descrow and mx row have a different property - mx uses nbsps?
				$headtxt .= '<span class="spc">&nbsp;</span>';
				
				$l=0;
				
				foreach my $blk (@wrkstems) {
					$headtxt .= '<span class="sb">';	
						
						$headtxt .= "$sindex{$blk}-$sindex{$stems{$blk}}";
						
						if (length "$sindex{$blk}-$sindex{$stems{$blk}}" < $collen[$l]) {
							$headtxt .= ('&nbsp;' x ($collen[$l] - length "$sindex{$blk}-$sindex{$stems{$blk}}"));
						}
						
					$headtxt .= '</span>';
					$l++;
				}	
				
				$headtxt .= '</div>'; #added
				$headtxt .= "\n";
			$headtxt_done = 1; 
		}
		
	
			$calc_headers = 1; # don't need to calculate header rows and column lengths the second time;

			#print "| @row1 | <br> \n";
			
			#print "| @row2 |\n";

			# header row
				if ($header == 15) { #actually print the header row
					print $headtxt;						
					$header = 0;
				}
								
			# row one
					print '<div class="srow">'; 
					
					print '<span class="taxa">';
						print $ters{$t};
					print '</span>';
					
					$l=0;
					foreach my  $blk (@row1) {

						print '<span class="sb">';
							print &css_colorRNAstr($blk);
							if (length $blk < $collen[$l]  ) {
								print '&nbsp;' x ($collen[$l]-length $blk);
							}
						print '</span>';

						$l++;
					}
					
					print '</div>';
					print "\n";	
					
			# row two 
					print '<div class="srow">'; 
					
					print '<span class="spc">';
						print '&nbsp;';
					print '</span>';

					$l=0;
					foreach my $blk (@row2) {
						print '<span class="sb">';
							print &css_colorRNAstr($blk);
							if (length $blk < $collen[$l]  ) {
								print '&nbsp;' x ($collen[$l]-length $blk)
							}
						print '</span>';
						
						$l++;
					}

					print '</div>';
					print "\n";	
		
			# bond row		
					print '<div class="brow">'; 
					print '<span class="spc">&nbsp;</span>';
	
					$l=0;
					
					foreach my  $blk (@wrkstems) {
						
						$curstate = $bond[$blk][0]; 

						print '<span class="sb">';	
							print '<span class="su';
								print "$bond[$blk][0]";
							print '">';
							
							for (my $p=0; $p< length $row1[$l]; $p++) { #was +1 - should loop through last index of bond[$blk]
							
								if ($bond[$blk][$p] == $curstate) {
									print $bond[$blk][$p]; 
								}

								else {
									print '</span>';
									print '<span class="su';
									
									print "$bond[$blk][$p]";
									print '">';
								
									print $bond[$blk][$p];

									$curstate = $bond[$blk][$p]; 
								}
							}
						
								
							print '</span>'; #end su#	
						
							if (length $blk < $collen[$l]  ) {
								print '&nbsp;' x ($collen[$l]-length $row1[$l]);
							}

							
						print '</span>'; #end sdq

						$l++;
					}

					print '</div>';
					print "\n";
			$header++;
		}  #end taxon loop 
}


sub descr_blkheader {
		local (*FH) = shift;
		my ($intlv) = @_;
		print FH '<div class="ilv_blk_header">'; 
		print FH "$deschead[$intlv][1]";
		print FH '</div>';
		print FH "\n";	
}

sub descr_blkstats () { # describe the composition of a given block
	# need to clean this up, seperate the parsing/stat generation from the reporting and make it more logical (use hash++ and library of chars, with remainders = other)
	#requires:
	local (*FH) = shift; # output handle
	my (
		$blk, 		# block index
		$slice,
		$interleave,
		$mode		# 0- to screen; 1 - css/html
	) = @_;

	my @str;
	my $c = 0;
	my $char;
	my @blcmp;
	my @statlength;
	
	local *STDOUT = *FH;
	
	#blcmp = block composition, in number of "nucleotides"
	foreach my $t ($slice->loop("Taxa")) {	
		# gather composition of basepairs for the current taxa      - need to declare as global; need to use str_rnafragmentcomp
		$blcmp[0][$t] = ($mx[$t][$blk][0] =~ tr/A//);
		$blcmp[1][$t] = ($mx[$t][$blk][0] =~ tr/U//);
		$blcmp[2][$t] = ($mx[$t][$blk][0] =~ tr/G//);
		$blcmp[3][$t] = ($mx[$t][$blk][0] =~ tr/C//);
		$blcmp[4][$t] = ($mx[$t][$blk][0] =~ tr/N//);
		$blcmp[5][$t] = ($mx[$t][$blk][0] =~ tr/-//);
		$blcmp[6][$t] = ($mx[$t][$blk][0] =~ tr/\?//);
		$blcmp[7][$t] = $nbd[$blk][0] - $blcmp[0][$t] - $blcmp[1][$t] - $blcmp[2][$t] - $blcmp[3][$t] - $blcmp[4][$t] - $blcmp[5][$t] - $blcmp[6][$t] ; # UNCHECKED - other

		for (my $z=0; $z<8; $z++) {
			 $blcmp[$z][$t] = 0 unless defined $blcmp[$z][$t] 
		}

		if ($nbd[$blk][3]== 1) { #don't calculate length for non bracketed regions
			push @statlength, ($nbd[$blk][0] - $blcmp[5][$t]);
		}
	}
		
	if ($mode == 0 ) {
		print "\n\n";
			
		print "$blk";	
		print " $nbd[$blk][4] - $nbd[$blk][5] ";
		print " (length: $nbd[$blk][0]; literal position: $nbd[$blk][1]-$nbd[$blk][2];";
	
		if ($nbd[$blk][3]==1) {
			print " bracketed ";
			@str = &blkstats(@statlength);
			print "@str";
		}
		
		else { print "paup position: $nbd[$blk][6]-$nbd[$blk][7] " 
		
		}
		
		print ")";
		print "\n";
		
		print "min-max\tmean\tstd_dev.\tvar\n";
		
		for (my $i= 0; $i<8; $i++) {
			print "$rna[$i]: ";
			
			@str = ();
			my @tmp;
			foreach my $tax ($slice->loop("Taxa")) {
				push @tmp, $blcmp[$i][$tax];
			}
			@str = &blkstats(@tmp);
			
			print "$str[0]-$str[1]\t";
			printf ("%0.2f", $str[2]);
			print "\t";
			printf ("%0.2f", $str[3]);
			print "\t";
			printf ("%0.2f", $str[4]);
			print "\t";
			
			print "\n";
		}
	} #end mode =0 

	elsif ($mode == 1) {
		my $curintlv = $interleave->{Interleaves}->{$blk}; 
		
		print '<div class="descblk">';
		print "\n";
		#header
			print '<div class="hrow1">';
				print '<span class="col1">';
				print $blk;
				print '</span>';
		
				print '<span class="hd1">';
                    		if (grep $_ eq $blk, %stems) { # tests to see if $blk is a key OR value in %stems, negate is: print 'Yep' if grep /^$foo$/, %hash;
					print "stem"
				}
				elsif (	$nbd[$blk][3]==1) {
					print "bracketed"
				}
				else {
					print "undefined"
				}
				print '</span>';

				print '<span class="hd2">';
					print $sindex{$blk};
				print '</span>';

			print '</div>';
			print "\n";
			
			print '<div class="hrow2">';
				print '<span class="col1">';
					print '&nbsp;';
				print '</span>';

				print '<span class="hd3">';
				print "length: $nbd[$blk][0];  literal position: $nbd[$blk][1]-$nbd[$blk][2];  ";
	
				if ($nbd[$blk][3]==1) {
					@str = &blkstats(@statlength);
					print " [min size-max size: $str[0]-$str[1], mean=";
					printf ("%0.2f", $str[2]);
					print ', std=';
					printf ("%0.2f", $str[3]);
					print ', var=';
					printf ("%0.2f", $str[4]);
					print '] ';	
				}

				else { 
					print " unexclusive position: $nbd[$blk][6]-$nbd[$blk][7] "; 
				}	
				print '</span>';
			print '</div>';

			print "\n";

			print '<div class="hrow3">';
				print '<span class="col1">';
				print '&nbsp;';
				print '</span>';

				print '<span class="col2">';
				print '&nbsp;';
				print '</span>';

				print '<span class="dcol">';
				print 'min-max';
				print '</span>';
				
				print '<span class="dcol">';
				print 'mean';
				print '</span>';
				
				print '<span class="dcol">';
				print 'std dev';
				print '</span>';
				
				print '<span class="dcol">';
				print 'var';
				print '</span>';
			print '</div>';
			print "\n";

		for (my $i= 0; $i<8; $i++) {
			print '<div class="drow">';

				print '<span class="col1">';
				print '&nbsp;';
				print '</span>';

				print '<span class="col2">';
				print $rna[$i];
				print '</span>';


			@str = ();
			my @tmp= ();
			foreach my $tax ($slice->loop("Taxa")) {
				push @tmp, $blcmp[$i][$tax];
			}
			@str = &blkstats(@tmp);

					print '<span class="dcol">';
					printf ("$str[0]-$str[1]");
					print '</span>';
				
					for (my $q = 2; $q <5; $q++) { # shorten the code
						print '<span class="dcol">';
						printf ("%0.3f", $str[$q]);
						print '</span>';
					}
						
			print '</div>';
			print "\n";
		}

		# spacer row
			print '<div>&nbsp;</div>';
			print "\n";

			print '</div>';
			print "\n";
		}
}

sub descr_seq_length () { # describes to STDOUT the combined sequence length for each taxon given the passed slice
		my $slice=shift;	
		my %outhash;
		my $seq;
		
		foreach my $t ($slice->loop("Taxa")) {
			foreach my $blk ($slice->loop("Blocks")) {
				$seq .= $mx[$t][$blk][0];
			}
			$seq =~ s/-//g;		#excludes dashes
			$seq =~ s/[0-9]//g;	#excludes numbers
			
			$outhash{$ters{$t}} = length ($seq);
			$seq = "";
		}
		my @bar = sort { $outhash{$a} <=> $outhash{$b} } keys %outhash;
		foreach my $foo (@bar) {print "$outhash{$foo} \t $foo \n" }
}

# series of subroutines to generate links between l,b,s,m-group html pages
sub html_blknav () {  # returns block navigation buttons
	my (	
		$interleave,	# interleave object
		$intlv,		# present interleave
		$pos,		# "postion" - a fileroot index
	) = @_;	

	my $link; # txt to return
	
	$link .= '<span class="bfbar">';
	
	#left
	if ( defined $interleave->left($intlv)  )  {
		$link .= "<span class=\"nfb\">";  #could simplify to an anonymous sub I suppose -> or a css spanned link
		$link .= "<a href=\"$fileroot[$pos]";
		$link .= $interleave->first_interleave;		#interleaves always start at zero  (nope they dont)
		$link .= ".html\">";
		$link .= '&lt;&lt;';
		$link .= "</a>";
		$link .= "</span>";

		$link .= "<span class=\"nfb\">";  #could simplify to an anonymous sub I suppose -> or a css spanned link
		$link .= "<a href=\"$fileroot[$pos]";
		$link .= $interleave->left($intlv);
		$link .= ".html\">";
		$link .= '&lt;';
		$link .= "</a>";
		$link .= "</span>";
	}

	#middle (home)	
	$link .= "<span class=\"nhome\">";
	$link .= '<a href="../interleave_index.html">';
	$link .= "home";
	$link .= "</a>";
	$link .= "</span>";

	#right
	if ( defined $interleave->right($intlv) ) {
		$link .= "<span class=\"nfb\">";
		$link .= "<a href=\"$fileroot[$pos]";
		$link .= $interleave->right($intlv);
		$link .= ".html\">";
		$link .= '&gt;';
		$link .= "</a>";
		$link .= "</span>";

		#rightmost fastforward
		$link .= "<span class=\"nfb\">";
		$link .= "<a href=\"$fileroot[$pos]";
		$link .=  $interleave->last_interleave;
		$link .= ".html\">";
		$link .= '&gt;&gt;';
		$link .= "</a>";
		$link .= "</span>";
	}

	$link .= '</span>';

	return $link;	
}

sub html_blknavlink () {  # returns block navigation buttons
	#requires
	my (
		$intlv,		# $intlv block to link to (current interleave block)
		$place,		# current "place" - 0- blkd; 1-mx, value not parsed if mode is 999
		$mode,		# mode - 0- inclusive (i.e. return link to current position); 1-exclusive (return links to all but current); 999 - all links
		$path,		# (root) path
		$length		# length - 0- long text names; 1- short text names;
	) = @_;
	
	my ($link, $tmppath);
	my %buttons;

	if ($length == 0) {
		 %buttons  = (0, "block composition", 1, "matrix", 2, "column composition", 3, "stems", 4, "base-pair frequency");
	}
	else {
		 %buttons  = (0, "bc", 1, "mx", 2, "cc", 3, "st", 4, "bpf");
	}
	
	$link .= '<span class="ivwnav">';
			$link .= '<span class="navdiv">';
			$link .= '|';
			$link .= '</span>';

	for (my $i=0; $i < $#fileroot+1; $i++) {	
		if (($mode == 0 and $place == $i) or ($mode ==1 and $place != $i) or ($mode == 999)) { 
			if ($path eq "root") {
				$tmppath = "$fileroot[$i]/$fileroot[$i]$intlv.html";
			}
			else {
				$tmppath = "../$fileroot[$i]/$fileroot[$i]$intlv.html";
			}

			$link .= "\<span class=\"iv$i\">";
			$link .= "\<a href=\"$tmppath\">"; # 		
			$link .= $buttons{$i};
			$link .= "</a>";
			$link .= "</span>";
			$link .= '<span class="navdiv">';
			$link .= '|';
			$link .= '</span>';
		}
	}
	
	$link .= '</span>';
	
	return $link;	
}

sub html_head () {
	# requires
	local (*FH) = shift;
	my $title = shift;
	my @css = @_; # can have multiple stylesheets
	
	my $sheet;
	
	print FH '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">',"\n";
	print FH "<html>\n";
	print FH "<head>\n";
	print FH "<title>$title</title>\n";

	foreach $sheet (@css) {
		print FH "<link rel=\"stylesheet\" href=\"$sheet\">\n";
	}

	print FH "</head>\n<body>\n";
}

sub html_foot () {
	local (*FH) = shift;
	print FH "</body>\n</html>";
}


sub css_colorRNAstr () { # css styles a string of acgunot-? characters 	
	my $txt = shift;
	# the order of translation is important because we don't want to inadvertantly translate the c,a,n in "span" and "class"

	$txt =~ s/(a+)/<spxx xxss =\"jj\">$&<\/spxx>/gi;
	$txt =~ s/(c+)/<spxx xxss =\"jp\">$&<\/spxx>/gi;
	$txt =~ s/(n+)/<spxx xxss =\"jl\">$&<\/spxx>/gi;

	$txt =~ s/(o+)/<span class =\"jo\">$&<\/span>/gi;

	$txt =~ s/(u+)/<span class =\"ju\">$&<\/span>/gi;
	$txt =~ s/(g+)/<span class =\"jg\">$&<\/span>/gi;
	$txt =~ s/(t+)/<span class =\"jt\">$&<\/span>/gi;
	$txt =~ s/(-+)/<span class =\"jd\">$&<\/span>/gi;
	$txt =~ s/(\?+)/<span class =\"jq\">$&<\/span>/gi;

	# some cleanup
	$txt =~ s/spxx/span/gi;
	$txt =~ s/xxss/class/gi;
	
	$txt =~ s/jj/ja/gi;
	$txt =~ s/jp/jc/gi;
	$txt =~ s/jl/jn/gi;

	return $txt;
}


sub str_colchars () { # return the a string containing the nth character of all elements 
	my ($pos) = shift;
    	my @array = @_;
	my $str;
	
	foreach my $i (@array) { 
		$str .=  substr($i, $pos, 1); 
	}	
	return $str;
}

sub str_totaluniquechars () { # pilfered from "http://www.newts.org/~troc/perl/uniqchar.perl"; returns the number of unique characters in a string *I
	my ($str) = shift;
	use integer;
	my %c;
	for (split //, $str) {$c{$_}=undef} scalar keys %c 
}

sub str_rnafragmentcomp () { # returns an array [0..7] containing the composition of the given string, for AUGCN-?[other]
	my $fragment = shift;
	$fragment = uc($fragment);
	
	my @tmparray;	
		my $i=0;
		foreach my $letter (@rna[0..6]) {
			$tmparray[$i] = ($fragment =~ s/([$letter])//g);
			$tmparray[$i] = 0 unless defined $tmparray[$i];
			$i++			
		}
		$tmparray[7]= length($fragment);
		$tmparray[8] = $fragment;
	return @tmparray;
}

sub str_unique_chars () { # code modified from O'Reliey Perl Cookbook, Chapter 1.5
	my %seen = ();
	my $str = shift;
	while ($str  =~ /(.)/g) {
    		$seen{$1}++;
	}
	return sort(keys %seen);
}

sub str_unique_chars_C () { # returns ord() instead of chars
	my %seen = ();		
	my $str = shift;
	while ($str =~ /(.)/g) {
    		$seen{ord($1)}++;
	}
	return sort(keys %seen);
}


sub str_out_mx_corrected () { # returns a string consiting of a legal alphabet of characters given a specific output type  *I
	my (	
		$str,	# input string
		$mode	# translation "mode"
	) = @_;
	
	if ($mode eq "mrbayes") {
		$str =~ s/X/N/g;	
	}

	elsif ($mode eq "phase") { 
		$str =~ s/(?)([^AUGCN\-\?])/N/g;
	}
	
	elsif ($mode eq "poy") {
	}

	return $str;
}









## ------------ unused probably unworking code ------------------- ##

# find invariant blocks (code currently cloned in INNASE OUT --------------------------------------fix to a call to here! NOT TESTED as of 03/04/04
sub unused_invariant_blocks () { # returns a list of invariant blocks for a range of data/taxa
	#requires
	# start taxon, end taxon, start block, end block 
	my ($starttaxa, $endtaxa, $startblk, $endblk) = @_;
	my ($invariant, $curchar, $invarchar, @invars);
		for (my $blk = $startblk; $blk < $endblk+1; $b++) {
			if ($nbd[$blk][3] == 1) {
				$invariant = 1;
				 for (my $p=0; $p < $nbd[$b][0]; $p++) { 		# loop through all positions in the string
					if ($invariant == 0 ) {	last }
					
					for (my $t = $starttaxa; $t < $endtaxa+1; $t++) {	 	# loop through the taxa - remember smaller subsets may have all dashes, while larger
						
						$mx[$t][$blk][0] =~ /.{$p}(.)/;
						$curchar = $1;
						
						if ($t == $starttaxa) {
							$invarchar  = $curchar; 
						}
						unless ($curchar eq $invarchar)  {$invariant = 0; last}	# length - 0- long text names; 1- short text names;har)  {$invariant = 0; last}
					}
				}
				if ($invariant == 1) {push @invars, $blk}
			}
		}
		return @invars
}

sub unused_useful_stemstr () { #returns a string containing only the bonding positions based on a character mask
	#requires
	my (
		$stem,	# stem string
	       	$mask	# mask string	
	) = @_;
	my $out;
	 
	if (length $stem != length $mask) {return "error";}
	#ERROR_CHK - check that mask consists of legal characters only
	
		for (my $i = 0 ; $i < length ($stem) + 1; $i++) {
			print $i;
			$mask =~ /.{$i}(.)/;
			unless  ($1 eq '.') {
				$stem =~ /.{$i}(.)/;
				$out .= $1;
			}	
		}	
	return $out;
}

sub continue () {
	print " data may have been exluded, if that's ok type 'yes' to continue ";
	my $answ;
	while (<STDIN>) {
		$answ = $_;
	}	
	die print "bye!" if not $answ eq "yes";
}

sub check() {  # NOT FUNCTIONING (YET)
	# checks your data verus the $fileroot.original file and reports errors
	
	my (	$input,
		$startblk,
		$endblk
	) = @_;
	my (%orig, @tmp, $name, $currow);
	
	# load "original" data
		open (INFILES, "data/$input->{Modelroot}.original") || die print "cannot open  - ";
			while (<INFILES>) {
				chomp;
				@tmp = split;
				$orig{$tmp[0]} = chomp($tmp[1]);
			}
		close (INFILES);
	
		my %revters = reverse %ters; # we can do this because the original is all nice indexed
				
		foreach $name (keys %orig) {
				$currow = "";
			for (my $blk = $startblk; $blk < $endblk+1; $blk++) {
				$currow .= $mx[$revters{$name}][$blk][0];
			}
			$currow =~ s/-//g;
			print $name, " ", $currow, "\n";
		}
}

1;

__END__
