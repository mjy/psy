
# das splitter - spilts to interleave a matrix based on the column ranges defined in split.txt 


#load the matrix into a matrix
open (infiles, "matrix.txt") || die print "cannot open datafile - ";
	@mtrx = <infiles>;
close (infiles);

#load the split data
#data must be in the format 0-123 125-111 etc... will die miserable if columns aren't set right, won't check for column problems

open (infiles, "split.txt") || die print "cannot open datafile - ";
	@splits = <infiles>;
close (infiles);

#load the original matrix into an array (might as well be associative
$j=0;

while (@mtrx[$j]) {
	@cr = (split (/\s+/, @mtrx[$j]));
	chomp($cr[1]);
	$mxhsh{$cr[0]}=$cr[1];
	$mxord{$j}=$cr[0];
	#print "MATRIX: $cr[0] $cr[1]";
	$j++;
}

#print "J : $j\n";
#

@spltz = split (/\s+/, @splits[0]); #takes the first line of data only, you can add whatever below a splits file

open(writehndl, ">split_matrix.txt") || die "couldn't open file to output too - ";

$b=0;


foreach $blk (@spltz) {
	chomp $blk;
	#print "$blk";
	@cols = split /\D/, $blk;
	print writehndl ("\t[--- BLOCK $b ($blk) of length ");
	print writehndl ($cols[1]-$cols[0]);
	print writehndl (" ---]\n");
	
	for ($k=0;$k<$j;$k++) {
		
		#very crude, don't really need the double hash bit but....
			#print "$mxord{$k}\t";
			#print substr($mxhsh{$mxord{$k}},$cols[0],($cols[1]-$cols[0])+1);
			#print "\n";
		
			print writehndl ("$mxord{$k}\t");
			print writehndl (substr($mxhsh{$mxord{$k}},$cols[0],($cols[1]-$cols[0])+1));
			print writehndl ("\n");
			
	}
		print writehndl "\n\n" ;
$b++;		
}


close (writehndl);
