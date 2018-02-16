
# re-order script ver 1.0 by Matt Yoder/2003
# a silly script to reorder a matrix of data based on an infile named "order.txt" and a matrix named "matrix.txt"
# probably could do this in one line if you really new Perl.

# load the order into a hash
open (infiles, "order.txt") || die print "cannot open datafile - ";
	@nameord =  <infiles>;
close (infiles);

$i=0;
#put all the names into a hash
while (@nameord[$i]) {
	chomp ($nameord[$i]);
	print $nameord[$i];
	$nhsh{$i}=$nameord[$i] ;
	$i++;
}

foreach $infile (@ARGV) {


#load the matrix into a matrix
open (infiles, $infile) || die print "cannot open $infile - ";
	@mtrx = <infiles>;
close (infiles);


$j=0;
#put all the data into another hash there is only names and data, two solid blocks allowed for each row!!!!
while (@mtrx[$j]) {
	@cr = (split (/\t+/, @mtrx[$j]));
	chomp($cr[1]);
	$mxhsh{$cr[0]}=$cr[1];
#	print "MATRIX: $cr[0] $cr[1]";
	$j++;
}

open(writehndl, ">$infile.reordered.txt") || die "couldn't open output file $infile.reordered.txt ";

for ($k=0;$k<$i;$k++) {
	#print writehndl" *$nhsh{$k}* ";
	print writehndl ($nhsh{$k}."\t".$mxhsh{$nhsh{$k}}."\n");
}
close (writehndl);
}
