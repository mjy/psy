# re-order script ver 1.0 by Matt Yoder/2003
# a silly script to reorder a matrix of data based on an infile named "order.txt" and a matrix named "matrix.txt"
# probably could do this in one line if you really knew Perl.

# load the order into a hash
open (infiles, "order.txt") || die print "cannot open ordering file called order.txt- ";
	@nameord =  <infiles>;
close (infiles);

#load the matrix into a matrix
open (infiles, "matrix.txt") || die print "cannot open matrix file matrix.txt- ";
	@mtrx = <infiles>;
close (infiles);

$i=0;

#put all the names into a hash
while (@nameord[$i]) {
	chomp ($nameord[$i]);
	print $nameord[$i];
	$nhsh{$i}=$nameord[$i] ;
	$i++;
}

$j=0;

#put all the data into another hash there is only names and data, two solid blocks allowed for each row!!!!
while (@mtrx[$j]) {
	@cr = (split (/\s+/, @mtrx[$j]));
	chomp($cr[1]);
	$mxhsh{$cr[0]}=$cr[1];
#	print "MATRIX: $cr[0] $cr[1]";
	$j++;
}

open(writehndl, ">reordered.txt") || die "couldn't open file to output too - ";

for ($k=0;$k<$i;$k++) {
	#print writehndl" *$nhsh{$k}* ";
	print writehndl ($nhsh{$k});
	print writehndl " " x (45 - length($nhsh{$k}));
	print writehndl ($mxhsh{$nhsh{$k}}."\n");
}
close (writehndl);

