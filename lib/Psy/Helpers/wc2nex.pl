# M.Yoder 2003
# 
# ugly use of Perl for sure, but it gets the job done
# 
# a perl script to convert winclada format .tre files (name taxa [either format], ordered)
# usage is perl wc2nex.pl <infile> <outfile> 
#	pass the full path to both infile and outfile

# open and read-in the infile
$infile = shift (@ARGV);
$outfile = shift (@ARGV);

print "$infile\n";

open (infiles, $infile) || die print "cannot open $infile - ";
	@inf = <infiles>;
close (infiles);


$l = 0; # the line number in the infile

# get the taxon names - place them in a hash
$taxanext=0;
$t=1;

	while (@inf[$l] !~ /;/) {
			
			#	chomp $inf[$l];
			
	
			if (@inf[$l] !~ /clados|torder/) {
				chomp ($inf[$l]);
				$taxa{$t}=$inf[$l];
				$t++;

				print "$inf[$l]\n";	
			}

			$l++;
	}
	
#print @taxa{5};

$l++; #escape the ; line and the clados- line
$l++; 

$tr=1; # tree counter

#get the trees - places them in a hash
	while (@inf[$l] !~ /proc \/;/) {
				
		print "\ngot here\n";
			#chomp $inf[$l];
			if ($inf[$l] =~ /tread 'trees from NONA'/) {
#				print $inf[$l];
#				split  /(tread 'trees from NONA')/ , $inf[$l];
				$inf[$l] =~ s/tread 'trees from NONA'//;
				$inf[$l] =~ s/\s\*\s//;
				$inf[$l] =~ s/;//;
				#print $inf[$l] ;
				chomp ($inf[$l]);
				
				$trees{$tr}=$inf[$l];
				print @trees{$tr};
				
				$tr++;
			}

			else {
				
				$inf[$l] =~ s/;//;
				$inf[$l] =~ s/\s\*\s//;
			
				chomp ($inf[$l]);
				
				$trees{$tr}=$inf[$l];
				#print $trees{$tr};
				
				$tr++;
			}

			$l++;
	}


	#converting a tree
#	$tmptree = $trees{5};
#	for ($j=1; $j<$t; $j++) { 
		#replace all the taxon names with their hash number
#		$tmptree =~ s/$taxa{$j}/$j/;

		
#		$tmptree =~ s/\s\s/,/;
#		$tmptree =~ s/\)\(/\),\(/;
#		$tmptree =~ s/\s\(/,\(/;
#		$tmptree =~ s/\s//;
#	}

#	print $tmptree;
	
	
open(whndl, ">$outfile") || die "couldn't open file to output too - ";

			# write the outfile
			
			print whndl "#Nexus\n\n";
			
			print whndl "Begin trees; \[Treefile converted from winclada\.tre format\] \n\n";
			print whndl "\t\tTranslate\n";
			
			for ($k=1; $k<$t; $k++) { # print the translate list
				print whndl "\t\t\t $k $taxa{$k},\n";
			}

			print whndl ";\n";
	#converting a tree
	
	
	for ($i=1; $i<$tr; $i++) { #loop through the trees
	
	$tmptree = $trees{$i}; #select the current tree
	
		for ($j=1; $j<$t; $j++) { 
			#replace all the taxon names with their hash number
			$tmptree =~ s/$taxa{$j}/$j/;
		
			#insert the needed commas
			$tmptree =~ s/\s\s/,/;
			$tmptree =~ s/\)\(/\),\(/;
			$tmptree =~ s/\s\(/,\(/;
			
			$tmptree =~ s/(\d\()/\$1,/;
			$tmptree =~ s/\s//;
			
		}

			print whndl "\ntree WNCL_$i = \[&U\] ";
			print whndl $tmptree;
			print whndl ";\n"
	}
	
	print whndl ("End;\n\n\n");
	
close (writehndl);

