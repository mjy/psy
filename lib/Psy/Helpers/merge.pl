# merge textfiles, takes first arguement as the filename to output too, remaining files as input

# meh- just use 'cat'

$outfile=shift (@ARGV);

die "file to merge to ($outfile) alread exists!\n" if (-z "$outfile" );

$usage = 'perl merge.pl <infile1> <infile2> ... <infilen>';

open(writehndl, ">$outfile") || die "couldn't open $outfile file to output too - \n".$usage;

	foreach $infile (@ARGV) { 

		open (infiles, $infile) || die print "cannot open $infile\n ";
	
			while (<infiles>) {
		
				print writehndl ($_);
			}
	
		print writehndl ("\n\n\n");

		close (infiles);
	}

close (writehndl);
