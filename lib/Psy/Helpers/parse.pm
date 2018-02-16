
# multiple packages for 's'equence 'p'arsing or 's'imple 'parsing' 

use strict;
use warnings;
use Data::Dumper;

use FindBin qw($Bin); 		# sets $Bin to the root directory
our $ROOT_DIR = $Bin;
our $VERSION = '0.02';

=head1 NAME

Raw seq parser

=head1 VERSION

Version 0.02

=head1 SYNOPSIS

Prolly not for public consumption.

Code to recursively parse a directory of .seq files and return some a file of nicely formatted concensus data that can be manually edited, then extracted.
Makes *a lot* of assumptions about how the seq files are named, but the general pattern is:
letter##_code.primer_lane.file_ending (e.g. g02_mjy234.D2_01.seq). 

You must include your primer pairs in DATA (bottom of this file), they are used in matches and can be abused to limit the output generated.

In its simplest form just do this:

	use sp;
	my $foo = sp->new('-dir' => '.');
	$foo->dumpAll;
	
After editing your output this returns a fasta formatted version of the data:
	
	use postsp;
	&postsp::extractToMx('out4.text');


=cut 


=head1 postsp

Post processing of sp output 

=cut


package postsp; 

sub extract {
	my $infile = shift;
	open (OUT, ">extracted.txt") or die "can't write to 'extracted.txt'";
	
	open (IN, $infile) or die "can't open $infile";
SEQ:		while (<IN>) {
			if (m/(^m:)\s(\S*)\s(\S*)\s/) {	
				my $t=$2;
				while (<IN>) {
					if (m/(^e:)\s(.*)/) {
						print OUT ">$t\n";
						print OUT $2, "\n\n";
						next SEQ;
					}
					elsif (m/^\n/) {
						next  SEQ;
					}
			}
		}
		elsif (m/^\n/) {
			next SEQ
		}
	}
	close (IN);
	close (OUT);
}


=head2 extractToMx

A hack to generate a table that is easily mappable to import to mx.  The primer hash values are gene_ids.

=cut


sub extractToMx {
	# otu gene notes attempt_complete concensus 
	my $infile = shift;

	my %gh = ( 
		'CF-CR' => '7',
		'D3F-D5R' => '8',
		'18SF-18SR' => '9',
		'18S8F-18S9R' => '10',
		'18SF1-18S4R' => '11',
		# '18SA2-18S9R' => '12', changed to 10 after import
		'28S1F-28SCR' => '13',
		'28S6F-28S7R' => '14',
		'28S8F1-28S8R' => '15',
		'28S8F2-28S10R' => '16',
		'D2QF-' =>		'17'
	);
	
	open (OUT, ">extracted4_mx.txt") or die "can't write to 'extracted.txt'";
	
	open (IN, $infile) or die "can't open $infile";
SEQ:		while (<IN>) {
			if (m/(^m:)\s([mM][Yy])(\d*)_(\d*_\w*)_(\S*)\s(\S*)\s/) {	
				my $t=$3;
				my $n=$4;
				my $g=$5;
				my $d;
				($_ =~ /DONE/) ? ($d = 'true') : ($d = 'false');
				while (<IN>) {
					if (m/(^e:)\s(.*)/) {
						print OUT "$t\t[series 1: $n]\t$gh{$g}\t$d\t";
						my $tmp = $2;
						$tmp =~ s/\W//g;
						print OUT $tmp , "\n";
						next SEQ;
					}
					elsif (m/^\n/) {
						next  SEQ;
					}
				}
			}
			elsif (m/^\n/) {
				next SEQ
			}
	}
	close (IN);
	close (OUT);
}

1;


=head1 sp

sp (sequence parser), for rapid pre-editing formatting of .seq files

=cut



package sp; # comes second so you can read __DATA__

# use Bio::SeqIO;
use File::Find;
use Data::Dumper;

my $INFILES; # required, need a global for file::find

sub new {
	my $type = shift;
	my %params = @_;	
	my $self = {};

	bless $self, $type;

	$self->_init(%params);

	return $self;                 
}

sub _init {
	my $self= shift;
	my %params = @_;

	while (<DATA>) {
		chomp;
		$_ eq "" and next;
		my @pr = split;
		$self->{'primers'}->{"$pr[0]"} = $pr[1];
	}

	$params{'-dir'} && $self->loadDir(%params);
	return 1;
}


=head2 loadDir

Loads (recursively) the seq files.

=cut


sub loadDir {
	my $self = shift;
	my %raw_params = @_;
	my %default_params = ('-format' => 'abi'); ## does nothing?
	my %params = (%raw_params, %default_params);

	&find(\&_wanted, $params{'-dir'}); 
	#print "\nfound: \n";
	#rint join "\n", @INFILES;


	my $l = 0; # number of files loaded	
	foreach my $file (keys %{$INFILES} ) {
		$l++;
		my @meta = $INFILES->{$file} =~ /(\w)(\d+)_(.+)\.(.+)_(\d+)/; # assumes - well- you figure it out

			# $meta[0] - 'A', 'B' etc.
			# $meta[1] - '02', '03' 
 			# $meta[2] - 'MY154', 'MY142'
			# $meta[3] - '18SA2'
			# $meta[3] - '016', '015', '003'

#		print join "\t", @meta; 
		
		print "reading: $file\n";
		
		if ((not defined  $meta[2]) or (not defined  $meta[1]) or (not defined  $meta[0])) {die "can't parse filename $file\n"}
			
		$meta[3] =~ /\./ && chop $meta[3]; # some files have periods- weird
		$meta[2] =~ tr/my/MY/;
		
		my $s = "$meta[2]_$meta[1]_$meta[0]_$meta[3]_$l";
		print "$s\n";
		
		$self->{seqs}->{$s} = {
						'match' => $meta[2].$self->threePrime($meta[3]), # .$meta[1] 
						'sample_pair' => $meta[0],
						'sample' => $meta[1], 
						'voucher_code' => $meta[2],
						'primer' => $meta[3],
						'filename' => $file
		};

		# print "match: ", $self->match("$meta[2]_$meta[1]_$meta[0]"), "\n\n";

		$self->{seqs}->{$s}->{seq} = &readSeq($file);
	}
	
	$self->{total_seqs} = $l;

	# when bioperl-ext work use this, for now grab the .seq files
	# foreach my $file (@INFILES) {
	#	$self->{_accession} = Bio::SeqIO->new(%params, -file => $file);
	#	print $self->{_accession}->id;
	#	$self->{_accession}++;
	# } 

	print "loaded $l seqs (stored:", $#{[keys %{$self->{seqs}}]} + 1, " , a missmatch implies overwritten data!)\n\n";

}

=head2 readSeq

Reads an individual sequence

## change this to the electropherogram once BioPerl-ext installs correctly

=cut


sub readSeq { 
	my $file = shift;
	open (IN, "$file") || die "can't open infile $file";
		my @seq = <IN>;
	close (IN);

		my $seq = join "", @seq;
		$seq =~ s/\s//gi; # strip whitespace and newlines
		$seq =~ s/\n//g; 

	return $seq
}

=head2 concensus

Returns a concensus sequence
- assumes str1 and str2 are equal lenght
- for columns opposite whitespace concensus returns lower case
- for non-matching columns "N" is entered
- otherwise the concensus letter is entered
	
=cut


sub concensus {
	my ($str1, $str2) = @_;
	(length ($str1) != length $str2) && return undef; # leave as is
	my $str;

	for (my $p=0; $p < length $str1; $p ++) {
		$str1 =~  /.{$p}(.)/;
		my $l1 = $1 ;
		$str2 =~  /.{$p}(.)/;
		my $l2 = $1;
		
		if ($l1 eq $l2) { $str .= uc($l1) }
		elsif ($l1 eq ' ') {$str .= lc($l2) }
		elsif ($l2 eq ' ') {$str .= lc($l1) }
		else {$str .= 'N'};
	}
	return $str;
}


=head2 concensusCue

Symbolizes a concensus string for easier visual-reference

=cut


sub concensusCue { 
	my $str = shift;
	$str =~ s/[^nN]/ /g;
	return $str;
}


=head2 aeStr

Align and equalize string: adds space offset (calc offset with &align) to begining of string 1, then caps then end of the then longer string with space so that both are equal length

=cut


sub aeStr { 
	# ! foo ? bar : stuff  requires the bracketing as below
	my ($offset, $str1, $str2) = @_;
	my $offsetstr = ' ' x abs(int($offset));

	# print "[ $str1 ] ", length ($str1), " \n [ $str2 ] ", , length ($str2), "\n";
	
	length($str1)<1 && die;
	length($str2)<1 && die;

	$offset < 0 ? $str2 = ($offsetstr.$str2) : ($str1 = $offsetstr.$str1);

	# print "[ $str1 ] ", length ($str1), " \n [ $str2 ] ", , length ($str2), "\n";
	
	# cap 
	my $dif = length($str1) - length($str2);
	
	# print "dif: $dif\n";

	($dif < 0) ? ($str1 .= ' ' x int(abs($dif))) : ($str2 .= ' ' x int(abs($dif)));

	# print "[ $str1 ] ", length ($str1), " \n [ $str2 ] ", , length ($str2), "\n";
	# if (length($str1) != length($str2)) {print "WARNING, not equal length!!!\n"}
	
	return ($str1, $str2);
}


=head2 align

A simple very rough "alignment" algorithm for near identical sequences of large length (e.g. forward and backwards reads). Returns the offest of seq1 to seq2; no gaps are inserted

=cut


sub align { 
	my %raw_params = @_;
	my %default_params = ( 
			'-min_inside' => 45,	# start a minimum of this number of characters in
			'-max_attempts' => 100,  # try to find a subseq with no N's and match it a maximum of this many times
			'-subseq_min' => 13,	# make the subseq a minimum of n nucs long
			'-subseq_max' => 35 	# make the subseq a max of m nucs long
		);

	my %params = (%default_params, %raw_params);	

	(length $params{'-seq1'} < $params{'-subseq_max'}) && return undef; # die "seq1 too short";
	(length $params{'-seq2'} < $params{'-subseq_max'}) && return undef; # die "seq2 too short";

	$params{'-min_inside'} > length $params{'-seq1'} && die "seq1 too short";

	# $params{'-min_inside'} > length $params{'-seq1'} && die "seq2 too short";

	#print $params{'-seq1'}, "\n", $params{'-seq2'} , "\n";

	my $found = 0;	
	my $ml; # max seq length
	my $i = 0;
	my $offset;
	length ($params{'-seq1'}) > length ($params{'-seq2'}) ? $ml = length ($params{'-seq2'})  : $ml = length ($params{'-seq1'}) ;

	do { # (min - max + 1) + max
		my $l = int( rand( $params{'-subseq_min'} - $params{'-subseq_max'} + 1 ) ) + $params{'-subseq_max'};
		
		# a length > min_inside and less than length min_inside and within range of the second sequence
		my $sp = int( rand( $params{'-min_inside'} - (length($params{'-seq1'}) - $params{'-min_inside'}) + 1 ) ) + (length($params{'-seq1'}) - $params{'-min_inside'}); # this used to be -seq2!
				
		# print "sp: $sp $l\n";
		my $subseq = substr($params{'-seq1'}, $sp, $l);

		print "\t aligning: $sp\t$l\t$subseq\t";
		
		if ($subseq =~ /[nN\?]/gi) {print " failed"; $i++; }
		else {
			$params{'-seq2'} =~ m/$subseq/g;
			if ($-[0]) {
				print " [";
				print $-[0], " ", $+[0];
				print "] ";
				$offset = ($-[0] - $sp);
				print " passed";
				$found = 1
			}
			else {
				print " failed";
			}
		}
		print "\n";
		$i++;
	} until ($found != 0 || $i> $params{'-max_attempts'});

	return $offset ;
}


=head2 pair

Pairs a forward and reverse sequence by elements of their filenames AND the <DATA> that describes primer pairs

=cut


sub pair { 
	# use a number of checks though it should be straightforward
	my $self = shift;
	my $seqin = shift; 

	foreach my $seq (keys %{$self->{seqs}}) {	
		if ( $self->match($seqin) eq $self->voucher_code($seq).$self->primer($seq)) { # MATCHES ON voucher + primer pair  
			print "matched: $seqin and $seq \n";		
			return $seq ; #$self->match($seq);
		}
	}
	return undef;
}

=head2  ### accessor methods ###

=head2 match

The match string, like 'MY154'

=cut

sub match {
	my $self = shift;
	my $seq = shift;
	if (@_) {$self->{seqs}->{$seq}->{match} = uc(shift)};
	if (defined $self->{seqs}->{$seq}->{match}) {return $self->{seqs}->{$seq}->{match}}
	else {return 0};
}


=head2 sample_pair

like 'A', 'B' etc.

=cut


sub sample_pair { # return the value for a sequence pair 
	my $self = shift;
	my $seq=shift;
	if (@_) {$self->{seqs}->{$seq}->{'sample_pair'} = uc(shift)};
	return $self->{seqs}->{$seq}->{'sample_pair'};
}

=head2 _pair_letter

Used to match. NO LONGER USED.

Assume its all nice ascii- this might not ultimately be the case (A/H wrapping).

=cut

sub _pair_letter { 
	my $l = shift;
	not defined $l and die;
	($l eq 'H') && (return 'A'); ## wrap around hack
	return chr(ord($l)+1) 
}

=head2 primer

Primer for a seq

=cut

sub primer {
	my $self = shift;
	my $seq=shift;
	if (@_) {$self->{seqs}->{$seq}->{primer} = uc(shift)};
	return $self->{seqs}->{$seq}->{primer};
}


=head2 threePrime

Accessor for the 3' for a 5'-3' primer pair, as read from <DATA>.

=cut

sub threePrime { 
	my $self = shift;
	my $pr = shift;
	return $self->{'primers'}->{$pr} || 'unmatched';
}

=head2 voucher_code

Voucher code 

=cut

sub voucher_code { 
	my $self = shift;
	my $seq = shift;
	if (@_) {$self->{seqs}->{$seq}->{voucher_code} = uc(shift)};
	return $self->{seqs}->{$seq}->{voucher_code};
}


=head2 seq

Seq (nucleotides) 

=cut


sub seq { 
	my $self = shift;
	my $seq = shift;
	if (@_) {$self->{seqs}->{$seq}->{seq} = shift};
	return $self->{seqs}->{$seq}->{seq};
}

sub filename { # filename accesor 
	my $self = shift;
	my $seqid = shift;
	return $self->{seqs}->{$seqid}->{filename};
}

### end accessor methods ###


sub rc { # reverse compliment a string - ## move to STRINGS library ultimately
	my $str = shift;	

	$str =~ s/[Aa]/1/ig;
	$str =~ s/[uU]/2/ig;
	$str =~ s/[gG]/3/ig;
	$str =~ s/[cC]/4/ig;
	$str =~ s/[tT]/5/ig;

	$str =~ s/1/T/ig;
	$str =~ s/2/A/ig;
	$str =~ s/3/C/ig;
	$str =~ s/4/G/ig;
	$str =~ s/5/A/ig;

	return scalar(reverse($str));
}

sub dumpOne {
my $self = shift;
	
}

sub dumpAll { 
	my $self = shift;
	my %raw_params = @_;
	my %default_params = ('-out_file' => 'out.text');
	my %params = (%default_params, %raw_params);

	open (OUT, ">".$params{'-out_file'} ) || die "can't open outfile!";

	my @done;
	my %unpaired;
	my %paired;
	my %nomatch;
	my @unalignable;
	
	my ($i, $j) = (0) x 2;
	my %foo;
	foreach my $seq (sort keys %{$self->{seqs}} ) {
			$foo{$self->primer($seq)} = '';	
		my $pair = $self->pair($seq);

		if ($pair) { # a pair for this seq is found
			
			$paired{$seq} = undef;
			$paired{$pair} = undef;
			
			# keep track of samples parsed
			push @done, $seq;
			push @done, $pair;
			
			my $seq1 = $self->seq($seq);
			my $seq2 = &rc($self->seq($pair));
			
			print OUT "m: ", $seq, "_", $self->primer($seq), "-" , $self->threePrime($self->primer($seq)), " ",  $self->filename($seq), " " , $self->filename($pair), "\n";
			
			#print "$seq\t$pair\n";

			my $offset;
			$offset = &align('-seq1' => $self->seq($seq), '-seq2' => $seq2);
	
			if (not defined $offset) {
				print OUT "UNALIGNABLE\n\n";
				push @unalignable, ("$seq-$pair  ".$self->primer($seq)."-".$self->primer($pair)) ;
				$j++;	
				next
			};

			# print "{ offset: $offset }  \n";
		
			my @aligned = &aeStr($offset, $seq1, $seq2); 
			# print "[ lengths: ", length($aligned[0]),  " " , length($aligned[1]) , " ]\n\n";	
		 	#	print join "\n", @aligned;

			my $con = &concensus($aligned[0], $aligned[1]);
			
			if ($con eq '') {
					print OUT "ERROR generating concensus seq\n\n";
					next
			}
			
			print OUT "c:  $con\n";		
			
			print OUT "1o: ", $aligned[0], "\n"; # original data
			print OUT "2o: ", $aligned[1], "\n";	
			
			print OUT "1e: ", $aligned[0], "\n"; 			# edit these lines
			print OUT "2e: ", $aligned[1], "\n";	
			print OUT "e:  $con\n";							# final edit
			print OUT "n:  ", &concensusCue($con),  "\n";	# visual cue 

			print OUT "\n\n";
			$i++;
		}
		else { # no pair
			$unpaired{$seq} = undef;
			$nomatch{$self->primer($seq)} = undef;
		}
	};

#	my %p = %unpaired;
	delete $unpaired{$_} for keys %paired;
	
	close OUT;

	print "\nprimers encountered:\n";
	print join "\n", keys %foo;
	
	print "\n\nfinished\n";
	print "wrote: $i pairs\n";
	print "unalignable: $j pairs\n";
	print "read but untreated (ESTIMATED): ", $self->{total_seqs} - ($i+$j) *2, "\n";

	print "\n\nUNPAIRED: (", $#{[keys %unpaired]}, ")\n";
	print "$_ ".$self->primer($_)."\n" for sort keys %unpaired;

	print "\n\nUNALIGNABLE:\n";
	print join "\n", @unalignable;
	
	return 1;
}


sub _wanted { # look at only those files that are abi/ab1 (now currently .seq)
	if ($_ =~ /\.seq/) { $INFILES->{$File::Find::name} = $_ }; 
}


1;





# primer pairs are currently stored in DATA, carefull not to include duplicate indices

__DATA__

18SF	18SR
18S8F	18S9R
D3F	D5R
18SF1	18S4R
18SA2	18S9R
CF	CR
28S1F	28SCR
28S6F	28S7R
28S8F1	28S8R
28S8F2	28S10R
D2QF	CR



