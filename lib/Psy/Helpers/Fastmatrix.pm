package Psy::Helpers::Fastmatrix;

use warnings;
use strict;

use Data::Dumper;
use Psy::Strings::Strings;
use Psy::Helpers::Gbname;  # required for only certain operations.

# use vars qw(@ISA @EXPORTER @EXPORT_OK);
# @EXPORT_OK = qw(new seq numCols loopTer);

=head1 NAME

Psy::Helpers::Fastmatrix 

=head1 VERSION

version 0.02 - Matt Yoder

=cut

our $VERSION = '0.02';

=head1 SYNOPSIS

Not used alone, but rather through its subclasses (e.g. fastfasta).
A class for loading and performing operations on simple sequence files (aligned or not).  
This is fast as in quickly derived code- definitely NOT fast as in speed.  

=cut

=head2 new

=cut


sub new {
	my $type = shift;
	my %params = @_;
  	my $self = {};
#	$self->{Interleaves} = {};	
	bless $self, $type;
	
	$self->_init(%params);

	# find the longest terminal here to buffer for 

	
    return $self; 
}

sub _init{
	# prototype for matrix loader, likeley pointless, but used in tests
	1;
}


=head2 seq

Accessor for seq given a terminal label.

=cut


sub seq { 
	my $self = shift;
	my $t = shift;
	if (@_) { $self->{mx}->{$t} = shift }
	return $self->{mx}->{$t}
}


=head2 seqLen

Return length of the seq.

=cut


sub seqLen { 
	my $self = shift;
	my $t = shift;
	my $s = $self->seq($t);
	$s =~ s/[-xX]//g; # assumes ? and nN are actually nucleotides, and that xX are not (hmmer goodness)
	return length($s)
}



=head2 numCols

NOT fast

=cut


sub numCols {
	my $self = shift;
	my @r = $self->loopTer;
	return length ( $self->seq($r[0])) 
}


=head2 allSeqs

Return an array of all the sequences

=cut


sub allSeqs {
	my $self = shift;
	my @r;
	foreach my $t ($self->loopTer) {
		push @r, $self->seq($t);
	}
	@r;	
}

=head2 aTerminal 

Returns (at "hash random" NOT TRUE RANDOM) a legal terminal label.

=cut

sub aTerminalLabel {
	my $self = shift;
	my @foo = $self->loopTer;
	$foo[0];
}


=head2 loopTer

Returns arrary of terminal labels

=cut


sub loopTer { #
	my $self = shift;
	return sort keys %{$self->{mx}};
}


=head2 loopBySeqLength

Returns arrary of terminal labels in order of sequence length

=cut


sub loopBySeqLength { 
	my $self = shift;
 	my %l; 
 	map { $l{$_} = $self->seqLen($_) } ($self->loopTer);
	return reverse sort { $l{$a} <=> $l{$b}  } (keys %l); # longest to shortest
}


=head2 loopGBnameHigher

Returns arrary of terminal labels in order of sequence length

=cut


sub loopGBnameHigher {
	my $self = shift;
 	my %n;
	my $gbname = Psy::Helpers::Gbname->new;
 	map { $n{$_} = $gbname->higher('-from_name' => 1, '-seq_obj' => $self->seqObj('-accessor' => $_)); 	 } ($self->loopTer);
	return  sort { $n{$a} cmp $n{$b}} (keys %n); # longest to shortest
}

=head2 loopByGBnameThenSeqLen

Returns arrary of terminal labels in order of GBname then sequence length

=cut


sub loopByGBnameThenSeqLen {
	my $self = shift;
	
 	my %l;
	my %n;
	
	my $gbname = Psy::Helpers::Gbname->new();
	map { $l{$_} = $self->seqLen($_) } ($self->loopTer);	
	map { $n{$_} = $gbname->higher('-from_name' => 1, '-seq_obj' => $self->seqObj('-accessor' => $_)); 	 } ($self->loopTer);

	# print Dumper($gbname->higher('-from_name' => 1, '-seq_obj' => $self->seqObj('-accessor' => $_));
	
	return reverse sort {$n{$b} cmp $n{$a} || $l{$a} <=> $l{$b}} ($self->loopTer); # by higher, then longest to shortest	
}


=head2 seqObj

Returns a pseudo-object for use in gbname (mimics the necessary bioperl seqobj structure), a cheat to use with gbname

=cut


sub seqObj { 
	my $self = shift;
	my %raw_params = @_;
	my %default_params = (
		'-classification' => ['unknown'],
		'-common_name' => 'no_common_name'
	);
	my %params = (%default_params, %raw_params);
	# require -accessor
	
	my $obj = {};
	$obj->{species}->{_classification} = $params{'-classification'};
	$obj->{primary_seq}->{display_id} =  $params{'-accessor'};
	$obj->{species}->{_common_name} =  $params{'-common_name'};

	bless $obj, "gbname"; # it just matters that its an object, not what kind
	return $obj;
}


=head2 longestTer

Returns an integer length of the longest terminal label

=cut


sub longestTer {
	my $self = shift;
	my $lt = 0;
	foreach my $t ($self->loopTer) {
		(length $t > $lt) && ($lt = length $t)
	}
	return $lt;
}


=head2 deleteFixedGapPositions

Deletes fixed positions of kind -char

=cut


sub deleteFixedGapPositions { 
	my $self = shift;
	my @seqs = $self->allSeqs;

	my @new_seqs = &delFixedPositions('-strings' => \@seqs, '-char' => '-');
	my $i=0;
	map {$self->{mx}->{$_} = $new_seqs[$i]; $i++} ($self->loopTer);
	
	1;
}


=head2 findOne

Searches all seqs for the passed regex, returning to STDOUT whether or not it was found.

=cut


sub findOne {
	my $self = shift;
	my $regex = shift;;
	my $i = 0;
	my $j = 0;
	foreach my $t ($self->loopTer) {
		if ($self->{mx}->{$t} =~ /$regex/ig) {
			print "found: $t\n";
			$i++;
		}
		else {
			print "notfound: $t\n";
			$j++;
		}
	}
	print "f: $i nf: $j\n";
}


=head2 findMany

Returns true if all members of -regexes are found in the seq for -terminal.

=cut


sub findMany {
	my $self = shift;
	my %params =  @_; # reference to a hash

	foreach my $r (@{$params{'-regexes'}}) {
		#	print "$r\n";
		if ($self->{mx}->{$params{'-terminal'}} =~ /$r/ig) { next }	
		else { return 0	}
	}
	return 1;
}


=head2 findManyOr

Returns true if any memeber of -regexes is found in the seq for  -terminal.

=cut


sub findManyOr {
	my $self = shift;
	my %params =  @_; # reference to a hash

	foreach my $r (@{$params{'-regexes'}}) {
		if ($self->{mx}->{$params{'-terminal'}} =~ /$r/ig) { return 1 }	
	}
	return 0;
}


=head2 boundOligo

Returns largest sequence bound by left most r1 and rightmost r2

=cut


sub boundOligo {
	my $self = shift;
	my ($t, $l, $r) = @_; # terminal, left and right regexes
	
	$self->{mx}->{$t} =~ /$l/i;
	my $left_pos = $-[0];
	
	$self->{mx}->{$t} =~ /$r/i;
	my $right_pos = $+[0];

	print "$left_pos, $right_pos len: ", $right_pos - $left_pos, "\n";
	return substr($self->{mx}->{$t}, $left_pos,  $right_pos - $left_pos);
	
}


=head2 trimOligo

As bound Oligo but trims left and right of r1,r2

=cut


sub trimOligo {
	my $self = shift;
	my ($t, $l, $r) = @_; # terminal, left and right regexes
	
	my $oligo = $self->{mx}->{$t};
	
	print "t: $t ";
	
	if ($self->{mx}->{$t} =~ /$l/i) {	
		print "l: $-[0] ";
		$oligo = substr($oligo, $-[0], length ($self->{mx}->{$t}) - $-[0] )
	}

	if ($oligo =~ /$r/i) {
		$oligo = substr($oligo, 0,  $+[0]) ; # +?
		print "r: $+[0]";
	}
	
	print "\n";
	
	return $oligo;
}


=head2 out_boundMx

Outputs a matrix that has only those sequences that are bound by -left_bound and -right_bound.

=cut


sub out_boundMx {
	my $self = shift;
	my %params = @_;
	my $i=0;

	my $length = 0;
	my $min = 999999;
	my $max = 0;
	my ($minter, $maxter);
	
	#	print $params{'-outfile'};
	my $of = $params{'-outfile'};

	open (OUT, ">$of" ) or die "can't open outputfile\n";

	my @re = ( $params{'-left_bound'}, $params{'-right_bound'} );
	print "\n";
	print join "\n", @re;
	print "\n";

	foreach my $t ($self->loopTer) {
		if ($self->findMany( '-terminal' => $t,	'-regexes' => \@re) == 1) {
			print OUT ">$t\n";
			my $s = $self->trimOligo($t, $params{'-left_bound'}, $params{'-right_bound'});
			my $l = length($s);
			 if ($l < $min) { $min = $l; $minter = $t}
			 if ($l > $max) { $max = $l; $maxter = $t}
			$length += $l;
			print OUT $s;
			print OUT "\n\n";
			$i++;
		}	
	}
	close OUT;
	print "wrote: $i sequences\n";
	print "avg. length: ", $length / $i, "\n";
	print "min: $min ($minter) \n";
	print "max: $max ($maxter) \n";
}


=head2 out_Fasta

Exports a fasta file.  Makes calls to genbank (but not by default)

=cut

sub out_Fasta {
	my $self = shift;
	my %params = @_;

		
	# print $params{'-outfile'};
	$params{'-outfile'} ||= "out.fas";
	$params{'-use_gbname'} ||= 0;
	$params{'-strip'} ||= 1; # a quick strip method for \? and -

	my $of = $params{'-outfile'};
	
	my $gbname = Psy::Helpers::Gbname->new;

	open (OUT, ">$of" ) or die "can't open outputfile $of\n";
	my $i=0;
	foreach my $t ($self->loopTer) {

		print OUT ">", ( $params{'-use_gbname'} == 1 ? $gbname->name('-in' => $t) : $t) , "\n";
		my $s = $self->seq($t);
		$s =~ s/[\?]/\-/g;
		print OUT $s  ;
		print OUT "\n\n";
		print ".";
		$i++;	
	}
	close OUT;
	print "wrote: $i sequences to Fasta matrix: $of \n";
}

=head2 out_splitRegex

Exports a oneLine formatted matrix with lines split and gapped to the found regexes.

=cut

sub out_splitRegex { 
	my $self = shift;
	my %params = @_;

	my $of = $params{'-outfile'} || "splitregex";
	$params{'-regexes'} || die "no regexes to split with";

	open (OUT, ">$of" ) or die "can't open outputfile\n";
	my $i=0;
	foreach my $t ($self->loopTer) {
		print OUT "$t\t\t\t", uc( &regexSplit('-str' => $self->seq($t), %params)) , "\n";
		$i++;	
	}
	close OUT;
	print "wrote: $i to one-line matrix\n";
}


=head2 out_splitRegexAligned

As out_splitRegex, but places found blocks on top of each other.

=cut


sub out_splitRegexAligned {
	my $self = shift;
	my %params = @_;

	my $of = $params{'-outfile'} || "splitRegexAligned";
	$params{'-regexes'} || die "no regexes to split with";
	$params{'-insert'} ||= '  ';
	
	my $mx; # $mx->{$ter}->{piece_index} = piece
	# use a 1b / 1 piece_index schema, where the b is the piece before 

	my @rgx = @{$params{'-regexes'}};	

	# build up a hash
	foreach my $t ($self->loopTer) {
		my $n = 0;
		my $s =	$self->seq($t);

		my $matched = 0;
		for (my $r = 0; $r < $#rgx +1; $r ++) {
			if ( $s =~ /  
						(.{$n})      # omit everything previously searched
						(.*?)        # allow for an 'insert', of shortest possible length prior to the regex in question
						($rgx[$r])   # find the regex
						(.*)		 # find everything after the regex
						/xi      
				) { 
				# print "$r [$1] [$2] [$3] [$4] $n \n";
			
				$mx->{$r}->{$t} = $3;
				$mx->{$r."b"}->{$t} = $2;
				$n = length ($1) + length($2) + length ($3);
				$matched = 1;
			}
		}

		#	($matched == 0) && (print OUT $self->seq($t));
		
		$mx->{'cap'}->{$t} = substr($s, $n, length($s)-$n);
		# print "\n";
	}

	my $lt = $self->longestTer;
	open (OUT, ">$of" ) or die "can't open outputfile\n";
	my $i=0;

	# print "rgx: ", $#rgx, "\n";
	# print $rgx[0];
	# print $rgx[1];

    # precalculate the max length of the regions
	my $ml;
	for (my $r = 0; $r < $#rgx + 1; $r++) {
		# find the longest $r AND $r.b ('b'efore)
		$ml->{$r.'b'} = &lengthLongestString( map {$mx->{$r.'b'}->{$_}} $self->loopTer );
		$ml->{$r} = &lengthLongestString( map {$mx->{$r}->{$_}} $self->loopTer );
	}
	
	# compose the gapped dequence strings  (this is ugly, but it works)
	my $data;
	foreach my $t ($self->loopTer) {
		my $hit = 0;
		my $str = '';	
		for (my $r = 0; $r < $#rgx + 1 ; $r++) {
			if ( $mx->{$r.'b'}->{$t} ) {
				$str .= &padleft($mx->{$r.'b'}->{$t}, " ", $ml->{$r.'b'}); 
				$hit = 1;
			}
			else {
				$str .=  &padleft('', " ", $ml->{$r.'b'}) ;
			}

			if ( $mx->{$r}->{$t}) {
				$str .= $params{'-insert'};
				$str .=  &padleft($mx->{$r}->{$t}, " ", $ml->{$r}); 
				$hit = 1;
			}
			else {
				$str .=  &padleft('', " ", $ml->{$r}) ;
				$str .= "  ";
			}	
		}

		if ($hit == 0) {
			$data->{$t} = $self->seq($t);
		}
		else {			
			$data->{$t} = $str;
			$data->{$t} .=  "  ".$mx->{'cap'}->{$t};
		}
	}

	
	# precalculate the max length of the overal strings (we want to pad them for nexus legalness
	my $maxlen =  &lengthLongestString(values %{$data});
	print "max length = $maxlen\n";
	
	# print Dumper($ml);

	# actually write the matrix
	foreach my $t ($self->loopTer) {
		print OUT &padright($t, " ", $lt + 2);
		$data->{$t} =~ s/\s/\-/gi;
		print OUT &padright($data->{$t}, "-", $maxlen );
		print OUT "\n";
	}
	
	close OUT;

	#print Dumper($mx);
	# print "wrote: $i to one-line matrix\n";
	return 1;
}


=head2 out_oneLine

Exports a oneLine formatted matrix

=cut


sub out_oneLine {
	my $self = shift;
	my %params = @_;

	my $of = $params{'-outfile'} || "oneLine.fas";

	open (OUT, ">$of" ) or die "can't open outputfile\n";
	my $i=0;
	foreach my $t ($self->loopTer) {
		print OUT "$t\t\t\t", uc( $self->seq($t)) , "\n";
		$i++;	
	}
	close OUT;
	print "wrote: $i to one-line matrix\n";
}


=head2 out_oneLineRegexCentered

Outputs a matrix that with sequence data centered (if found) over -regex.

=cut


sub out_oneLineRegexCentered { # realy quick and dirty
	my $self = shift;
	my %params = @_;

	my ($ll, $lt, $ls) = (0) x 3; # longest left, longest terminal, longest full sequence (pre- gap + seq len)

	$lt = $self->longestTer;
	
	$params{'-regex'} ||= 'GTTGAATCT';

	my @found;
	
	# pre-loop to find those that pass, and determine length longest
	foreach my $t ($self->loopTer) {
		if ($self->seq($t) =~ /$params{'-regex'}/i) {
			($-[0] > $ll) && ($ll = $-[0]);

			( ( $ll + length $self->seq($t) ) > $ls) && ($ls = $ll + length $self->seq($t));
			
			print "$ll\t$ls\n";
			unshift @found, ($t); # dump all found first
		}
		else {	
			(length ( $self->seq($t) ) > $ls) && ( $ls = length ($self->seq($t)));
			push @found, $t; 	
		}
	}
	
	my $of = $params{'-outfile'} || "oneLine.fas";

	open (OUT, ">$of" ) or die "can't open outputfile\n";
	my $i=0;

	foreach my $t (@found) {
		print OUT &padright($t, " ", $lt + 1); 
		
		my $fp = 0; # found point
		if ($self->seq($t) =~ /$params{'-regex'}/i) { # gap the front half
			$fp = $-[0];
			print OUT "-" x ($ll - $fp); # gap or space
		}	
	
		my $seq = $self->seq($t);
		$seq =~ s/[tT]/U/gi;
		
		print OUT uc( $seq);
		print "$ll, $fp, $ls,", length($self->seq($t)) ,"\n";
		print OUT "-" x ($ls - length($self->seq($t)) - ($ll - $fp) ), "\n";
		
		$i++;	
	}
	close OUT;
	print "wrote: $i to one-line matrix\n";
}



=head2 unfound

Compares a file containing a list of genbank Accessions to the current dataset, and writes a file of missing sequences

=cut



sub unfound { 
	my $self = shift;
	my %params = @_;
	
	my @cmp; # terminals to find
	
	# read in, doesn't assume anything about taxa, just appends all seq data to the preceeding '>taxon'
	open (IN, $params{'-compare_file'} ) || die print "cannot open $params{'-compare_file'} \n";
		while (<IN>) {
			my $row = $_;
			chomp $row;
			
			$row =~ s/^(\s+)//; # strip space from begining of the row, Justin Case
			# print "[$row] ";
			
			if ($row =~ /_/) { 
					push @cmp, substr($row, 0, $+[0]-1);	
			}
			else {
				push @cmp, $row;
			}		
		}
	close IN;

	# print join "|\n", @cmp;
	
	my $i = 0;	
	my $gbname = Psy::Helpers::Gbname->new;

	open (OUT, ">$params{'-outfile'}") || die print "cannot open $params{'-outfile'} to write too\n";
		foreach my $t ($self->loopTer) {
			if (not (grep $_ eq $t, @cmp) ) { 
				print OUT ">";
				print OUT $gbname->name('-in' => $t);
				print OUT "\n";
				print OUT $self->seq($t);
				print OUT "\n\n";
				$i++;			
				print "not found: $t\n";
			}
			else {
				print "found: $t\n";
			}
		}
	close OUT;
	print "wrote: $i unfound sequences\n";
}


=head2 clustalAlign

Aligns the matrix, then reloads it to memory.  Must have clustalw/x in your path as 'clustalw'

Usage:
	$foo->clustalAlign(
		'-clustal_params' => ''
	)
	
Tested on clustalW 1.81.

=cut

sub clustalAlign {
	my $self = shift;
	my %raw_params = @_;
	my %default_params = ('-clustal_params' => '');
	my %params = (%default_params, %raw_params);

	# output the file
	$self->out_Fasta('-outfile' => 'clustal.in');
	
	# run the alignment
	
	# hack for differences in output type from PC to Mac - *nix systems might be different yes
	#
	# still some weirdness here, perhaps based on clustal/clustalw??
	
	my $end = 'msf';  # other possibility is 'aln'
	$end = 'msf' if $^O eq 'darwin'; 

	print "aligning...";
    `clustalw clustal.in -output=gcg $params{'-clustal_params'} -align `; # used to be -infile=
	print "done\n";
	
	use Psy::Helpers::Fastmsf;
	# read it back in
	my $in = Psy::Helpers::Fastmsf->new('-infile' => "clustal.$end"); # curiously -output=gcg returns two different file endings clustal.msf vs clustal.aln for PC- ARG
	
	$self->merge($in);
	# clean up
	unlink(qw/clustal.dnd clustal.in/);
	unlink("clustal.$end");
}

=head2 muscleAlign


Aligns the matrix, then reloads it to memory.  Must have muscle in your path as 'muscle'

Usage:
	$foo->clustalAlign(
		'-clustal_params' => ''
	)

=cut

sub muscleAlign {
	my $self = shift;
	my %raw_params = @_;
	my %default_params = ('-muscle_params' => '');
	my %params = (%default_params, %raw_params);

	# output the file
	$self->out_Fasta('-outfile' => 'muscle.in');
	
	print "aligning...";
    `muscle -in muscle.in -out muscle.out $params{'-muscle_params'} -stable`;
	print "done\n";
	
	# read it back in
	use Psy::Helpers::Fastfasta;
	my $in = Psy::Helpers::Fastfasta->new('-infile' => "muscle.out");

	#	print Dumper($in);
	
	$self->merge($in);

	# Muscle eliminates null rows, while clustalw apparently keeps them in.  This creates problems for all '-' or '?' only lines.
	# DANGEROUS (but necessary) assume all lines of data are aligned to the same size, and trim (this is tricky) all others (i.e. ones muscle didn't touch) to the same size.
	# If I'm right all others should be completely ? or - or other characters that muscle doesn't reconize
	
	my $l = length($in->seq($in->aTerminalLabel));
	for my $t ($self->loopTer) {
		my $s = $self->seq($t);
		
		if ($l < length($s)) { # remainder too long
			$self->seq($t, substr($s, 0, $l));
		}
		elsif ($l > length($s)) { # aligned too short
			$self->seq($t, $s.('-' x ($l - length($s)) )    )
		}
		
		#	Sequence (original): ----, new:----, Taxon index: 10, aligned length: 5.
		
		die "Sequence (original): $s, new:", $self->seq($t), ", Taxon index: $t, aligned length: $l.\n" if (length($self->seq($t)) != $l);
	}
	
	# clean up
	unlink(qw/muscle.out muscle.in/);
}







=head2 merge

Merges two 'Fast' based matrices.  Terminal data is OVERWRITTEN by matrix two's data if a terminal is present in both matrices.

=cut

sub merge {
	my $self = shift;
	my $merge = shift;
	%{$self->{mx}} = (%{$self->{mx}}, %{$merge->{mx}});
}


__DATA__
[C][CT][T][G][ACT][GT][AG][AG][A][CT][CT][CGT][ACG][AC][AG][A][AGT]?[AGT][ACT]?[C][G][AG][AGT][ACGT][AG][T]?[G]?[AG][ACG][AG][AG][AGT][A][ACT]?[ACT]?[AC][AC][AGT]?[CT][AGT][T]?[ACGT][ACT][AGT][AG]
d2_5_prime	CTYTGAAKAGAGAGTYMAANAGTRCGTGAAACYKYYYRGRG
d2_3_prime	DCCCGTCTTGAAACACGGACCAAGRAGT
d3_5_prime	GCGTACACGTTGGGACCCGAAAGATGGTGAACTATGCCTGGT
d2_internal1	CCYGAGAAACCYRWWWGNTCGAA
d2_foo	RRRRGWYYCWDCGKCWRYD
d2_sir_foo	RRRRRWYYYWDYRKYWRYD
d2_sir_foo2	RRRRRWYYYWDYRKYWRY
d2_sir_foo	RRRRRWYYYWDYRKY
d2_sir_foo	RRRRRWYYYWDYR
little_foo	WGRRGAGATTCAWCGTTCARYR
foo_foo	CCGTCTTGAAACACGG
d3_3_prime2	GCGTACACGTTGG
blk_156_d3_3p_core [AG][ACT][AG][AG][ACGT]?[AGT][ACGT][AG][G][]?[CT][G][C][AC][AC][ACT]

=head1 AUTHOR

'Matt, C<< <m{j}yoder@{tee}[aye](em)(you).domain4unis> >>

=head1 BUGS

Please report any bugs or feature requests to
<bug-psy-helpers-fastmatrix@rt.cpan.org>, or through the web interface at
<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Psy>.
I will be notified, and then you'll automatically be notified of progress on
your bug as I make changes.

=head1 ACKNOWLEDGEMENTS

=head1 COPYRIGHT & LICENSE

Copyright 2005 'Matt, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1; # End of Psy::Helpers::Fastmatrix
