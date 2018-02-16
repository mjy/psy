package matrix;

use vars qw(@ISA @EXPORTER @EXPORT_OK $VERSION);

use strict;
use warnings;
use Data::Dumper;
use Carp;

use Psy::Psy qw(generatedByHeader version rootpath projectname starttime $ROOT_DIR $PSY_BASE_DIR $PSY_ROOT_DIR);
 @ISA = qw(Psy);

 use Psy::Matrix::Format;
 use Psy::Matrix::Interleave;
 use Psy::Matrix::Block;
 use Psy::Matrix::Column;
 use Psy::Matrix::Slice;
 use Psy::Matrix::Terminal;
 use Psy::Matrix::Structure;
 use Psy::Io::Io;
 use Psy::Dna::Alphabet;
 use Psy::Strings::Strings;

=head1 NAME

Psy::Psy

=head1 VERSION

Version 0.01

=cut


$VERSION = '0.01';


=head1 SYNOPSIS


my $foo = matrix->new(
						'-filename' => 'Evanioidea.0.04.Nex',
						'-path' => '../../data/evan/',
						'-number_terminals' => '59'  # n-1
					);

#print "\n";
#print Dumper($foo);
#print "\n";
#print $foo->{blocks}->version();

use Psy::Matrix;

my $foo = matrix->new();
    ...

=head1 OBJECT

=cut


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
	
	$self->starttime;

	%params = $self->_verify(%params);
	croak 'error initializing matrix object' if $params{'-error'}; # ERROR HANDLING (does nothing)

	$self->{blocks} = {};
	$self->{interleaves} = {};

	$self->_interleaveAdd("original");
	
	$self->nters($params{'-number_terminals'});
	$self->{mx} = [];
	
	$params{'-matrix_label'} ? $self->label($params{'-matrix_label'}) : $self->label('no_label'); 
	$self->{terminals} = {};
	$self->{structures} = {};
	
	$self->_structureAdd(%params, '-structure_name' => "original");
	 
	$self->{_terminals_accession} = 0;
	$self->{_block_accession} = 0; # the index of the next new block
	
	$self->_load(%params) if $params{'-matrix_file'};

	$self->_initOrigSlice(); # keep here 	

	$self->_calcBlockPositions; 

    # run some error checking, this essentially unparsed blocks
    my @mx_check = $self->sliceData();
    for (my $i = 0; $i < $#mx_check + 1; $i++) {
        die "!! Uninitialized 'cell', died on piece $i bounded with these data: ", $mx_check[$i - 1], ", missing data, ", $mx_check[$i + 1], "\n" if not $mx_check[$i]
    }

    $self->{alphabet} = Psy::Dna::Alphabet->new(
										'-type' => 'custom',
									   	'-alphabet' => join "",  &uniqueChars(join "", $self->sliceData)  ## this is weird
								   	);	
								
    print "Parsed alphabet: ", $self->{alphabet}->alphabet, "\n";

	$self->version($VERSION);
	return 1;
}
	
sub _verify {
	my $self=shift;
	my %params = @_;
		#	foreach my $p (keys %params) { # check matrix parameters here
		#	print "$p ";
		# }
	return %params;
}

=head2 _load

Internal.  Used to load Psy/jRNA formats and Stockholm format.  

Needs an over-ride ultimately, likely to set to a single filetype, or better yet to a parser/lexer.

=cut


sub _load {
	my $self = shift;
	my %raw_params = @_;
	my %default_params = (
		'-data_offset' => 17 # should be the same as out_mx:outcolstart - used to grab (dataoffset) charcharters from description blocks
	);
	
	my %params = (%default_params, %raw_params);
	
	open (IN, ($params{"-path"}.$params{"-matrix_file"})) || die ("can't open ", $params{'-path'}.$params{'-matrix_file'} );  # ERROR HANDLING
		my @mtrx = <IN>; # slurp the file, might be bad for very large matrices
	close (IN);

	print "\n\n\nfile loaded: ", $params{"-path"}.$params{"-matrix_file"}, " \n";
	my $ff = fileformat->new();
	
	if ($params{'-fileformat'}) { 
		$self->inputformat($params{'-fileformat'})
	}
	else {
		$self->inputformat($ff->detect(join "",@mtrx)); 
	}
			
	die ("can't detect format of ", $params{"-filename"}, ", you might want to explictly define it with -fileformat.\n") if not $self->inputformat();  # ERROR HANDLING

	print "reading input as type: ", $self->inputformat() ,"\n";

	# escape to non Stockholm/jRNA formats here.  A kludge at present, but not many multi-paritioned formats are envisioned/out there.  See also ARC format for a likely candidate.
	if ( $self->inputformat eq 'PHASE_simulate') { $self->_load_PHASE_simulate(%params, '-mx' => \@mtrx); return 1 }; 

	my (
		$i,           # input file line position
		$j,           # pre-data input file line position
		$cur_intlv,   # current interleave working within
		$cur_ter,     # index of the current terminal being worked on				
		$z_blk_pos,   # zeroed block position at the start of a given interleave (keep Psy) 
		$blk_tot,	  # sum total of blocks at the start of processing an interleave, increments on completion of interleave processing	
	) = (0) x 6;      # current total number of blocks
	
	my ($tmpline1, $tmpline2, $tmpline3);
	my $piece;
	my $brackets;	
	my (%descblk1, %descblk2, %descblk3);
	
	my $id = $ff->data_start($self->inputformat()); # simplify to regex below?
	print "filetype ID: $id\n";

	# $self->_blockAdd; # assume there is at least one datablock ## this might be wrong
		
	while ($mtrx[$j] !~ /\s*$id\s*/i ) { # scan to find the line that will identify the file format 
		$j++;
		 $j > 1000 && die "\nCan't guess fileformat from the first 1000 lines.\n"
	}
	
		# error checking - move to detect?
		#print "$j\n";
		#	if ( ($mtrx[$j] =~ /ntax/i) and ( $self->inputformat() eq "nexus") ) {
			#		print "HERE";
			#$mtrx[$j] =~ /(ntax)(\s*=\s*)(\d*)/i;
			#	die "ntax=$3 does not match your set number of taxa [", $self->nters()  ,"]\n" unless ($self->nters() == $3-1); 
			# } 	

	print "found $id at line: $j\n";

	$i=$j+1; #mark the first line of data
	
	my $meta;
	my @cur_blks;
	
	# some position/state flags 
	my ($havedesc, $havedata, $haverowmeta) = (0) x 3;

	my $end = $ff->data_end($self->inputformat());
	print "Data end at: $end\n";

    print "Parsing interleave: ";
	while ($mtrx[$i] !~ /$end/) { # find the first line after which we can assume matrix data
      # match a data line - this line has a taxon/sequence in it		
		if ($mtrx[$i] !~ /(^(\[)|^(#=)|^(\s*#=)|^(\s+\[)|^([\n\r])|^(\s+[\n\r]))/)  { 				
			# print "data line: ", substr( $mtrx[$i], 0, 5), "\n";
		
			# grab the column/start/end/length etc info for each new block
			if ($haverowmeta == 0) { 
				## there are simpler ways to grab sizes (split), but we also want to grab equally sized chunks in equal positions in other rows!	
				my $tmprow = $mtrx[$i];	 # needs this scope?			
				chomp $tmprow;
				$meta = {};
				$meta = lineWordMeta('-line' => $tmprow); # simple - but remember the taxon is grabbed too so subtract its length when doing _addBlock
		
				# word_start 1 is the line the data actually starts on	
				
				#print Dumper($meta);
				
				@cur_blks = sort { $a <=> $b} keys %{$meta} ; ## change this scope and use it to gather descriptive data
				shift @cur_blks; # the first block is a taxon label
	
				# print "keys :", keys %{$meta}, ": \n";
				
				foreach my $k (@cur_blks) {
				
					# print $meta->{$k}->{'fc'}; 	
					# print Dumper ($meta->{$k});					
					
					my $brk;
					if ($meta->{$k}->{'fc'} eq "[") {
						$brk = 1;  		
					}
				   	else {
						$brk = 0;
					}
					
					my $cb = $self->_blockAdd( # this is just the meta data, not the actual matrix
									  '-bracketed'	    => $brk,
									  '-blk_length'     => $brk == 1 ? $meta->{$k}->{word_length} - 2 : $meta->{$k}->{word_length}  ## kludge - we need a global validation for block length etc					  
					   				);
					
					# debug
					# print  "blk_tot $blk_tot [$cb]\n";
					# $self->{blocks}->{$cb}->debugDescribe;						
					$self->{interleaves}->{'original'}->assign($cur_intlv, $cb); # assign the block to the original interleave
				}
				$haverowmeta = 1;	# have the data for this intlv block
			} 

			# grab and store the actual data
			my @cr = split (/\s+/, $mtrx[$i]); # split the blocks into @cr
			
			$cur_intlv == 0 ? $self->_terminalAdd('-label' => shift @cr) :	shift @cr ; # grab or trash a taxon

			my $piece2;
			$z_blk_pos = 0; # must retain larger scope
			foreach my $dblk (@cr) {
				# store the raw data only  
				if ($dblk =~ m/\[/) {	# bracketed data
					$dblk =~ /(^\[)(.+)(\])/;
					$piece2 = $2;
					die "ERROR reading input matrix at line ",  ($i + 1) ,". Check for whitespace where its not supposed to be." if not $piece2;
					chomp($piece2) 
				}
				else {		# not bracketed data
					$piece2 = $dblk;
					chomp($piece2)
				}

				# kludge check for a couple characters that may have been causing problems
                if ($piece2 =~ /[\,\.\!\:\;]/) {
                    die "!! within interleave position $z_blk_pos ($piece2) has non-alphabetical data in it at or around line ", ( $i + $j ) , "\n";
                }


				$self->cell('-ter'=> $cur_ter, '-blk'=> ($blk_tot + $z_blk_pos), '-data' => $piece2 );  	
				$z_blk_pos++;
			}	
		$cur_ter++;	 
	} # end parse a data line
	
	# a non-data line
	elsif ($havedesc == 0) { # grab the meta data in for the given interleave
		my $temp =  quotemeta ($ff->mask($self->inputformat()));

		# Nexus format - ## needs changing to a more tokenized version, but only after careful review
		if ( $mtrx[$i] =~ /^(\[Block)|^(\s+\[Block)/ ) { # test to see if this is the [block ** REWORK THIS LOGIC TO MAKE IT UNIVERSAL 	 
		   # should ultimately have variable number of meta rows, right now output functions are closely tied to these description blocks
			
			$tmpline1 = $mtrx[$i];
			$tmpline2 = $mtrx[($i+1)];
			$tmpline3 = $mtrx[($i+2)];
				
			if (length($tmpline1 or $tmpline2 or $tmpline3) < $params{'-data_offset'}) {die "problem with input file [Block associated regions"};
			
			chomp $tmpline1;
			chomp $tmpline2;
			chomp $tmpline3;
				
			$descblk1{$cur_intlv} = $tmpline1;
			$descblk2{$cur_intlv} = $tmpline2; # assume the second and third blocks are descriptive
			$descblk3{$cur_intlv} = $tmpline3;

			my $lb1 = substr ($tmpline1, 1, $params{'-data_offset'});
				$lb1 =~ tr/ //s; 
			my $lb2 = substr ($tmpline2, 1, $params{'-data_offset'});
				$lb2 =~ tr/ //s; 
			my $lb3 = substr ($tmpline3, 1, $params{'-data_offset'});
				$lb1 =~ tr/ //s; 
			
			#print Dumper ($self->{interleaves}->{'original'});
			$self->{interleaves}->{'original'}->label(
														'-interleave' => $cur_intlv,	
														'-label' => '0',
														'-value'=> $lb1
													);

			$self->{interleaves}->{'original'}->label(
														'-interleave' => $cur_intlv,	
														'-label' => '1',
														'-value'=> $lb2 
													);

			$self->{interleaves}->{'original'}->label(
														'-interleave' => $cur_intlv,	
														'-label' => '2',
														'-value'=> $lb3 
													);
			$havedesc = 1;
		}
			
		# Stockholm format
		elsif ( $mtrx[$i] =~ /$temp/) { # is this the mask block?
			$descblk3{$cur_intlv} = $mtrx[$i];
			#	print " $cur_intlv, $descblk3{$cur_intlv} \n ";		
			$havedesc = 1;
		}
		# don't put havedesc =1 here- it needs to fail too	
	} # end grab meta for interleaves

	# do some checking for presence of data and metadata
	if ($cur_ter > $self->nters() ) {
		$havedata = 1;
	}

	# capture "descriptions" and mask
	if ( ($havedata == 1) and ($havedesc == 1)) { 
	    print " $cur_intlv ";

		my $q = 1; # start minus the taxon row!

		foreach my $blk ($self->{interleaves}->{'original'}->blocks($cur_intlv)) {
            # print "$descblk1{$cur_intlv}\n";
			my @d1 = &wordsFromLineMeta($descblk1{$cur_intlv}, $meta); # bracketed blocks
            my @d2 = &wordsFromLineMeta($descblk2{$cur_intlv}, $meta);
			my @d3 = &wordsFromLineMeta($descblk3{$cur_intlv}, $meta);

			#print "-- @d3 --\n\n ";
			# some warnings
			# (not @d1 or not @d2 or not @d3) && (die "\nDied on parsing matrix.  You may have incorrectly identified the number of taxa in your matrix, or your matrix is not formatted correctly.\n\n");
			
            # print "\n";
            # print "d1: ", join " ] [", @d1;
			
            # print join " ", @d1;
            # print "\n";
			
			# don't edit whitespace in 4 and 5!!
            # print "mask: $blk $q $d1[$q] $d3[$q] \n";;

			if ( $self->inputformat() =~ /nexus/i) {

				$self->blk($blk)->label('-label_name' => 'label1', '-label' => $d1[$q]);
				$self->blk($blk)->label('-label_name' => 'label2', '-label' => $d2[$q]);
				
				$self->{structures}->{'original'}->label('-blk'=>$blk,'-label_name'=>'label1', '-label' => $d1[$q] ); #  substr ($descblk1{$cur_intlv}, $meta[$blk][1], ($meta[$blk][2]-$meta[$blk][1])+1) );
				$self->{structures}->{'original'}->label('-blk'=>$blk,'-label_name'=>'label2', '-label' => $d2[$q] ); #  substr ($descblk2{$cur_intlv}, $meta[$blk][1], ($meta[$blk][2]-$meta[$blk][1])+1) );
			}


			if ($self->{blocks}->{$blk}->bracketed == 0) {		
				$self->{structures}->{'original'}->mask('-blk'=> $blk, '-mask' => $d3[$q] ); # substr ($descblk3{$cur_intlv}, $meta[$blk][1], ($meta[$blk][2]-$meta[$blk][1])+1) );	
			}
			else { # assume its not bracketed (there are no other options) - in this case have to slice brackets from mask
				$self->{structures}->{'original'}->mask('-blk'=> $blk, '-mask' => substr( $d3[$q],1,-1) ); 
			}
			
			# error check		
			if (($self->{structures}->{'original'}->mask('-blk' => $blk) =~ /\s/ ) and ($self->{blocks}->{$blk}->bracketed == 0) ){
				print "\n!! WARNING: mask '", $self->{structures}->{'original'}->mask('-blk' => $blk), "' of block $blk contains whitespace (at line ", ($i - $j + 1), ")\n";
                print "!! If this does not appear to be the case then an error exists before this point !! \n";
                print "!! If the error is mask related sequence data for that mask include: '", $self->cell('-ter'=> 1, '-blk'=> $blk ), "'";
                print "\n\n"
             }
				
			$q++;	
		}
		$havedesc = 1;
		$meta = ();
	}

	if ((($havedata == 1) and ($havedesc == 1)) 
		or (($havedata ==1) and ($cur_ter > $self->nters()))) {	
		$cur_ter = 0;
		$haverowmeta = 0;
		$havedesc = 0;
		$havedata = 0;
		$blk_tot += $z_blk_pos;	
		$cur_intlv++;
	}		

	$i++; #next line in the matrix	
	
	# print "  blktot: $blk_tot\n"; 
								
	if ($i > $#mtrx) {die "incorrectly terminating infile?\n";} 	# some basic error checking for incorrectly terminating nexus files

	} # done parsing matrix
	
	print "\n\nfile read\n";
	print "block range: 0 - ", $self->totalBlocks-1, "\n";
	
	# close MXERR;

	# if (-z "err/mx_err.txt") {
	#	print "No detectable errors in matrix input parsing.\n";
	# }
	# else {
	# die "ERRORS PRESENT! See: err/mx_err.txt";
	# }
		
	return 1;
}


=head2 _load_PHASE_simulate

NOT WORKING!

Loads a file generated by PHASE 2.0(beta)'s simulate package.  

At this point it handles only a MIXED model, where the second model is RNA with a repeating () structure (the first parition is essentially ignored. 

You must include a -struct_size of 2 or larger.  For example if your simulate structure is () then struct_size = 2. (()) = 4 etc. 

=cut

sub _load_PHASE_simulate {
	my $self = shift;
	my %params = @_;

	my ($tot_tax, $tot_nucs,  $l) = (0) x 3;
	my $t = 1;
	my @mtrx = @{$params{'-mx'}};
	
	(not $params{'-struct_size'}) && die 'no -struct_size passed to _load_PHASE_simulate'; 
	 $params{'-struct_size'} < 2 && die '-struct_size to small';
	
	my $hl = $params{'-struct_size'} / 2; # block length
	my ($tot_blks, $tot_len) = (0) x 2;
	while ( $#mtrx > $l) {
			my $row = $mtrx[$l];
			chomp $row;
	
			if ($row eq '#structure') {
				print "STRUCTURE!\n"; $l++;	next; 
			}
		 	if ( $row =~ /^#/) {$l++; next;} # skip # lines
			
			$row =~ s/^(\s+)//; # strip space from begining of the row, Justin Case
			
			if ($row =~ /^(\d+)\s(\d+)/ && ($tot_tax == 0)) { # grab the totals
				$tot_tax = $1; $tot_nucs = $2;
				#	print "$l $1 $2\n";
			}	

			if ($row =~ /^([os]\d+)/ && ($tot_tax != 0)) { # seq data	
				$self->_terminalAdd('-label' => $t); # no interleaved data here 	
				# print "TAXA\n";
				my $r1 = $mtrx[$l+1];
				chomp $r1;
				$self->cell('-ter' => $t, '-blk' => 0, '-data' => $r1);
				my $r2 = $mtrx[$l+2];
				chomp $r2;
				# print length "$r2\n\n";
			 if ( $tot_blks == 0) {	 # bad to put this here, as its hit lots, but allows for more automation
					$tot_blks = length($r2) / $hl;
					$tot_len = length($r2);
					}
				map {
					# print "b ", ($_ + 1), ," ", $_ * $hl, " " , substr($r2, $_ * $hl , $hl), "\n"; 
					$self->cell('-ter' => $t, '-blk' => $_ + 1, '-data' => substr($r2, $_ * $hl , $hl));
				} (0..(length($r2) / $hl)); 		
							
				$t++;
				$l += 2;
			}
			$l++;
	}

	# error check to see if we have a valid read, could be done earlier but then would be continuous rechecking?
	
	die "mask size ($hl) for blocks ($tot_blks) not compatible with seq length ($tot_len)\n" if (($tot_blks ) * $hl) != $tot_len; 
	
	# make structure
	$self->_structureAdd('-structure_name' => 'original');	

	$self->_blockAdd( '-bracketed'  => 1); # the very first block is the non-RNA data


	
	map { 
		$self->{structures}->{'original'}->mask('-blk'=> $_, '-mask' => "(" x $hl  );	
		$self->{structures}->{'original'}->mask('-blk'=> ($_ + 1), '-mask' => ')' x $hl);
		$self->{structures}->{'original'}->helix($_, $_ + 1);
		$self->_blockAdd( '-bracketed'  => 0,'-blk_length'     => $hl);
		$self->_blockAdd( '-bracketed'  => 0,'-blk_length'     => $hl);
	} map ({$_ * 2} 1..$tot_blks/2);

	# @odd=map {$_ * 2 +1} 0..($n-1)/2

	# print Dumper($self->{structures}->{'original'});
	return 1
}


### Methods that return objects (references to)


=head2 ter

Returns the terminal _object_ for the passed index.

	Smx->ter(2);

=cut

sub ter { 
	my $self = shift;
	my $t = shift;
	return  $self->{terminals}->{$t};	
}


=head2 blk 

Returns the block meta-data _object_ (not index) for a given block ## COMPLETE to full accession? (low priority).
	Smx->blk(2);

=cut


sub blk { # 
	my $self = shift;
	my $blk = shift;
	return $self->{blocks}->{$blk}
}

### Methods that return sequence data as strings

=head2 cell

Returns the sequence data for a -ter, -col pair.

	Smx->cell('-ter' => 1, '-col' => 2);


=cut


sub cell {
	my $self = shift;
	my %params = @_; 
	
	die 'no -ter to cell' if not defined $params{'-ter'};
	die 'no -blk to cell' if not defined $params{'-blk'};
	
	if (defined $params{'-data'}) {
		$self->{mx}->[$params{'-ter'}][$params{'-blk'}] = $params{'-data'};
		return 1;
	}
	else {
		return $self->{mx}->[$params{'-ter'}][$params{'-blk'}]; 
	}
	return 0;
}


=head2 rawCell

Returns sequence data minus the regex -gapchar for a -ter, -col pair.


	Smx->cell('-ter' => 1, '-col' => 2, '-gapchar' => '-'); OR
	Smx->cell('-ter' => 1, '-col' => 2, '-gapchar' => '[ACT]');


=cut

sub rawCell { 
	my $self = shift;
	my %params = @_;

	croak "no -tax passed to rawCell\n" if not defined $params{'-ter'};
	croak "no -blk passed to rawCell\n" if not defined $params{'-blk'};

	$params{-gapchar} = '-' if not defined $params {-gapchar};
	
	my $seq = $self->{mx}->[$params{'-ter'}][$params{'-blk'}];
	$seq =~ s/$params{'-gapchar'}//g;	

	return $seq;
}


=head2 yx

Returns sequence data in a yx accessor for matrix, only returns values! 

	Smx->yx(1,2);


=cut

sub yx { #
	my $self = shift;
	if (@_) {
		my ($y, $x) = @_;
		return $self->{mx}->[$y][$x];
	}
	return 0;
}


=head2 rawRowData

As in rawCell, but across a given row- overwritten by a -slice.

	$mx->rawRowData(
		'-gapchar' => '-',
		'-u2t' => 0,
		'-slice' => $my_slice,
		'-legal_alphabet' => $alphabet
	);


=cut

sub rawRowData { 
	my $self = shift;
	my %raw_params = @_;
	my %default_params = (
			'-gapchar' => '',
		   	'-u2t' => 0,
			'-blk_spacer' => '',
			'-bracket_blocks' => 0
		);
	my %params = (%default_params, %raw_params);

	croak "no -ter passed to rawRowData\n" if not defined $params{'-ter'};	
	$params{'-slice'} = $self->origSlice() if not defined $params{'-slice'};
	
	croak '-slice not a slice in rawRowData' if not ref($params{'-slice'}) eq 'slice';
    $params{'-legal_alphabet'}->add($params{'-blk_spacer'}) unless $params{'-blk_spacer'} eq '';
    $params{'-legal_alphabet'}->add("[]") if $params{'-bracket_blocks'} == 1;
	my $seq;

	foreach my $blk ($params{'-slice'}->loop('Blocks')) {
		if (($self->blk($blk)->bracketed == 1) and $params{"-bracket_blocks"} == 1) {
			$seq .= "[@{$self->{mx}->[$params{'-ter'}]}[$blk]]$params{'-blk_spacer'}"
		}
		else {
			$seq .= @{$self->{mx}->[$params{'-ter'}]}[$blk];
		  	$seq .=	$params{'-blk_spacer'}
		}
	}

	#	$seq .= join $params{'-blk_spacer'}, @{$self->{mx}->[$params{'-ter'}] } [ $params{'-slice'}->loop('Blocks') ];

	$seq =~ s/$params{'-gapchar'}//g; # strips the gapchar from the raw row data ## use a alphabet->clean here?

	if ($params{'-u2t'} == 1) { $seq =~ s/[uU]/T/gi }; ## use clean
	if ($params{'-legal_alphabet'}) { $seq = $params{'-legal_alphabet'}->clean('-str' => $seq, %params) } 
	
	return $seq;
}



=head2 column

Returns a column vector (string).

Requires either -col OR -pos and -blk.

If -blk is passed then -col is reference to position in block (0 is first position).

Usage: 
	Smx->column(
		'-blk' => 0,
		'-pos' => 1,
	);
	OR
	Smx->column(
		'-col' => 32
	);



See also columns package.  No bracketed columns can be returned.


=cut

sub column { 
	my $self = shift;
	my %raw_params = @_;
	my %def_params = ('-mode' => 'excludebracketed', '-blk' => $self->colInBlk(%raw_params), '-slice' => $self->origSlice );
	my %params = (%def_params, %raw_params);
	# print Dumper(%params);
	
	$params{'-pos'} ||= $params{'-col'} - $self->blk($params{'-blk'})->excludedStart  ;
	# print "pos: ", $params{'-pos'}, "\n";
	
	croak "no -blk in column \n" if $params{'-blk'} < 0;
	
	return -1 if $self->{blocks}->{$params{'-blk'}}->blkLength < $params{'-pos'};

	my $col;		
	foreach my $t ($params{'-slice'}->loop('Taxa')) {
		my $blk = $self->yx( $t, $params{'-blk'});
		$blk =~ /.{$params{'-pos'}}(.)/;
		$col .= $1	
	}
	return $col;
}


=head2 rowAsArray

-returns the sequence data for a given taxa for a given -slice as an array of blocks of length -max_seq_length
-strips the -gapchar regex if provided prior to calculating the length
-cleans the data if -alphabet and -clean set

=cut


sub rowAsArray {

	my $self = shift;
	my %raw_params = @_;
	
	die 'no -ter to matrix::rowAsArray' if not defined $raw_params{'-ter'};
	
	my %default_params = (
		'-max_seq_length' => 60,
		'-gapchar' => "",
		'-strip_gaps' => 0,
		'-clean' => 0,  ## added perhaps incorrectly, does it break out_fasta?
		'-slice' => $self->origSlice
	);
	my %params = (%default_params, %raw_params);
	my $seq = $self->rawRowData(%params);

	if ($params{'-strip_gaps'} == 1) { $seq =~ s/-//g };
	
	if ($params{'-clean'} == 1) { ## easier to implement here because of length requirements
		$seq = $params{'-legal_alphabet'}->clean(
											%params,
											'-str' => $seq,
										);	
	}
	return ( $seq =~ /.{1,$params{'-max_seq_length'}}/g ); # NO SPAACES IN THIS REGEX!!
}

=head2 rowData

-returns all blocks of data for a given row, overwritten by slice if provided
-strips the -gapchar regex if provided prior to calculating the length
-cleans the data if -alphabet and -clean set

=cut



sub rowData {
	my $self = shift;
	my %params = @_;
	croak 'no -tax passed to row' if not defined $params{'-tax'};
	# check 
	$params{'-slice'} = $self->origSlice() if not defined $params{'-slice'};
	croak '-slice not a slice in row_slice()' if not ref($params{'-slice'}) eq 'slice';

	return @{$self->{mx}->[$params{'-tax'}]} [ $params{'-slice'}->loop('Blocks') ];
}


=head2 blockData

Returns an array of all blocks of data for a given -blk, takes -slice too.

	Smx->blockData(
		'-blk' => 0,
		'-slice' => $my_slice,
	);


=cut


sub blockData {
	my $self = shift;
	my %params = @_;
	
	croak 'no -blk passed to block_column' if not defined $params{'-blk'};

	$params{'-slice'} = $self->origSlice if not defined $params{'-slice'};
	croak '-slice not a slice in blockData()' if not ref($params{'-slice'}) eq 'slice';

	my @results;
	foreach my $tax ($params{'-slice'}->loop('Taxa')) {
		push @results, $self->{mx}->[$tax][$params{'-blk'}]  
	}
	return (@results); # returns as list -> perhaps not correct
}


=head2 sliceData

Returns an array (DON'T DEPEND ON ORDER) containing all data in the passed slice, if no slice is passed the original slice (all data) is returned.

	Smx->sliceData(
		'-slice' => $my_slice,
	);


=cut

sub sliceData {
	my $self = shift;
	my %params = @_;
	
	$params{'-slice'} = $self->origSlice() if not defined $params{'-slice'};
	croak '-slice not a slice in data_slice()' if not ref($params{'-slice'}) eq 'slice';

	my @results;
	foreach my $tax ($params{'-slice'}->loop('Taxa')) {
		push @results, @{$self->{mx}->[ $tax] } [ $params{'-slice'}->loop('Blocks') ]
	}
	return (@results); # returns as list 
}


### Methods that return meta data derived from the matrix


sub alignedMask { # returns a Stockholm GC RF mask (i.e. marks aligned positions with an x, unaligned with .
	my $self = shift;
	my %params = @_;
	my $str = '';
	
	foreach my $blk ($params{'-slice'}->loop("Blocks")) {
		if ( $self->blk($blk)->bracketed == 1 ) {
	   		$str .= "." x $self->blk($blk)->blkLength }
		else { 
			$str .= "x" x $self->blk($blk)->blkLength
		}	
	}
	return $str; 
}

sub blkLengths { # returns an array of lengths of each aligned (unbracketed) block by default (overridden by -mode for different slice types)
	my $self = shift;
	my %raw_params = @_;
	my %default_params = ('-mode' => 'unbracketed', '-slice' => $self->origSlice);
	my %params = (%default_params, %raw_params);

	my @r;
	
	map { push @r, $self->blk($_)->blkLength  } ($self->loopBlocks(%params) );
	return @r;
}




=head2 numBracketedBlocks

Returns the count of the total number of bracketed blocks, overides with -slice

=cut

sub numBracketedBlks {
my $self = shift;
	my %raw_params = @_;
	my %default_params = (
	);
	my %params = (%default_params, %raw_params);

	$params{'-slice'} = $self->origSlice() if not defined $params{'-slice'};

    my $count = 0;

	foreach my $blk ($params{'-slice'}->loop('Blocks')) {
		$count++ if $self->blk($blk)->bracketed == 1
	}

	return $count
}



=head2 colInBlk 

Returns the block index to which -col (the column) can be found within.

Two modes: -mode => includebracketed | excludebracketed, defaults to later.


=cut


sub colInBlk { 
	
	## could change the algorithim to be more bubble-sort like to increase speed
	## could add -slice option, with rebuilding column indices for the slice
	
	my $self = shift;
	my %params = @_;

	$params{'-mode'} ||= 'excludebracketed';

	foreach my $i ($self->loopBlocks('-mode'=> 'all')) {
	 #	print "$i\n";
		## crude- could use better bubble-sort type approach	
		if ($params{'-mode'} eq 'includebracketed' ) {
			return $i if (($params{'-col'} >= $self->{blocks}->{$i}->literalStart) and ($params{'-col'} <= $self->{blocks}->{$i}->literalEnd) ); 
	 	}
		if ($params{'-mode'} eq 'excludebracketed' ) {
			if (defined	$self->{blocks}->{$i}->excludedStart) {
				return $i if (($params{'-col'} >= $self->{blocks}->{$i}->excludedStart) and ($params{'-col'} <= $self->{blocks}->{$i}->excludedEnd) ); 
			}
		}
	}	
	return -1;
}


=head2 loopBlocks

Returns an ordered array of block indices for the given slice of particular type, similar functions availible for slice objects.

Hmm.  Likely need to merge with slice object.

Legal values for 
	-mode => < unbracketed | bracketed | [all] | fiveprime > 

=cut


sub loopBlocks { 	
	my $self = shift;
	my %raw_params = @_;
	my %default_params = ('-mode' => 'all', '-slice' => $self->origSlice);
	my %params = (%default_params, %raw_params);

	my @blks;
	
	foreach my $blk ($params{'-slice'}->loop('Blocks')) { # keys %{$self->{blocks}} ) {	
		if ($params{'-mode'} eq 'all') {push @blks, $blk; next} 
		elsif ($params{'-mode'} eq 'bracketed' and $self->blk($blk)->bracketed == 1) {push @blks, $blk; next}
		
		elsif ($params{'-mode'} eq 'unbracketed' and $self->blk($blk)->bracketed == 0) {push @blks, $blk; next}

		elsif ($params{'-mode'} eq 'fiveprime' and  $self->blk($blk)->{original}->helix($blk)) {push @blks, $blk; next}
	}
	return (sort { $a <=> $b} @blks);
}


=head2 plan

Returns a basic 'plan', which is essentially a directive to fuse or not fuse paritions.
	'-mode' => < [basic] (recoded bracketed blocks) | all (recode all blocks) | all_fused (returns 1 big partition = the full slice) >
	'-slice' => '$self->origSlice'
	
=cut



sub plan {
	my $self = shift;
	my %raw_params = @_;
	my %default_params = ('-mode' => 'basic', '-slice' => $self->origSlice);
	my %params = (%default_params, %raw_params);

	my $p;
		
	if ( $params{'-mode'} eq 'all_fused') {
		$p->{1}->{blks} = ([$params{'-slice'}->loop('Blocks')]);
		$p->{1}->{kword} = 1;
		$p->{1}->{type} = 'trans'
	}

	else {	
		for my $blk ($params{'-slice'}->loop('Blocks')) {	

			if ($self->blk($blk)->bracketed == 1) {
				$p->{$blk}->{blks} = ([$blk]);
				$p->{$blk}->{kword} = 1; 
				$p->{$blk}->{type} = 'trans' # its always translated
			}
			else {
				if ($params{'-mode'} eq 'all') {
					$p->{$blk}->{blks} = ([$blk]);
					$p->{$blk}->{kword} = 1;
					$p->{$blk}->{type} = 'trans'
				}
				else {
					$p->{$blk}->{blks} = ([$blk]);
					$p->{$blk}->{kword} = 0; 
					$p->{$blk}->{type} = 'orig';
				}			
			}
		}

	}	
	return $p;
}



=head2 fragmentHash
	
Returns a hash of unique fragments, and their total. 

=cut


sub fragmentHash { 
	my $self = shift;
	# requires -blk, defaults to origSlice if not -slice
	my %params = @_;

	$params{'-slice'} ||= $self->origSlice;

	my %hash;
	foreach my $t ($params{'-slice'}->loop('Taxa')) {
		$hash{$self->yx($t, $params{'-blk'})}++;
	}
	return %hash;	
}


=head2 lastBlockIndex

Returns the index of the last block.

=cut


sub lastBlockIndex {
	my $self = shift;
	my @blocks = $self->loopBlocks;
	return pop @blocks;
}


### Methods that manipulate/change the matrix

=head2 insertBlock

Not tested/finished.

=cut

sub insertBlock {
	# inserts meta data for a block 
	# requires -blk

	my $self = shift;
	my %raw_params = @_;
	my %default_params	= ('-position' => 0, '-blk_length' => 0, '-bracketed' => 0, '-interleave' => 'original');
	my %params = (%default_params, %raw_params);
	
	# bump blocks to right of position

	my $blocks = $self->{blocks}; # store the old blocks;

	# delete everything post position
	for (my $i = $params{'-position'} + 1; $i < $self->lastBlockIndex +1; $i++) {
		$self->_deleteBlock('-blk' => $i);
	}
	
	# add new block
	my $cb = $self->_blockAdd(%params);

	## deleting blocks reduces accession !!! should work fine to reduce index then add it back up
	
	# add the block to the interleave to its right;																				
	$self->{interleaves}->{$params{'-interleave'}}->assign($self->{interleaves}->{$params{'-interleave'}}->blockFrom($params{'-position'}), $cb);	

	# add other blocks back	
	my $blk = $params{'-position'} + 1;
	until (not defined $self->{blocks}->{$blk}) {
		$self->{blocks}->{FIXME} =  $self->{blocks}->{$blk} ;
		$blk++;
	}
	
	# if $my_obj is a reference just say $foo{2} = $my_obj. Now they both (1, 2) point to the same reference
	# <blueberryCoffee> if $my_obj is not a reference then say $foo{1} = \$my_obj, and repeat for 2
	
	return 1;
}

sub removeBlock {
	# nukes a block completely, collapses indicies etc., the overal algorithm is crude, but it works
	my $self= shift;
	my %raw_params = @_;
	my %default_params = ('-structure' => 'original' );
	my %params = (%default_params, %raw_params) ;

	# delete data/meta data/structure refernces from the matrix
	$self->_deleteBlock(%params);

	# rebuild origSlice
	## might be wrong
	$self->{_origSlice}->remove('Blocks', $params{'-blk'});
	
	# collapse indicies
	my $mx = $self->{mx}; # store the old data
	$self->{mx} = undef; # wipe the old data
	
	my (@front_range, @back_range);
	unless ($params{'-blk'} == 0) { @front_range = ( 0..$params{'-blk'}-1 ) };
	unless ($params{'-blk'} == $self->lastBlockIndex) { @back_range = ($params{'-blk'} +1 .. $self->lastBlockIndex )};

	foreach my $ter ($self->loopTers) {
		# write front range
		@{$self->{mx}->[$ter]}[@front_range] = @{$mx->[$ter]}[@front_range];
		# write the back range
		my @foo = map $_ - 1, @back_range;
		@{$self->{mx}->[$ter]}[@foo] = @{$mx->[$ter]}[@back_range];
	}

	# collapse accession
	$self->{_block_accession}--;
	
	# collapse structure indicies	
	$self->{structures}->{$params{'-structure'}}->_collapse(%params);
	
	# adjust original slice
	$self->_initOrigSlice;
	
	# adjust original interleave
	$self->{interleaves}->{'original'}->remove_blocks($params{'-blk'});
	
	# recalculate meta
	$self->_calcBlockPositions;	
	return 1;
}

sub explodeBlocks { ## DO THIS
	# takes mode, or hash with splits
	my $self= shift;
	my %raw_params = @_;
	my %default_params = ( );
	my %params = (%default_params, %raw_params) ;

	# store the block to explode
	
	# delete the block to explode
	
	# insert the new blocks
	
	# insert the new data
	
	# rebuild the origSlice

}


sub generateSafenames {
	# generates an additional safeName for each taxon in the matrix   (no slice override)
	my $self = shift;
	foreach my $t ( $self->loopTers) { $self->{terminals}->{$t}->safeName($t) }
   	return 1;	
}


=head2 terLabel2Index 

Returns the index of a given terminal -label

=cut 


sub terLabel2Index {
	my $self = shift;
	my %params = @_; # requires -label

	die 'no label passed to terLabel2Index' if not defined $params{'-label'};
	
	foreach my $ter ($self->loopTers) {
		return $ter if $self->{terminals}->{$ter}->label eq $params{'-label'};
	}
	return undef;
}


=head2 totalChars

Returns the total number of characters for a given slice (defaults to whole matrix) for the given -slicemode
	$mx->totalChars(
		'-slicecharmode' => 'all'
	);

Options:
	'-slicecharmode' => < all | [unbracketed] | bracketed >

=cut


sub totalChars {
	my $self = shift;
	my %params = @_;
	
	$params{'-slicecharmode'} = 'unbracketed' if not defined $params{'-slicecharmode'};
	$params{'-slice'} = $self->origSlice() if not defined $params{'-slice'};
	croak '-slice not a slice in total_chars()' if not ref($params{'-slice'}) eq 'slice';

	my $total = 0;
	
	foreach my $blk ($params{'-slice'}->loop("Blocks")) {
		if ( $self->{blocks}->{$blk}->{bracketed} ) {
			if ( ($params{'-slicecharmode'} eq 'all') or ($params{'-slicecharmode'} eq 'bracketed'  ) ) {
				$total += $self->{blocks}->{$blk}->blkLength;
			}
		}
		else {
			if ( ($params{'-slicecharmode'} eq 'all') or ($params{'-slicecharmode'} eq 'unbracketed'  ) ) {
				$total += $self->{blocks}->{$blk}->blkLength;
			}
		}
	}
	return $total;
}

=head2 loopTers

Returns a sorted array of _indicies_ to the taxa in the matrix (no slice override)

=cut


sub loopTers {
	my $self = shift;
	return (sort {$a <=> $b} ( keys %{$self->{terminals}} ))	
}

=head2 lengthLongestTer

Returns the length of the longest taxon name in a slice

=cut


sub lengthLongestTer{ 
	my $self = shift;
	my %params = @_;
	
	$params{'-slice'} ||= $self->origSlice;	

	my $longest = 0;
	foreach my $t ($params{'-slice'}->loop("Taxa")) {
        #  print $self->ter($t)->label;
        # print "-- $t -- \n";
		$longest = length($self->ter($t)->label) if length($self->ter($t)->label) > $longest;	
	}
	return $longest;
}

=head2 numberer

Numbers blocks or columns. Includes one big string with line endings. Deprecated reference to a particular interleave

		-mode' => 'unbracketed',
		-numbering_mode' => 'col',
		-blk_spacer' => '',
		-gap_chr' => ' ', # inserted b/w blocks
		-number_bracketed' => 0, # number the bracketed blocks, or not
		-start_col' => 1,
		-original_indexing' => 0 
	

=cut

sub numberer () {
	my $self = shift;
	my %raw_params = @_;
	my %default_params	= (
		'-mode' => 'unbracketed',
		'-numbering_mode' => 'col',
		'-blk_spacer' => ' ',
		'-number_bracketed' => 0, 
		'-start_col' => 1,
		'-original_indexing' => 0 
	);	
		
	my %params = (%default_params, %raw_params);
	$params{'-slice'} = $self->origSlice unless defined $params{'-slice'};
	croak "Can't number bracketed blocks if no bracketed blocks are passed!" if $params{'-number_bracketed'} == 1 and $params{'-mode'} eq 'unbracketed';
	
	my ($numrows, $longest) = (0) x 2;
	my ($t, $outstring);
	my $llt ;
    $llt=	$self->lengthLongestTer(%params);

	# find the longest item (number of digits) 	
	if ($params{'-numbering_mode'} eq 'col' ) {	# sum chararcters by mode
		$longest = $self->totalChars(%params, '-slice_chr_mode' => $params{'-mode'});
	}
	elsif ($params{'-numbering_mode'} eq 'blk' ) { 
		$longest = $params{"-slice"}->total('Blocks'); 
	}
	else { croak "illegal mode [$params{'-numbering_mode'}]\n"}
	
	$numrows = length("$longest");

	for (my $i=0; $i < $numrows; $i++) {  # print this many lines 
		$outstring .= "[";
		$outstring .= ' ' x $llt;
		#print &padright("["," ", $startcol);  #put $i. to check number of rows
	
		my ($bi,$pos);
		$bi = 0; # reindex the blocks here		
		$pos = $params{'-start_col'};

		foreach my $blk ($params{'-slice'}->loop("Blocks")) {
			if ($params{'-numbering_mode'} eq 'col') {
				my ($start, $end); # for the present block

				# calculate start/end positions
				if ($params{'-original_indexing'} == 1) {
					if ($params{'-number_bracketed'} == 0) { # use literal 
						$start = $self->blk($blk)->literalStart;	
						$end = $self->blk($blk)->literalEnd;
					}
					else {
						$start = $self->blk($blk)->unbracketedStart;	
						$end = $self->blk($blk)->unbracketedEnd;
					} 
				}
				else {
					$start = $pos;
					$end = $pos + $self->blk($blk)->blkLength -1;	
				}	

				# numbering or not
				if ($params{'-number_bracketed'} == 0 and $self->blk($blk)->bracketed == 1) { # gap it
					$outstring .= '[' if ($params{'-bracket_blocks'} == 1);	
					$outstring = $outstring.&padright("", "-", $self->blk($blk)->blkLength); #+1 is correct!!!
					$outstring .= ']' if ($params{'-bracket_blocks'} == 1);	
				}
				else { # number it
					$outstring .= '[' if (($self->blk($blk)->bracketed == 1) and ($params{'-bracket_blocks'} == 1));	
					for (my $k = $start; $k < $end + 1; $k++) { # using unbracketed blocks only now, number individual positions	
											$t = &padleft($k,"0",$numrows);
						$t =~ /.{$i}(.)/; # there's likely a better way to choose either the ith character or print a zero.
						$outstring = $outstring.$1;
					}
					$outstring .= ']' if (($self->blk($blk)->bracketed == 1) and ($params{'-bracket_blocks'} == 1));	
				}

				# between blocks
				if (not $params{"-blk_spacer"} eq '') {
					$outstring = $outstring.$params{"-blk_spacer"}
				}

				if (($params{'-number_bracketed'} == 1 and $self->blk($blk)->bracketed == 1) or ($self->blk($blk)->bracketed == 0)) {
					$pos += $self->blk($blk)->blkLength;
				}
			}

			elsif ($params{'-numbering_mode'} eq 'blk') { # number blocks by their original position
				if ($params{'-original_indexing'} == 1 ) { # number using the position as originally read
					$t = &padleft($blk, "0" , $numrows);
				} 
				else { # number sequentially
					$t = &padleft($bi, "0" , $numrows);
				}	
			
				$t =~ /.{$i}(.)/; 	
				my $t1 = &padleft("$1"," ", $self->blk($blk)->blkLength);
				
				if ($self->blk($blk)->bracketed == 0) {
					$outstring = $outstring.$t1;
					if (not $params{"-blk_spacer"} eq '') {
						$outstring = $outstring.$params{"-blk_spacer"}
					}
				}

				elsif (($self->blk($blk)->bracketed == 1)) { # and ($params{'-mode'} eq 'bracketed')) { # also bracket this block to indicate original bracketing


					$outstring .= '[' if ($params{'-bracket_blocks'} == 1);	
					$outstring = $outstring."$t1";
					$outstring .= ']' if ($params{'-bracket_blocks'} == 1);	
					
					if (not $params{"-blk_spacer"} eq '') {
						$outstring = $outstring.$params{"-blk_spacer"}
					}
				}	
			$bi++;	
			}
	}
	 		$outstring = $outstring."]";
		$outstring .= "\n";
}
	return $outstring;
}






### Internal methods- not to be called outside

=head2 _calcBlockPositions

(Re)calculate the excluded/unexcluded starts/ends

=cut


sub _calcBlockPositions {
	my $self = shift;
	my ($bl, $ubl) = (0) x 2;
	
	foreach my $blk ($self->loopBlocks) {
		if ($self->{blocks}->{$blk}->bracketed == 1) {
			$self->{blocks}->{$blk}->literalStart($bl);	
			$self->{blocks}->{$blk}->literalEnd($bl + $self->{blocks}->{$blk}->blkLength -1);	
		  }
		else {
			$self->{blocks}->{$blk}->excludedStart($ubl);
			$self->{blocks}->{$blk}->excludedEnd( $ubl + $self->{blocks}->{$blk}->blkLength - 1);
			$self->{blocks}->{$blk}->literalStart($bl);	
			$self->{blocks}->{$blk}->literalEnd($bl + $self->{blocks}->{$blk}->blkLength -1);	
 	
			$ubl += $self->{blocks}->{$blk}->blkLength;
		}
		$bl += 	$self->{blocks}->{$blk}->blkLength;
		  
		#print "blk: $blk\n";
		#print "ubl: $ubl\n";	
		#print "bl: $ubl\n\n";
	}
	return 1;
}

sub _deleteBlock { 
	# internal function to remove the data/metadata/structure refs for all terminals for a given -blk  (does not collapse etc -> see removeBlock for that) 
	my $self= shift;
	my %raw_params = @_;
	my %default_params = ('-structure' => 'original' );
	my %params = (%default_params, %raw_params) ;

	die 'no -blk passed to _deleteBlock' if not defined $params{'-blk'};
	
	# -mode => [collapse|preserve] - collapses the metadata and block indices of the downstream blocks, or not
	
	# delete the block metadata
	delete $self->{blocks}->{$params{'-blk'}};
	
	# delete the block data
	foreach my $ter ($self->loopTers) {
		delete $self->{mx}->[$ter][$params{'-blk'}];
	}
		
	# delete structure references to block
	$self->{structures}->{$params{'-structure'}}->_deleteBlock(%params);	

	return 1
}


sub _blockAdd { 
	# adds a block (meta data) to the matrix, returns the index of the *created* block
	my $self = shift;
	my %raw_params = @_;
	my %default_params = ('-blk_index' => $self->{_block_accession});
	my %params = (%default_params, %raw_params);
	
	$self->{blocks}->{$params{'-blk_index'}} = block->new(%params);
	#print $self->{_block_accession};
	
	if ( $self->{_block_accession} == $params{'-blk_index'} ) { # we added a block in sequential order
		#print " $self->{_block_accession} ";
		$self->{_block_accession}++;	
		return $self->{_block_accession} -1;
	}
	return $params{'-blk_index'};
}

sub _interleaveAdd { ## used?
	my $self = shift;
	if (@_) {
			my $name = shift;
			$self->{interleaves}->{$name} = interleave->new();
			#$self->{_interleaves_accession}++;
		return 1;
	}
	return '';
}

=head2 _terminalAdd

Makes a terminal by setting the label for it.  Not for external use.

=cut


sub _terminalAdd {
	my $self = shift;	
	if (@_) {
		my %params = @_; 
		$self->{terminals}->{$self->{_terminals_accession}} = terminal->new(%params);
			$self->{_terminals_accession}++;
		return 1;	
	}
	return 0;
}

sub _structureAdd {
	my $self = shift;
	if (@_) {
		my %params = @_; 
			$self->{structures}->{$params{'-structure_name'}} = structure->new();
			$self->{_structure_accession}++;
		return 1;	
	}
	return '';
}

sub _initOrigSlice {
	my $self = shift;
		$self->{_origSlice} = slice->new();
		$self->{_origSlice}->blocks(0..$self->{_block_accession}-1); 
		$self->{_origSlice}->taxa(0..$self->{_terminals_accession}-1);
	return 1;
	
}

=head1 Accessors

=cut

sub nters { # accessor for # of terminals
	my $self = shift;
    if (@_) {$self->{nters} = shift}
	return $self->{nters};
}

sub inputformat { # accessor for the inputformat
	my $self = shift;
	if (@_) {$self->{inputformat} = shift;}
	return $self->{inputformat} || undef;
}

sub interleaves { ## THIS SHOULD RETURN CLONE accessor for the interleave object(s) 
	my $self = shift;
	if (@_) {$self->{interleave} = shift;}
	return $self->{interleave};
}


=head2 label

Accessor for the matrix label

=cut


sub label {
	my $self = shift;
	if (@_) {$self->{label} = shift;}
	return $self->{label};
}


=head2 totalBlocks

Returns the total blocks in the matrix (_block_accession)

=cut


sub totalBlocks {
	my $self = shift;
	return $self->{_block_accession}  ;  # total is 0+1 
}


=head2 origSlice

Returns a slice object that is the clone of _origSlice.  

=cut


sub origSlice {
	my $self = shift;
	return $self->{_origSlice}->clone;
}


=head2 origStructure

Returns reference to the ORIGINAL structure (not a clone-> yet)

=cut


sub origStructure {
	my $self = shift;
	return $self->{'structures'}->{'original'};
}


=head2 origInterleave

Returns the original interleave (not a clone).


=cut


sub origInterleave { 
	my $self = shift;
	#print Dumper ($self->{_origSlice});
	return ($self->{interleaves}->{'original'});
}


### Alignment

=head2 alignBlock 

Was clustalAlignBlock.  Align an individual block using an algorithmic aligner.  After alignment removes all columns with only gaps. If no -blk is passed aligns all bracketed blocks.

Needs work with paths, don't create an Output object prior to using this at present (Jan/16/06).

Required:
	-blk => @foo  (an list of blocks to align)
	
Options:
	-unbracket => <[1] | 0>  when on bracketed blocks that are aligned will be turned into unbracketed blocks
	-taxa => @bar (align only for these taxa)
	-align_method => < [clustalw] | muscle >
	-clustal_params => '' (additional parameters to be passed to clustal)
	-muscle_params => ''
	
=cut


sub alignBlock {
	my $self = shift;
	my %raw_params = @_;
	
	my %default_params = ('-align_method' => 'clustalw', '-clustal_params' => '', '-muscle_params' => '', '-unbracket' => 1);
	my %params = (%default_params, %raw_params);

	$params{'-blk'} || die "no -blk to alignblock\n";
	$params{'-taxa'} ||= [ $self->loopTers ];
	
	my $s = slice->new;
	$s->blocks($params{'-blk'});
	$s->taxa($params{'-taxa'});
	
	if ($params{'-unbracket'} == 1) {
		print "unbracketing $params{'-blk'}\n";
		$self->blk($params{'-blk'})->bracketed(0);
	}
	
	my $fmx = $self->blk2Fastmatrix('-slice' => $s);

	print "\naligning block: ", $params{'-blk'}, "\n";

	if ( $params{'-align_method'} eq 'clustalw') {
		$fmx->clustalAlign(%params);
	}
	elsif ($params{'-align_method'} eq 'muscle') {
		$fmx->muscleAlign(%params);
	}
	else { die 'Not a legal align_method to mx->alignBLock.'}
	
	$fmx->deleteFixedGapPositions;
	$self->fastmatrix2Blk(%params, '-fastmatrix' => $fmx);

	# might have changed the length of a block
	$self->blk($params{'-blk'})->blkLength( length($self->yx(0, $params{'-blk'})) );
	
	print "\n";
	return 1;
}


=head2 alignSliceBlocks

Aligns (<foo>AlignBlock) all the blocks in the -slice, defaults to all bracketed blocks if no slice is passed.

Options:
	-align_method => < [clustalw] | muscle >
	-slice
	-mode (sensu slice)
	
=cut

sub alignSliceBlocks {
	my $self = shift;
	my %raw_params = @_;
	
	my %default_params = ('-align_method' => 'clustalw');
	my %params = (%default_params, %raw_params);
	
	if (not defined $params{'-slice'}) {
		my $s = slice->new;
		$s->blocks($self->loopBlocks('-mode' => 'bracketed'));
		$params{'-slice'} = $s;
	}
		
	for my $blk ($params{'-slice'}->loop('Blocks')) {
		$self->alignBlock('-blk' => $blk, %params)
	}
	return 1;
}

=head2 blk2Fastmatrix 

Translates a slice to a Fastmatrix object, concatenating all blocks.  
Note that rawRowData returns origSlice if no -slice is passed.
Returns the Fastmatrix object. 

=cut

sub blk2Fastmatrix {
	my $self = shift;
	
	my %raw_params = @_;
	my %default_params = ();
	my %params = (%default_params, %raw_params);

	use Psy::Helpers::Fastmatrix;
	my $fm = Psy::Helpers::Fastmatrix->new;
	
	map {my $d = $self->rawRowData('-ter' => $_, %params); $fm->{mx}->{$_} = $d} $self->loopTers(%params);
	return $fm;
}

=head2 fastmatrix2Blk 

Reverse of blk2Fastmatrix.  Merges a fastmatrix to a block.  

=cut

sub fastmatrix2Blk {
	my $self = shift;
	my %params = @_;
	$params{'-labels_are_indicies'} ||= 1;
	
	defined $params{'-blk'} || die 'no -blk to fastmatrix2blk';
	$params{'-fastmatrix'} || die 'no -fastmatrix to fastmatrix2blk';
	my $fm = $params{'-fastmatrix'}; # simplify
	if ($params{'-labels_are_indicies'} == 1) {
		map { $self->cell(%params, '-ter' => $_, '-data' => $fm->seq($_)) }	$fm->loopTer
	}
	else {
		map { $self->cell(%params, '-ter' => $self->terLabel2Index('-label' => $_), '-data' => $fm->seq($_)) } $fm->loopTer
	}
	
	return 1;
}



### Misc

=head2 origSliceHash

For some functions (fusing blocks) it is useful to have a hash of label => @array.  This returns a default for those functions, with one block per array.

=cut


sub origSliceHash {
	my $self = shift;
	my %r;
	map {$r{$_} = ([$_]) } $self->origSlice->loop('Blocks'); 
	return \%r;
}

### Debugging and tests

sub _debugBlockRanges {
	my $self = shift;
	foreach my $blk ($self->loopBlocks) {
		print $blk, "\t";
		print $self->{blocks}->{$blk}->bracketed, "\t";
		print $self->{blocks}->{$blk}->blkLength, "\t";	
		print $self->{blocks}->{$blk}->literalStart, "\t";
		print $self->{blocks}->{$blk}->literalEnd, "\t";
		print $self->{blocks}->{$blk}->excludedStart, "\t";
		print $self->{blocks}->{$blk}->excludedEnd, "\t";
		print "\n";
	}
	print "\n";
	return "WHWHAHHTT $self\n\n\n";
}

sub _debug_blockLengths {
	my $self = shift;
	print "debuging block Lengths (data from first row)";	
	foreach my $blk ($self->{_origSlice}->loop("Blocks")) {
		print $self->blk($blk)->blkLength, " ", $self->yx(0,$blk), " ", length $self->yx(0,$blk) ==  $self->blk($blk)->blkLength ? "passed" : "failed",  "\n";
	}
}


1;
