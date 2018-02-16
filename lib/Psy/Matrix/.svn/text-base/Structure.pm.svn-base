package structure;

use vars qw(@ISA @EXPORTER @EXPORT_OK $VERSION);

$VERSION = '0.01';

use strict;
use warnings;
use Data::Dumper;
use Carp;

use Psy::Psy qw(generatedByHeader version rootpath projectname starttime $ROOT_DIR $PSY_BASE_DIR $PSY_ROOT_DIR);
@ISA = qw(Psy);

use Psy::Io::Io;


=head1 NAME

Psy::Matrix::Matrix

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Builds a structure object, at present based on a helix-index file (but see also Psy::Matrix::_load_PHASE_simulate). 

    use Psy::Matrix:Structure;	
	my $fs = structure->new(
						'-helix_index' => 'Evaniodea.index.0.02.txt',
						'-path' => '../../data/evan/',		
					);


	print Dumper($s);
	print $s->{blocks}->version();

=cut

=head2 new

Create a new object, not for direct reference/use.

=cut

sub new {
	my $type = shift;
	my %params = @_;	
  	my $self = {};
	
	bless $self, $type;

	$self->_init(%params);

	#print "\n[", ref ($self) ,"]\n";;
    return $self;                 
}

sub _init {
	my $self= shift;
	my %params = @_;

	%params = $self->_verify(%params); # probably a bad way to do this
	
	# $self->structure_name(%params);
	
	$self->{block_masks} = {}; 		# hash with block index pointing to mask
	$self->{helices} = {}; 			# hash with block index (5') pointing to block index (3') 
	$self->{block_labels} = {};		# hash with block index pointing to label

	$self->_load(%params) if defined $params{'-helix_index'};

	$self->version($VERSION);
	return 1;
}
	
sub _verify {
	my $self=shift;
	my %params = @_;

	#	print "\nSTRUCTURE PARAMS: ";
	
	# $foreach my $p (keys %params) { # check matrix parameters here
	#		print "$p, ";
	# }
	# print "\n";
	return %params;
}

sub _collapse { # collapses all the indices to the right of -blk (assumes -blk has been delete prior)
	my $self = shift;
	my %params = @_; # requires '-blk' - collapse if > than

	my $blkmasks = $self->{block_masks};
	my $helices = $self->{helices};
	my $blocklabels = $self->{block_labels};

	# collapse {block_masks}
	foreach my $blk (keys %{$self->{block_masks}}) {
		if ($blk > $params{'-blk'}) { # only delete/rebuild if in range
			my $tmp = $self->{block_masks}->{$blk}; 
			$self->{block_masks}->{$blk} = undef;
	
			$self->{block_masks}->{$blk-1} = $tmp; 
		}
	}
	
	# collapse {helices}
	foreach my $blk (keys %{$self->{helices}}) {
		if ( ($blk > $params{'-blk'}) or ($self->threePrimePair(%params) > $blk) ) { # at least one side needs to be updated
			my $three = $self->threePrimePair(%params);
			my $five = $blk;

			if ($three > $blk) {$three--};
			if ($five > $blk) {$five--};
			$self->deleteHelix(%params); # trash old
			$self->helix($five, $three); # write new
		}
	}

	# collapse {block_labels}
	foreach my $blk (keys %{$self->{block_helices}}) {
		if ($blk > $params{'-blk'}) {
			my $label = $self->label('-blk' =>$blk);
			delete $self->{block_labels}->{$blk};
			$self->label('-blk' => $blk-1, '-label' => $label);	
		}
	}

	return 1	
}

# sub structure_name {
#	my $self=shift;
#	if (@_) {$self->{structure_name} = shift};
#	return $self->{structure_name};
	#
#}




=head2 helix

An accessor, if two blocks passed makes a helix with block 1 = 5', block 2 = 3'; with one block passed returns 3'.

=cut

sub helix {  
	my $self = shift;
	if (@_) {
		my $s1 = shift;
		if (@_) {
			my $s2 = shift;
			$self->{helices}->{$s1} = $s2;
			$self->{_helix_accession}++;
			return 1;
		}
		else {
			return $self->{helices}->{$s1};  # returns 3' end !!
		}
	}
	return 0; # don't change!, required to allow helix to detect for paired
}


##### PORT specific block label to Block
=head2 label

Accessor for a -label of -blk.
	$s->label('-blk' => 1, '-label_name' => 'foo');
	OR 
	$s->label('-blk' => 1);	
=cut


sub label {
	my $self=shift;
	my %params = @_;
	$params{'-blk'} || return ;
	$params{'-label_name'} || return;
	if ($params{'-label'}) { 
		$self->{block_labels}->{$params{'-blk'}}->{$params{'-label_name'}} = $params{'-label'};
	}
	return 	$self->{'block_labels'}->{$params{'-blk'}}->{$params{'-label_name'}} ;
}


=head2 stockholm1RFLine

Returns a stockholm legal RF line

=cut


sub stockholm1RFLine {
	my $self = shift;
	my %raw_params = @_;
	# needs -slice, $mx
	my %default_params = ('-max_seq_length' => 80, '-no_mask_char' => '.');
	my %params = (%default_params, %raw_params);

	my $mask;
	foreach my $bi ($params{'-slice'}->loop("Blocks")) {
		my $blk =  $params{'-mx'}->block($bi);
		if ($blk->bracketed == 1) {
			$mask .= '.' x $blk->length;
		}	
		else {
		}
		
	}

}


=head2 phase2Mask

Return a formated, PHASE 2.0 legal mask in one full string.

=cut


sub phase2Mask { 
	my $self = shift;
	my %raw_params = @_;
	# needs -slice
	my %default_params = ('-max_seq_length' => 80, '-no_mask_char' => '.');
	my %params = (%default_params, %raw_params);
	my %mask_chars = $self->phansyMaskChars(%params);
	my $full_mask;
	
	foreach my $blk ($params{'-slice'}->loop("Blocks")) {
		my $blk_mask = $self->mask('-blk' => $blk);
		if (not defined $mask_chars{$blk}) { $blk_mask =~ s/./\./g}
		else {$blk_mask =~ s/[\(\)]/$mask_chars{$blk}/g};
		$full_mask .= $blk_mask;
	}

	my @full = $full_mask =~ /.{1,$params{'-max_seq_length'}}/g; # NO SPACES IN THIS REGEX!!
	
	return join "\n", @full;
}


=head2 phase2Class

Returns a semi-complete class line for use in PHASE 2.0 analyses

=cut


sub phase2Class { 
	my $self = shift;
	my %raw_params = @_;
	# needs -slice
	my %default_params = ('-max_seq_length' => 80);
	my %params = (%default_params, %raw_params);

	my $full_mask;
	
	foreach my $blk ($params{'-slice'}->loop("Blocks")) {
		my $blk_mask = $self->mask('-blk' => $blk);
		$blk_mask =~ s/[^\(\)]/1/g;
		$blk_mask =~ s/[\(\)]/2/g;
		$full_mask .= $blk_mask;
	}

	my @full = $full_mask =~ /.{1,$params{'-max_seq_length'}}/g; # NO SPAACES IN THIS REGEX!!
	
	return join "\n", @full;
}


=head2 phansyMaskChars

Returns a hash of block->mask-letter. Used for translating () masks to <> or aA masks.

Algorithm could be very likely be greatly simplified by recursion, but it *appears* to work as is (its ugly!).

=cut


sub phansyMaskChars { # 

	my $self = shift;
	my %params = @_; # requires a -slice
	my %mask; # the hash to return

	my $i = 0; # index to mask character pair
	
	# map the possible mask pairs to a hash of hashes, can add more letters if needed to %mchars
	my %mc; 
   	my %mchars = qw(( ) [ ] { } < > A a B b C c D d E e F f G g H h I i J j K k L l M m N n O o P p Q q R r S s T t U u V v W w X x Y y Z z);
	map {$mc{$i} = ([$_, $mchars{$_}]); $i++} (sort keys %mchars);
	
	# print Dumper(%mc);
	
	# build an ordered hash of arrays containing helices, so we can do look-aheads
	my @helices = $self->loopHelices(%params);
	my $hel;
	my $h = 0;
	map { $hel->{$h}= [$_, $self->helix($_)]; $h++} (@helices);
	
	# set some flags
	my $reset = 0;
	my $in_nested = 0;
	 
	$i = 0;
	# loop through the helices and assign mask pairs
	foreach my $k (sort {$a <=> $b} keys %{$hel}) {
		# print $self->helix($k);
		# print " ";
		# print $k, " ", @{$hel->{$k}}[0], " ",  @{$hel->{$k}}[1], ": ";	
		my $nested = 1;
		my $cl = @{$hel->{$k}}[1]; 
		my $cur_reset = 0;
			
		for (my $c = $k+1; $c < $#helices + 1; $c++) {	
			# print " $cl";  
	        $cur_reset++; 

			#(@{$hel->{$k}}[0] == @{$hel->{$k}}[1] + 1) && last; # they are side by side, must be nested

			if (@{$hel->{$c}}[1] > @{$hel->{$k}}[1]) { 
				if (@{$hel->{$k}}[1] > @{$hel->{$c}}[0] ) { $nested =0 }  			
				last;
			}
						
			#	(@{$hel->{$c}}[1] > @{$hel->{$k}}[0]) && ($nested = 0);  # added

			if ( $cl <  @{$hel->{$c}}[1]) {
				$nested = 0; 
			}

			#if ( @{$hel->{$k}}[1] > @{$hel->{$c}}[0] ) {
				#$nested = 0;
				#}	

			$cl = @{$hel->{$c}}[1]; 	
		}	
		
		($reset < $cur_reset) ? ($reset = $cur_reset -1) : $reset--;

		if ($nested == 1) {
		   	if ($in_nested == 0) {
				$in_nested = 1;
			};
			# print "\t nested |$i| ";

			# set the mask
			$mask{ @{$hel->{$k}}[0] } = @{$mc{$i}}[0];
			$mask{ @{$hel->{$k}}[1] } = @{$mc{$i}}[1];
		}
		else {
			if ($in_nested ==1) {
				$in_nested = 0;
				$i++;
			}
			# print  "\t |$i| ";

			# set the mask
			$mask{ @{$hel->{$k}}[0] } = @{$mc{$i}}[0];
			$mask{ @{$hel->{$k}}[1] } = @{$mc{$i}}[1];

			$i++; 
		};

		($reset == 0)  && ($i = 0);
		# print	"[ $cur_reset $reset ] \n";	
	 }
	 # print Dumper(%mask);
	return %mask;
}


=head2 sliceMask

Returns the mask for a slice as a string, substitutes {bracket_chars} if passed, inserts -blk_spacer character between blocks if passed. If a -mx and -bracketed_blocks == 1 then mask will be spaced accordingly.

=cut


sub sliceMask {
	my $self = shift;
	my %raw_params = @_;
	my %default_params = (
			'-blk_spacer' => '',
			'-bracket_blocks' => 0
		);
	my %params = (%default_params, %raw_params);
	my $str;

	if ($params{'-bracket_blocks'} == 1) { 	# require a matrix
		croak '-mx not a matrix in matrix' if not ref($params{'-mx'}) eq 'matrix';
	}	
		
	foreach my $blk ($params{'-slice'}->loop("Blocks")) {
		if (($params{'-bracket_blocks'} == 1) and ($params{'-mx'}->blk($blk)->bracketed == 1 )) {
			$str .= " ".$self->mask('-blk' => $blk)." ";	
		}
		else {
			$str .= $self->mask('-blk' => $blk);
		}
			$str .= $params{'-blk_spacer'};
	}
	
	# replace the backeted chars if replacements are there
	if ($params{'-bracket_chars'}) {
		#my $c = quotemeta($params{'-bracket_chars'}->{left});
		$str =~ s/\(/$params{'-bracket_chars'}->{left}/g;	
		$str =~ s/\)/$params{'-bracket_chars'}->{right}/g;		
	}
	
	# $str =~ s/\s/./g; # proof for whitespace in mask # CAN'T DO THIS
	return $str;
}





=head2 mask

Returns the mask for a block.

=cut


sub mask {
	my $self=shift;
	my %params = @_;
	if ($params{'-mask'}) {
	   	croak 'no -blk to ->mask in structure.pm\n' if not defined $params{'-blk'};	
		$self->{block_masks}->{$params{'-blk'}} = $params{'-mask'};
	}
	return $self->{block_masks}->{$params{'-blk'}} || ''; 
}


=head2 maskPair

Returns hash with paired masks under '5prime' and '3prime' indicies.

=cut


sub maskPair {
	my $self= shift;
	my %params = @_;

	croak "no block passed to maskPair" if not defined $params{'-blk'};
	
	my %hash;
	
	$hash{'5prime'} = $self->mask(%params);
   	$hash{'3prime'} = $self->mask('-blk' => $self->helix($params{'-blk'}) );	

	return %hash;
}


=head2 bracketsInMask

Returns a count of the number of parens for a given block.

=cut


sub bracketsInMask {  
	my $self = shift;
	my %params = @_;
	croak "no block passed\n" if not defined $params{'-blk'};
	my $temp = $self->{block_masks}->{$params{'-blk'}};
	return $temp =~ tr/[\(\)]//; 
}


=head2 nucPairs

Returns a hash with keys {pairs) pointing to 5':3' pair positions (MrBayes 'pairs') and {nucs} a list of all positions (column indicies).

To count from 0 set zero_beers = 1

=cut



sub nucPairs { 
	my $self = shift;
	my %params = @_;

	$params{'-zero_beers'} ||= 0; # set to 1 to count from 0

	my $hash = {};
	$hash->{nucs} = [];
	$hash->{pairs} = [];
	$hash->{nucs_ns} = [];
	
	die "no matrix passed" if not defined $params{'-mx'};	
	croak '-mx not a matrix in matrix' if not ref($params{'-mx'}) eq 'matrix';

	$params{'-slice'} ||= $params{'-mx'}->origSlice; 

	my $mx = $params{'-mx'}; 	# easier to work with
	my ($p, $val, @r, @l, $tmp, $rc, $lc, @cs, @pairs);

	foreach my $helix ($self->loopHelices(%params)) { # passes the slice
		my @lefts;
		my @rights;
		my $p = 0;
		@lefts = map { $p++; $p + $params{'-slice'}->excludedBlkStart('-mx' => $mx, '-blk_index' => $helix)  - $params{'-zero_beers'} if $_ eq '('} ( split(//, $self->mask('-blk' => $helix) ));
		$p = 0;
		@rights = map { $p++; $params{'-slice'}->excludedBlkEnd('-mx' => $mx, '-blk_index' => $self->helix($helix)) + 2 - $params{'-zero_beers'}  - $p if $_ eq ')'} ( split(//, reverse ( $self->mask('-blk' => $self->helix($helix) ) )));	
		
		foreach my $val (@lefts) {
			push (@l, $val) if length $val > 0
		} # strip nulls from map, should be an easier way
		foreach $val (@rights) {
			push (@r, $val) if length $val > 0
		}

		# debug to screen
		# print "[helix $helix-", $self->helix($helix) ,":\t@l | @r \n";
		# if ($#l != $#r) {print "<- unequal"};
		
		while (@l) {
			$tmp = "";
			$lc = shift @l;
			push (@cs, $lc); #make list of stem columns
			$tmp .= $lc;
			$tmp .= ":";
			$rc = shift @r;
			$tmp .= $rc;
			push (@{$hash->{nucs}}, $lc);
			push (@{$hash->{nucs}}, $rc);
			push (@{$hash->{pairs}}, $tmp); #make list of pairs
		}
		
		# get the non-stem indecies (hmm- returns array of letters, not integers)
		my %foo = map {$_ => undef} ((1-$params{'-zero_beers'})..$params{'-mx'}->totalChars('-slice' => $params{'-slice'})- $params{'-zero_beers'});
		map {delete $foo{$_}} @{$hash->{nucs}};
		
		@{$hash->{'nucs_ns'}} = sort {$a <=> $b} (map {int($_)} keys %foo);
	}
	return $hash
}


=head2 basePairs

Returns a hash of basepairs encountered, strictly adhering to $alphabet;

Not done.

=cut


sub basePairs { 
	my $self = shift;
	my %params = @_;
	# use column objects?
	
	die "no matrix passed" if not defined $params{'-mx'};	
	croak '-mx not a matrix in matrix' if not ref($params{'-mx'}) eq 'matrix';
}



=head2 alignedHelices

Returns a hash of keys (mask1, mask2, stem1, stem2) with the stems/masks aligned

=cut


sub alignedHelices {
	# requires -mx, -ter, -blk
	
	my $self = shift;
	my %params = @_;

	my %hash;
	
	croak "no matrix passed" if not defined $params{-mx};
    croak "not a matrix in alignedStems\n" if not ref($params{'-mx'}) eq 'matrix';
	croak "no terminal passed" if not defined $params{'-ter'};
	croak "no helix 5' end passed" if not defined $params{'-blk'};

	croak "-blk is not a helix member in alignedHelices" if not defined $self->helix($params{-blk});
	
	$params{'-structure'} = 'original' if not defined $params{'-structure'}; # can override to another structure if you want
	
	my %masks = $self->maskPair(%params);
	
	$hash{stem1} = $params{-mx}->cell(%params);
	$hash{stem2} = scalar reverse ( $params{'-mx'}->cell( %params, '-blk' => $self->helix($params{'-blk'}) )); 
	
	$hash{mask1} = $masks{'5prime'};
	$hash{mask2} = scalar reverse ($masks{'3prime'});
	
	my $pair = "";	
	
	my $p = 0;
	
	do {	
		$pair ="";
		$hash{mask1}  =~ /.{$p}(.)/;
		#print "[1: $1] ";
	        
		$pair .= $1 if defined $1; # this and line below work to delimit problems with non def $1 - but not sure why there should be such!
		$hash{mask2} =~ /.{$p}(.)/;
		$pair .= $1 if defined  $1;  
	
		#print "$p - $pair ";
		
		if ($pair eq '.)') {
			substr ($hash{stem2}, $p ,0 ) = '-'; 
			substr ($hash{mask2}, $p, 0 ) = '.';
			#print " 1";
		}
			
		elsif  ($pair eq '(.') {
			#print " 2";
			substr ($hash{stem1}, $p ,0 ) = '-'; 
			substr ($hash{mask1} , $p, 0 ) = '.';			
		} 
		$p++;

	} until (length $pair == 0); 

	return %hash; 
}



=head2 block_helix_label

A label accessor.

=cut

sub block_helix_label {
	my $self = shift;
	if (@_) {	
		my $s1 = shift; 
		if (@_) { 
			my $s2 = shift;
			$self->{block_helix_labels}->{$s1} = $s2;
			return 1
		}
		else {
			return $self->{block_helix_labels}->{$s1};
		}
	}
	return '';
}

=head2 deleteHelix
 
Deletes a helix.  Unproven.

=cut

sub deleteHelix { 
	my $self = shift;
	my %params = @_;
	delete $self->{helices}->{$params{'-blk'}};
	$self->{_helix_accession}--;
	return 1
}


=head2 threePrimePair
 
Given a 3' helix index return the 5' complement (missing the faster way to reverse hashes- fix).

=cut


sub threePrimePair {
	my $self = shift;
	my %params = @_;
	
	croak "no 3' block passed to 3primePair\n" if not defined $params{'-blk'};
	
	foreach my $five (values %{$self->{helices}}) {
		if  (defined $self->helix($five) ) {
			return $self->helix($five) if $self->helix($five) == $params{'-blk'};
		}
	}
	return 0;
}


sub _deleteBlock {
	# deletes all references to a given block from a structure object
	my $self = shift;
	my %params = @_; # needs -blk

	die "no block passed to _deleteBlock" if not defined $params{'-blk'};
	
 	delete $self->{'block_masks'}->{$params{'-blk'}};
	delete $self->{'helices'}->{$params{'-blk'}};
	delete $self->{'block_labels'}->{$params{'-blk'}};

	# also delete references to the 3' helices
	delete $self->{'helices'}->{$self->threePrimePair(%params)} ;
		
	return 1;	
}


sub _build_helices {
	my $self = shift;
	my %params = @_; # needs -mx, -helix_slice
	
	print "initializing helicies with matrix: ", $params{'-mx'}->label() , "\n";

	# some verification code here
	croak "not a matrix in _build_helices" if not ref($params{'-mx'}) eq "matrix"; ## don't check ref

	$params{'-helix_slice'} = $params{'-mx'}->origSlice if not defined $params{'-helix_slice'};
	
	croak "not a slice in _build_helices" if not ref($params{'-helix_slice'}) eq "slice"; ## don't check ref
	## can bracketed blocks be made to helices?
	 if ($params{'-helix_slice'} ) {
		my ($err1, $err2) = (0) x 2;
		foreach my $blk ($params{'-helix_slice'}->loop("Blocks")) {
			if ($params{'-mx'}->{blocks}->{$blk}->bracketed == 0) {

				foreach my $p ($params{'-helix_slice'}->loop("Blocks")) {
					if ($params{'-mx'}->{blocks}->{$p}->bracketed == 1) { next }
					
					print "[$blk $p] range error - you probably have whitespace at end of stem index?\n" if not defined ($self->block_helix_label($p));

					if ($self->block_helix_label($blk)."'" eq $self->block_helix_label($p)) {
						$self->helix($blk, $p);
	
						# print "[$blk, $p] ";
						# print $self->mask('-blk' => $blk)  , $self->mask('-blk' => $p); 
						# print $self->bracketsInMask('-blk' => $blk), $self->bracketsInMask('-blk' => $p);  
						
						# error check 
						$err1 =  $self->bracketsInMask('-blk' => $blk);
						$err2 =  $self->bracketsInMask('-blk' => $p);
						
						if ($err1 != $err2 ) {
							print "!! ERROR in stem pairing, mismatched masks: "; # See err\/stem_err.txt\n
							print " block $blk has $err1 parens [", $self->mask('-blk' => $blk), "] and block $p has $err2 parens [", , $self->mask('-blk' => $p), "].\n";	
							#print STEMERR "Uneven-masks error in blocks [$deschead[$orig_interleave->{Interleaves}->{$blk}][0]/$deschead[$orig_interleave->{Interleaves}->{$stems{$blk}}][0]] stem mask pair [$blk-$stems{/$blk}] - masks are: $nbd[$blk][8] [$err1 brackets]) and $nbd[$stems{$blk}][8] [$err2 brackets])\n"; 
						}
						if ($err1 == 0 and $err2 == 0) {
							print "ERROR in stem pairing!!!!! Defined helix has no stems. See err\/stem_err.txt\n";
							#print STEMERR "Defined helix has no stems. Interleave: [$deschead[$orig_interleave->{Interleaves}->{$blk}][0]/$deschead[$orig_interleave->{Interleaves}->{$stems{$blk}}][0]] Block pair:[$blk ($sindex{$blk})-$stems{$blk} ($sindex{$stems{$blk}})] - [$err1 brackets]) and [$err2 brackets])\n"; 
						}		
						#$stemtot++;
						last
					}
				}
			}
		}
	 }
	 # check that regions not defined as stems don't have brackets in them and vice versa
}


=head2 loopHelices
 
Returns an array of 5' helix sides within a slice 

!! needs to limit those to fully (3'/5') bound by slice only !!  
## reduntant with mx->loopblocks (keep this one) ##

# slice dependant, it's not a helix unless the passed slice has both blocks?

=cut


sub loopHelices { 

	my $self = shift;
	my %params = @_;


	
	my @helices;
	foreach my $blk ($params{'-slice'}->loop('Blocks')) {
		if ( defined($self->helix($blk))) {
			# print "{    $blk ", $self->helix($blk), "\n";
			push (@helices , $blk) if $params{'-slice'}->containsBlk($self->helix($blk)); 
		}
	}
	return @helices;
}


=head2 loopStems
 
Returns an array of all block indices involved in pairing given the particular -slice

=cut


sub loopStems {
	my $self = shift;
	my %params = @_;
	my @stems;
	
	my @h = $self->loopHelices(%params);
	# print @h;

	foreach my $blk (@h) {
		if (defined $self->helix($blk)) {
			push (@stems, $blk);
	   		push (@stems, $self->helix($blk))
		} 
	}
	return sort @stems;	
}


=head2 loopBasePairs
 
Returns an array of basepairs for a given 5' block index, assumes mask is made of '()'

Requires -mx, -blk, -ter

=cut


sub loopBasePairs { 
	my $self = shift;
	my %params = @_;
	
	my @basepairs;
	return if not $self->helix($params{'-blk'}); # if the -blk isn't a helix there are no basepairs 
	
	my %data = $self->alignedHelices(%params);
	
	for (my $i=0; $i < length $data{mask1}; $i++) {
		if (substr($data{mask1},$i,1) eq '(' ) {
			my $c = substr($data{stem1},$i,1);
			$c .=  substr($data{stem2},$i,1);
			push @basepairs, $c;
		}
	}
	return @basepairs;
}


=head2 loopNonBasePairs 

NOT IMPLEMENTED
Returns an array of the positions of non paring columns 

=cut

sub loopNonBasePairs {

}


=head2 pseudoKnotPairs 

A hack.

Takes an array representing a hash of indices represent helix pairs.

Returns an array pointing to the 5' blocks which, when removed along with their 3' pair, will leave remain blocks psuedoknot free.

Note that there are multiple solutions possible, and that the result returned here is very inefficient (i.e removes the fewest number of helices making the overall structure pseudoknot free).
Someone with more experience will need to solve the most efficient problem.

Use with caution, and possibly in combination with manually selected helicies.

=cut

sub pseudoKnotPairs {
    my $self = shift;
    my %p = @_;
    my @r;

	# print "\n\n";
    
    my $hel;
    my $h = 0;
    my @keys = sort {$a <=> $b} keys %p;

    map { $hel->{$h} = [$_, $p{$_}]; $h++} (@keys);   
    foreach my $k (sort {$a <=> $b} keys %{$hel}) {
		# print $k, " ", @{$hel->{$k}}[0], " ",  @{$hel->{$k}}[1], ": ";
        my $nested = 1;
        my $cl = @{$hel->{$k}}[1];
            
        for (my $c = $k+1; $c < $#keys + 1; $c++) {    
			#   print " $cl";  
            if (@{$hel->{$c}}[1] > @{$hel->{$k}}[1]) { 
                if (@{$hel->{$k}}[1] > @{$hel->{$c}}[0] ) { $nested =0}              
                last;
            }
            if ( $cl <  @{$hel->{$c}}[1]) {
                $nested = 0; 
            }
            $cl = @{$hel->{$c}}[1];     
        }

        if ($nested == 1) {
			# print "\tok!"; 
        }
        else {
			push @r, @{$hel->{$k}}[0];
			# print "\tremove me!";
        }
		#   print "\n";
    }
	# print "\n\n";

    return @r;
}

	
=head2 loopPairs
 
Returns a array of 1-n where n is the number of brackets.

Requires -mx, -blk 

=cut


sub loopPairs { 
	my $self = shift;
	my %params = @_;
	my $t = $self->bracketsInMask(%params);
	return [1..$t]
}


=head2 pairPos

Returns a two value array of column indices (unbracketed position) for the -pos position of the helix defined by the 5' -blk.
Requires -blk, -mx, -pos.
Similar to nucPairs but notably counts from ZERO.  Primarily for internal use.


=cut

sub pairPos {
	my $self = shift;
	my %params = @_;

	return if not $self->helix($params{'-blk'}); # if the -blk isn't a helix there are no basepairs 
	$params{'-pos'} || die "no -pos to pairPos\n";
	
	($self->bracketsInMask(%params) < $params{'-pos'}) && die 'requested pair greater than pairs present';
	
	my ( $val, @r, @l  );
	my $mx = $params{'-mx'};

	my @lefts;
	my @rights;
	my $p = 0;
	@lefts = map { $p++; $p + $mx->{blocks}->{$params{'-blk'}}->excludedStart - 1  if $_ eq '('} ( split(//, $self->mask(%params) ));
	$p = 0;
	@rights = map { $p++; $mx->{blocks}->{ $self->helix($params{'-blk'}) }->excludedEnd  + 1  - $p if $_ eq ')'} ( split(//, reverse ( $self->mask('-blk' => $self->helix($params{'-blk'}) ) )));	

	foreach my $val (@lefts) {
		push (@l, $val) if length $val > 0
	} # strip nulls from map, should be an easier way
	foreach $val (@rights) {
		push (@r, $val) if length $val > 0
	}

	return ($l[$params{'-pos'}-1], $r[$params{'-pos'}-1] );  
}


sub _load {
	my $self = shift;
	my %params = @_;

	#$slice '-slice'
	# -warnings (global)
	
	#my ( 
	#	$input,
	#    $slice,
	#	$warnings, # 0- don't die on warnings; 1- die on warnings
	#) = @_;

	my $tm;

	# chdir($rootdir);
	# &io_confirmdir("err");
	# chdir($rootdir);
	# open (STEMERR, ">err/stem_err.txt") || die print "can't open stem_err file";

	#read/error check stem index
	my %uniques_counter;
	
	croak "no index in _load" if not defined $params{'-helix_index'};
	
	open (IN, ($params{'-path'}.$params{'-helix_index'})) || die ("can't open ", $params{'-path'}.$params{'-helix_index'} );# ERROR HANDLING

	#print "\nopened stem index ", $params{"-path"}.$params{"-helix_index"}, "\n";
	
	while (<IN>) {
		chomp;
		my ($key, $val) = split;
		foreach ($key, $val) { die "undefined pair in stem index" if not defined};
			if (not $val eq '?') { $uniques_counter{$val}++; die "\nduplicate name $val in stem index, you should use probably use '?' instead.\n\n" if $uniques_counter{$val} > 1;} # allows only '?' to be duplicated- possibly extend
			#print $key, $val;
			$self->block_helix_label($key, $val);
		}
	close (IN);
	# close STEMERR;

	return 1;
	
	# if (-z "err/stem_err.txt" ) { print "No detectable errors in stem-index parsing.\n"; }
	# else { die "ERRORS PRESENT!  See: err/stem_err.txt" if $warnings==1; }
}


1;


