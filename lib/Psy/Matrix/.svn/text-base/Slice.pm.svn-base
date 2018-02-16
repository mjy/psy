
package slice;

use vars qw(@ISA @EXPORTER @EXPORT_OK $VERSION);

use strict;
use warnings;
use Data::Dumper;
use Carp; 

use Psy::Psy;  #qw(version); ## used?
use Psy::Strings::Strings;
@ISA = qw(Psy);

=head1 NAME

Psy::Matrix::Slice

=head1 VERSION

Version 0.01

=cut


our $VERSION = '0.01';


=head1 SYNOPSIS

A slice object.  Slices are essentially character and/or terminal partitions (sets).  Characters are used in the sense of blocks, not columns of data (though a block can be length 1).  
The object uses hashes rather than arrays so that set operations are easy.

	use Psy::Matrix::Slice;	
	my $s = slice->new();

=cut

=head2 new

=cut


sub new {  
	my $type = shift;	
  	my $self = {};
	my %params = @_;
	
	bless $self, ref($type) || $type ;
	
	$self->_init(%params);
    return $self;                 
}

sub _init {
	my $self= shift;
	my %params = @_;
		
  	$self->{Kind} = "";	# "include" or "exclude"
	$self->{Blocks} = {};
	$self->{Taxa} = {};	 	

	$self->version($VERSION);
	return 1;
}
	
sub _verify {
	my $self=shift;
	my %params = @_;
		foreach my $p (keys %params) { # check matrix parameters here
	}
		
	return %params;
}


=head2 clone

Returns a clone of the object

=cut

sub clone {
	my $slice = shift;
	my $self = $slice->new(%$slice);

	$self->blocks( keys %{$slice->blocks} );
	$self->taxa( keys %{$slice->taxa} );
	
	return $self;
}


	# not implemented
	sub kind { 
		my $self = shift;
		if (@_) {my $val = shift; if (($val eq "include") || ($val eq "exclude")) {$self->{Kind} = $val } else { } }
		return $self->{Kind};
	}
	

=head2 blocks

Adds/accesses block has (indicies, not objects).  Returns a hash, not an array.
Usage:	
	$s->blocks(0, 1, 2, 45..65, 3);
	@foo = $s->blocks;


=cut


sub blocks {
	my $self = shift;
	if (@_) { foreach (@_) { $self->{Blocks}->{$_} = undef}}
	return $self->{Blocks};
}

=head2 taxa
	
As blocks but for taxa.

=cut


sub taxa {
	my $self = shift;
	if (@_) { foreach (@_) { $self->{Taxa}->{$_} = undef}}
	return $self->{Taxa};
}


=head2 remove
	
Removes blocks OR taxa.

Usage:

	my @remove_me = (3, 9, 10, 11, 13, 14..25);
	$s->remove('Blocks', @remove_me);	
OR

	$s->remove('Taxa', 3, 9, 10, 2..4);  # note overlaps (e.g. 3) don't fail, and order doesn't matter


=cut


sub remove {
	my $self = shift;
	my $kind = shift;
	if (($kind eq "Taxa") || ($kind eq "Blocks")) {
		if (@_) { foreach (@_) { delete $self->{$kind}->{$_}} }
	}
	else { die "can't remove $kind from slice object"}
	1;
}

=head2 loop
	
Returns ordered list of taxa OR blocks as array.

Usage:

	my @foo = $s->loop('Taxa');
	OR
	for my $blk ($s->loop('Blocks') { ... }


=cut


sub loop {
	my ($self, $kind) = @_;
	my @array = keys %{$self->{$kind}};
	return (sort {$a <=> $b} @array);
}
	
=head2 union
	
Returns a union of two hashes, that are referenced like $first->{$hash}, $second->{$hash}.

Used internally like $s1->{'Blocks'}, $s2->{'Blocks'}

=cut


sub union  { ## make internal?
	my $first = shift;
	my $second = shift; 
	my $hash = shift;
	my %union = (%{$first->{$hash}}, %{$second->{$hash}});
	return %union;
}

=head2 difference
	
As union, but returns difference.

=cut


sub difference {
	my $first = shift;
	my $second = shift; 
	my $hash = shift;
	
	my %result = %{$first->{$hash}};
	delete $result{$_} for keys %{$second->{$hash}};
	return %result;
}

=head2 total
	
Returns total blocks/taxa.
Usage: 	

	print "total blocks: ", $s->total('Blocks'); 


Not zero beers.

=cut


sub total { # total blocks/taxa, not zero beers
	my $self = shift;
	my $kind = shift;
	my @array = keys %{$self->{$kind}};
	return $#array + 1;
}

=head2 first

Return the first block OR taxon in the slice (as ordered by integer index).

=cut


sub first {
	my $self = shift;
	my $kind = shift;
	my @array = sort {$a <=> $b} (keys %{$self->{$kind}});
	return (shift @array)
}

=head2 final

Return the final ('last' is reserved word) block OR taxon in the slice (as ordered by integer index).

=cut


sub final {
	my $self = shift;
	my $kind = shift;
	my @array = sort {$a <=> $b} (keys %{$self->{$kind}});
	return (pop @array)
}


=head2 lengthLastBlk

Returns the string length of the last block (e.g. 2000 has length 4). 

=cut


sub lengthLastBlk { ## edited from length_last_block
	my $self =shift;
	my @array = sort {$a <=> $b} (keys %{$self->{Blocks}});
	my $blk = pop @array;
	return length("$blk");
}	my %raw_params = @_;


=head2 describe

Returns a text description of the slice, formatted for human reading.

=cut


sub describe {
	my $self = shift;
	print "\nslice\n";
	print "total tax: ", $self->total("Taxa"), "\ttotal blocks: ", $self->total("Blocks"), "\n";
	print "taxa:\t", join (" ", sort {$a <=> $b} keys %{$self->taxa}), "\n";
	print "blocks:\t", join (" ", sort {$a <=> $b} keys %{$self->blocks}),"\n"	
}



=head2 decribeCoded

Returns a text description of the slice, formatted for eval() or use in batch-file.

=cut


sub decribeCoded { # was describe_coded
	my $self = shift;
	my $s = 'my $slice = slice->new()';
	$s .= "\n";
	$s .= '$slice->blocks';
	$s .= &arrayAsRange($self->loop('Blocks'));
	$s .= "\n";
	$s .= '$slice->taxa';
	$s .= &arrayAsRange($self->loop('Taxa'));
	$s .= "\n";
	return $s
}
	
=head2 random

Returns an index to a randomly selected taxa OR block.

=cut

sub random {
	my $self = shift;
	my $type = shift;	

	my @keys = keys %{$self->$type};
	return $keys[int rand(@keys)]; 
}

=head2 newRandomSlice

Returns a reference to a new slice object with randomly selected blocks and taxa.

Usage:

	$foo = $s->new_random_slice();
	OR 
	$foo = $s->new_random_slice("numtaxamax" => 5, "numtaxamin"=> 5, "numblocksmax" => "all" );


=cut


sub newRandomSlice {		
	my $self = shift;
	my %params = @_; # num[taxa|blocks][min|max] => #|'all';
	my $err;
	
	my $new_slice = slice->new;

	# check for initialization
	foreach my $p (qw/numtaxa numtaxamin numblocksmax numblocksmin/) {
		$params{$p} = 'all' if not defined $params{$p};
	}	

	print "\n";
	
	if (not ($params{"numtaxamax"} eq 'all')) {
		if ( $params{"numtaxamax"} > $self->total("Taxa")  ) { $err .= ' number of Taxa greater than total in matrix'; }
	}
	
	if (not ($params{"numblocksmin"} eq 'all')) {	
		if ( $params{"numblocksmax"} > $self->total("Blocks")) { $err .= ' number of Blocks greater than total in matrix'; }
	}
	
	die "$err - usage: slice->new_random_slice( [ num[taxa|blocks][min|max] ]=>#|'all')" if defined $err;
	
	my ($numtaxa, $numblocks ) = (0) x 2;


	# generate the slice
	if ($params{numtaxamax} eq "all" ) {
		$new_slice->taxa(keys %{$self->taxa});
	}
	else {
		$numtaxa = int(rand ($params{numtaxamax} - $params{numtaxamin} +1)) + $params{numtaxamin};
		while ($numtaxa > $new_slice->total("Taxa") ) {
			$new_slice->taxa( $self->random("taxa") );
		}	
	}

	if ($params{numblocksmax} eq "all") {
		$new_slice->blocks(keys %{$self->blocks});
	}
	else {
		$numblocks = int(rand ($params{numblocksmax} - $params{numblocksmin} +1)) + $params{numblocksmin};		
		while ($numblocks > $new_slice->total("Blocks") ) {
			$new_slice->blocks( $self->random("blocks") );
		}	
	}	
	return $new_slice;
}


=head2 pruneUninformativeTaxa

Somewhat missnamed, removes taxa with less than -cutoff % (0.00 - 1.00 data) from $self.

Must pass a matrix object as -mx

=cut 

sub pruneUninformativeTaxa {  # not proofed yet
	my $self = shift;
	my %default_params = ('-cutoff' => 0); # % data from 0-1, taxa with <= cutoff are excluded
	my %raw_params = @_;

	my %params = (%default_params, %raw_params);
	
	croak '-mx not a matrix in matrix' if not ref($params{'-mx'}) eq 'matrix';
	my $tmpdata;
	print "\npruning uniformative taxa, cutoff: ", $params{'-cutoff'}, "\n", $self->total("Taxa"), " starting taxa\n";
		
	foreach my $t ($self->loop("Taxa")) {
			$tmpdata ="";
			foreach my $blk ($self->loop("Blocks")) {
				$tmpdata .= $params{'-mx'}->yx($t, $blk);
			}
			my $l = length($tmpdata);
			
			$tmpdata =~ s/(?)[\-N\?]//gi; ## perhaps alphabet->clean ultimately
			my $d =  (1- ($l - length ($tmpdata) ) / $l);
			if (  $d  <= $params{'-cutoff'} ) {
				$self->remove("Taxa", $t);
				print "nuking taxon ", $params{'-mx'}->ter($t)->label, " [%";
				printf("%.2f", $d * 100);
				print "]\n";
			}	
		}
	print $self->total("Taxa"), " remaining taxa.\n\n";
	1;
}


=head2 collapse

Collapses slices based on adjoining bracketed or unbracketed blocks.  Returns a HOH that looks like this:

	'0' => {
			 'type' => 0,
			 'slice' => bless( {
								 'Blocks' => {
											   '0' => undef
											 },
								 'VERSION' => '0.01',
								 'Taxa' => {},
								 'Kind' => ''
							   }, 'slice' )
		   },
	'1' => {
			 'type' => 1,
			 'slice' => bless( {
								 'Blocks' => {
											   '1' => undef,
											   '2' => undef
											 },
								 'VERSION' => '0.01',
								 'Taxa' => {},
								 'Kind' => ''
							   }, 'slice' )
		   }
		   ...


=cut 


sub collapse {
	my $self= shift;
	my %raw_params = @_;
	my %default_params = ( );
	
	my %params = (%default_params, %raw_params);
	(my $mx = $params{'-mx'}) || croak('no -mx passed to slice->collapse'); # simplify things
		
	my $r; # the result hash to return
	my @blks = $self->loop('Blocks');
	my $i = 0;
	
	$r->{$i}->{'slice'} = slice->new;
	$r->{$i}->{'slice'}->blocks($blks[0]);
	$r->{$i}->{'type'} = $mx->blk(0)->bracketed;	
	shift @blks; # ditch block 0
	
	foreach my $blk (@blks) {	
		if ( $r->{$i}->{'type'} == $mx->blk($blk)->bracketed) {
			$r->{$i}->{'slice'}->blocks($blk);
		}
		else { # initiate a new partition
			$i++;
			$r->{$i}->{'slice'} = slice->new;
			$r->{$i}->{'slice'}->blocks($blk);
			$r->{$i}->{'type'} = $mx->blk($blk)->bracketed;
		}
	}
	return $r
}

=head2 describeCollapsed

Not really an object method, but called as such for now.  Describes the breakdown of the results of $s->collapse 

=cut


sub describeCollapsed {	
	my $self = shift;
	my $c = shift; # the results of $self->collapse
	my ($bracketed, $unbracketed) = (0) x 2;
	foreach my $k (sort keys %{$c}) {
		$c->{$k}->{type} == 0  ?  $unbracketed++ : $bracketed++;
		my $h = $c->{$k}->{slice}->blocks;
		print "partition $k: ", join " ", keys( %{$h}), "\n";
	}
	print $unbracketed + $bracketed , " partitions. $unbracketed unbracketed, $bracketed bracketed\n";

}

=head2 prune

Removes bracketed or unbracketed blocks from the slice.  
	
	-modes => < [bracketed] | unbracketed >

=cut


sub prune { 
	my $self = shift;
	my %raw_params = @_;
	my %default_params = ('-mode' => 'bracketed'); # make compatible with mx->loopBlocks

	my %params = (%default_params, %raw_params);

	$self->remove('Blocks', $params{'-mx'}->loopBlocks(%params));
}


=head2 containsBlk

Returns true if it contains a block at the passed index (not block object)
	
=cut


sub containsBlk {
	my $self = shift;
	my $blk = shift;
	die 'no block index to slice->containsBlk' if not defined $blk; # 0 must test true
	(grep /$blk/, keys %{$self->{Blocks}}) && return 1;
	0
}

=head2 unexcludedBlkStart

The following 4 methods are hacks to take a matrix, and a block index and return a start/end value.  See structure->nucPairs for usage.

Requires:
	-mx
	-blk_index

=cut

sub excludedBlkStart {
	my $self = shift;
	my %params = @_;

	return undef if $self->containsBlk( $params{'-blk_index'} ) != 1;
	
	my $p = 0;
	for my $blk ($self->loop('Blocks')) {
		if (($params{'-mx'}->blk($blk)->bracketed == 0) and (not $blk == $params{'-blk_index'})) { # kludge
			($p += $params{'-mx'}->blk($blk)->blkLength ) ;
		}
		last if $blk == $params{'-blk_index'}; # leave early 
	}
	$p;
}

sub excludedBlkEnd {
	my $self = shift;
	my %params = @_;

	return undef if $self->containsBlk($params{'-blk_index'}) != 1;
	my $p = 0;
	for my $blk ($self->loop('Blocks')) {
		if (($params{'-mx'}->blk($blk)->bracketed == 0)) { # kludge
			($p += $params{'-mx'}->blk($blk)->blkLength ) ;
		}
		last if $blk == $params{'-blk_index'}; # leave early 
	}
	$p - 1;
}

sub unexcludedBlkStart {
	my $self = shift;
	my %params = @_;

	return undef if $self->containsBlk($params{'-blk_index'}) != 1;
	my $p = 0;
	for my $blk ($self->loop('Blocks')) {
		if (not $blk == $params{'-blk_index'}) { # kludge
			($p += $params{'-mx'}->blk($blk)->blkLength ) ;
		}
		last if $blk == $params{'-blk_index'}; # leave early 
	}
	$p;
}

sub unexcludedBlkEnd {
	my $self = shift;
	my %params = @_;

	return undef if $self->containsBlk($params{'-blk_index'}) != 1;
	my $p = 0;
	for my $blk ($self->loop('Blocks')) {
		$p += $params{'-mx'}->blk($blk)->blkLength;	
		last if $blk == $params{'-blk_index'}; # leave early 
	}
	$p -1;
}

1;


__END__

