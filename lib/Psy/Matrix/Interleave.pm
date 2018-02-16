
package interleave; 

use vars qw(@ISA @EXPORTER @EXPORT_OK $VERSION);

$VERSION = '0.01';

use strict;
use warnings;
use Data::Dumper;
use Carp;

use Psy::Psy qw(generatedByHeader version rootpath projectname starttime $ROOT_DIR $PSY_BASE_DIR $PSY_ROOT_DIR);
@ISA = qw(Psy); ## ??
  
=head1 NAME

Psy::Matrix::Interleave

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

	Interleave objects.  Stores information regardin the distribution of blocks across interleaves.

    use Psy::Matrix::Interleave;	
	my $foo = interleave->new();

=cut


sub new {
	my $type = shift;	
  	my $self = {};
	$self->{Interleaves} = {};	
	bless $self, $type;

	$self->_init();
	
    return $self; 
}

sub _init {
	my $self= shift;
	my %params = @_;
	
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

sub label {
	my $self = shift;
	my %params = @_;
	if ($params{'-value'}) {$self->{labels}->{$params{'-interleave'}}->{$params{'-label'}} = $params{'-value'} }
	else { $self->{labels}->{$params{'-interleave'}}->{$params{'-label'}} };
}


sub assign {
	my $self = shift;
	my $interleave = shift;
	foreach (@_) { $self->{Interleaves}->{$_} = $interleave};
}

sub blockFrom { # return the interleave the block is placed in
	my $self = shift;
	my %params = @_;
	foreach my $i ($self->interleaves) {
		if (($params{'-blk'} <= $self->{interleaves}->{$i}->last_block_in_interleave) and ($params{'-blk'} >= $self->{interleaves}->{$i}->first_block_in_interleave)) {
			return $i;
		};
	}
	return -1; # can't find the block index
}

sub remove_blocks {
	my $self = shift;
	my @removals = @_;
	foreach my $blk (keys %{$self->{Interleaves}}) {
		delete $self->{Interleaves}->{$blk} unless grep ($_ == $blk, @removals);
	 };
}

sub interleaves { 
	my $self = shift;
	my %uniques;
	foreach my $vals (values %{$self->{Interleaves}}) { $uniques{$vals} = undef};
	return (sort {$a <=> $b} (keys %uniques));
}

sub total_interleaves {
	my $self = shift;
	my %uniques;
	foreach my $vals (values %{$self->{Interleaves}}) { $uniques{$vals} = undef};
	return int(keys %uniques);
}

sub last_block_in_interleave {
	my $self = shift;
	my $intlv = shift;
	my @blks;
	map {push (@blks, $_) if $self->{Interleaves}->{$_} == $intlv } (keys %{$self->{Interleaves}}); 
	my @sorted = sort {$a <=> $b} @blks;
	return pop @sorted;
}

sub first_block_in_interleave {
	my $self = shift;
	my $intlv = shift;
	my @blks;
	map {push (@blks, $_) if $self->{Interleaves}->{$_} == $intlv } (keys %{$self->{Interleaves}}); 
	my @sorted = sort {$a <=> $b} @blks;
	return shift @sorted;
}

sub right {
	my $self = shift;
	my $intlv = shift;
	my @i = $self->interleaves;
	my $p=0;
	foreach my $c (@i) {
		if ($intlv == $c) {
			if (defined $i[$p+1]) { return $i[$p+1] }
			else {return undef}
		}
		$p++;
	}
}

sub left {
	my $self = shift;
	my $intlv = shift;
	my @i = $self->interleaves;
	my $p=0;
	foreach my $c (@i) {
		if ($intlv == $c) {
			if ($self->first_interleave == $c) { return undef }
			elsif (defined $i[$p-1]) { return $i[$p-1] }
			else { return undef }
		}
		$p++;
	}
}

sub blocks {
	my $self = shift;
	my $intlv = shift;
	my @blks;
	map {push (@blks, $_) if $self->{Interleaves}->{$_} == $intlv } (keys %{$self->{Interleaves}}); 
	return (sort {$a <=> $b} @blks);
}

sub first_interleave {
	my $self = shift;
	my @array = $self->interleaves;
	return shift @array;
}


sub last_interleave {
	my $self = shift;
	my @array = $self->interleaves;
	return pop @array;
}

sub describe {
	my $self = shift;
	print "\ninterleave(s)\n";
	print "interleaves: ", join " ", $self->interleaves,"\n";
	foreach my $i ($self->interleaves) {
		print "interleave: $i \t";
		foreach my $blk ($self->blocks($i)) {
			print " $blk ";
		}
		print "\n";
	}
}


# STRIKE THIS - use intlv_build

	#sub intlv_working { # returns an interleave based on user input
		#my $self = shift;
		#my %params = @_;
				
		# my $working_interleave = interleave->new;
		#compute interleave
	#	if ($mode == 2) { # return an interleaved results with interleave size $outintlvsize
	#		$working_interleave->intlv_build(	
	#											%params,
	#											'-asfasf' => ' '
	#												$slice, 0, $size 
	#										 );
	#	}
			
	#	elsif ($mode == 1) { # return format as in infile
	#		print "WARNING!  You are using/modifying the original interleave format, and may be excluding characters! Check the output.\n";
	#		$model_intlv->remove_blocks(keys %{$slice->{Blocks}});
	#		$working_interleave = $model_intlv;
	#	}		
		
	#	else { # default to not interleaved
	#		$working_interleave->intlv_build($slice, 1, 1);
	#	}
		
	#	return $working_interleave;
	#}

	# NOT WORKING YET - is a merger of build/working 
	sub intlv_build () { # populate an interleave object, ultimately pass it a matrix too (nbd part) # required:	
		my $self = shift; 	#interleave object	
		my %params = @_;

		croak "not a slice in intlv_build\n" if not ref($params{'-slice'}) eq 'slice';
		croak "not a mx in intlv_build\n" if not ref($params{'-mx'}) eq 'mx';
		croak "not a slice in intlv_working\n" if not ref($params{'-slice'}) eq 'slice';

		$params{'-build_mode'} = 'columns' if not defined $params{'-build_mode'};	# valid values = 'columns' | 'divisor'	
		$params{'-mode'} = 'single' if not defined $params{'-mode'}; 				# options are 'single', 'original', 'sized' 
		
		if ( not defined $params{'-size'} ) {	$params{'-size'} = 140; print "no size passed to working interleave, setting to 140\n";}

		$params{'-divisor'} = 140 if not defined $params{'-build_mode'};
		
		# my $working_interleave = interleave->new;

	#	my (	
	#		$slice,		# reference to a slice object
	#		$mode,		# 0 - columns; 1 - divisor
	#		$divisor	# mode == 0 - max number of columns; mode == 1 - number of interleaves to divide dataset into
	#	) = @_;
		
		my ($size);
		
		if ($params{'-build_mode'} eq 'columns') { # ($mode == 0) { 
			$size = $params{'-divisor'}; 
			}
		elsif ($params{'-build_mode'} eq 'divisor' ) { # $mode == 1) {
			my $totlen = 0;
			foreach my $blk ($params{'-slice'}->loop("Blocks")) {
				$totlen += $params{'-mx'}->{blocks}->{$blk}->blk_length;   # $totlen += $nbd[$blk][0];
			}
			print $totlen;
			$size = int( $totlen / $params{'-divisor'} ) + 1
		}	
		
		my $curintlv = 0;
		my $curlen = 0;

		foreach my $blk ($params{'-slice'}->loop("Blocks")) {
			if (($params{'-mx'}->{blocks}->{$blk}->blk_length + $curlen) > $size) { 
				$curintlv++; 
				$self->assign($curintlv, $blk);
				$curlen = $params{'-mx'}->{blocks}->{$blk}->blk_length;
			}
			else {
				$self->assign($curintlv, $blk);
				$curlen += $params{'-mx'}->{blocks}->{$blk}->blk_length
			}				
		}
	}


1;

__END__
