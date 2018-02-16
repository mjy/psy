package block;

use strict;
use warnings;
use Data::Dumper;
use Carp;

use vars qw(@ISA @EXPORTER @EXPORT_OK $VERSION);

use Psy::Psy qw(generatedByHeader version rootpath projectname starttime $ROOT_DIR $PSY_BASE_DIR $PSY_ROOT_DIR);
@ISA = qw(Psy); ## hmm...

=head1 NAME

Psy::Matrix::Block

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Contains meta-data for blocks.  Internal use only for all intents and purposes.
	
    use Psy::Matrix::Block
    my $foo = block->new; # see _init


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

	$self->{literal_start} = $params{'-literal_start'}; 
	$self->{literal_end} =  $params{'-literal_end'};
	$self->{excluded_start} =  $params{'-excluded_start'};  
	$self->{excluded_end} =  $params{'_excluded_end'}; 
	$self->{bracketed} =  $params{'-bracketed'}; 
	$self->{blk_length} =  $params{'-blk_length'};
	$self->{labels} = ();
	
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




=head2 label

Accessor for a -label of -blk.

=cut


sub label {
	my $self=shift;
	my %params = @_;
	if ($params{'-label'}) { 
		$self->{labels}->{$params{'-label_name'}} = $params{'-label'};
	}
	return 	$self->{labels}->{$params{'-label_name'}} ;
}




=head2 literalStart

Accessor for the starting (nucleotide) position of this block, ** as originally read from the matrix **.
All blocks are included.

=cut



sub literalStart {
	my $self = shift;
    if (@_) { $self->{literal_start} = shift ;}
	return $self->{literal_start};
}


=head2 literalEnd

As literalStart.
All blocks are included.

=cut


sub literalEnd {
	my $self = shift;
    if (@_) { $self->{literal_end} = shift ;}
	return $self->{literal_end};
}


=head2 excludedStart

Accessor for the starting position of this block, as originally read from the matrix.
Bracketed blocks are _excluded_ in calculating this position.

=cut


sub excludedStart {
	my $self = shift;
    if (@_) { $self->{excluded_start} = shift ;}
	defined $self->{excluded_start} ? return $self->{excluded_start} : return -1  ;
}


=head2 excludedEnd

Ad excludedStart.
Bracketed blocks are _excluded_ in calculating this position.

=cut


sub excludedEnd {
	my $self = shift;
    if (@_) { $self->{excluded_end} = shift ;}
	defined $self->{excluded_end} ? return $self->{excluded_end} : return -1;
}


=head2 bracketed

Accessor for bracketed property.  1 == bracketed

=cut


sub bracketed { # 0/1
	my $self = shift;
    if (@_) { $self->{bracketed} = shift;}
	return $self->{bracketed}  
}



=head2 blkLength

Length including gaps.

=cut


sub blkLength {
	my $self = shift;
    if (@_) { $self->{blk_length} = shift ;}
	return $self->{blk_length};
}




=head2 debugDescribe

Prints a line containing the various metadata.

=cut


sub debugDescribe {
	my $self = shift;
	print "[", join " ", ( $self->literalStart, $self->literalEnd, $self->excludedStart, $self->excludedEnd, $self->blkLength, $self->bracketed), "]" ;
   	print "\n";
}

1;


