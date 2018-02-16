package fileformat;

use vars qw(@ISA @EXPORTER @EXPORT_OK $VERSION);

$VERSION = '0.01';

use strict;
use warnings;
use Data::Dumper;

use Psy::Psy;

# @ISA = qw(Psy);
 
#my $foo = type->new();
#print Dumper($foo);
#print $foo->formats();

use Exporter;
@EXPORT_OK = qw ();

=head1 NAME

Psy::Matrix::Format

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

	A matrix-format object.  A hack used for detection of filetypes, storing of file-specific meta-data.

    use Psy::Matrix::Slice;	
	my $foo = format->new();

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
		$self->{"nexus"} -> {"data_start"} = 'matrix';
		$self->{"nexus"} -> {"identifier"} = 'nexus';
		$self->{"nexus"} -> {"data_end"} = ';';
		$self->{"nexus"} -> {"mask"} = "[mask";
		$self->{"nexus"} -> {"descriptor1"} = '[Block_';
		
		$self->{"stockholm"} -> {"identifier"} = 'stockholm';
		$self->{"stockholm"} -> {"data_start"} = 'stockholm';
		$self->{"stockholm"} -> {"data_end"} = '//';
		$self->{"stockholm"} -> {"mask"} = '#=GC SS_cons';
	return 1;
}
	
sub _formats {
	my $self= shift;
	return keys %{$self};
}

sub _identifiers {
	my $self= shift;
	return keys %{$self};
}

sub identifier {
	my $self = shift; 
	my $format = shift;
	return $self->{$format}->{'identifier'};
}

sub mask {
	my $self = shift; 
	my $format = shift;
	return $self->{$format}->{'mask'};
}

sub detect {
	my $self = shift;
	my $allfile = shift;
	foreach my $f ($self->_formats) {
		my $id = $self->identifier($f);
		return $f if $allfile =~ /$id/gi;	# simplifiable to 1 line?
	}
	return '';
}

sub data_end {
	my $self = shift;
	my $format = shift;
    if (@_) { $self->{$format}->{data_end} = shift ;}
	return $self->{$format}->{data_end};
}

sub data_start {
	my $self = shift;
	my $format = shift;
    if (@_) { $self->{$format}->{data_start} = shift ;}
	return $self->{$format}->{data_start};
}


## SUBS BELOW HERE ARE NOT IMPLEMENTED AT PRESENT

sub format {
	my $self = shift;
    if (@_) { $self->{format}->{$_} = {} ;}
}

sub start {
	my $self = shift;
    if (@_) { $self->{start}->{$_} = {} ;}
}

sub descriptor1 {
	my $self = shift;
    if (@_) { $self->{format}->{$_} = {} ;}
}

sub _verify {
	my $self=shift;
	my %params = @_;
		foreach my $p (keys %params) { # check matrix parameters here
	}
		
	return %params;
}

1;


