package wrapper;

use vars qw(@ISA @EXPORTER @EXPORT_OK $VERSION);

$VERSION = '0.01';

use strict;
use warnings;
use Data::Dumper;
use Carp;


=head1 NAME

Psy::Output::Wrapper

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.02';

=head1 SYNOPSIS

NOT PRESENTLY USED ANYWHERE

Referenced internaly only. Used to wrap text output with. 
	
=cut

sub new {
	my $type = shift;
	my %params = @_;	
  	my $self = {};
	
	bless $self, $type;

	$self->_init;

	#print "\n[", ref ($self) ,"]\n";;
    return $self;                 
}

sub _init {
	my $self = shift;
	$self->{'row'}->{'pre'} = "";
	$self->{'row'}->{'post'} = "";
	$self->{'taxa'}->{'pre'} = "";
	$self->{'taxa'}->{'post'} = "";
	$self->{'block'}->{'pre'} = "";
	$self->{'block'}->{'post'} = "";
	$self->{'bracketed_block'}->{'pre'} = "";
	$self->{'bracketed_block'}->{'post'} = "";

	return 1;
}

sub row {
	my $self = shift;
	croak 'no position passed' if not (@_);
	my $pos = shift;
	if (@_) { $self->{row}->{$pos} = shift; }

	return $self->{row}->{$pos};
}

sub taxa {
	my $self = shift;
	croak 'no position passed' if not (@_);
	my $pos = shift;
	if (@_) { $self->{taxa}->{$pos} = shift; }

	return $self->{taxa}->{$pos};
}

sub block {
	my $self = shift;
	croak 'no position passed' if not (@_);
	my $pos = shift;
	if (@_) { $self->{block}->{$pos} = shift; }

	return $self->{block}->{$pos};
}

sub bracket_block {
	my $self = shift;
	croak 'no position passed' if not (@_);
	my $pos = shift;
	if (@_) { $self->{bracketed_block}->{$pos} = shift; }

	return $self->{bracketed_block}->{$pos};
}

1;

__END__
