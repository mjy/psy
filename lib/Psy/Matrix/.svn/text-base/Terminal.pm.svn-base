package terminal;

use vars qw(@ISA @EXPORTER @EXPORT_OK $VERSION);

use strict;
use warnings;

use Carp;
use Data::Dumper;

use Psy::Strings::Strings;
use Psy::Psy;
@ISA = qw(Psy);

=head1 NAME

Psy::Matrix::Column

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

	A terminal object.

    use Psy::Matrix:Terminal;	
	my $foo = terminal->new();

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

	$self->_verify(%params);	
	
	$self->{label} = $params{'-label'};
	
	$self->version($VERSION);
	return 1;
}
	
sub label {
	my $self = shift;
	if (@_) {
		$self->{label} = shift;
	}
	return $self->{label};
}

sub safeName { # takes an arguement as base of name
	my $self = shift;
	
	if (not defined $self->{safename} and @_) {
		my $prefix = shift;
		my $safename = $self->label;
		
		$safename =~ s/(?)[^\w]//g; #strip non-alphanumerics
		$safename =~ s/_//g; # strip underscore (POY doesn't like in trees)
	
		$safename = substr ($safename,0,7); #grab the first seven chars
		$safename .= &strings::str_padleft($prefix ,"0" , 3);
	
		$self->{safename} = $safename; 		
	}
	elsif (not defined $self->{safename} and not @_) {croak "safename not previously generated!"}
	
	return $self->{safename};
}

sub _verify {
	my $self=shift;
	my %params = @_;
		foreach my $p (keys %params) { # check matrix parameters here
	}
		
	return %params;
}

# write sub ter2genbank_name {}

1;

__END__


####
#

#    sub AUTOLOAD {
#               my $self = shift;
#               my $type = ref($self)
#                           or croak "$self is not an object";
#
#               my $name = $AUTOLOAD;
#               $name =~ s/.*://;   # strip fully-qualified portion
#
#               unless (exists $self->{_permitted}->{$name} ) {
#                   croak "Can't access `$name' field in class $type";
#               }
#
#               if (@_) {
#                   return $self->{$name} = shift;
#               } else {
#                   return $self->{$name};
#               }
#           }

