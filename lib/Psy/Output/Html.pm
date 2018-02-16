=head1 NAME

Psy::Output::Html


=head1 VERSION

Version 0.01

=cut

$VERSION = '0.01';

=head1 SYNOPSIS

Html output.  

=cut

package Html;

use strict;
use warnings;
use Data::Dumper;
use Carp;
 use diagnostics;

use vars qw(@ISA);
use Psy::Dna::Alphabet;

@ISA = qw(output);

our $OUT_DIR = "html/";
our $LEGAL_CHARS =  Psy::Dna::Alphabet->new(); # defaults to the RNA alphabet
my $path = $OUT_DIR;

sub process {
	my $self = shift;
	my %raw_params = @_;
	my %default_params = (
		'-reports' => [qw/block_patterns pair_stats/],
	);
	# requires -mx
		
	$default_params{'-slice'} ||= $raw_params{'-mx'}->origSlice;
	my %params = (%default_params, %raw_params);

 	print "PATH: $path\n";
	
	# gather needed data
	my $data = &output::mxData(%params); # get a basic data object
	my $tt = Template->new($self->ttConfig) || die $tt::ERROR, "\n";  #params{'-config'}); 

	# - Call the various reports, note we want to keep the original output object ($self), yet use methods in Html
	# $r is the name of a template to use in Output/template/html
	foreach my $r (@{$params{'-reports'}}) {
		print " generating: $r\n";

		$params{'-file_name'} = "test_$r.html"; 
		$tt->process($r, $data, "$path/$params{'-file_name'}")|| croak $tt->error(), " died trying to process Html template";
	}
	
	1;
}


__DATA__

