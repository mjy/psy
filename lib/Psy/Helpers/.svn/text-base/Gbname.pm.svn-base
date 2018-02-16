package Psy::Helpers::Gbname;

use warnings;
use strict;

use Data::Dumper;
use BIO::DB::Genbank;
use Scalar::Util qw(blessed);
use Carp;

=head1 NAME

Psy::Helpers::Gbname

=head1 VERSION

# version 0.02 - Matt Yoder
# requires Scalar::Util, hopefully a common package

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Object for generating taxon labels for genbank sequence objects or strings containing a genbank accession number.
Can call genbank directly or generate from passed object.

#my $foo = gbname->new;
#print $foo->name('-in' => 'AY764386_foogurt');

A seq obj looks like this:

	$seq->{species}->{_classification} ... 

=head1 OBJECT

=cut

=head2 new

Make a new gbname object.

=cut

sub new {
	my $type = shift;
	my %params = @_;
  	my $self = {};

	bless $self, $type;

	$self->{gb} =new Bio::DB::GenBank; # (-retrievaltype => 'tempfile' ); # use with lookups in case needed
#	$self->seqObjLookup('U39957');
	$self->_init(); # at present just loads from DATA

    return $self; 
}

sub _init {
	my $self = shift;
	my %params = @_;
	while (<DATA>) {
		chomp;
		next if $_ eq '';
		my @row = split / /, $_;

		# print join "|", @row;
		# print "\n";
		my $key = shift @row;
		$self->{namesets}->{$key} = \@row;
		1
	}
}

=head2 nameSet 

Returns an array given the set name (sets found in __DATA__)

=cut


sub nameSet { 
	my $self = shift;
	my $set = shift;
	return @{$self->{namesets}->{$set}};
}


=head2  nameIntersect

Returns the FIRST (may produce undesired results, be careful) intersection of a namset array and a supplied array (e.g. Genbank classification)
 
=cut

sub nameIntersect { # 
	my $self = shift;
	my %params = @_;
	# -names is an array reference
	my @ns = $self->nameSet($params{'-higher_name_set'}); 

	foreach my $n ( @{$params{'-names'}} ) { 
		if (grep $_ eq $n, @ns) { return $n };
	}
	return ''
}

sub _seqClassification { # returns the _classification array for a seq object ## bioperl redundant? 
	my $seq = shift;

	(not defined $seq) && return ();
	return @{$seq->{species}->{_classification}};
}

sub _seqCommonName { # returns the _classification array for a seq object ## bioperl redundant?
	my $seq = shift;
	return $seq->{species}->{_common_name} || "";
}

sub _seqDisplayId { # returns the _classification array for a seq object ## bioperl redundant?
	my $seq = shift;
	return $seq->{primary_seq}->{display_id};
}


=head2 seqObjLookup

Returns a seq object for a given name, legal formats are ( accession || accession_foo_barBlort)
 
=cut

sub seqObjLookup { 
	my $self=shift;
	my $id = shift;
	$id eq '' and return 0;
	($id =~ /_/) && ( $id = substr($id,0,$+[0]-1)); # grab the accessor
	# print "finding: $id\n";
	return $self->{gb}->get_Seq_by_id($id) || croak "no such genbankID $id\n";	
}


=head2 seqFamilyGroupNames

Returns (subfamily, family) names given a seq object
 
=cut

sub seqFamilyGroupNames { 
	my $self = shift;
	my $seq = shift;	# bioperl seq object OR object from fastgenbank
	my @cl = &_seqClassification($seq); 
	#print @{$seq->{species}->{_classification}};
	my %names;
	my @foo;
	push @foo, grep $_ =~ m/idae$/gi, @cl;
		$names{family} = shift @foo || 'undef';
	push @foo, grep $_ =~ m/inae$/gi, @cl;
		$names{subfamily} = shift @foo || 'undef';
	push @foo, grep $_ =~ m/oidea$/gi, @cl;
		$names{superfamily} = shift @foo || 'undef';
	return %names; 
}


=head2 summaryHigher
 
=cut

sub summaryHigher {
	my $self = shift;
	my @names = sort (keys %{$self->{summary}->{higher}});
	print "\n\n";
	foreach my $h (@names) {
		print "$h ->", $self->{summary}->{higher}->{$h}, "\n"
	}
	return 1;
}

=head2 summaryNoHigher


=cut

sub summaryNoHigher {
	my $self = shift;
	$self->{summary}->{nohigher} and (print join "\n", @{$self->{summary}->{nohigher}})
}


=head2 name

Return a formated name, lots of options here

=cut

sub name { 
	my $self = shift;
	my %raw_params = @_;
	my %default_params = (
		'-accession_label' => 1, 				# on or off
		'-accession_label_position' => 'pre', 	# pre or post (POST NOT BACKWARDS LEGAL)
		'-species_label' => 'gb_common_name', 			# legal: original, gb_common_name, none
		'-higher_label' => 'family', 			# legal: set, family, subfamily, superfamily
		'-higher_name_set' => 'Higher',			# see __DATA__ sets
		'-cap_higher' => 1, 					# capitilize the higher portion
		'-truncate_higher' => 0, 				# 0 - no trucate, other truncate to that number
		'-undefined_label' => '', 				# appended to undefined names
		'-summarize' => 1 						# store a summary of names generated if object is called > 1 times
	);	
	# requires -in, where -in == a blessed seqobject or text  containing a genbank accession in the two legal forms (accession || accession_foo)
		
	my %params = (%default_params, %raw_params);
	my $seqObj;
	
	# need to set three things, accession, common name, higher classification
	
	# allow -in to be either a seq object or a accession name, a little kludge
	if 	(&blessed($params{'-in'})) { # blessed so its a seq object
		$seqObj = $params{'-in'};
	}
	else {
		$seqObj = $self->seqObjLookup($params{'-in'})
 	};
		
	my $name = '';

	# accession	 pre
	if ($params{'-accession_label'} == 1 && ($params{'-accession_label_position'} eq 'pre')) { $name .= &_seqDisplayId($seqObj)};

	# species
	if ($params{'-species_label'} eq 'original') {
		$params{'-in'} =~ /_/;
		(length $name != 0) && ($name .= "_");
		$name .= substr($params{'-in'}, $+[0], length($params{'-in'}) - $+[0]) 
	}
	elsif ($params{'-species_label'} eq 'gb_common_name') {
	 	(length $name != 0) && ($name .= "_");		
		$name .= &_seqCommonName($seqObj)		
	}

	# higher
	my $tmp_higher = $self->higher(%params, '-seq_obj' => $seqObj);

	if (length $tmp_higher != 0) { # there is a higher tax to add
		(length $name != 0) && ($name .= "_");	
		if 	($params{'-truncate_higher'} != 0) {
			$tmp_higher = substr($tmp_higher, 0, $params{'-truncate_higher'})
		}
		if 	($params{'-cap_higher'} == 1) {
			$tmp_higher = uc($tmp_higher)
		}	
		$name .= $tmp_higher;
	}

	# accession post
	if ($params{'-accession_label'} == 1 && ($params{'-accession_label_position'} eq 'post')) { $name .= "_".&_seqDisplayId($seqObj)};
	
	# no label available
	if (length $name == 0) { 
		$name = 'undefined'.$params{'-undefined_label'}
	}
	
	# final processing 
	$name =~ s/\W/_/gi; # strip any whitespace just in case
	
	if ($params{'-summarize'}) {
		push @{$self->{summary}->{namelist}}, $name;
		if (length $tmp_higher > 0) {
			$self->{summary}->{higher}->{$tmp_higher}++
		}
		else {
			push @{$self->{summary}->{nohigher}}, $seqObj->_seqDisplayId; 
		}
	}	
	return $name;
}


=head2 higher

Split from name so it can be accessed alone

=cut


sub higher { 
	my $self = shift;
	my %raw_params = @_;
	my %default_params = (
		'-higher_label' => 'family', 			# legal: set, family, subfamily, superfamily
		'-higher_name_set' => 'Higher',			# see __DATA__ sets
		'-from_name' => 0						# boolean, when 1 then assume its from the higher name
	);
	# require -seq_obj => $seqObj 
	
	my %params = (%default_params, %raw_params);

	my $seqObj = $params{'-seq_obj'};

	# assume the last part of the string (after _) is the higher- yak (yet another kludge)
	if ( $params{'-from_name'} == 1) {
		$seqObj->{'primary_seq'}->{'display_id'} =~ m/(.+)_(.*)\z/;
		if ($2) {
			return $2}
		else {
			return 'undefined'
		}
	}

	my $tmp_higher = '';
	if ($params{'-higher_label'} =~ /fam/) { # matches fam/subfam/superfam options
		my %fams = $self->seqFamilyGroupNames($seqObj);
		$tmp_higher = $fams{$params{'-higher_label'}};
	}
	elsif ($params{'-higher_label'} eq 'set') {
		my @cl = &_seqClassification($seqObj); # should point to this anonymously I guess
		$tmp_higher = $self->nameIntersect(	%params, '-names' => \@cl );
	}

	if ( $tmp_higher eq '') {return  'undefined'};	
	return  $tmp_higher
}

=head1 AUTHOR

'Matt, C<< <m{j}yoder@{tee}[aye](em)(you).domain4unis> >>

=head1 BUGS

Please report any bugs or feature requests to
C<bug-psy-helpers-gbname@rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Psy>.
I will be notified, and then you'll automatically be notified of progress on
your bug as I make changes.

=head1 ACKNOWLEDGEMENTS

=head1 COPYRIGHT & LICENSE

Copyright 2005 'Matt, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1; # End of Psy::Helpers::Gbname


=head2 DATA section

Data contains name-sets, first value is key, remaining are set values, add as you want.  

No lines below the namesets!

=cut

__DATA__

Arthropoda Hexapoda Crustacea Pauropoda Diplopoda Chilopoda Symphyta Arachnida Xiphosura Pycnogonida
Higher Protura Diplura Collembola Thysanura Archaeognatha Ephemeroptera Odonata Plecoptera Embiidina Phasmida Orthoptera Matophasmatodea Zoraptera Isoptera Mantodea Blattaria Dermaptera Grylloblattodea Grylloblatta Psocodea Thysanoptera Hemiptera Megaloptera Raphidioptera Neuroptera Coleoptera Strepsiptera Diptera Mecoptera Siphonaptera Phthiraptera Embioptera Psocoptera Phasmatodea Trichoptera Lepidoptera Hymenoptera Anostraca Notostraca Spinicaudata Laevicaudata Cladocera Nectiopoda Cephalocarida Copepoda Myodocopa Podocopa Branchiura Tantulocarida Thecostraca Mystacocarida Malacostraca Xiphosura Paligradi Araneae Amblypygi Thelyphonida Palpigradi Schizomida Ricinulei Acari Opiliones Scorpiones Pseudoscorpiones Solifugae Symphyta Symphyla Chilopoda Diplopoda Diplostraca Maxillopoda Pycnogonida Scutigerella Pauropoda Pentastomida Nemertea Tardigrada
Myriapoda Pauropoda Diplopoda Symphyta
Chelicerata Arachnida Xiphosura
Endopterygota Megaloptera Rhapidioptera Neuroptera Coleoptera Strepsiptera Diptera Mecoptera Siphonaptera Trichoptera Lepidoptera Hymenoptera
Hexapoda Collembola Protura Diplura Insecta

