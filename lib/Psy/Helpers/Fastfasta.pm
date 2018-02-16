package Psy::Helpers::Fastfasta;

use warnings;
use strict;
use Data::Dumper;

use Psy::Strings::Strings;
use Psy::Helpers::Fastmatrix;
# use Psy::Helpers::Gbname; ## might want to remove this

 our @ISA = ('Psy::Helpers::Fastmatrix');

=head1 NAME

Psy::Helpers::Fastfasta - The great new Psy::Helpers::Fastfasta!

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

A fastmatrix subclass for fasta

my $bar= fastfasta->new('-infile' => 'ichneumonoidea_28S_trimmed.fas');
$bar->out_unBoundMx;

$bar->unfound('-compare_file' => '18Shaves.txt', '-outfile' => '18Shavenots.txt');

my $foo = fastfasta->new('-infile' => 'hymenoptera_18s.all');
$foo->out_oneLineRegexCentered('-outfile' => 'hymenoptera_18s.aln');
 $foo->findOne(''); 

# [C][CT][T][G][ACT][GT][AG][AG][A][CT][CT][CGT][ACG][AC][AG][A][AGT]?[AGT][ACT]?[C][G][AG][AGT][ACGT][AG][T]?[G]?[AG][ACG][AG][AG][AGT][A][ACT]?[ACT]?[AC][AC][AGT]?[CT][AGT][T]?[ACGT][ACT][AGT][AG]
# 
#$foo->out_unBoundMx('-outfile' => "foo.txt");
#$foo->findOne('[AG][ACT][AG][AG][ACGT]?[AGT][ACGT][AG][G][CT][G][C][AC][AC][ACT]');
#$foo->boundOligo('Ichneumonidae_EMA302909',
#	'[C][CT][T][G][ACT][GT][AG][AG][A][CT][CT][CGT][ACG][AC][AG][A][AGT]?[AGT][ACT]?[C][G][AG][AGT][ACGT][AG][T]?[G]?[AG][ACG][AG][AG][AGT][A][ACT]?[ACT]?[AC][AC][AGT]?[CT][AGT][T]?[ACGT][ACT][AGT][AG]',
#		'[G][C][AG][CT][A][C][A][CT][G][T][T][G][G]'

#);

#	 $foo->out_boundMx( # 28S D2-D3
#		'-outfile' => 'ichs_28S_trimmed.fas',
#		'-left_bound' => '[C][CT][T][G][ACT][GT][AG][AG][A][CT][CT][CGT][ACG][AC][AG][A][AGT]?[AGT][ACT]?[C][G][AG][AGT][ACGT][AG][T]?[G]?[AG][ACG][AG][AG][AGT][A][ACT]?[ACT]?[AC][AC][AGT]?[CT][AGT][T]?[ACGT][ACT][AGT][AG]', # finds 681
#		'-right_bound' => '[A]?[AGT][AG][CT][C][CT][A][A][A][G][G][C][AG][ACT][A][A][GT][AG][A][A][A]' # finds 562
#	 );


 # $foo->out_boundMx( # 28S D1-D3
 #	'-outfile' => 'arth_28S_D1-D3_trimmed.fas',
 #	'-left_bound' => '[C][C][CT][CG][C][CT][AG][A][A][CT][T][T][A][A][G][C][A][T][A][T][ACT][AC][ACGT][T][ACT][A][G][CG][AG][G][AC][G][G][A][A][AG][A][G][AG][A][A]',
 #	'-right_bound' => '[AG][ACT][AG][AG][ACGT]?[AGT][ACGT][AG][G][CT][G][C][AC][AC][ACT]'
 # );


#$foo->out_boundMx( # 18S Large (these are good bounds)
#	'-outfile' => 'arthropoda_18S_large_trimmed.fas',
#	'-left_bound' => '[AG][ACGT][A][ACT][AG][G][C][T][C][AC][AGT][T][A][ACGT]?[AT]', #'[G]?[ACG][CT]?[T]?[C]?[A][A]?[AG]?[C]?[C]?[A]?[T][G][C][AG][ACT]?[AG][CT][CGT]?[CGT][ACGT][ACT]?[CG][T][ACGT][ACG][ACGT][ACGT]?[G]?[C]?[G]?[T]?[AGT][CG]?[ACGT]',
#	'-right_bound' => 'GTTTCCGTAGGTGAACCTGCGGAAGGATCATTA' # '[G][A][CT][T][AG][ACG][G][T][C][C][C][T][G][CT][C][CT][T][T][GT][CG][CT]?[T]?[ACG]?[C]?[A]?' # good bound
#);

#$foo->out_boundMx( # 18S Mystery (based on small)
#	'-outfile' => 'arthropoda_18S_small_trimmed_mystery.fas',
#	'-left_bound' => '[ACT]?[AGT]?[GT]?[ACT]?[T]?[AG]?[A]?[AG]?[AG]?[AC]?[G]?[CT]?[CT]?[C]?[G]?[C]?[AGT]?[ACGT]?[GT]?[CT]?[CT]?[AGT]?',
#	'-right_bound' => '[T]?[T]?[AC]?[AGT]?[AGT]?[AG]?[AG]?[A]?[CT]?[CGT]?[ACG]?[CGT]?[ACGT]?[ACGT]?[CGT]?[ACG]?[ACGT]?[G]?[AG]?[C]?[ACGT]?[A]?[ACGT]?[ACGT]?[CGT]?[ACGT]?[AGT]?[ACGT]?[AGT]?[G]?[C]?[T]?[ACG]?[A]?[A]?[ACGT]?[T]?[A]?[A]?[ACG]?[ACGT]?[AGT]?'
#);

#$foo->out_boundMx( # 18S small
#	'-outfile' => 'arthropoda_18S_small_trimmed.fas',
#	'-left_bound' => '[ACT][AGT][GT][ACT][T][AG][A][AG][AG][AC]?[G][CT][CT][C][G][C]?[AGT][ACGT][GT][CT][CT][AGT]',
#	'-right_bound' => '[T][T][AC][AGT][AGT][AG][AG]?[A][CT][CGT][ACG][CGT][ACGT]?[ACGT]?[CGT][ACG][ACGT][G][AG][C]?[ACGT]?[A]?[ACGT]?[ACGT][CGT][ACGT][AGT][ACGT][AGT][G][C]?[T]?[ACG][A][A][ACGT][T]?[A]?[A]?[ACG]?[ACGT][AGT]?'
#);


#$foo->findMany('[G][C][AG][CT][A][C][A][CT][G][T][T][G][G]');


#$foo->query( 
#	'-query' => '(AB036200)' # ("Hymenoptera"[Organism] AND "28S"[All Fields]) || ("Hymenoptera"[Organism] AND "28S"[All Fields])
#);

# $foo->writeFasta('-outfile' => 'arth_18S.fas');

=cut

=head2 new

=cut

sub new {
	my $type = shift;
	my %params = @_;
  	my $self = {};
# 	$self->{Interleaves} = {};	
	bless $self, $type;
	
	$self->_init(%params);

    return $self; 
}

sub _init {
	my $self = shift;
	my %params = @_;
	
	my $current_ter;
	
	# read in, doesn't assume anything about taxa, just appends all seq data to the preceeding '>taxon'
	open ( IN, $params{'-infile'} ) || die print "cannot open $params{'-infile'} \n ";
		my $i = 0;
		while (<IN>) {	
			my $row = $_;
			chomp $row;
			
			$row =~ s/^(\s+)//; # strip space from begining of the row, Justin Case
			
			# >gi|66865690|gb |DQ026302.1| Ca
			
			if ($row =~ /\>/) { 
				if ($row =~ /(>gi\|\d+\|)(\w*)\|(\w*)(\.\d\|)/) { # genbank style
					$current_ter = $3;	
				}
				elsif (	$row =~ /(>)(\w*)/) {
					$current_ter = $2;	
				}

				$i++;
			}
			elsif (length $row != 0) {
				$self->{mx}->{$current_ter} .= $row;
			}
			else {
			}
		}	
		print "read $i terminals\n";
	close IN;
}

1;

__DATA__
d2_5_prime	CTYTGAAKAGAGAGTYMAANAGTRCGTGAAACYKYYYRGRG
d2_3_prime	DCCCGTCTTGAAACACGGACCAAGRAGT

=head1 AUTHOR

'Matt, C<< <m{j}yoder@{tee}[aye](em)(you).domain4unis> >>

=head1 BUGS

Please report any bugs or feature requests to
<bug-psy-helpers-fastfasta@rt.cpan.org>, or through the web interface at
<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Psy>.
I will be notified, and then you'll automatically be notified of progress on
your bug as I make changes.

=head1 ACKNOWLEDGEMENTS

=head1 COPYRIGHT & LICENSE

Copyright 2005 'Matt, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1; # End of Psy::Helpers::Fastfasta
