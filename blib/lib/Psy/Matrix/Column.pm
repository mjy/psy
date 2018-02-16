package column;

use Psy::Strings::Strings;
use Psy::Dna::Iupac;

use vars qw(@ISA @EXPORTER @EXPORT_OK $VERSION);

$VERSION = '0.02';

use strict;
use warnings;
use Data::Dumper;

use Psy::Psy;
@ISA = qw(Psy);




=head1 NAME

Psy::Matrix::Column

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

-This object does not reference terminals.
-alphabet is passed throughout as \@

	Builds a column object.

    use Psy::Matrix:Column;	
	my $col = column->new(@my_array);
	

=cut


sub new {
	my $type = shift;
	my @seqs = @_; # reference to array of seqs, column to return (starts at 1)
	my $self = {};
	
	bless $self, $type;

	$self->_init(@seqs);
	
	return $self; 
}

sub _init {
	my $self = shift;
	my @seqs = @_;
	$self->{total_seqs} = $#seqs + 1; # zero beers
	$self->version($VERSION);

	##  check blocks for whitespace - not legal?!

	foreach my $seq (@seqs) {	
		# change this to allow for column vectors to be added at will, with indicies derived from the matrix
		my @cols = split //, $seq;
		for (my $i=0; $i < length($seq); $i++) {
			$self->{colvec}->{$i} .= $cols[$i];		# the column vector
			$self->{cols}->{$i}->{$cols[$i]}++;		# a count of each unique element (also a method, but useful enough to make "global")
		}
	}
}

sub _verify {
	my $self=shift;
	my %params = @_;
	foreach my $p (keys %params) { # check matrix parameters here
	
	}
	return %params;
}


=head2 returnColHash

Returns the hash count of the column

	$col->returnColHash(1);


=cut


sub returnColHash {	
	my $self = shift;
	my $col = shift;
	die "column not defined" if not defined $col;		
	return $self->{cols}->{$col};
}


=head2 returnColVec

Returns a string representing the column.

	$col->returnColVec(1);
	

=cut


sub returnColVec {	
	my $self = shift;
	my $col = shift;
	die "column not defined" if not defined $col;		
	return $self->{colvec}->{$col};
}

=head2 totalSeqs

returns the total number of original sequences passed

=cut


sub totalSeqs { 
	## make this true accessor
	my $self = shift;
	return $self->{total_seqs};
}

sub _column_alphabet { # returns reference to array of all unique characters in -column
	my $self=shift;
	my %params = @_;
	return "no col" if not defined $params{'-column'};
	my @alpha = keys %{$self->{cols}->{$params{'-column'}}};
	@alpha = (sort {$a cmp  $b} ( @alpha) ); # could probably merge with return, but  \(sort...) doesn't seem to work?
	return \@alpha;
}

sub _alphabet { # returns the whole alphabet for a columns object, probably not the most efficient way of doing this
	my $self = shift;
	my %keys;
	my @alpha;
	foreach my $v (values %{$self->{cols}} ) {
		%keys =  (%keys, %{$v});
	}
	@alpha = (sort {$a cmp $b} ( keys %keys) ); # could probably merge with return, but  \(sort...) doesn't seem to work?
	return \@alpha;
}
	
sub _maxInCol {
	my $self = shift;
	my %params = @_;
	#$params{'-alphabet'} =  $self->_alphabet() if not defined $params{'-alphabet'};
	die "no column" if not defined $params{'-column'};		
	my @empty = \();
	$params{'-exclude'} = \@empty if not defined $params{'-exclude'};
	
	my $excluded = $self->totalSeqs() - $self->_row_count(%params) ;
	# print "[$excluded]\t";
	my $max = 0;
	my $char = ("X");
	
	foreach my $l ( @{$self->_column_alphabet(%params)} ) {
		if ( not  (grep $_ eq $l, @{$params{'-exclude'}}) ) {
			if ($self->{cols}->{$params{'-column'}}->{$l} > $max) {
				$char = $l;
				$max = $self->{cols}->{$params{'-column'}}->{$l};
			}
			elsif  ($self->{cols}->{$params{'-column'}}->{$l} == $max) {
				$char .= $l;
			}
		}
	}

	#print "{$max}";
	
	if ( ($self->totalSeqs() - $excluded) == 0) {
		return (
			'char' => 'X',
			'val' => 0,
			'nrows' => 0 
		);
	}
	else {
		return (
			'char'=> $char,
			'val'=> $max / ($self->totalSeqs() - $excluded), # - $excluded
			'nrows' => $self->totalSeqs() - $excluded
		); 
	}
}

sub _row_count { # returns count of '-alphabet' minus '-excluded' rows
	my $self = shift;
	my %params = @_;

	return "undef" if not defined $params{'-column'};

	if (not defined $params{'-alphabet'}) { print "no alphabet to _row_count"; return -1 }
	my @empty = ();
	$params{'-exclude'} = \@empty if not defined $params{'-exclude'};

	my $count = 0;
	foreach my $char ( @{$self->_column_alphabet(%params)}) {
		 if ( not (grep $_ eq $char, @{$params{'-exclude'}}) ) {
			$count += $self->{cols}->{$params{'-column'}}->{$char};
		}; 
	}
	return $count;
}
	
sub _col_by_percent {
	my $self = shift;
	my %params = @_;

	return {'no column' => 0 } 	if not defined $params{'-column'};
	$params{'-exclude'} = () if not defined $params{'-exclude'};

	my $excluded = $self->totalSeqs - $self->_row_count(%params);
	my %result;

	foreach my $l (keys %{$self->{cols}->{$params{'-column'}}} ) {
		if ($self->totalSeqs() - $excluded == 0) {
			$result{$l} = -1;										# figure this out!!!!!
		}
		else {
			$result{$l} = $self->{cols}->{$params{'-column'}}->{$l} / ($self->totalSeqs() - $excluded);
		}
	}		
	return {%result};
}


=head2 pluralConcensus

Returns the non-iupac plural concensus for a given column

	$col->pluralConcensus(
		'-column' => 1
	);


=cut


sub pluralConcensus { 
	my $self = shift;
	my %params = @_;

	die "no column passed to pluralConcensus" if not defined $params{'-column'};
	$params{'-alphabet'} = $self->_alphabet() if not defined $params{'-alphabet'};
	$params{'-concensuscutoff'} = .1 if not defined $params{'-concensuscutoff'};
	$params{'-concensusminrows'} = 0 if not defined $params{'-concensusminrows'};
		
	my $rowcount = $self->_row_count(%params);

	if ( $rowcount / $self->totalSeqs() < $params{'-concensusminrows'}) { # check versus percent of rows of data present overall
		my $s = int( $rowcount - ($params{'-concensusminrows'} * $self->totalSeqs() )); # return the number of additional rows that would allow a pass of the test
		$s = "{".$s."}";	   
		return $s;
	}	

	die "-concensuscutoff out of range: ".$params{'-concensuscutoff'} if (( $params{'-concensuscutoff'} > 1) or (  $params{'-concensuscutoff'} <  0));	
	my $chars = $self->_col_by_percent(%params);
	my $pcon;
	
	foreach my $key (keys %{$chars}) {
		if ( not  (grep $_ eq $key, @{$params{'-exclude'}}) ) {
			$pcon .= $key if (${$chars}{$key} >= $params{'-concensuscutoff'});
		}
	}
	
	return &alphabetize($pcon) if $pcon;
	return 'X';	
}

sub _columns {
	my $self = shift;
	my @array = keys %{$self->{cols}};
	return (sort {$a <=> $b} @array)
}

sub _score { 
	# a column score, reflects 1) the number of rows with data present; 2) the level of ambiguity in the column
	# this is off the top of my head, I have no clue how relevent it is, but it seems to make *some* sense
	# could be futher refined to reflect scewed ratios by scoring some measure of equivalence, i.e. if AT present, is it 90%A/10%T or 50%A/50%T
	my $self = shift;
	my %params = @_;
	my $score = 0;
	my $rc = $self->_row_count(%params);
	return $score if $rc == 0;					# should always pair with 'X' in 'combined' tables
	
	my $con = $self->pluralConcensus(%params, '-concensustype' => 'residue');
	if ($con eq 'X') { # failed the _pluarl_concensus for cutoff reasons
		$score = 0;
	}
	elsif ($con =~ /\{/) { # another check for failure as above (this one is the one caught I think always
		$score = 0;
	}
	else {
		$score = ( $self->workingAlphabetLength(%params) + 1 - length($con) ) / $self->workingAlphabetLength(%params); 
	}
	$score = $score * $rc / $self->totalSeqs;
	return $score;
}



=head2 workingAlphabetLength

returns the total number of unique characters in the input

=cut


sub workingAlphabetLength {
	my $self = shift;
	my %params = @_;
	return @{$params{'-alphabet'}} - @{$params{'-exclude'}} + 1;
}


=head2 primerRegex

Returns the full primer regex (PERL) for a given columns object.
	-reverse_comp = 1 will return the reverse compliment

=cut


sub primerRegex { 
	my $self = shift;
	my %params = @_;
	$params{'-reverse_comp'} ||= 0;
	my $regex = '';
	foreach my $c ($self->_columns) {
		my $c = $self->regexChar('-column' => $c);
		$c && ($regex .=  $c); 
	 }
	if ($params{'-reverse_comp'} == 1) { return &reverseComp($regex) }	
	else { return $regex;} 
}


=head2 vimRegex

returns a Vim legal regex 

=cut


sub vimRegex { ## not complete
	my $self = shift;
	my %params = @_;
	my $regex = '';
	 foreach my $c ($self->_columns) {
		my $c = $self->regexChar('-column' => $c);
		$c && ($regex .=  $c); #  $self->regexChar('-column' => $c)
	 }
	return $regex;
}
	
sub return_col_vec {	# returns the vector
	my $self = shift;
	my $col = shift;
	die "column not defined" if not defined $col;		
	return $self->{colvec}->{$col};
}



=head2 regexChar

Returns a regex formated position for a given '-column', A MAJOR KLUDGE WITH LOTS OF REDUNDANCY
!!! CUTOFFS NOT ENFORCED !!!

=cut


sub regexChar { 

	my $self = shift;
	my %params = @_;
	
	$params{'-type'} ||= 'vim';
	$params{'-u2t'} ||= 0;
	
	my $out = '';
	
	my $in = $self->return_col_vec($params{'-column'});

	$in =~ s/n//ig; # don't need to examine these
	$in =~ s/\?//ig; # or these
	
	# $in =~ s/-//ig;

	($params{'-u2t'} == 1) && ($in =~ s/[uU]/T/ig); # REQUIRED, don't remove on tweak
	
	# print " {";	
	# print @{$self->_column_alphabet(%params)};
	# print "}";
	
	my $iu = Psy::Dna::Iupac->new;
	
	my $c='';
	return if length $in == 0;
	# not the fastest algorithm, but perhaps better than looping through every character, particularly in long sequences
	
	# replace Iupac codes with coresponding chars;
	foreach my $ichar ($iu->loopIupac) {
		my $tr = $iu->nuc($ichar);
		$in =~ s/$ichar/$tr/gi;
	}

	foreach my $c (split //, 'ACGT') {
		if ($in =~ /$c/i) { $out .= $c};  
	}	

	#	print "OUT: $out";
	#	$out =  $iu->iupac($out);
	
	($out eq '') && (return $out);
	
	length $out == 1 && ( return $out ); # doesn't need []
	
	# make it ? (gapped)
	if ($in =~ /-/) {
		if ($params{'-type'} eq 'vim') {
			$out = "[$out]\\\?" 
		}
		else {
			$out = "[$out]\?" # perl style
		}
	}
	else {
		$out = "[$out]"
	}
	
	return $out;
}



=head2 table

Generate a text table with various reports.

=cut



sub table {
	my $self = shift;
	my %params = @_;

	# paramter checks
	my @empty = ();

	$params{'-type'} = 'percent' if not defined $params{'-type'};
	if ( not grep ($params{'-type'} eq $_, ('tstv', 'count', 'percent', 'max', 'concensus', 'combined', 'score', 'figure', 'primer'))) {  	# check for legal table type
		print "WARNING: type '", $params{'-type'}, "' not recognized, reverting to 'count'\n";
		$params{'-type'} = 'count';
	}
	
	$params{'-offset'} ||= 0; #  if not defined $params{'-offset'};
	$params{'-alphabet'} ||=  $self->_alphabet(); # if not defined $params{'-alphabet'};
	$params{'-header'} ||= 'false'; # if not $params{'-header'} eq ('true');
	$params{'-blockgaps'} ||= 'false'; # if not defined	$params{'-blockgaps'};
	$params{'-exclude'} = \@empty if not defined $params{'-exclude'};
	$params{'-concensuscutoff'} = .1 if not defined $params{'-concensuscutoff'};
	$params{'-concensustype'} = 'iupac' if not defined $params{'-concensustype'};
	$params{'-concensusminrows'} = 0 if not defined $params{'-concensusminrows'};

	if (($params{'-concensusminrows'} < 0) or (	$params{'-concensusminrows'} > 1 )) { die " -concensusminrows out of range: ",$params{'-concensusminrows'} }
	
	#if ($params{'-concensustype'} eq 'iupac') {
	   my $myiupac = iupac->new();
	 #}	   
	
	# print column headers if requested
	if ($params{'-header'} eq 'true') {
		print "[alphabet: ", @{$params{'-alphabet'}}, "; excluded: ", @{$params{'-exclude'}}, "; table type: ", $params{'-type'};

		print (" (% cutoff: ", $params{'-concensuscutoff'}, "; % min rows: ", $params{'-concensusminrows'}, "; type: ", $params{'-concensustype'}, ")") if $params{'-type'} eq 'concensus';	
		
		print (" (% cutoff: ", $params{'-concensuscutoff'}, "; % min rows: ", $params{'-concensusminrows'}, ")") if (($params{'-type'} eq 'tstv') or ($params{'-type'} eq 'combined') or ($params{'-type'} eq 'figure') ) ;	

		print "]\n";
		
		# column headers
		unless ($params{'-type'} eq 'primer') {	print "col#"}; # most tables have this

		# max
		if ($params{'-type'} eq 'max') {print "\tchr(s)\tmax %\tnrows";}
		
		# concensus 
		elsif ($params{'-type'} eq 'concensus') {print "\tconc.";}
		
		# combined
		elsif ($params{'-type'} eq 'combined') {print "\tALL     Ncon    NconRC  Ncon%   NconRC% Icon    IconRC  Icon%  IconRC%  tstv   tstvRC   tstv%   tstvRC% score" ;} # fixed column width
		
		# figure
		elsif ($params{'-type'} eq 'figure') {print "\tcon\ttstv\tgapped"; }

		# score
		elsif ($params{'-type'} eq 'score') {print "\tcon\tscore"; }

		# primer
		elsif ($params{'-type'} eq 'primer') {print "\nWARING: Primer does not use the cutoffs listed above!!!\n"; }

		# percent / count
		else {
			foreach my $char ( @{$params{'-alphabet'}} ) {
				if ( not  (grep $_ eq $char, @{$params{'-exclude'}}) ) {
					print "\t$char";
				}
			}
		}
		print "\n";
	 } 
		
	foreach my $col ($self->_columns()) {
		unless ($params{'-type'} eq 'primer') {print $col+$params{'-offset'}};
	
		# max	
		if  ($params{'-type'} eq 'max') {
			my %max = $self->_max_in_col(%params, '-column' => $col);
			print "\t", $max{'char'};
			printf("\t%.5f", ( 100*$max{'val'} ));
			print "\t", $max{'nrows'};
		}
	
		# concensus
		elsif ( $params{'-type'} eq 'concensus') {
			print "\t"; 
			if ($params{'-concensustype'} eq 'iupac') {
				print $myiupac->iupac( $self->pluralConcensus(%params, '-column' => $col));
			} 
			elsif ($params{'-concensustype'} eq 'residue') {
				print $self->pluralConcensus(%params, '-column' => $col);
			}
			else { print "blOooRff" }
		}
			
		# combined
		elsif ($params{'-type'} eq 'combined') {
				print "\t";
				# col, all, all %, all % minrow (rows); iupac all, iupac %, iupac % minrow (rows), tstv all, tstv %, tstv minrow (rows)

				print &Jrna::str_padright ( join ("", @{ $self->_column_alphabet('-column' => $col) } ), " ", 8);
				
				print &Jrna::str_padright ($self->pluralConcensus(%params, '-column' => $col, '-concensusminrows' => 0, '-concensuscutoff' => 0), " ", 8);
				print &Jrna::str_padright ($self->pluralConcensus(%params, '-column' => $col, '-concensuscutoff' => 0), " ", 8);
				print &Jrna::str_padright ($self->pluralConcensus(%params, '-column' => $col, '-concensusminrows' => 0), " ", 8);
				print &Jrna::str_padright ($self->pluralConcensus(%params, '-column' => $col), " ", 8);
				
				print &Jrna::str_padright ( $myiupac->iupac($self->pluralConcensus(%params, '-column' => $col, '-concensusminrows' => 0, '-concensuscutoff' => 0)), " ", 8);
				print &Jrna::str_padright ( $myiupac->iupac($self->pluralConcensus(%params, '-column' => $col, '-concensuscutoff' => 0)), " ", 8);
				print &Jrna::str_padright ( $myiupac->iupac($self->pluralConcensus(%params, '-column' => $col, '-concensusminrows' => 0)), " ", 8);
				print &Jrna::str_padright ( $myiupac->iupac($self->pluralConcensus(%params, '-column' => $col)), " ", 8);

				print &Jrna::str_padright ( $myiupac->tstv($self->pluralConcensus(%params, '-column' => $col, '-concensusminrows' => 0, '-concensuscutoff' => 0)), " ", 8);
				print &Jrna::str_padright ( $myiupac->tstv($self->pluralConcensus(%params, '-column' => $col, '-concensuscutoff' => 0)), " ", 8);
				print &Jrna::str_padright ( $myiupac->tstv($self->pluralConcensus(%params, '-column' => $col, '-concensusminrows' => 0)), " ", 8);
				print &Jrna::str_padright ( $myiupac->tstv($self->pluralConcensus(%params, '-column' => $col)), " ", 8);

				printf ("%.3f", 100 * $self->_score(%params, '-column' => $col)); # score
		}	

		# tstv
		elsif ($params{'-type'} eq 'tstv') {
			print "\t";
			print $myiupac->tstv($self->pluralConcensus(%params, '-column' => $col));
		}

		# figure (alla Joe and Larry)
		elsif ($params{'-type'} eq 'figure') {
			print "\t";
			print $myiupac->iupac($self->pluralConcensus(%params, '-column' => $col));

			print "\t";
			print $myiupac->tstv($self->pluralConcensus(%params, '-column' => $col, '-concensusminrows' => 0, '-concensuscutoff' => 0 ));

			print "\t";
			print $self->_has('-col' => $col, '-flag' => "gapped", '-has' => '-') || "";
		}
		
		# score
		elsif ($params{'-type'} eq 'score') {
			print "\t";
			print $myiupac->iupac($self->pluralConcensus(%params, '-column' => $col));
			print "\t";
			print $self->_score(%params, '-column' => $col);
		}

		# primer
		elsif ($params{'-type'} eq 'primer') {
			print $self->regexChar('-column' => $col);
		}

		# count/percent
		else {	
			foreach my $char (@{$params{'-alphabet'}} ) {
				if ( not  (grep $_ eq $char, @{$params{'-exclude'}}) ) {	
					if ($params{'-type'} eq 'count') {
						print  "\t", ($self->{cols}->{$col}->{$char} || 0);
					}
					elsif ($params{'-type'} eq 'percent') {
						printf("\t%.3f", ( 100 * ( ${$self->_col_by_percent(%params, '-column' => $col)}{$char} || 0 )  )); # "$char";
					}
				}
			}
		}	

		unless ($params{'-type'} eq 'primer') {	print "\n"};
	}	
}

1;

__END__




