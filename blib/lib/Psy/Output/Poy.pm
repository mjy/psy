package Poy;

use strict;
use warnings;

use Data::Dumper;
use Psy::Dna::Alphabet;
use Psy::Strings::Strings;

use vars qw(@ISA);
@ISA = qw(output);

=head1 NAME

Psy::Output::Poy

=head1 VERSION

Version 0.01

=cut


our $VERSION = '0.01';


=head1 SYNOPSIS

use Psy::Psy;


Returns the the total number of blocks exported when successfull.

=cut

=head1 OPTIONS

Defaults are in (), valid options are in [ ]

	-mode (do_all) [do_all fso_all do_unbrack_fso_brack do_brack_fso_unbrack fso_unbrack do_unbrack do_brack fso_brack]
		legal values are:
		'do_all'	-	DO on all blocks
		'fso_all'	-	FSO on all blocks
	
		'do_unbrack_fso_brack'	- DO on unbracketed, FSO on bracketed
		'do_brack_fso_unbrack' 	- DO on bracketed, FSO on unbracketed
	
		'fso_unbrack'	-	FSO on unbracketed, bracketed excluded
		'do_unbrack'	-	DO on unbracketed blocks, bracketed excluded
	
	 	'do_brack'	-	DO on bracketed, unbracketed excluded
		'fso_brack'	-	FSO on bracketed, unbracketed excluded
	
	-min_terminals (80) [0 to 100]
		see POY manual
		
	-use_prealigned (0) [0 1]
		detect blocks that are prealigable and indicate them as such
		
	-strip_bracketed_gaps [0 1]
	
	-use_deletegapsfrominput (0) [0 1]
	
	-legal_alphabet (IUPAC + '#*? ') [any POY legal character]
		the legal POY alphabet, is set by default and shouldn't likely be used
		
	-u2t (0) [0 1]
		translates U to T when set to 1
		
	-clean (1) [0 1]
		
	-replace ('N') 
	
	-heuristics (contents of <DATA> ) [any legal POY options, as a string]
		overwrites the heuristics in <DATA>, included in the .cmd file

Note on usage of ?, -, and N, we assume:
  ? is missing data (positions are stripped)
  N is ambiguous/missing data whose positions should NOT be stripped (however paritions of only N will be stripped)
  - is true gap (always stripped for POY).

=cut



our $OUT_DIR = "analyses/poy";
our $LEGAL_CHARS = Psy::Dna::Alphabet->new('-type' => 'iupac'); 
$LEGAL_CHARS->add('#*? ');  #  '\#\*\? ';

sub process {
	my $self = shift;
	my %raw_params = @_;
	my %default_params = (
		'-mode' => 'do_all',	
		'-min_terminals' => 80, 		 # NOTE NOT AS whole number here!!
		'-use_prealigned' => 0,
		'-strip_bracketed_gaps' => 1,	 # strip gaps from bracketed blocks
		'-use_deletegapsfrominput' => 0, # place a '-deletegapsfrominput in front of *all* blocks
		'-legal_alphabet' => $LEGAL_CHARS, 
		'-collapse' => 0,
		'-u2t' => 0, 			# translate U to T? 
		'-clean' => 1, 			# translate non legal_alphabet chars to N 
		'-replace' => 'N',  	# should go with the alphabet object ultimately
		'-outgroup_index' => 0,
	);
	# requires -mx
	
	$default_params{'-slice'} = $raw_params{'-mx'}->origSlice if not defined $default_params{'-slice'};

	my %params = (%default_params, %raw_params);
	my $path = $OUT_DIR;

	# $slice->remove("Blocks", $slice->Jrna::slice_stems($mode));
#	$params{'-path'} ||= $OUT_DIR;
	$params{'-file_name'} ||= 'poy.cmd';
 	print "output path: $path\n";
	
	my $fsflag;	# flag to mark fixesstates toggle, simplifying commandfile
	my $paflag;	# flag to mark prealigned toggle, simplifying commandfile
	my $chkdata;
	my $tmpdata;	
	my @excludedblks; 
	my @prealigned;

	my @fs_partitions; # not used in collapse kludge
	my @do_partitions; # not used in collapse kludge

	my $i=0; # track number of datafiles exported

	my $mx = $params{'-mx'} || croak('no matrix passed'); # simplify things

	print "\n\ngenerating POY output\n\n";

	print "|----------------------------\n";
	print "| WARNING! WARNING! WARNING! \n";
	print "| * check your .cmd file, fix it manually, the heuristic search options within are only included as examples \n";
	print "| * delete all old script generated files before generating new ones\n";
	print "| * please check your input and results carefully, automated analyses still require that you\n"; 
	print "| try and understanding what you're doing\n";	
	print "| * check your random seed to see that it is indeed different (or not)\n";
	print "|----------------------------\n\n";
	
	# print "\nlegend: DO- direct optimization; FSO- fixed states optimization.\n\n";
	
	if ($params{'-mode'} eq 'do_all') { print "generating commandfile with DO on all blocks" }
	elsif ($params{'-mode'} eq 'fso_all') { print "generating commandfile with FSO on all blocks" }
	
	elsif ($params{'-mode'} eq 'do_unbrack_fso_brack') { print "generating commandfile with DO on unbracketed, FSO on bracketed" }
	elsif ($params{'-mode'} eq 'do_brack_fso_unbrack') { print "generating commandfile with DO on bracketed, FSO on unbracketed" }
	
	elsif ($params{'-mode'} eq 'fso_unbrack') { print "generating commandfile with FSO on unbracketed, bracketed excluded" }
	elsif ($params{'-mode'} eq 'do_unbrack') { print "generating commandfile with DO on unbracketed blocks, bracketed excluded" }
	
	elsif ($params{'-mode'} eq 'do_brack') { print "generating commandfile with DO on bracketed, unbracketed excluded " }
	elsif ($params{'-mode'} eq 'fso_brack') { print "generating commandfile with FSO on bracketed, unbracketed excluded " }

	print "\n";
	

	# need to setup a results directory, we should be at the right path at this point by calls from Output::Output
 	&Psy::Io::Io::confirmDir("$path/data/"); 
	chdir('../../../'); # this is a hack, and will break with passing custom paths
	&Psy::Io::Io::confirmDir("$path/results/"); 
	chdir('../../../'); # this is a hack, and will break with passing custom paths

	# for filenames
	my $llblk = $mx->origSlice->lengthLastBlk;

	# print "FOOOOO", $Psy::PSY_BASE_DIR;
	#	chdir($Psy::PSY_BASE_DIR);
	#	&Psy::Io::Io::confirmDir($self->basePath);
	# 	&Psy::Io::Io::confirmDir("$path/data/"); 
	
	
	# not used in POY 4?
	# ($params{'-use_deletegapsfrominput'} == 1) && print PCMD " -deletegapsfrominput";	# could code it to leave them in, but POY seems to not be very happy
	
	# gather needed data
	my $data;					   
	$data = &output::mxData(%params); # get a basic data object
	$data->{randseed} = int(rand(999));

	# two output modes now, collapse, or not,
	# collapse collapses consecutive bracketed or unbracketed blocks together, treating them as one

	# a KLUDGE separate output mode
	if ($params{'-collapse'} == 1) {
		# create a new slice
		my $s = $params{'-slice'}->collapse ('-mx' => $mx);
		$i++;
		foreach my $k (keys %{$s}) {
		
			# DO/FSO FLAGS
			if ($s->{$k}->{type} == 0) { # unbracketed
				if ( grep $_ eq $params{'-mode'}, qw(fso_all do_brack_fso_unbrack fso_unbrack) ) { 
					push @fs_partitions, $k;
					#	print PCMD " -fixedstates";
				}
				else {
					push @do_partitions, $k;
					#	print PCMD " -nofixedstates";
				}
			}
			else {
				if ( grep $_ eq $params{'-mode'}, qw(do_unbrack_fso_brack fso_brack) ) { 
					push @fs_partitions, $k;
					# print PCMD " -fixedstates";
				}
				else {
					push @do_partitions, $k;	
				 #	print PCMD " -nofixedstates";
				}
			}
			
			# open a file to write out to
			my $tmpname = "pd";
			$tmpname .= &padleft($k, "0", $mx->origSlice->lengthLastBlk);
			open (PDATA, ">$path/data/$tmpname") || die "couldn't open $path/data/$tmpname to output POY to - \n";

			#		print PCMD " data/$tmpname";
	
			for my $t ($params{'-slice'}->loop("Taxa")) {
				$tmpdata = '';
				foreach my $blk ($s->{$k}->{slice}->loop('Blocks')) {	
					$tmpdata .=  $mx->yx($t,$blk);
				}
				
				$tmpdata =~ s/\?/N/g; 	# - POY barfs on "?"
				$tmpdata =~ s/-//g; 	# all data should have gaps stripped (-deletegapsfrominput)
				

					# know for sure:
					# three options for ? - one leave in/2 translate/3-strip
					# stripping/translating some characters -> should it really be done?
				# $q = quotemeta('?');
				
					#$tmpdata =~ s/-//g;

					# $tmpdata =~ s/U/T/ig;			# translate U to T
					# $tmpdata =~ s/(?)([^ATGUC\?])/N/gi;  	# translate everything that's left to N
				
					print PDATA ">", $mx->ter($t)->label, "\n";
					print PDATA uc($tmpdata); # just "in case" (hehe)
					print PDATA "\n\n";				
			}	
		}	

		# print the script file here
		# return $i; # exit here 
	}

	else { # we're not collapsing (should merge with shared code above)
		
		# based on mode some blocks are excluded beforehand
		if (($params{'-mode'} eq 'do_brack') or ($params{'-mode'} eq 'fso_brack')) {
			$params{'-slice'}->prune('-mode' => 'unbracketed', '-mx' => $mx);
		}
		elsif (($params{'-mode'} eq 'do_unbrack') or ($params{'-mode'} eq 'fso_unbrack')) {
			$params{'-slice'}->prune('-mode' => 'bracketed', '-mx' => $mx);
		}

		print "\n\n";
		$params{'-slice'}->describe;
		print "\n";	
	
		# check for bracketed blocks with no data once dashes OR all N blocks are removed
		for my $blk ($params{'-slice'}->loop("Blocks")) { 
			# next if ((grep($params{'-mode'} eq $_, qw(do_unbrack fso_unbrack))) and ($mx->blk($blk)->bracketed == 1)); # should be removed by default

			my $xcld = 1;
			my $pre = 1;

			for my $t ($params{'-slice'}->loop("Taxa")) {
				$chkdata = $mx->yx($t,$blk);
				$chkdata = uc($chkdata); # make sure everything is Kosher
				
				if ($LEGAL_CHARS->strOk($chkdata)) { # check for illegal characters in this block
					$xcld=2; last;
				}
					
				if ( $mx->blk($blk)->bracketed == 0) { # run the pre-aligned check (all unbracketed blocks with no gaps only) NOTE: This has yet to makes sense and at present isn't used.
					if ($chkdata =~ /-/) { $pre = 0;};
				}			

				#check for illegal length in this block # rewrite with alphabet ultimately
				$chkdata =~ s/-|N|\?//gi;  ## possibly also 'X'	
				if (length $chkdata != 0) { # if true there is at least 1 block with data
					$xcld = 0;	
				}
			}
			
			if ($xcld != 0) {push  @excludedblks, $blk;
				print "excluded block $blk for consiting of -,n,or \? only\n" if $xcld == 1; 
				print "excluded block $blk for POY illegal and untranslatable characters for fasta format: ($chkdata)\n" if $xcld == 2;
			}

		
			if  ($mx->blk($blk)->blkLength  == 1) {	
				if ($mx->blk($blk)->bracketed == 1) {
					if (grep $_ eq $params{'-mode'}, qw(do_all do_brack do_brack_fso_unbrack) ) {
						print "excluded block $blk for illegal length 1 of a bracketed block\n";
						push @excludedblks, $blk;
					}	
				}
				else {
					if (grep $_ eq $params{'-mode'}, qw(do_all do_unbrack_fso_brack do_unbrack) ) {
						print "excluded block $blk for illegal length 1 of a unbracketed block\n";
						push @excludedblks, $blk;
					}	
				} 
				
			} # end taxa loop
			if ( ($pre == 1) && ($mx->blk($blk)->bracketed == 0)) {push (@prealigned, $blk)};
		} # end block loop
		
		# write command and data files in non collapse mode

		for my $blk ($params{'-slice'}->loop("Blocks")) {
			if (not grep($_ == $blk, @excludedblks)) { 
				next if	( (grep ($params{'-mode'} eq $_, qw(fso_unbrack do_unbrack))) and ($mx->blk($blk)->bracketed == 1)) ;
				
				# open a file to write out to
				my $tmpname = "pd";
				$tmpname .= &padleft($blk, "0", $params{'-slice'}->lengthLastBlk);
				open (PDATA, ">$path/data/$tmpname") || die "couldn't open $path/data/$tmpname to output POY too - \n";
			
				# prealigned flag, only used on all blocks if -use_prealigned set to 1
				if ($params{'-use_prealigned'} == 1) {
					if ( grep $_ == $blk, @prealigned) { # huh? need to grok this again
						push @prealigned, $blk;
					} 
				}
				
				# DO/FSO FLAGS
				if ($mx->blk($blk)->bracketed == 0) {
					if ( grep $_ eq $params{'-mode'}, qw(fso_all do_brack_fso_unbrack fso_unbrack) ) { 
						push @fs_partitions, $blk;
					}
					else {
						push @do_partitions, $blk;
						# print PCMD " -nofixedstates";
					}
				}
				else {
					if ( grep $_ eq $params{'-mode'}, qw(do_unbrack_fso_brack fso_brack) ) { 
						push @fs_partitions, $blk;
						# print PCMD " -fixedstates";
					}
					else {
						push @do_partitions, $blk;
						# print PCMD " -nofixedstates";
					}
				}
					
				# datafile 

				for my $t ($params{'-slice'}->loop("Taxa")) {		
					$tmpdata =  $mx->yx($t,$blk);
		
					$tmpdata =~ s/\?//g; 	
				
					# should do this through Alphabet ultimately	
					$tmpdata =~ s/(?)([^ATGUC\-])/N/gi;  # translate ambiguities to Ns (including ? but NOT -)
					
					# at this point you should have ATGUC or -			

					# if the block is prealigned and you're in prealigned mode leave the gap
	    			if ( grep $_ eq $blk, @prealigned ) { 
						# do nothing at present	
					}
					else {
						$tmpdata =~ s/\-//g; 	
					}
				
					# $q = quotemeta('?');
				
					if ($params{'-u2t'} == 1) {
					 $tmpdata =~ s/U/T/ig;			# translate U to T
					}	

					# if there is data left to print, print it
					if (length($tmpdata) > 0) {
						print PDATA ">", $mx->ter($t)->label, "\n";
						print PDATA uc($tmpdata); # just in case (hehe)
						print PDATA "\n\n";		
					}		
			}
				$i++;
				close PDATA;
			}
		}


	} 	# end non-collapse

	# write the command file (need to get this into a .tt

	if ($#fs_partitions > 0) {
		$data->{'fs'} = 1;	
		$data->{'fs_partitions'} = \@fs_partitions;
		$data->{'fs_file_str'} =  join(", ", map "\"pd".&padleft($_, "0", $llblk)."\"", @fs_partitions)
	}	

	if ($#do_partitions > 0) {
		$data->{'do'} = 1;	
		$data->{'do_partitions'} = \@do_partitions;
		$data->{'do_file_str'} =  join(", ", map "\"pd".&padleft($_, "0", $llblk)."\"", @do_partitions) 	
	}

	# this isn't used yet...
	if ($#prealigned > 0) {
		$data->{'pa'} = 1;	
		$data->{'prealigned'} = \@prealigned;
		$data->{'pa_file_str'} =  join(", ", map "\"pd".&padleft($_, "0", $llblk)."\"", @prealigned) 	
		#
	}

	# print some debug output
	print "all fs blocks (", $#fs_partitions +1, "): @fs_partitions \n";
	print "all do blocks (", $#do_partitions +1, "): @do_partitions \n";
	print "all excluded blocks (", $#excludedblks +1, "): @excludedblks \n";
	print "prealigned blocks  (",$#prealigned+1, "): @prealigned\n";
	
	$data->{'filename'} = sub {my $blk = shift; return "pd".&padleft($blk, "0", $llblk)};

	$data->{outgroup} = $params{'-mx'}->{'terminals'}->{$params{'-outgroup_index'}}->label;

	my $tt = Template->new($self->ttConfig) || die $Template::ERROR, "\n";  	
	$tt->process('poy_cmd.tt', $data, "$path/poy.cmd") || die 'died trying to process poy cmd template', $tt->error(); 
	print "$i blocks exported\n";
	
	return $i;
}


# default search heuristics are listed here, are overwritten by '-heuristics'
__DATA__

__END__

