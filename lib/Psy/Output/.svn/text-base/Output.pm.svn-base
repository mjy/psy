package output;

use strict;
use warnings;
use Data::Dumper;
use Carp;
use Cwd;

use Template;

use vars qw(@ISA @EXPORTER @EXPORT_OK $VERSION);
use Psy::Psy qw($PSY_BASE_DIR);

our $AUTOLOAD;

use Psy::Io::Io qw(confirmdir);
use Psy::Output::Wrapper;
use Psy::Dna::Alphabet;
use Psy::Strings::Strings;
use Psy::Analysis::Seqpair;

=head1 NAME

Psy::Output::Phase

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Base class for Psy output


=cut


## output the default loaded matrix if none passed (export ok these functions?)
## need a global list of legal output object types to compare autoload with!


sub new {
	my $type = shift;
	my %raw_params = @_;
	my $self = {};

	bless $self, $type;
	
	$self->_init(%raw_params);

    return $self;     
}


=head2 _init

Called internally only.  Does several things, primarily setting paths up and providing default parameters for Template::Toolkit output.

=cut


sub _init {
	#print "\n\nconfig:";
	#print Dumper (%params);
	#print "FOO: ", $params{'-config'}->{'OUTPUT_PATH'};
	# print "ROOT_DIR: $Psy::ROOT_DIR\n";
	
	# $chdir($Psy::ROOT_DIR);
	
	my $self= shift;
	my %raw_params = @_;
	my %default_params = (
							'-project_name' => 'foo',						
						);	

	my %params = (%default_params, %raw_params); # these are stored by the default params method with the object and later accesses by the AUTLOADED methods
	
	my $basepath;	
	if  (defined $params{'-path'}) {$basepath = $params{'-path'} }
	else { $basepath = 'out/'.$params{'-project_name'}} #  $Psy::PSY_BASE_DIR.'/

	$self->basePath($basepath);

	print "output to: ", $self->basePath($basepath), "\n";
	
	# $self->defParams(\%params); # pass by reference ##  (change this to object methods)

	$self->projName($params{'-project_name'});  ## need this, but try and get more directly from Psy object ultimately
	
	# file names are set in their resepective module by default
	# the output path and file are set as concat strings, *not* via the config hash sent to the Template object
	# the INCLUDE_PATH for the config is used for non-matrix files (e.g. control files) file use
	
	my $cpath = $self->basePath;
 	my $cwd = getcwd(); # current working directory
	
	# defaults for the Template::Toolkit configuration, note the searchpath for templates ugliness
	
	my $template_libs;
	
	map {$template_libs .= "$_/Psy/Output/templates/html;$_/Psy/Output/templates;"}  @INC; # point to the templates, should work across platforms

	print "basepath: ", $self->basePath, "\n";

	# print "O ", $cwd, "\n";
	# print "template libs: $template_libs\n";
	
	my %config = ( # a default configuration hash used in Template::Toolkit
					'EVAL_PERL' => 1,
				   	'TRIM' => 1,
					'DELIMITER' => ';', # because we want to use with DOS too (which uses ':')
					'INCLUDE_PATH' => "$template_libs;$cwd/lib/Psy/Output/templates;$cwd/lib/Psy/Output/templates/html;", 
					'POST_CHOMP' => 0
				); 
				
	$self->ttConfig(\%config);

	# print Dumper(%config), "\n\n";

	# check for base path creation
	chdir($Psy::PSY_BASE_DIR); # must exist
	print &Psy::Io::Io::confirmDir($self->basePath); # confirm that the root output directory is set and changedir there
	return 1;
}


=head2 AUTOLOAD


Used to allow for different calls on an output object, each of which refers to a seperate module in Output::

Example

	my $o = output->new;
	$o->Poy( ... );
	$o->Nexus( ... );

Since no 'Poy' or 'Nexus' methods exists AUTOLOAD is called, the module is loaded, and the function 'process' in that module is called.

=cut


sub AUTOLOAD {  # creates a specific output type which inherits from output
	my $self = shift;
	my $type = ref($self) or croak "$self is really not an object";

	my %params = @_;
	$params{'-psy_ver'} = $VERSION;
	
	# my %raw_params = @_;	
	# my %default_params = $self->defParams;
	# my %params = (%default_params, %raw_params ); #	%raw_params

	#	$print "raw_params: ";
	#	$print (%raw_params);	
	#	print "\nAutoload defParams: ", $self->defParams;
	
	print "\nAUTOLOAD: $AUTOLOAD\n";
	
	my $name = $AUTOLOAD;
	$name =~ s/.*://;   # strip fully-qualified portion

	return if $AUTOLOAD =~ /::DESTROY$/;
	
	print "name: $name \n";

	require "Psy/Output/$name.pm" || croak "output module not found";

	## check permissible functions here
	
	# call Psy::Output::<foo>::process
	my $result = $name->can('process')->($self, %params); ## this dies if the package name is not correct - see if that can be captured! 

	# return to the home path in case AUTOLOAD is called again
	chdir($Psy::PSY_BASE_DIR); # must exist
	chdir($self->basePath);
	return $result;	
}

sub DESTROY { ## make redundant in AUTOLOAD
	my $self = shift;
	# warn "destroying a ", Dumper($self), "\n";
};
		   

=head2 ttConfig

Accessor for the default Template::Toolkit configuration hash

=cut


sub ttConfig { # TemplateToolkit configuration
	my $self = shift;
	if (@_) {
		$self->{tt_config} = shift;
	}
	return %{$self->{tt_config}};	
}


sub projName { ##  this should inherit from Psy??
	my $self = shift;
	if (@_) {
		$self->{project_name} = shift;
	}
	return $self->{project_name};	
}


sub defParams {
	my $self = shift;
	if (@_) {
		$self->{default_params} = shift;
	}
	return %{$self->{default_params}};
}


=head2 basePath

Accessor for the base path

=cut


sub basePath {
	my $self = shift;
	if (@_) {
		$self->{base_path} = shift;
	}
	return $self->{base_path};
}

sub _verify {
	my %params = @_;
	print "OUTPUT _verify: ";
	foreach my $p (keys %params) { # check matrix parameters here
		print "$p\t";
	}		

	# verify directory to output too, defaults to analyses/output_type/project_name/
	# get a filehandle 
	# get a character library?
	
	print "\n";
	return %params;
}


=head2 matrix

Base method for rendering a matrix.  Not used at present.

=cut


sub matrix { 
	# prints a simple matrix using -slice and -wrapper
	my %params = @_;

	die "no matrix passed" if not defined $params{'-mx'};	
	croak '-mx not a matrix in matrix' if not ref($params{'-mx'}) eq 'matrix';

	$params{'-slice'} = $params{'-mx'}->origSlice() if not defined $params{'-slice'};

	$params{'-wrapper'} = wrapper->new() if not defined $params{'-wrapper'};
	
	my $wr = $params{'-wrapper'}; # short form for code use
	my $mx = $params{'-mx'};
	
	$params{'-block_gap_char'} = " " if not defined $params{'-block_gap_char'};
	#$params{'-unbracketed_block'};
	#$params{'-taxon_gap'};

	foreach my $ter ($params{'-slice'}->loop('Taxa')) {
		print $wr->row('pre');
		
		print $wr->taxa('pre');
			print $mx->{terminals}->{$ter}->label;
		print $wr->taxa('post');

		foreach my $blk ($params{'-slice'}->loop('Blocks')) {
			print $wr->block('pre');
				print $mx->yx($ter, $blk); 
				print $params{-block_gap_char};			
			print $wr->block('post');
		}

		print $wr->row('post');
		print "\n";
	}
}


=head2 mxData

Collect a large number of references to subs etc. in one place, the returned hash is used in many Output modules, tweak with caution.

=cut


sub mxData { 
	# returns a hash primarily for  use in Template objects- needs to be overridden in some cases
	my %params = @_;
	die "no matrix passed" if not defined $params{'-mx'};	
 
	my $data;
	
	$params{'-slice'} = $params{'-mx'}->origSlice() if not defined $params{'-slice'};
	$params{'-wrapper'} ||= wrapper->new(); # not used
	$params{'-u2t'} ||= 1; ## MIGHT NOT BE CORRECT
	
	my $wr = $params{'-wrapper'}; # short form for code use
	my $mx = $params{'-mx'};
	my $llt = $mx->lengthLongestTer(%params);

	# looping (array)
	$data->{terminals} = [$params{'-slice'}->loop('Taxa') ];
	$data->{blocks} = [$params{'-slice'}->loop('Blocks') ] ;
	$data->{fiveprimeblocks} = sub { $mx->loopBlocks('-mode' => 'fiveprime', %params) } ; ## not working
	$data->{helices} = sub {$mx->origStructure->loopHelices('-slice' => $mx->origSlice)};
	$data->{bracketed} = sub {$mx->loopBlocks('-mode' => 'bracketed', %params)};
	$data->{unbracketed} = sub {$mx->loopBlocks('-mode' => 'unbracketed', %params)};
	$data->{pairs_in_block} = sub {my $blk = shift;  $mx->origStructure->loopPairs('-blk' => $blk)}; ## origStructure 
	
	# seq data (string)
	$data->{cell} = sub {my ($y, $x) = @_; $mx->yx($y,$x) };
	$data->{row} = sub {my ($t) = shift; $mx->rawRowData('-ter' => $t, %params)}; # includes bracketed and unbracketed by default
	$data->{blk_mask} = sub {my ($blk) = @_; $mx->origStructure->mask('-blk' => $blk) };
	
	# composition/columns (string)
	$data->{block_regex} = sub {my $blk = shift; my $c = column->new($mx->blockData('-blk' => $blk)); $c->primerRegex};
	$data->{block_revcomp_regex} =  sub {my $blk = shift; my $c = column->new($mx->blockData('-blk' => $blk)); $c->primerRegex('-reverse_comp' => 1)};

	# below is an example for other types of mapping, see Pairmeta.pm
	$data->{mi} = sub { # slow, only works for blks who are 5' helices!
		my ($p, $blk) = @_;
		my @cols = $mx->origStructure->pairPos('-pos' => $p, '-blk' => $blk, '-mx' => $mx);
		my $l = $mx->column('-col' => $cols[0]);
	   	my $r = $mx->column('-col' => $cols[1]);
	   	my $sp = Psy::Analysis::Seqpair->new($l, $r); 
		return $sp->mi;
	};
	
	# labels (string)
	$data->{rowlabel} = sub {my ($t) = shift; my $ter = $mx->{terminals}->{$t}->label};
	$data->{justifiedrowlabel} = sub {my ($t) = shift; my $tr = $mx->ter($t); &padright($tr->label, " ", $llt+1)};
	$data->{data_head_aligner} = ' ' x ($llt) ;
	$data->{mask_head_aligner} = "[mask".(" " x ($llt - 4)); # doesn't have closing ]
	$data->{column_numbering} = sub {$mx->numberer(%params, '-numbering_mode' => 'col')};
	$data->{block_numbering} =  sub {$mx->numberer(%params, '-numbering_mode' => 'blk')};
	
# -- temporary, hopefully replaced with Stockholm syntax in future, meh also hacked and redundant (pass label_name

	$data->{blocklabel1} = sub {my $blk = shift; my $str = ""; $str = $mx->blk($blk)->label('-label_name' => 'label1'); $str.$params{'-blk_spacer'} if $params{'-blk_spacer'} ;};
	$data->{blocklabel2} = sub {my $blk = shift; my $str = ""; $str = $mx->blk($blk)->label('-label_name' => 'label2'); $str.$params{'-blk_spacer'} if $params{'-blk_spacer'} ;};
	
	# sums (string)
	$data->{total_chars} = $params{'-mx'}->totalChars(%params);
	$data->{total_ters} = $params{'-slice'}->total('Taxa');
   
	# only if there is a structure object passed (need to modify for multiple structures)
	
	$data->{nucpairs} = join ", ", @{ $mx->origStructure->nucPairs('-mx' => $mx, %params)->{'pairs'} } ;
	
	$data->{paired_positions} = join " ", @{$mx->origStructure->nucPairs('-mx' => $mx, %params)->{'nucs'}} ;
	
	$data->{unpaired_positions} = join " ", @{$mx->origStructure->nucPairs('-mx' => $mx, %params)->{'nucs_ns'}};
	
	$data->{rowmask} = $mx->{structures}->{'original'}->sliceMask(%params);  # apparently working 
	
	$data->{wrapper} = $wr; # does nothing right now AFAIKT 

  # misc
  $data->{randseed} = int(rand(10000));

	## set below to pointer to sub!
	foreach my $ter ($params{'-slice'}->loop('Taxa')) {	
		my @raw = $params{'-mx'}->rowAsArray(%params, '-ter' => $ter);
		if ($params{'-u2t'} == 1) {  map $_ =~ s/[uU]/T/gi , @raw;}   ## should do this thing with alphabet objects, sub 2rna 	
		$data->{multi_row_seq}{$ter} = join "\n", @raw;
	}	

	$data->{output_ver} = $VERSION;
	
	return $data;
}









1;
 
__END__


!! OLD REFERENCE CODE IGNORE !!

sub matrix {

	my %params = @_;

	# required:
	#local (*FH) = shift; 	# filehandle	
	#my (
		#$outslice,			# reference to a slice object
		#$outmode,			# output mode: 0 - bare matrix; 1- nexus; 2- phylip (phymil); 3- css; 4- xml [not currently implemented]; 5- tnt; 6 - Mr. Bayes; 
		#$outintlv,			# interleave form:  0- one big block; 1- original format; 2- blocksize
		#$outintlvsize,		# if outintlv == 2 then =interleaved block size; if outintlv == 3 then ==interleave block to output
		#$outbrackets,		# output bracketed blocks (1) or not (0)
		#$outgapped,			# gap (1) between blocks, or not (0)  ---- may be some problems with NOT plus intlv_numseq (YES THERE ARE!!!)
		#$outgapchar,		# gap character (a string, can be length > 1, e.g. "  ") 
		#$outunbracket,		# unbracket (1) or not (0) - unbracket the characters, include them in count etc **** figure out depenedancies here ****
		#$outdescblks,		# include (1) the [Block from the input file
		#$outnumcols,		# number columns: 0-don't; 1- on top; 2- on bottom; 3- on both; 
		#$outcolstart,		# first column of matrix data position (carefull, must be longer than the longest taxon name - perhaps compute this in mtrx stats
		#$outnumblocks,		# number blocks: 0-don't; 1- on top; 2- on bottom; 3- on both
		#$outdeschead,		# include description block "headers" for interleave blocks 0-no; 1-yes
		#$stemmode,			# for modes see slice_stems
		#$outbracketedAsN	# output bracketed characters as "N"
		#) = @_; # 15

		# RENDER A MATRIX
		# - no headers - just wrappers?  Tricky for multiple interleaves!!!
		
	# check parameters (? perhaps globally )
	die "no matrix passed" if not defined $params{'-mx'};	
	croak '-mx not a matrix in matrix' if not ref($params{'-mx'}) eq 'matrix';

	$params{'-slice'} = $self->orig_slice();

	$params{'-interleave'} = $self->orig_interleave() if not defined $params{'-interleave'}; 
	croak '-interleave not an interleave in matrix' if not ref($params{'-interleave'}) eq 'interleave';

	$params{'-fh'} = *STDOUT if not defined $params{'-fh'}; # use the local version of this

	# local $params{'-fh'}; # use this
	
	# possible parameters
	# -fh
	# -slice
	# -interleave
	# -bracketedchar
	# -outputbracketed
	# -blockgapchar
	# -unbracketblocks
	# 
	
	$params{'-bracketedchar'} = if not defined $params{'-bracketedchar'};
	 	
 	 
	# set the filehandle	

	
	#$outbracketedAsN = 0 unless defined $outbracketedAsN;
	
	#die "outmode not defined in out_mx" unless defined($outmode);
	
	# *MXOF = *FH;
		
	#my ($totchars, $intlvstart, $blk) = (0) x 2;
	
	#$outslice->remove("Blocks", $outslice->Jrna::slice_stems($stemmode) );

	# my $working_interleave = $orig_interleave->Jrna::intlv_working($outslice, $outintlv, $outintlvsize);

	#$outcolstart = (&tax_longest_name($outslice) + 1); # if ($outmode == 3); ## ---------- REMOVE FROM PARAMETER LIST --------------
	
#	print "\n\n--- generating a matrix with the following:---";
	$working_interleave->describe;
	#$outslice->describe;
		
	#compute total characters (ultimately a slice method)
	#foreach $blk (keys %{$outslice->blocks}) {
		#	if ($outunbracket == 1) {
			#		$totchars += $nbd[$blk][0];
			#	}
#		else {
#			if ($nbd[$blk][3] == 0) { 
#				$totchars += $nbd[$blk][0];
#			} 
#		}
#	}

	#end pre-parse
	
	# headers - remove to individual calls ultimately!!!
#	if (($outmode == 1) or ($outmode ==  6)) { #make it a paup file
	
#		print MXOF "#Nexus\n\n"; # leave this as first line, thus supposedly readable in TNT -> but its not right now
#		
#		print MXOF &str_file_header();
#		print MXOF "\n\n";
	
#		print MXOF "Begin data;\n\n";
#		print MXOF "\tDimensions\n";
#		print MXOF "\t\tntax = ", $outslice->total("Taxa"), "\n";
#		print MXOF "\t\tnchar = $totchars;\n";
#		print MXOF "\tFormat \n";		
#		print MXOF "\t\tdata = RNA\n" unless $outmode == 5;  # needed for (6) Bayes, but not tnt (5)?
#		print MXOF "\t\tmissing = ?\n";
#		print MXOF "\t\tgap = - \n"; 
#		if ($outintlv == 1 or $outintlv == 2 ) { print MXOF "\t\tinterleave\n"; }    
#		print MXOF "\t\t";
#		print MXOF 'symbols="0 1 2 3 4 5 6 7 8 9 A B C D E F G H I J K L M N P Q R S U V W X Y Z !";'."\n" unless $outmode == 6;
#		if ($outmode == 6) {print MXOF ";\n"};
		
#		print MXOF "\n\tMatrix\n\n";		
#	}
		
#	elsif ($outmode == 2)  { # make it a phylip file
#		print MXOF int (keys %{$outslice->taxa}); 
#		print MXOF " $totchars\n";
#	}

#	elsif ($outmode == 3) { # make it html
#		# no special computation needed as of yet
#	}

#	elsif ($outmode == 5) { #make it TNT		
#		print MXOF "nstates dna;\n"; # JRNA DOES DNA, if you want change it for morphology by leaving it out
#		print MXOF "xread\n";
		
		#note order switched relative to PHLYIP, silly competitive silliness
#		print MXOF "$totchars ";
#		print MXOF int (keys %{$outslice->taxa}); 
#		print MXOF "\n";
#	}

#	my @outdesclist = (4 , 5); # crude way of doing this, 4 and 5 are the $nbd references for a given block, for each new desc. block added will have to increment this too

#	my ($tmpblkfocus, $dblk); # dblk is odd

	foreach my $intlv ($params{'-interleave'}->interleaves) {  #  $working_interleave->interleaves) {
	



	}
	
	#headers
		#should make these subroutines, more clearly called and re-orderd for different matrices
		
		#include description header on top? - only really makes sense when including block numbering as well ************* make this a subroutine!!! **********
		if ($outdeschead == 1) {
			$tmpblkfocus = $nbd[$working_interleave->first_block_in_interleave($intlv)][9];
			print MXOF '<div class="descrow">';
			print MXOF "[ $deschead[$nbd[$tmpblkfocus][9]][0] $deschead[$nbd[$tmpblkfocus][9]][1] - (";
			
			foreach my $blk ($working_interleave->blocks($intlv)) {	
				if ($nbd[$dblk][9] == $tmpblkfocus) {
					print MXOF " $dblk";
				}
				else {
					print MXOF ")  <|  |> $deschead[$nbd[$tmpblkfocus][9]][0] $deschead[$nbd[$tmpblkfocus][9]][1] - ($dblk";
					$tmpblkfocus = $nbd[$dblk][9];
				}			
			}
			print MXOF ")  ]";
			print MXOF "</div>\n";
		}

		# description of individual blocks
		if ($outdescblks == 1) {
			for (my $desc=0; $desc < $#outdesclist+1; $desc ++) { #loop through the total list of interleaved descriptions
		
				print MXOF &css_tagswitch($outmode, 0, 0, "descrow");
				print MXOF &css_tagswitch($outmode, 1, 0, "rightholder");
				
				print MXOF &padright("["," ", $outcolstart);
				
				print MXOF &css_tagswitch($outmode, 1, 1, "rightholder");
					
				print MXOF &css_tagswitch($outmode, 1, 0, "desctext");
			
				foreach my $blk ($working_interleave->blocks($intlv)) {
					if ($nbd[$blk][3]==1 and $outbrackets == 1) { 
				
						print MXOF $nbd[$blk][$outdesclist[$desc]];		#may turn this into a series of descriptive blocks, one after another. (i.e. 4+5 etc)
						print MXOF &str_gapchar($outgapped, $outgapchar);	
					}
	
					elsif ($nbd[$blk][3]==1 and $outbrackets == 0) {} # do nothing

					else {
						print MXOF $nbd[$blk][$outdesclist[$desc]];		#may turn this into a series of descriptive blocks, one after another. (i.e. 4+5 etc)
						print MXOF &str_gapchar($outgapped, $outgapchar);
					}
				}
				
				print MXOF &css_tagswitch($outmode, 1, 1, "desctext");
		
			print MXOF "]";
			
				print MXOF &css_tagswitch($outmode, 0, 1, "descrow");
				print MXOF "\n";
			}
		}

		# mask on top -> fix for gap mode!!
		if ($outmode == 1 or ($outmode ==6) ) { # output the mask for nexus files
			print MXOF &padright("["," ",$outcolstart);
		
			foreach my $z ($outslice->loop("Blocks")) {
				if ($nbd[$z][3]==0) {print MXOF $nbd[$z][8]}
				#else { print MXOF" " x $nbd[$z][0]};
			}
			print MXOF "]\n";
		}
		
		# block numbered top?			
		if (($outnumblocks == 1) or ($outnumblocks == 3)) {
			print MXOF $working_interleave->Jrna::intlv_numseq($intlv, 1, $outmode, 1, $outgapped, $outgapchar, $outbrackets, $outcolstart);
		}
	
		# column numbered top?
		if (($outnumcols == 1) or ($outnumcols == 3)) {
			print MXOF $working_interleave->Jrna::intlv_numseq($intlv, 0, $outmode, 1, $outgapped, $outgapchar, $outbrackets, $outcolstart);
		}
			
		# matrix		
			foreach my $t ($outslice->loop("Taxa")) { 						# loop through taxa
				print MXOF &css_tagswitch($outmode, 0, 0, "mxrow");
				print MXOF &css_tagswitch($outmode, 1, 0, "taxa");
				print MXOF &padright($ters{$t}, ' ', $outcolstart);		# print taxon	
				print MXOF &css_tagswitch($outmode, 1, 1, "taxa");
				print MXOF &css_tagswitch($outmode, 1, 0, "seq");
				
				foreach my $blk ($working_interleave->blocks($intlv)) {		# loop through blocks			
					if (($nbd[$blk][3]== 1) and ($outbrackets == 1)) { 		# bracketed block and print brackets
		
						print MXOF "[" if $outunbracket == 0;

						if ($outmode == 3) {
							print MXOF &css_colorRNAstr($mx[$t][$blk][0]);
						}
	
						else {	
							if ($outbracketedAsN == 1) {
								print MXOF "N" x $nbd[$blk][0];
							}
							else {
								print MXOF $mx[$t][$blk][0];	
							}
						}		
					
						print MXOF "]" if $outunbracket == 0;;
						print MXOF &str_gapchar($outgapped, $outgapchar);
					}


					elsif (($nbd[$blk][3] ==1) and ($outbrackets == 0)) {			# bracketed block and don't print brackets
					}

					else  {	# ( 0 and 0 ) 											# non-bracketed block 
						if ($outmode == 3) {
							print MXOF &css_colorRNAstr($mx[$t][$blk][0]);
						}
						
						elsif ($outmode == 6) {
							print MXOF &str_out_mx_corrected($mx[$t][$blk][0], "mrbayes")
						}
						
						else {
							print MXOF $mx[$t][$blk][0];
						}
						
					print MXOF &str_gapchar($outgapped, $outgapchar);
					}
				}

				print MXOF &css_tagswitch($outmode, 1, 1, "seq");
				print MXOF &css_tagswitch($outmode, 0, 1, "mxrow");

			print MXOF "\n";
			
		}	# end taxon loop
	
		# blocks numbered at bottom
		if ($outnumblocks == 2 or ($outnumblocks ==3)) {
			print MXOF $working_interleave->Jrna::intlv_numseq($intlv, 1, $outmode, 1, $outgapped, $outgapchar, $outbrackets, $outcolstart);
		}
	
		# column numbered at bottom
		if ($outnumcols == 2 or ($outnumcols == 3)) {
			print MXOF $working_interleave->Jrna::intlv_numseq($intlv, 0, $outmode, 1, $outgapped, $outgapchar, $outbrackets, $outcolstart);
		}
		
		print MXOF &css_txtswitch($outmode, "\n\n<div class=\"spacer\"></div>\n");
	}
	 	
	if (($outmode == 1) or ($outmode == 6)) { print MXOF ";\nend;"; }

	if ($outmode == 5) { print MXOF ";\n"; } #close the matrix on tnt output
	
	return ($totchars);  # need to output end position for numbering here

	
	# requires -mx

}

