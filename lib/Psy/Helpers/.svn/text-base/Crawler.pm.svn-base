package Crawler;

use strict;
use vars qw(@ISA);


use warnings;
use Data::Dumper;	
use lib qw(lib/);

use Psy::Psy;
use Psy::Output::Output;

use File::Spec;
use File::Util;
use Cwd;

=head1 Psy::Helpers::Crawler

=cut

=head1 SYNOPSIS

Crawls through folders, finding data and trees, translates data, cues analyses, summarizes results (tree files).

Tree files must end in .tre
e,g, 
	foo.tre
	bar.tre

Datafiles must be matrices (PAUP legal is safe, but others might work) in the format "data"_<type>_<misc>*
e,g,
	data_rna18S.nex
	data_rna28s.fas
	data_rna28s_a.fas

<type> acts as a group name whereby results from similar types are indicated as such in the output and results reporting
<misc> is not required, allows multiple types in the same folder, which should be rare

You need a file called meta.txt in the -rootfolder, it has two columns, folder_name/dataset_name <tab> number of taxa.

At present doesn't do anything with structure, but could be implemented. 



Required parameters:
	

Optional parameters:

=cut

our $OUT_DIR = "foo";
our $LEGAL_CHARS =  Psy::Dna::Alphabet->new(); 



use FindBin qw($Bin);		# sets $Bin to the root directory
our $BASE_DIR = $Bin;	# this directory is the directory the script is executed from??
print "BASEDIR: $BASE_DIR;\n";

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

	$self->{'input_folder'} = $params{"-root_folder"};
	$self->{'datafile_index'} = 0; # each read file is uniquely indexed
	$self->{'analysis_index'} = 0; # as is each analysis
	
	$self->_crawl(%params);   # crawl through the folders and gather a report
 
	# load meta data about the inputfiles, so far just number taxa
	open (IN, "$self->{'input_folder'}/meta.txt") || die "can't open the meta file";
		my	@contents = <IN>;
	close IN;
	foreach my $l (@contents) {
		my @f = split /\s/, $l;
		$self->{'meta'}->{$f[0]} = $f[1];
	}

	# load a list of analyses to run
	open (IN, "$self->{'input_folder'}/analyses.txt") || die "can't open the meta file";
		@contents = <IN>;
	close IN;
	my $i = 0;
	foreach my $l (@contents) {
		my @f = split /\s/, $l;
		$self->{'analysis_loop'}->{$i} = \@f; # each line has a list of kword sizes
		$i++;
	}

	return 1;
}

sub run {
	my $self = shift;
	
	my @params = @_;	

	grep '-convert', @params && $self->_convert(); # convert the data to the required input
	grep '-analyze', @params && $self->_analyze(); # run some analyses

	#	$params{'-collect'} && $self->_collect(); # collect the results

}

sub _crawl {
	# builds the catalog of files
	my $self = shift;
	
	my %params = @_;	
	$params{'-root_folder'} ||= File::Spec->curdir();
	
	my $f = File::Util->new();
	print  join "\n", $f->list_dir( File::Spec->catdir( $params{'-root_folder'} ),'--dirs-only');
	
	$self->datasets( $f->list_dir( File::Spec->catdir( $params{'-root_folder'} ),'--dirs-only') );
	
	print "\n\n";
	
	foreach my $d ($self->datasets) {
		my @files = $f->list_dir(  File::Spec->catdir($params{'-root_folder'}, $d) , '--files-only');
		$self->_indexFiles(
			'-folder' =>  $d ,
		    '-files' => \@files
		);
		print $d, "\t",  File::Spec->catdir( $params{'-root_folder'}, $d) , "\n";
		print  join "\n", $f->list_dir( File::Spec->catdir( $params{'-root_folder'}, $d) , '--files-only');	
		print "\n\n";
	}
	# loop through the folder...
	# {folder}->{tree}
	# {folder}->{data}->{file1}
	# {folder}->{data}->{file2}
	# folder acts as the root	
}


sub analysis_loop {
	my $self = shift;
	return sort keys %{$self->{'analysis_loop'}}
}

sub analysis_params {
	my $self = shift;
	my $f = shift;
	return $self->{'analysis_loop'}->{$f}
}

sub _convert {
	my $self = shift;
	my %params = @_;

	my @analysis = <DATA>; # hack to get this appended to the end
	
	my $dir = getcwd;
	
	# loop through the input
	
	for my $f ($self->datasets) {
		for my $i ($self->datafileIndexes('-folder' => $f) ) {

			# loop through the analyses
			for my $analysis ($self->analysis_loop) { # a loop through integers
						
	#			print "$f\t$i\n";
				
				# load the data
				my $file = $self->datafile('-folder' => $f, '-index' => $i);
				my $foo = "$dir\\$self->{'input_folder'}/$f/";
				my $p = Psy->new(
					'-path' => File::Spec->catdir( $foo)."/", # hmm catdir needs to be integrated to Psy
					'-matrix_file' => $file, # so a hack around
					'-matrix_label' => 'temp',
					'-number_terminals' => $self->{'meta'}->{"$f/$file"}
					);
			
				# print Dumper($p);
					
				# export the translation (MODULARIZE LATER)
				my $o = output->new();
				my $output = $i."_".$self->datafile_index('-folder' => $f, '-index' => $i)."_$analysis";

				$o->Kword(
					'-mx' => $p->mx('-matrix_label' => 'temp'),
					'-translation_mode' => 'all', 
					'-out_format' => 'tnt',
					'-file_name' => $output,
				 	'-kword_size' => $self->analysis_params($analysis)
				
				);
				
				open (APPEND, ">>analyses/kword/$output.tnt") || die "can't append the analysis analyses/kword/$output";
					print APPEND @analysis; 
				close APPEND;
				
				# cue the file to analysis list ...
				$self->add_analysis("$output\.tnt");
							
				chdir("$dir");
			}
			# contat the required analysis to the end of the tnt file here-

		}

	}
	
	# create a new output object
	# load the data
	# generate the output
	# destroy the data (free memory)
	# log the 
	1;
}

sub add_analysis {
	my $self = shift;

	if (@_) { $self->{'analyses'}->{$self->{'analysis_index'}} = shift; 	
		$self->{'analysis_index'} += 1;
	}
}

sub return_analysis {
	my $self = shift;
	return $self->{'analyses'}->{shift @_}
}


sub analyses {
	my $self = shift;
	my @foo;
	map {push @foo, $self->return_analysis($_) } (sort keys %{$self->{'analyses'}});
	return @foo
}

sub results {
	my $self = shift;
	return sort keys %{$self->{'result'}};
}


sub _analyze {
	my $self = shift;	
	foreach my $r ($self->analyses) {
		print "R:::::::::::", $r;
		`tntc mxram 384 run $BASE_DIR/out/foo/analyses/kword/$r;`;

		# concatonate the result file 
		open (RESULT, ">>results.tre") || die "can't append the analysis analyses/kword/$r";

				# open a file that we can append results to
				 
				open (TREES, "result") || die "can't open a tree results file";
					my @r = <TREES>;
				close TREES;
					# TNT specific
					# don't want the first or last lines in the tnt output
					shift @r;
					pop @r;

				# log the total results in a result array

				$self->{'result'}->{$r} = $#{r} + 1; # the number of trees
					
						# steps to convert STUPID TNT tree format to Newick (hint to Pablo- \s and , are both characters)

						map {$_ =~ s/ /\,/g } @r;
						map {$_ =~ s/\,\)/\)/g } @r;
						map {$_ =~ s/\)\(/\)\,\(/g } @r;
						map {$_ =~ s/\*/\;/g } @r;
					print RESULT @r;
		close RESULT;
		unlink(qw/result/);
		
		open (RESULT, ">>tree_index.txt") || die "can't open tree log";
		for my $l ($self->results) {
			print RESULT "$l\t";
			print RESULT $self->{'result'}->{$l};
			print RESULT "\n";
		}
		close RESULT;
			
		
	}
	



}


# accessors

=head2 _indexFiles

Indexes the -files array from a given -folder to either tree or data.

=cut


sub _indexFiles {
	my $self = shift;
	my %params = @_;
	print "indexed: ", $params{"-files"};
	foreach my $f (@{$params{'-files'}})  {
	
		if ($f =~ /\.tre$/) {
			push @{$self->{'datafiles'}->{$params{'-folder'}}->{'trees'}}, $f; # an array
		}
		elsif ($f =~ /^data/) {
		 	$self->{'datafiles'}->{$params{'-folder'}}->{'data'}->{ $self->{'datafile_index'} } = $f ; # a hash
			$self->{'datafile_index'} += 1;
		}
	}
	1;
}


=head2 datasets

List the datasets as an array

=cut

sub datasets {
	my $self = shift;
	if (@_) { foreach (@_) { ($self->{'datasets'}->{  $_   } = undef ) unless $_ =~ /\..*/  } 
		}
	return (sort keys %{$self->{'datasets'}})
}




=head2 datafile_index

Returns a unique index reference for the given -file in -folder.

=cut

sub datafile_index {
	my $self = shift;
	my %params = @_;
	$self->{'datafiles'}->{$params{'-folder'}}->{'data'}->{$params{'-index'}} 
}
	



sub datafile {
	my $self = shift;
	my %params = @_;
	$self->{'datafiles'}->{$params{'-folder'}}->{'data'}->{$params{'-index'}} 
}
	



=head2 datafiles

Return the datafiles for -folder

=cut

sub datafiles {
	my $self = shift;
	my %params = @_;
	return (sort values %{ $self->{'datafiles'}->{$params{'-folder'}}->{'data'}});
}

sub datafileIndexes {
	my $self = shift;
	my %params = @_;
	return (sort keys %{ $self->{'datafiles'}->{$params{'-folder'}}->{'data'}});
}



=head2 treefiles

Return the treefile names for -folder.  'trees' is an array, not has hash

=cut

sub treefiles {
	my $self = shift;
	my %params = @_;
	return sort  $self->{'datafiles'}->{ $params{'-folder'}  }->{'trees'}  ;
}



# this is the analysis that gets added to the tnt files at present

__DATA__


taxname=;
hold 1000;
tsave *result; 
sectsch: xbuf;
xmult /
	hits 5
	replications 10
	rss
	css
	xss
	fuse 1
	drift 1
	ratchet 2
	autoconst 3
	level 5
	chklevel 3; 
bbreak;
collapse;
save;
tsave /;
quit;



