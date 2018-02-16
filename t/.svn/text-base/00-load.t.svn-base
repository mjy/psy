#!perl -T

use Test::More ('tests' => 40);


BEGIN {
	# there are a number of CPAN required packages, only Statistics::Descriptive , these are listed first
	use_ok( 'Statistics::Descriptive');
	use_ok( 'Statistics::Distributions');
	use_ok( 'Template' );
	use_ok( 'Data::Dumper' );
	use_ok( 'BIO::DB::Genbank' );
	use_ok( 'Scalar::Util' );

	# the remaining modules are packaged with Psy
	use_ok( 'Psy::Psy' );
	#use_ok( 'Psy::Analysis::Blossum' );
	#use_ok( 'Psy::Analysis::Blossum2' );
	use_ok( 'Psy::Analysis::Seqpair' );
	use_ok( 'Psy::Dna::Alphabet' );
	use_ok( 'Psy::Dna::Iupac' );
	use_ok( 'Psy::Helpers::Bayes2gelr' );
	use_ok( 'Psy::Helpers::Fastfasta' );
	use_ok( 'Psy::Helpers::Fastgenbank' );
	use_ok( 'Psy::Helpers::Fasthmmer' );
	use_ok( 'Psy::Helpers::Fastmatrix' );
	use_ok( 'Psy::Helpers::Gbname' );
	use_ok( 'Psy::Io::Io' );
	use_ok( 'Psy::Matrix::Block' );
	use_ok( 'Psy::Matrix::Column' );
	use_ok( 'Psy::Matrix::Format' );
	use_ok( 'Psy::Matrix::Interleave' );
	use_ok( 'Psy::Matrix::Slice' );
	use_ok( 'Psy::Matrix::Structure' );
	use_ok( 'Psy::Matrix::Template' );
	use_ok( 'Psy::Matrix::Terminal' );
	use_ok( 'Psy::Misc::Colour' );
	use_ok( 'Psy::Output::Wrapper' );
	use_ok( 'Psy::Output::Vimsyntax' );
	use_ok( 'Psy::Output::Tnt' );
	use_ok( 'Psy::Output::Stockholm' );
	use_ok( 'Psy::Output::Poy' );
	use_ok( 'Psy::Output::Phase' );
	use_ok( 'Psy::Output::Pairedparsimony' );
	use_ok( 'Psy::Output::Output' );
	use_ok( 'Psy::Output::Nexus' );
	use_ok( 'Psy::Output::Mitable' );
	use_ok( 'Psy::Output::Fasta' );
	use_ok( 'Psy::Output::Ct' );
	use_ok( 'Psy::Output::Column_template' );
	use_ok( 'Psy::Output::Column' );
	use_ok( 'Psy::Output::Blockmeta' );
	use_ok( 'Psy::Output::Blkmapped' );
	use_ok( 'Psy::Strings::Composition' );
	use_ok( 'Psy::Strings::Strings' );
}

diag( "Testing Psy::Psy $Psy::VERSION, Perl $], $^X" ); #  
