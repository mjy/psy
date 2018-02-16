#!/usr/local/bin/perl

# ver 0.02

# a crude multi/find replace, using two columns, a key->value
# usage: 'replace infile <filename> replacefile <filename>
# the outfile is optional, if not included

use File::Copy;

my $usage = "\nusage: 'perl replace.pl infile <filename> replacefile <filename>'";

my %p = @ARGV;

# because multiple replaces (-i) destroy infile, the real backup is copied to .original first

copy($p{infile}, $p{infile}.".original");

open (REP, $p{replacefile}) || die "can't open ", $p{replacefile}, $usage, "\n";
	while (<REP>) {
		
		chomp;
		#	next if $_ = "";
		@r = split /\s+/, $_;
		`perl -pi.bak -e "s/$r[0]/$r[1]/g" $p{'infile'}`
	}
close REP;

