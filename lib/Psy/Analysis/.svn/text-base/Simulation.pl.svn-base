

package simulation;


# need to: un-hard-code paths/filenames, add a paramter array etc. etc.
# requires the hmmer suite of tools to be in the $PATH as numerous calls to the shell are made

sub test_hmmer () {	
	my $my_data = shift;		# a slice object

	my @insize = qw/1 2 5 10 20 50 100 200/;
	my $builds;
	
	my $root_file = "$data->{Modelroot}.stk";
	print "root: $root_file\n";

	my $new_file;
	my $new_slice = slice->new;
	
	my $models;
	
	# run the all case just once
		$new_slice = $my_data->new_random_slice( 
														"numtaxamax" => 'all', 
														"numtaxamin" => 'all',
														"numblocksmin" => "all",
														"numblocksmax" => "all"
												);
		$data->out_stockholm($new_slice);
		`rename $root_file ht_290_rep_0`;
		`hmmbuild -F --nucleic -n ht_290_rep_0 -A models.hmm ht_290_rep_0`; # 
		$models->{"ht_290_rep_0"} = defined;				
		
	chdir 'd:\lib\jrna\analyses\stockholm\ichs2' or die "can't change to dir"; # gotta run it from the root directory at present

		foreach my $size (@insize) {
			for (my $i=1; $i<11; $i++) {
				$new_file = "ht_$size"."_rep_$i";	
			
				$models->{$new_file} = defined;
				
				do {
					$new_slice = $my_data->new_random_slice( 
														"numtaxamax" => $size, 
														"numtaxamin" => $size,
														"numblocksmax" => "all"
													);	
													
					&Jrna::slice_prune_uninformative_taxa($new_slice); # strip zeroed taxa-> will force a regeneration

				} until ($new_slice->total("Taxa") == $size);
				
				$data->out_stockholm($new_slice);	
				
				`rename $root_file $new_file`;
				`hmmbuild -F --nucleic -n $new_file -A models.hmm $new_file`; # 
			}
	}

	# all datasets are calibrated
	`hmmcalibrate --mean 850 --sd 300 models.hmm`;
	
	`hmmindex models.hmm`;
	
	foreach my $mdl (keys %{$models}) {
		`hmmfetch models.hmm $mdl > current_model`;
		`hmmsearch current_model ichs_seqs_28S_query_a.txt > results\\$mdl.out`
	}

	#`cd results`;
	#`copy /b *.out all_output`;

	print Dumper($models);
}

1;


__END__
