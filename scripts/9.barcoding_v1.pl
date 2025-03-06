#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Find;
use Data::Dumper;
$|=1;
#array of hash plus transposed matrix(two-dimensional array)

my $global_options=&argument();
my $input_directory=&default("input","input");
my $pattern=&default("9","pattern");
my $output_directory=&default("output","output");

my $osname=$^O;
if ($osname eq "MSWin32") {
	system("del/f/s/q $output_directory") if (-e $output_directory);
}elsif ($osname eq "cygwin") {
	system("rm -rf $output_directory") if (-e $output_directory);
}elsif ($osname eq "linux") {
	system("rm -rf $output_directory") if (-e $output_directory);
}elsif ($osname eq "darwin") {
	system("rm -rf $output_directory") if (-e $output_directory);
}
mkdir ($output_directory) if (!-e $output_directory);


my $pattern1=".fasta";
my $pattern2=".fas";
my $pattern3=".fa";
my $pattern4=".fsa";
my @filenames;
find(\&target,$input_directory);
sub target{
	if (/$pattern1/ or /$pattern2/ or /$pattern3/ or /$pattern4/){
		push @filenames,"$File::Find::name";
	}
	return;
}


my (%percentage_variable,%percentage_pis,%number_variable,%number_pis);
while (@filenames) {
	my $filename_fasta=shift @filenames;
	my $target_name=substr($filename_fasta,0,rindex($filename_fasta,"\."));
	$target_name=~ s/\s/_/g;
	my $filename=substr($target_name,rindex($target_name,"\/")+1);
	my $output_filename_temp="$input_directory/$filename\_temp";

	open(my $input_temp,"<",$filename_fasta);
	open(my $output_temp,">",$output_filename_temp);
	my $row=<$input_temp>;
	print $output_temp $row;
	while (<$input_temp>){
		$_=~ s/\r|\n//g;
		if ($_=~ /^>/) {
			print $output_temp "\n".$_."\n";
		}else{
			print $output_temp $_;
		}
	}
	print $output_temp "\n";
	close $input_temp;
	close $output_temp;


	open(my $input,"<",$output_filename_temp);
	my ($header,$sequence,$length,@id,@array);
	while (defined ($header=<$input>) && defined ($sequence=<$input>)) {
		$header=~ s/\r|\n//g;
		$sequence=~ s/\r|\n//g;
		$sequence=uc $sequence;
		$length=length $sequence;
		push @id,$header;

		for (0..$length-1){
			push @{$array[$_]},substr($sequence,$_,1);
		}
	}
	#print Dumper \@array;
	close $input;
	unlink("$output_filename_temp");


	my $total_length=@array;
	my (@variable_array,@variable_array2);#variable/polymorphic sites, first array_ref is fisrt column
	my (@invariable_array,@invariable_array2);#invariable/monomorphic sites, first array_ref is fisrt column
	my (@gap_array,@gap_array2);#gap sites, first array_ref is fisrt column
	my (@variable_sites_2nucleotides_array,@variable_sites_2nucleotides_array2);
	my (@variable_sites_3nucleotides_array,@variable_sites_3nucleotides_array2);
	my (@variable_sites_4nucleotides_array,@variable_sites_4nucleotides_array2);
	my (@parsimony_informative_sites_array,@parsimony_informative_sites_array2);
	my (@non_parsimony_informative_sites_array,@non_parsimony_informative_sites_array2);
	for my $column (0..$#array){#1 of 36 columns
		for my $nucleotide ($array[$column]){#1 of 4 nt in each column
			if ((grep {$_ ne $nucleotide->[0]} @$nucleotide) and ((!grep {$_ eq "-"} @$nucleotide) and (!grep {$_ eq "N"} @$nucleotide) and (!grep {$_ eq "?"} @$nucleotide)) and (@$nucleotide >= 2)) {
				push @variable_array,[@$nucleotide];
				push @variable_array2,$column;
				#delete $array[$column];#undefined
				my %nt;
				foreach my $key (@$nucleotide) {
					$nt{$key}++;
				}
				my @pis_array=values %nt;
				if ((keys %nt) == 2) {#A T T T
					push @variable_sites_2nucleotides_array,[@$nucleotide];
					push @variable_sites_2nucleotides_array2,$column;
				}elsif (((keys %nt) == 3) and (@$nucleotide >= 3)) {#A A T G
					push @variable_sites_3nucleotides_array,[@$nucleotide];
					push @variable_sites_3nucleotides_array2,$column;
				}elsif (((keys %nt) == 4) and (@$nucleotide >= 4)) {#A T G C
					push @variable_sites_4nucleotides_array,[@$nucleotide];
					push @variable_sites_4nucleotides_array2,$column;
				}
				if (((keys %nt) >= 2) and (!grep {$_ == 1} @pis_array) and (@$nucleotide >= 4)) {#A A T T
					push @parsimony_informative_sites_array,[@$nucleotide];
					push @parsimony_informative_sites_array2,$column;
				}else {#A T T T
					push @non_parsimony_informative_sites_array,[@$nucleotide];
					push @non_parsimony_informative_sites_array2,$column;
				}
			}elsif ((!grep {$_ ne $nucleotide->[0]} @$nucleotide) and ((!grep {$_ eq "-"} @$nucleotide) and (!grep {$_ eq "N"} @$nucleotide) and (!grep {$_ eq "?"} @$nucleotide) ) and (@$nucleotide >= 2)) {#A A A A, not N N N N, not ? ? ? ?
				push @invariable_array,[@$nucleotide];
				push @invariable_array2,$column;
				#delete $array[$column];#undefined
			}elsif (((grep {$_ eq "-"} @$nucleotide) or (grep {$_ eq "N"} @$nucleotide) or (grep {$_ eq "?"} @$nucleotide)) and (@$nucleotide >= 2)) {#- A T T or N A T T or ? A T T
				push @gap_array,[@$nucleotide];
				push @gap_array2,$column;
				#delete $array[$column];#undefined
			}
		}
	}
	##@array=grep {defined $_} @array;#constant sites,remove undefined elements




	#variable/polymorphic sites with all nucleotides number and position in alignment matrix
	if (($pattern == 1) or ($pattern == 9) or ($pattern == 10)) {
		open(my $output1,">","$output_directory/$filename\_p1_variable_sites_allnucleotides.txt");#Variable (polymorphic) sites (>= 2 species)
		open(my $output1_2,">","$output_directory/$filename\_p1_variable_sites_allnucleotides.fasta");#Variable (polymorphic) sites (>= 2 species)
		print $output1 "Total sequence length of alignment matrix $filename:\n";
		print $output1 $total_length."\n\n";
		print $output1 "Number of variable/polymorphic site(s) with all nucleotides in alignment matrix $filename:\n";
		if (@variable_array != 0) {
			my $variable_number=@variable_array;
			print $output1 $variable_number."\n\n";
			$number_variable{$filename}=$variable_number."\t".$total_length;
			print $output1 "Percentage of variable/polymorphic site(s) with all nucleotides in alignment matrix $filename:\n";
			print $output1 sprintf("%.4f",$variable_number/$total_length)."\n\n";
			$percentage_variable{$filename}=sprintf("%.4f",$variable_number/$total_length);
			my @variable_position;
			foreach (@variable_array2) {
				push @variable_position,$_+1;
			}
			my $variable_position=join " ",@variable_position;
			print $output1 "Position of variable/polymorphic site(s) with all nucleotides in alignment matrix $filename:\n";
			print $output1 $variable_position."\n\n";
			print $output1 "Variable/polymorphic site(s) with all nucleotides in alignment matrix $filename:\n";
		}

		#variable/polymorphic sites with all nucleotides, first array_ref is first row
		if (@variable_array != 0) {
			my @transposed_variable_array=map {my $x=$_;[map {$variable_array[$_][$x]} (0..$#variable_array)]} (0..$#{$variable_array[0]});
			for (0..(@id-1)){
				print $output1 $id[$_],"\n",@{$transposed_variable_array[$_]},"\n";
				print $output1_2 $id[$_],"\n",@{$transposed_variable_array[$_]},"\n";
			}
			close $output1;
			close $output1_2;
		}else {
			close $output1;
			close $output1_2;
			unlink("$output_directory/$filename\_p1_variable_sites_allnucleotides.txt");
			unlink("$output_directory/$filename\_p1_variable_sites_allnucleotides.fasta");
		}
	}




	#variable/polymorphic sites with 2 nucleotides number and position in alignment matrix
	if (($pattern == 2) or ($pattern == 10)) {
		open(my $output2,">","$output_directory/$filename\_p2_variable_sites_2nucleotides.txt");#Variable (polymorphic) sites with 2 type (e.g., A,T) of nucleotides (>= 2 species)
		open(my $output2_2,">","$output_directory/$filename\_p2_variable_sites_2nucleotides.fasta");#Variable (polymorphic) sites with 2 type (e.g., A,T) of nucleotides (>= 2 species)
		print $output2 "Total sequence length of alignment matrix $filename:\n";
		print $output2 $total_length."\n\n";
		print $output2 "Number of variable/polymorphic site(s) with 2 nucleotides in alignment matrix $filename:\n";
		if (@variable_sites_2nucleotides_array != 0) {
			my $variable_sites_2nucleotides_number=@variable_sites_2nucleotides_array;
			print $output2 $variable_sites_2nucleotides_number."\n\n";
			my @variable_sites_2nucleotides_position;
			foreach (@variable_sites_2nucleotides_array2) {
				push @variable_sites_2nucleotides_position,$_+1;
			}
			my $variable_sites_2nucleotides_position=join " ",@variable_sites_2nucleotides_position;
			print $output2 "Position of variable/polymorphic site(s) with 2 nucleotides in alignment matrix $filename:\n";
			print $output2 $variable_sites_2nucleotides_position."\n\n";
			print $output2 "Variable/polymorphic site(s) with 2 nucleotides in alignment matrix $filename:\n";
		}

		#variable/polymorphic sites with 2 nucleotides, first array_ref is first row
		if (@variable_sites_2nucleotides_array != 0) {
			my @transposed_variable_sites_2nucleotides_array=map {my $x=$_;[map {$variable_sites_2nucleotides_array[$_][$x]} (0..$#variable_sites_2nucleotides_array)]} (0..$#{$variable_sites_2nucleotides_array[0]});
			for (0..(@id-1)){
				print $output2 $id[$_],"\n",@{$transposed_variable_sites_2nucleotides_array[$_]},"\n";
				print $output2_2 $id[$_],"\n",@{$transposed_variable_sites_2nucleotides_array[$_]},"\n";
			}
			close $output2;
			close $output2_2;
		}else {
			close $output2;
			close $output2_2;
			unlink("$output_directory/$filename\_p2_variable_sites_2nucleotides.txt");
			unlink("$output_directory/$filename\_p2_variable_sites_2nucleotides.fasta");
		}
	}




	#variable/polymorphic sites with 3 nucleotides number and position in alignment matrix
	if (($pattern == 3) or ($pattern == 10)) {
		open(my $output3,">","$output_directory/$filename\_p3_variable_sites_3nucleotides.txt");#Variable (polymorphic) sites with 3 type (e.g., A,T,G) of nucleotides (>= 3 species)
		open(my $output3_2,">","$output_directory/$filename\_p3_variable_sites_3nucleotides.fasta");#Variable (polymorphic) sites with 3 type (e.g., A,T,G) of nucleotides (>= 3 species)
		print $output3 "Total sequence length of alignment matrix $filename:\n";
		print $output3 $total_length."\n\n";
		print $output3 "Number of variable/polymorphic site(s) with 3 nucleotides in alignment matrix $filename:\n";
		if (@variable_sites_3nucleotides_array != 0) {
			my $variable_sites_3nucleotides_number=@variable_sites_3nucleotides_array;
			print $output3 $variable_sites_3nucleotides_number."\n\n";
			my @variable_sites_3nucleotides_position;
			foreach (@variable_sites_3nucleotides_array2) {
				push @variable_sites_3nucleotides_position,$_+1;
			}
			my $variable_sites_3nucleotides_position=join " ",@variable_sites_3nucleotides_position;
			print $output3 "Position of variable/polymorphic site(s) with 3 nucleotides in alignment matrix $filename:\n";
			print $output3 $variable_sites_3nucleotides_position."\n\n";
			print $output3 "Variable/polymorphic site(s) with 3 nucleotides in alignment matrix $filename:\n";
		}

		#variable/polymorphic sites with 3 nucleotides, first array_ref is first row
		if (@variable_sites_3nucleotides_array != 0) {
			my @transposed_variable_sites_3nucleotides_array=map {my $x=$_;[map {$variable_sites_3nucleotides_array[$_][$x]} (0..$#variable_sites_3nucleotides_array)]} (0..$#{$variable_sites_3nucleotides_array[0]});
			for (0..(@id-1)){
				print $output3 $id[$_],"\n",@{$transposed_variable_sites_3nucleotides_array[$_]},"\n";
				print $output3_2 $id[$_],"\n",@{$transposed_variable_sites_3nucleotides_array[$_]},"\n";
			}
			close $output3;
			close $output3_2;
		}else {
			close $output3;
			close $output3_2;
			unlink("$output_directory/$filename\_p3_variable_sites_3nucleotides.txt");
			unlink("$output_directory/$filename\_p3_variable_sites_3nucleotides.fasta");
		}
	}




	#variable/polymorphic sites with 4 nucleotides number and position in alignment matrix
	if (($pattern == 4) or ($pattern == 10)) {
		open(my $output4,">","$output_directory/$filename\_p4_variable_sites_4nucleotides.txt");#Variable (polymorphic) sites with 4 type (e.g., A,T,G,C) of nucleotides (>= 4 species)
		open(my $output4_2,">","$output_directory/$filename\_p4_variable_sites_4nucleotides.fasta");#Variable (polymorphic) sites with 4 type (e.g., A,T,G,C) of nucleotides (>= 4 species)
		print $output4 "Total sequence length of alignment matrix $filename:\n";
		print $output4 $total_length."\n\n";
		print $output4 "Number of variable/polymorphic site(s) with 4 nucleotides in alignment matrix $filename:\n";
		if (@variable_sites_4nucleotides_array != 0) {
			my $variable_sites_4nucleotides_number=@variable_sites_4nucleotides_array;
			print $output4 $variable_sites_4nucleotides_number."\n\n";
			my @variable_sites_4nucleotides_position;
			foreach (@variable_sites_4nucleotides_array2) {
				push @variable_sites_4nucleotides_position,$_+1;
			}
			my $variable_sites_4nucleotides_position=join " ",@variable_sites_4nucleotides_position;
			print $output4 "Position of variable/polymorphic site(s) with 4 nucleotides in alignment matrix $filename:\n";
			print $output4 $variable_sites_4nucleotides_position."\n\n";
			print $output4 "Variable/polymorphic site(s) with 4 nucleotides in alignment matrix $filename:\n";
		}

		#variable/polymorphic sites with 4 nucleotides, first array_ref is first row
		if (@variable_sites_4nucleotides_array != 0) {
			my @transposed_variable_sites_4nucleotides_array=map {my $x=$_;[map {$variable_sites_4nucleotides_array[$_][$x]} (0..$#variable_sites_4nucleotides_array)]} (0..$#{$variable_sites_4nucleotides_array[0]});
			for (0..(@id-1)){
				print $output4 $id[$_],"\n",@{$transposed_variable_sites_4nucleotides_array[$_]},"\n";
				print $output4_2 $id[$_],"\n",@{$transposed_variable_sites_4nucleotides_array[$_]},"\n";
			}
			close $output4;
			close $output4_2;
		}else {
			close $output4;
			close $output4_2;
			unlink("$output_directory/$filename\_p4_variable_sites_4nucleotides.txt");
			unlink("$output_directory/$filename\_p4_variable_sites_4nucleotides.fasta");
		}
	}




	#parsimony informative sites number and position in alignment matrix
	if (($pattern == 5) or ($pattern == 9) or ($pattern == 10)) {
		open(my $output5,">","$output_directory/$filename\_p5_parsimony_informative_sites.txt");#Parsimony informative sites (>= 4 species)
		open(my $output5_2,">","$output_directory/$filename\_p5_parsimony_informative_sites.fasta");#Parsimony informative sites (>= 4 species)
		print $output5 "Total sequence length of alignment matrix $filename:\n";
		print $output5 $total_length."\n\n";
		print $output5 "Number of parsimony informative site(s) in alignment matrix $filename:\n";
		if (@parsimony_informative_sites_array != 0) {
			my $parsimony_informative_sites_number=@parsimony_informative_sites_array;
			print $output5 $parsimony_informative_sites_number."\n\n";
			$number_pis{$filename}=$parsimony_informative_sites_number."\t".$total_length;
			print $output5 "Percentage of variable/polymorphic site(s) with all nucleotides in alignment matrix $filename:\n";
			print $output5 sprintf("%.4f",$parsimony_informative_sites_number/$total_length)."\n\n";
			$percentage_pis{$filename}=sprintf("%.4f",$parsimony_informative_sites_number/$total_length);
			my @parsimony_informative_sites_position;
			foreach (@parsimony_informative_sites_array2) {
				push @parsimony_informative_sites_position,$_+1;
			}
			my $parsimony_informative_sites_position=join " ",@parsimony_informative_sites_position;
			print $output5 "Position of parsimony informative site(s) in alignment matrix $filename:\n";
			print $output5 $parsimony_informative_sites_position."\n\n";
			print $output5 "Parsimony informative site(s) in alignment matrix $filename:\n";
		}

		#parsimony informative sites, first array_ref is first row
		if (@parsimony_informative_sites_array != 0) {
			my @transposed_parsimony_informative_sites_array=map {my $x=$_;[map {$parsimony_informative_sites_array[$_][$x]} (0..$#parsimony_informative_sites_array)]} (0..$#{$parsimony_informative_sites_array[0]});
			for (0..(@id-1)){
				print $output5 $id[$_],"\n",@{$transposed_parsimony_informative_sites_array[$_]},"\n";
				print $output5_2 $id[$_],"\n",@{$transposed_parsimony_informative_sites_array[$_]},"\n";
			}
			close $output5;
			close $output5_2;
		}else {
			close $output5;
			close $output5_2;
			unlink("$output_directory/$filename\_p5_parsimony_informative_sites.txt");
			unlink("$output_directory/$filename\_p5_parsimony_informative_sites.fasta");
		}
	}




	#non-parsimony informative sites number and position in alignment matrix
	if (($pattern == 6) or ($pattern == 10)) {
		open(my $output6,">","$output_directory/$filename\_p6_non-parsimony_informative_sites.txt");#Non-parsimony informative sites (>= 2 species)
		open(my $output6_2,">","$output_directory/$filename\_p6_non-parsimony_informative_sites.fasta");#Non-parsimony informative sites (>= 2 species)
		print $output6 "Total sequence length of alignment matrix $filename:\n";
		print $output6 $total_length."\n\n";
		print $output6 "Number of non-parsimony informative site(s) in alignment matrix $filename:\n";
		if (@non_parsimony_informative_sites_array != 0) {
			my $non_parsimony_informative_sites_number=@non_parsimony_informative_sites_array;
			print $output6 $non_parsimony_informative_sites_number."\n\n";
			my @non_parsimony_informative_sites_position;
			foreach (@non_parsimony_informative_sites_array2) {
				push @non_parsimony_informative_sites_position,$_+1;
			}
			my $non_parsimony_informative_sites_position=join " ",@non_parsimony_informative_sites_position;
			print $output6 "Position of non-parsimony informative site(s) in alignment matrix $filename:\n";
			print $output6 $non_parsimony_informative_sites_position."\n\n";
			print $output6 "Non-parsimony informative site(s) in alignment matrix $filename:\n";
		}

		#non-parsimony informative sites, first array_ref is first row
		if (@non_parsimony_informative_sites_array != 0) {
			my @transposed_non_parsimony_informative_sites_array=map {my $x=$_;[map {$non_parsimony_informative_sites_array[$_][$x]} (0..$#non_parsimony_informative_sites_array)]} (0..$#{$non_parsimony_informative_sites_array[0]});
			for (0..(@id-1)){
				print $output6 $id[$_],"\n",@{$transposed_non_parsimony_informative_sites_array[$_]},"\n";
				print $output6_2 $id[$_],"\n",@{$transposed_non_parsimony_informative_sites_array[$_]},"\n";
			}
			close $output6;
			close $output6_2;
		}else {
			close $output6;
			close $output6_2;
			unlink("$output_directory/$filename\_p6_non-parsimony_informative_sites.txt");
			unlink("$output_directory/$filename\_p6_non-parsimony_informative_sites.fasta");
		}
	}




	#invariable site number and position in alignment matrix
	if (($pattern == 7) or ($pattern == 10)) {
		open(my $output7,">","$output_directory/$filename\_p7_invariable_sites.txt");#Invariable (monomorphic) sites (>= 2 species)
		open(my $output7_2,">","$output_directory/$filename\_p7_invariable_sites.fasta");#Invariable (monomorphic) sites (>= 2 species)
		print $output7 "Total sequence length of alignment matrix $filename:\n";
		print $output7 $total_length."\n\n";
		print $output7 "Number of invariable/monomorphic site(s) in alignment matrix $filename:\n";
		if (@invariable_array != 0) {
			my $invariable_number=@invariable_array;
			print $output7 $invariable_number."\n\n";
			my @invariable_position;
			foreach (@invariable_array2) {
				push @invariable_position,$_+1;
			}
			my $invariable_position=join " ",@invariable_position;
			print $output7 "Position of invariable/monomorphic site(s) in alignment matrix $filename:\n";
			print $output7 $invariable_position."\n\n";
			print $output7 "Invariable/monomorphic site(s) in alignment matrix $filename:\n";
		}

		#invariable sites,first array_ref is first row
		if (@invariable_array != 0) {
			my @transposed_invariable_array=map {my $x=$_;[map {$invariable_array[$_][$x]} (0..$#invariable_array)]} (0..$#{$invariable_array[0]});
			for (0..(@id-1)){
				print $output7 $id[$_],"\n",@{$transposed_invariable_array[$_]},"\n";
				print $output7_2 $id[$_],"\n",@{$transposed_invariable_array[$_]},"\n";
			}
			close $output7;
			close $output7_2;
		}else {
			close $output7;
			close $output7_2;
			unlink("$output_directory/$filename\_p7_invariable_sites.txt");
			unlink("$output_directory/$filename\_p7_invariable_sites.fasta");
		}
	}




	#gap site number and position in alignment matrix
	if (($pattern == 8) or ($pattern == 10)) {
		open(my $output8,">","$output_directory/$filename\_p8_gap_sites.txt");#Sites with gaps/missing data (>= 2 species)
		open(my $output8_2,">","$output_directory/$filename\_p8_gap_sites.fasta");#Sites with gaps/missing data (>= 2 species)
		print $output8 "Total sequence length of alignment matrix $filename:\n";
		print $output8 $total_length."\n\n";
		print $output8 "Number of gap site(s) in alignment matrix $filename:\n";
		if (@gap_array != 0) {
			my $gap_number=@gap_array;
			print $output8 $gap_number."\n\n";
			my @gap_position;
			foreach (@gap_array2) {
				push @gap_position,$_+1;
			}
			my $gap_position=join " ",@gap_position;
			print $output8 "Position of gap site(s) in alignment matrix $filename:\n";
			print $output8 $gap_position."\n\n";
			print $output8 "Gap site(s) in alignment matrix $filename:\n";
		}

		#gap sites,first array_ref is first row
		if (@gap_array != 0) {
			my @transposed_gap_array=map {my $x=$_;[map {$gap_array[$_][$x]} (0..$#gap_array)]} (0..$#{$gap_array[0]});
			for (0..(@id-1)){
				print $output8 $id[$_],"\n",@{$transposed_gap_array[$_]},"\n";
				print $output8_2 $id[$_],"\n",@{$transposed_gap_array[$_]},"\n";
			}
			close $output8;
			close $output8_2;
		}else {
			close $output8;
			close $output8_2;
			unlink("$output_directory/$filename\_p8_gap_sites.txt");
			unlink("$output_directory/$filename\_p8_gap_sites.fasta");
		}
	}
}


if (($pattern == 1) or ($pattern == 9) or ($pattern == 10)) {
	if (%percentage_variable) {
		open(my $output9,">","$output_directory/markers_by_variable_sites.txt");
		print $output9 "Sequence_matrix\tpercentage_of_variable_sites\tnumber_of_variable_sites\ttotal_length_of_matrix\n";
		foreach my $key (sort {$percentage_variable{$b} <=> $percentage_variable{$a}} keys %percentage_variable) {
			print $output9 $key."\t".$percentage_variable{$key}."\t".$number_variable{$key}."\n";
		}
		close $output9;
	}
}
if (($pattern == 5) or ($pattern == 9) or ($pattern == 10)) {
	if (%percentage_pis) {
		open(my $output10,">","$output_directory/markers_by_parsimony_informative_sites.txt");
		print $output10 "Sequence_matrix\tpercentage_of_parsimony_informative_sites\tnumber_of_parsimony_informative_sites\ttotal_length_of_matrix\n";
		foreach my $key (sort {$percentage_pis{$b} <=> $percentage_pis{$a}} keys %percentage_pis) {
			print $output10 $key."\t".$percentage_pis{$key}."\t".$number_pis{$key}."\n";
		}
		close $output10;
	}
}








##function
sub argument{
	my @options=("help|h","input|i:s","pattern|p:i","output|o:s");
	my %options;
	GetOptions(\%options,@options);
	exec ("pod2usage $0") if ((keys %options)==0 or $options{'h'} or $options{'help'});
	if(!exists $options{'input'}){
		print "***ERROR: No input directory is assigned!!!\n";
		exec ("pod2usage $0");
	}elsif(!exists $options{'pattern'}){
		print "***ERROR: No statistic pattern is assigned!!!\n";
		exec ("pod2usage $0");
	}elsif(!exists $options{'output'}){
		print "***ERROR: No output directory is assigned!!!\n";
		exec ("pod2usage $0");
	}
	return \%options;
}

sub default{
	my ($default_value,$option)=@_;
	if(exists $global_options->{$option}){
		return $global_options->{$option};
	}
	return $default_value;
}


__DATA__

=head1 NAME

    barcoding_v1.pl version 1.0

=head1 COPYRIGHT

    copyright (C) 2025 Xiao-Jian Qu

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 DESCRIPTION

    calculate and extract variable sites, invariable sites, and gap sites

=head1 SYNOPSIS

    barcoding_v1.pl -i -p -o
    Copyright (C) 2025 Xiao-Jian Qu
    Please contact <quxiaojian@sdnu.edu.cn>, if you have any bugs or questions.

    [-h -help]    help information.
    [-i -input]   required: (default: input) input directory name containing fasta alignment matrices.
    [-p -pattern] required: (default: 9) 1/2/3/4/5/6/7/8/9/10, statistic pattern.
                  1=variable sites all nucleotides, 2=variable sites 2 nucleotides,
                  3=variable sites 3 nucleotides, 4=variable sites 4 nucleotides,
                  5=parsimony informative sites, 6=non-parsimony informative sites,
                  7=invariable sites, 8=gap sites, 9=1+5, 10=1+2+3+4+5+6+7+8.
    [-o -output]  required: (default: output) output directory name.

=cut
