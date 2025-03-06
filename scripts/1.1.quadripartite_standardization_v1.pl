#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Find;
no warnings "uninitialized";
use Data::Dumper;
$|=1;

my $global_options=&argument();
my $input_directory=&default("input","input");
my $reverse_complement=&default("N","rc");
my $same=&default("N","same");
my $length=&default("1000","length");
my $output_directory=&default("output","output");


#generate output directory
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


#generate random number
my $maxLenth=14;
my @a = (0..9,'%','a'..'z','A'..'Z','-','+','_');
my $random = "%".join ('',map ($a[int rand @a],0..($maxLenth-1)))."%";

my $pattern1=".fasta";
my $pattern2=".fas";
my $pattern3=".fa";
my $pattern4=".fsa";
my @sequence_filenames;
find(\&target,$input_directory);
sub target{
	if (/$pattern1/ or /$pattern2/ or /$pattern3/ or /$pattern4/){
		push @sequence_filenames,"$File::Find::name";
	}
	return;
}


#generate first row of quadripartite_structure_coordinate.txt
open(my $output_coordinate1,">","$output_directory/quadripartite_structure_coordinate.txt");
print $output_coordinate1 "FastaFileNames\tAdjustmentForms\tLSC(JLA-JLB)\tIRb(JLB-JSB)\tSSC(JSB-JSA)\tIRa(JSA-JLA)\n";
close $output_coordinate1;


#generate screen.log file for debug
my $screenlog = $output_directory."/"."screen.log";
my $screenperf = $output_directory."/"."screen.perf";


#batch process each fasta file
my $j=0;
@sequence_filenames=sort @sequence_filenames;
while (@sequence_filenames) {
	$j++;
	my $filename_fasta=shift @sequence_filenames;
	my $target_name=$filename_fasta;#input/Amborella_trichopoda.fasta
	$target_name=substr($target_name,0,rindex($target_name,"\."));#input/Amborella_trichopoda
	$target_name=~ s/\s/_/g;
	my $filename=substr($target_name,rindex($target_name,"\/")+1);#Amborella_trichopoda
	my $target_path=$target_name;
	$target_path=substr($target_name,0,rindex($target_name,"\/"));#input
	my $target_name_random = $target_path."/".$random."_".$filename;
	my $IR_temp_random = $target_path."/".$random."_".$filename."_IR_temp";
	my $blast_IR_temp_random = $target_path."/".$random."_".$filename."_IR_temp_blast";


	#from interleaved to non-interleaved
	open(my $input_ag,"<",$filename_fasta);
	open(my $output_ag,">",$target_name_random);
	my $row_ag=<$input_ag>;
	print $output_ag $row_ag;
	while ($row_ag=<$input_ag>){
		$row_ag=~ s/\r|\n//g;
		if ($row_ag=~ /^>/) {
			print $output_ag "\n".$row_ag."\n";
		}else{
			print $output_ag $row_ag;
		}
	}
	print $output_ag "\n";
	close $input_ag;
	close $output_ag;


	#generate sequence and length_cp of fasta file, for following locate quadripartite structure LSC-IRb-SSC-IRa
	open (my $in_fasta,"<",$target_name_random);
	my @fasta;
	while (<$in_fasta>){
		$_=~ s/\r|\n//g;
		push @fasta,$_;
	}
	close $in_fasta;

	my ($header,$sequence,$length_cp,$rc_sequence);
	while (@fasta){
		my $head=shift @fasta;
		$head=~ s/(\s)+/_/g;
		$head=~ s/_$//g;
		$header=$1 if ($head=~ /^>(.+)$/);
		$header=~ s/\\/_/g;
		$header=~ s/\//_/g;
		$header=~ s/\:/_/g;
		$header=~ s/\*/_/g;
		$header=~ s/\?/_/g;
		$header=~ s/\"/_/g;
		$header=~ s/\</_/g;
		$header=~ s/\>/_/g;
		$header=~ s/\|/_/g;

		$sequence=shift @fasta;
		$sequence= uc $sequence;
		$length_cp=length $sequence;

		$rc_sequence=reverse $sequence;
		$rc_sequence=~ tr/ACGTacgt/TGCAtgca/;

		if ($reverse_complement eq "N") {
			$sequence=$sequence;
		}elsif ($reverse_complement eq "Y") {
			$sequence=$rc_sequence;
		}
	}


	#fasta_sequence_with_length_cp*2
	open (my $in_fasta_2,"<",$target_name_random);
	open (my $out_fasta_2,">",$IR_temp_random);
	while (<$in_fasta_2>) {
		$_=~ s/\r|\n//g;
		if ($_=~ /^>/) {
			print $out_fasta_2 "$_\n";
		}elsif ($_!~ /^>/) {
			$_= uc $_;
			if ($reverse_complement eq "N") {
				print $out_fasta_2 $_.$_."\n";
			}elsif ($reverse_complement eq "Y") {
				my $rev_coms=reverse $_;
				$rev_coms=~ tr/ACGTacgt/TGCAtgca/;
				$_=$rev_coms;
				print $out_fasta_2 $_.$_."\n";
			}
		}
	}
	close $in_fasta_2;
	close $out_fasta_2;
	unlink("$target_name_random");


	#blast search
	if ($same eq "N") {
		if ($osname eq "MSWin32") {
			#IRb and IRa
			system ("makeblastdb.exe -in $IR_temp_random -hash_index -dbtype nucl -blastdb_version 4 -logfile $screenlog");
			system ("blastn.exe -task blastn -query $IR_temp_random -db $IR_temp_random -outfmt 6 -perc_identity 99 -out $blast_IR_temp_random");
		}elsif ($osname eq "cygwin") {
			#IRb and IRa
			system ("makeblastdb -in $IR_temp_random -hash_index -dbtype nucl -blastdb_version 4 -logfile $screenlog");
			system ("blastn -task blastn -query $IR_temp_random -db $IR_temp_random -outfmt 6 -perc_identity 99 -out $blast_IR_temp_random");
		}elsif ($osname eq "linux") {
			#IRb and IRa
			system ("makeblastdb -in $IR_temp_random -hash_index -dbtype nucl -blastdb_version 4 -logfile $screenlog");
			system ("blastn -task blastn -query $IR_temp_random -db $IR_temp_random -outfmt 6 -perc_identity 99 -out $blast_IR_temp_random");
		}elsif ($osname eq "darwin") {
			#IRb and IRa
			system ("makeblastdb -in $IR_temp_random -hash_index -dbtype nucl -blastdb_version 4 -logfile $screenlog");
			system ("blastn -task blastn -query $IR_temp_random -db $IR_temp_random -outfmt 6 -perc_identity 99 -out $blast_IR_temp_random");
		}
	}elsif ($same eq "Y") {
		if ($osname eq "MSWin32") {
			#IRb and IRa
			system ("makeblastdb.exe -in $IR_temp_random -hash_index -dbtype nucl -blastdb_version 4 -logfile $screenlog");
			system ("blastn.exe -task blastn -query $IR_temp_random -db $IR_temp_random -outfmt 6 -perc_identity 100 -out $blast_IR_temp_random");
		}elsif ($osname eq "cygwin") {
			#IRb and IRa
			system ("makeblastdb -in $IR_temp_random -hash_index -dbtype nucl -blastdb_version 4 -logfile $screenlog");
			system ("blastn -task blastn -query $IR_temp_random -db $IR_temp_random -outfmt 6 -perc_identity 100 -out $blast_IR_temp_random");
		}elsif ($osname eq "linux") {
			#IRb and IRa
			system ("makeblastdb -in $IR_temp_random -hash_index -dbtype nucl -blastdb_version 4 -logfile $screenlog");
			system ("blastn -task blastn -query $IR_temp_random -db $IR_temp_random -outfmt 6 -perc_identity 100 -out $blast_IR_temp_random");
		}elsif ($osname eq "darwin") {
			#IRb and IRa
			system ("makeblastdb -in $IR_temp_random -hash_index -dbtype nucl -blastdb_version 4 -logfile $screenlog");
			system ("blastn -task blastn -query $IR_temp_random -db $IR_temp_random -outfmt 6 -perc_identity 100 -out $blast_IR_temp_random");
		}
	}

	unlink ("$IR_temp_random.nhd");
	unlink ("$IR_temp_random.nhi");
	unlink ("$IR_temp_random.nhr");
	unlink ("$IR_temp_random.nin");
	unlink ("$IR_temp_random.nog");
	unlink ("$IR_temp_random.nsd");
	unlink ("$IR_temp_random.nsi");
	unlink ("$IR_temp_random.nsq");
	unlink("$screenperf");


	#locate quadripartite structure LSC-IRb-SSC-IRa
	my %IR;
	open (my $input_IR,"<",$blast_IR_temp_random);
	while (<$input_IR>) {
		$_=~ s/\r|\n//g;
		my ($c1,$c2,$c3,$ir_length,$c5,$c6,$qs,$qe,$ss,$se,$c11,$c12)=split /\t/,$_;
		if (($ir_length != 2 * $length_cp) and ($ir_length != $length_cp) and ($ir_length >= $length) and ($qe <= $length_cp) and ($se <= $length_cp) and ($qs < $ss)) {
			$IR{$ir_length}=$qs."\t".$qe."\t".$ss."\t".$se;
		}
	}
	close $input_IR;
	unlink ("$blast_IR_temp_random");


	#if same=Y is not working for plastomes with length and nucleotide of IRb and IRa not be exactly the same,
	#the script will jump to same=N, equivalent to running blast search again.
	if (scalar keys(%IR) == 0) {
		if ($osname eq "MSWin32") {
			#IRb and IRa
			system ("makeblastdb.exe -in $IR_temp_random -hash_index -dbtype nucl -blastdb_version 4 -logfile $screenlog");
			system ("blastn.exe -task blastn -query $IR_temp_random -db $IR_temp_random -outfmt 6 -perc_identity 99 -out $blast_IR_temp_random");
		}elsif ($osname eq "cygwin") {
			#IRb and IRa
			system ("makeblastdb -in $IR_temp_random -hash_index -dbtype nucl -blastdb_version 4 -logfile $screenlog");
			system ("blastn -task blastn -query $IR_temp_random -db $IR_temp_random -outfmt 6 -perc_identity 99 -out $blast_IR_temp_random");
		}elsif ($osname eq "linux") {
			#IRb and IRa
			system ("makeblastdb -in $IR_temp_random -hash_index -dbtype nucl -blastdb_version 4 -logfile $screenlog");
			system ("blastn -task blastn -query $IR_temp_random -db $IR_temp_random -outfmt 6 -perc_identity 99 -out $blast_IR_temp_random");
		}elsif ($osname eq "darwin") {
			#IRb and IRa
			system ("makeblastdb -in $IR_temp_random -hash_index -dbtype nucl -blastdb_version 4 -logfile $screenlog");
			system ("blastn -task blastn -query $IR_temp_random -db $IR_temp_random -outfmt 6 -perc_identity 99 -out $blast_IR_temp_random");
		}

		unlink ("$IR_temp_random.nhd");
		unlink ("$IR_temp_random.nhi");
		unlink ("$IR_temp_random.nhr");
		unlink ("$IR_temp_random.nin");
		unlink ("$IR_temp_random.nog");
		unlink ("$IR_temp_random.nsd");
		unlink ("$IR_temp_random.nsi");
		unlink ("$IR_temp_random.nsq");
		unlink("$screenperf");

		open(my $warning,">>","$output_directory/warning.txt");
		print $warning "$filename need to be checked for following reason(s):\n";
		open (my $input_IR,"<",$blast_IR_temp_random);
		while (<$input_IR>) {
			$_=~ s/\r|\n//g;
			my ($c1,$c2,$c3,$ir_length,$mismatch,$gapopen,$qs,$qe,$ss,$se,$c11,$c12)=split /\t/,$_;
			if (($ir_length != 2 * $length_cp) and ($ir_length != $length_cp) and ($ir_length >= $length) and ($qe <= $length_cp) and ($se <= $length_cp) and ($qs < $ss)) {
				$IR{$ir_length}=$qs."\t".$qe."\t".$ss."\t".$se;

				if (($mismatch != 0) and ($gapopen == 0)) {
					if ($ss < $se) {
						print $warning "Nucleotide of IRb and IRa may not be exactly the same!\n";
					}elsif ($ss > $se) {
						print $warning "Nucleotide of IRb and IRa may not be exactly the same!\n\n";
					}
				}elsif (($mismatch == 0) and ($gapopen != 0)) {
					if ($ss < $se) {
						print $warning "Length of IRb and IRa may not be exactly the same!\n";
					}elsif ($ss > $se) {
						print $warning "Length of IRb and IRa may not be exactly the same!\n\n";
					}
				}elsif (($mismatch != 0) and ($gapopen != 0)) {
					if ($ss < $se) {
						print $warning "Nucleotide and length of IRb and IRa may not be exactly the same!\n";
					}elsif ($ss > $se) {
						print $warning "Nucleotide and length of IRb and IRa may not be exactly the same!\n\n";
					}
				}

				if ($ss < $se) {
					print $warning "IRb and IRa may be forward repeat rather than inverted repeat!\n\n";
				}
			}
		}
		close $input_IR;
		unlink ("$blast_IR_temp_random");

		if (scalar keys(%IR) == 0) {
			print $warning "Minimum allowed IR length $length may be larger than actual maximum IR length!\n\n";
		}
		close $warning;
	}
	unlink ("$IR_temp_random");


	#generate four junctions for IRb (JLB and JSB) and IRa (JSA and JLA)
	if (scalar keys(%IR) != 0) {
		my (@IR_length,@boundary);
		foreach my $key (sort {$b <=> $a} keys %IR) {
			push @IR_length,$key;
			push @boundary,$IR{$key};
		}
		my ($AA,$BB,$CC,$DD)=split /\t/,$boundary[0];
		my ($JLB,$JSB,$JLA,$JSA);
		if ($CC <= $length_cp) {
			if ($CC > $DD) {
				$JLB=$AA;
				$JSB=$BB;
				$JLA=$CC;
				$JSA=$DD;
			}elsif ($CC < $DD) {
				$JLB=$AA;
				$JSB=$BB;
				$JLA=$DD;
				$JSA=$CC;
			}
		}elsif ($CC > $length_cp) {
			$JLB=$AA;
			$JSB=$BB;
			$JLA=$CC-$length_cp;
			$JSA=$DD;
		}
		#print "$JLB\t$JSB\t$JSA\t$JLA\n";

		#if (defined $JLB) {
		#	print "     repeat_region   "."$JLB..$JSB"."\n";
		#	print "                     /note=\"inverted repeat B\""."\n";
		#	print "                     /rpt_type=\"inverted\""."\n";
		#	print "     repeat_region   "."complement($JSA..$JLA)"."\n";
		#	print "                     /note=\"inverted repeat A\""."\n";
		#	print "                     /rpt_type=\"inverted\""."\n";
		#}


		#generate pre-/post-adjustment coordinate, generate post-adjustment complete sequence,
		#and generate post-adjustment LSC, IRb, SSC, and IRa sequences
		open(my $output_sequence,">","$output_directory/$filename.fasta");
		open(my $output_sequence2,">","$output_directory/$filename\_LSC_IRb_SSC_IRa.fasta");
		open(my $output_coordinate,">>","$output_directory/quadripartite_structure_coordinate.txt");
		my ($LSC,$IRb,$SSC,$IRa);
		my ($JLB_adjusted,$JSB_adjusted,$JSA_adjusted,$JLA_adjusted);
		if ($JLA == $length_cp) {
			if ((($JLB-1)-1+1) >= (($JSA-1)-($JSB+1)+1)) {
				$LSC=substr($sequence,0,(($JLB-1)-1+1));
				$IRb=substr($sequence,($JLB-1),($JSB-$JLB+1));
				$SSC=substr($sequence,(($JSB+1)-1),(($JSA-1)-($JSB+1)+1));
				$IRa=substr($sequence,($JSA-1),($JLA-$JSA+1));
				print $output_sequence ">$header\n";
				print $output_sequence $LSC.$IRb.$SSC.$IRa."\n";

				print $output_sequence2 ">$header\_LSC\n";
				print $output_sequence2 $LSC."\n";
				print $output_sequence2 ">$header\_IRb\n";
				print $output_sequence2 $IRb."\n";
				print $output_sequence2 ">$header\_SSC\n";
				print $output_sequence2 $SSC."\n";
				print $output_sequence2 ">$header\_IRa\n";
				print $output_sequence2 $IRa."\n";

				$JLB_adjusted=$JLB;
				$JSB_adjusted=$JSB;
				$JSA_adjusted=$JSA;
				$JLA_adjusted=$JLA;

				#print $output_coordinate "$filename\tPre-adjustment\t$JLB\t$JSB\t$JSA\t$JLA\n";
				#print $output_coordinate "$filename\tPost-adjustment\t$JLB_adjusted\t$JSB_adjusted\t$JSA_adjusted\t$JLA_adjusted\n";

				my $JLB_temp=$JLB-1;
				my $JSB_temp=$JSB+1;
				my $JSA_temp=$JSA-1;
				my $JLA_temp=$JLA+1;
				my $JLB_adjusted_temp=$JLB_adjusted-1;
				my $JSB_adjusted_temp=$JSB_adjusted+1;
				my $JSA_adjusted_temp=$JSA_adjusted-1;
				my $JLA_adjusted_temp=$JLA_adjusted+1;

				print $output_coordinate "$filename\tPre-adjustment\t1-$JLB_temp\t$JLB-$JSB\t$JSB_temp-$JSA_temp\t$JSA-$JLA\n";
				print $output_coordinate "$filename\tPost-adjustment\t1-$JLB_adjusted_temp\t$JLB_adjusted-$JSB_adjusted\t$JSB_adjusted_temp-$JSA_adjusted_temp\t$JSA_adjusted-$JLA_adjusted\n";
			}elsif ((($JLB-1)-1+1) < (($JSA-1)-($JSB+1)+1)) {
				$SSC=substr($sequence,0,(($JLB-1)-1+1));
				$IRa=substr($sequence,($JLB-1),($JSB-$JLB+1));
				$LSC=substr($sequence,($JSB+1-1),(($JSA-1)-($JSB+1)+1));
				$IRb=substr($sequence,($JSA-1),($JLA-$JSA+1));
				print $output_sequence ">$header\n";
				print $output_sequence $LSC.$IRb.$SSC.$IRa."\n";

				print $output_sequence2 ">$header\_LSC\n";
				print $output_sequence2 $LSC."\n";
				print $output_sequence2 ">$header\_IRb\n";
				print $output_sequence2 $IRb."\n";
				print $output_sequence2 ">$header\_SSC\n";
				print $output_sequence2 $SSC."\n";
				print $output_sequence2 ">$header\_IRa\n";
				print $output_sequence2 $IRa."\n";

				$JLB_adjusted=(($JSA-1)-($JSB+1)+1)+1;
				$JSB_adjusted=(($JSA-1)-($JSB+1)+1)+($JLA-$JSA+1);
				$JSA_adjusted=$length_cp-($JSB-$JLB+1)+1;
				$JLA_adjusted=$length_cp;

				#print $output_coordinate "$filename\tPre-adjustment\t$JSA\t$JLA\t$JLB\t$JSB\n";
				#print $output_coordinate "$filename\tPost-adjustment\t$JLB_adjusted\t$JSB_adjusted\t$JSA_adjusted\t$JLA_adjusted\n";

				my $JLB_temp=$JLB-1;
				my $JSB_temp=$JSB+1;
				my $JSA_temp=$JSA-1;
				my $JLA_temp=$JLA+1;
				my $JLB_adjusted_temp=$JLB_adjusted-1;
				my $JSB_adjusted_temp=$JSB_adjusted+1;
				my $JSA_adjusted_temp=$JSA_adjusted-1;
				my $JLA_adjusted_temp=$JLA_adjusted+1;

				print $output_coordinate "$filename\tPre-adjustment\t$JSB_temp-$JSA_temp\t$JSA-$JLA\t1-$JLB_temp\t$JLB-$JSB\n";
				print $output_coordinate "$filename\tPost-adjustment\t1-$JLB_adjusted_temp\t$JLB_adjusted-$JSB_adjusted\t$JSB_adjusted_temp-$JSA_adjusted_temp\t$JSA_adjusted-$JLA_adjusted\n";
			}
		}elsif (($JLB-1+$length_cp) == $length_cp) {
			if ((($JSA-1)-($JSB+1)+1) >= ($length_cp-($JLA+1)+1)) {
				$IRa=substr($sequence,($JLB-1),($JSB-$JLB+1));
				$LSC=substr($sequence,(($JSB+1)-1),(($JSA-1)-($JSB+1)+1));
				$IRb=substr($sequence,($JSA-1),($JLA-$JSA+1));
				$SSC=substr($sequence,(($JLA+1)-1),($length_cp-($JLA+1)+1));
				print $output_sequence ">$header\n";
				print $output_sequence $LSC.$IRb.$SSC.$IRa."\n";

				print $output_sequence2 ">$header\_LSC\n";
				print $output_sequence2 $LSC."\n";
				print $output_sequence2 ">$header\_IRb\n";
				print $output_sequence2 $IRb."\n";
				print $output_sequence2 ">$header\_SSC\n";
				print $output_sequence2 $SSC."\n";
				print $output_sequence2 ">$header\_IRa\n";
				print $output_sequence2 $IRa."\n";

				$JLB_adjusted=(($JSA-1)-($JSB+1)+1)+1;
				$JSB_adjusted=$JLA-($JSB+1)+1;
				$JSA_adjusted=$length_cp-($JSB-$JLB+1)+1;
				$JLA_adjusted=$length_cp;

				#print $output_coordinate "$filename\tPre-adjustment\t$JSA\t$JLA\t$JLB\t$JSB\n";
				#print $output_coordinate "$filename\tPost-adjustment\t$JLB_adjusted\t$JSB_adjusted\t$JSA_adjusted\t$JLA_adjusted\n";

				my $JLB_temp=$JLB-1;
				my $JSB_temp=$JSB+1;
				my $JSA_temp=$JSA-1;
				my $JLA_temp=$JLA+1;
				my $JLB_adjusted_temp=$JLB_adjusted-1;
				my $JSB_adjusted_temp=$JSB_adjusted+1;
				my $JSA_adjusted_temp=$JSA_adjusted-1;
				my $JLA_adjusted_temp=$JLA_adjusted+1;

				print $output_coordinate "$filename\tPre-adjustment\t$JSB_temp-$JSA_temp\t$JSA-$JLA\t$JLA_temp-$length_cp\t$JLB-$JSB\n";
				print $output_coordinate "$filename\tPost-adjustment\t1-$JLB_adjusted_temp\t$JLB_adjusted-$JSB_adjusted\t$JSB_adjusted_temp-$JSA_adjusted_temp\t$JSA_adjusted-$JLA_adjusted\n";
			}elsif ((($JSA-1)-($JSB+1)+1) < ($length_cp-($JLA+1)+1)) {
				$IRb=substr($sequence,($JLB-1),($JSB-$JLB+1));
				$SSC=substr($sequence,(($JSB+1)-1),(($JSA-1)-($JSB+1)+1));
				$IRa=substr($sequence,($JSA-1),($JLA-$JSA+1));
				$LSC=substr($sequence,(($JLA+1)-1),($length_cp-($JLA+1)+1));
				print $output_sequence ">$header\n";
				print $output_sequence $LSC.$IRb.$SSC.$IRa."\n";

				print $output_sequence2 ">$header\_LSC\n";
				print $output_sequence2 $LSC."\n";
				print $output_sequence2 ">$header\_IRb\n";
				print $output_sequence2 $IRb."\n";
				print $output_sequence2 ">$header\_SSC\n";
				print $output_sequence2 $SSC."\n";
				print $output_sequence2 ">$header\_IRa\n";
				print $output_sequence2 $IRa."\n";

				$JLB_adjusted=($length_cp-($JLA+1)+1)+1;
				$JSB_adjusted=($length_cp-($JLA+1)+1)+($JSB-$JLB+1);
				$JSA_adjusted=$length_cp-($JLA-$JSA+1)+1;
				$JLA_adjusted=$length_cp;

				#print $output_coordinate "$filename\tPre-adjustment\t$JLB\t$JSB\t$JSA\t$JLA\n";
				#print $output_coordinate "$filename\tPost-adjustment\t$JLB_adjusted\t$JSB_adjusted\t$JSA_adjusted\t$JLA_adjusted\n";

				my $JLB_temp=$JLB-1;
				my $JSB_temp=$JSB+1;
				my $JSA_temp=$JSA-1;
				my $JLA_temp=$JLA+1;
				my $JLB_adjusted_temp=$JLB_adjusted-1;
				my $JSB_adjusted_temp=$JSB_adjusted+1;
				my $JSA_adjusted_temp=$JSA_adjusted-1;
				my $JLA_adjusted_temp=$JLA_adjusted+1;

				print $output_coordinate "$filename\tPre-adjustment\t$JLA_temp-$length_cp\t$JLB-$JSB\t$JSB_temp-$JSA_temp\t$JSA-$JLA\n";
				print $output_coordinate "$filename\tPost-adjustment\t1-$JLB_adjusted_temp\t$JLB_adjusted-$JSB_adjusted\t$JSB_adjusted_temp-$JSA_adjusted_temp\t$JSA_adjusted-$JLA_adjusted\n";
			}
		}elsif (($JLA != $length_cp) and (($JLB-1+$length_cp) != $length_cp)) {
			if (($JLA < $JLB) and ($JLB < $JSB) and ($JSB < $JSA)) {
				if ((($JLB-1)-$JLA+1) >= (($JSA-1)-$JSB+1)) {
					$LSC=substr($sequence,($JLA+1-1),(($JLB-1)-($JLA+1)+1));
					$IRb=substr($sequence,($JLB-1),($JSB-$JLB+1));
					$SSC=substr($sequence,($JSB+1-1),(($JSA-1)-($JSB+1)+1));
					$IRa=substr($sequence,($JSA-1),($length_cp-$JSA+1)).substr($sequence,0,($JLA-1+1));
					print $output_sequence ">$header\n";
					print $output_sequence $LSC.$IRb.$SSC.$IRa."\n";

					print $output_sequence2 ">$header\_LSC\n";
					print $output_sequence2 $LSC."\n";
					print $output_sequence2 ">$header\_IRb\n";
					print $output_sequence2 $IRb."\n";
					print $output_sequence2 ">$header\_SSC\n";
					print $output_sequence2 $SSC."\n";
					print $output_sequence2 ">$header\_IRa\n";
					print $output_sequence2 $IRa."\n";

					$JLB_adjusted=(($JLB-1)-($JLA+1)+1)+1;
					$JSB_adjusted=$JSB-($JLA+1)+1;
					$JSA_adjusted=(($JSA-1)-($JLA+1)+1)+1;
					$JLA_adjusted=$length_cp;

					#print $output_coordinate "$filename\tPre-adjustment\t$JLB\t$JSB\t$JSA\t$JLA\n";
					#print $output_coordinate "$filename\tPost-adjustment\t$JLB_adjusted\t$JSB_adjusted\t$JSA_adjusted\t$JLA_adjusted\n";

					my $JLB_temp=$JLB-1;
					my $JSB_temp=$JSB+1;
					my $JSA_temp=$JSA-1;
					my $JLA_temp=$JLA+1;
					my $JLB_adjusted_temp=$JLB_adjusted-1;
					my $JSB_adjusted_temp=$JSB_adjusted+1;
					my $JSA_adjusted_temp=$JSA_adjusted-1;
					my $JLA_adjusted_temp=$JLA_adjusted+1;

					print $output_coordinate "$filename\tPre-adjustment\t$JLA_temp-$JLB_temp\t$JLB-$JSB\t$JSB_temp-$JSA_temp\t$JSA-$JLA\n";
					print $output_coordinate "$filename\tPost-adjustment\t1-$JLB_adjusted_temp\t$JLB_adjusted-$JSB_adjusted\t$JSB_adjusted_temp-$JSA_adjusted_temp\t$JSA_adjusted-$JLA_adjusted\n";
				}elsif ((($JLB-1)-$JLA+1) < (($JSA-1)-$JSB+1)) {
					$SSC=substr($sequence,($JLA+1-1),(($JLB-1)-($JLA+1)+1));
					$IRa=substr($sequence,($JLB-1),($JSB-$JLB+1));
					$LSC=substr($sequence,($JSB+1-1),(($JSA-1)-($JSB+1)+1));
					$IRb=substr($sequence,($JSA-1),($length_cp-$JSA+1)).substr($sequence,0,($JLA-1+1));
					print $output_sequence ">$header\n";
					print $output_sequence $LSC.$IRb.$SSC.$IRa."\n";

					print $output_sequence2 ">$header\_LSC\n";
					print $output_sequence2 $LSC."\n";
					print $output_sequence2 ">$header\_IRb\n";
					print $output_sequence2 $IRb."\n";
					print $output_sequence2 ">$header\_SSC\n";
					print $output_sequence2 $SSC."\n";
					print $output_sequence2 ">$header\_IRa\n";
					print $output_sequence2 $IRa."\n";

					$JLB_adjusted=(($JSA-1)-($JSB+1)+1)+1;
					$JSB_adjusted=(($JSA-1)-($JSB+1)+1)+(($length_cp-$JSA+1)+($JLA-1+1));
					$JSA_adjusted=$length_cp-($JSB-$JLB+1)+1;
					$JLA_adjusted=$length_cp;

					#print $output_coordinate "$filename\tPre-adjustment\t$JSA\t$JLA\t$JLB\t$JSB\n";
					#print $output_coordinate "$filename\tPost-adjustment\t$JLB_adjusted\t$JSB_adjusted\t$JSA_adjusted\t$JLA_adjusted\n";

					my $JLB_temp=$JLB-1;
					my $JSB_temp=$JSB+1;
					my $JSA_temp=$JSA-1;
					my $JLA_temp=$JLA+1;
					my $JLB_adjusted_temp=$JLB_adjusted-1;
					my $JSB_adjusted_temp=$JSB_adjusted+1;
					my $JSA_adjusted_temp=$JSA_adjusted-1;
					my $JLA_adjusted_temp=$JLA_adjusted+1;

					print $output_coordinate "$filename\tPre-adjustment\t$JSB_temp-$JSA_temp\t$JSA-$JLA\t$JLA_temp-$JLB_temp\t$JLB-$JSB\n";
					print $output_coordinate "$filename\tPost-adjustment\t1-$JLB_adjusted_temp\t$JLB_adjusted-$JSB_adjusted\t$JSB_adjusted_temp-$JSA_adjusted_temp\t$JSA_adjusted-$JLA_adjusted\n";
				}
			}elsif (($JLB < $JSB) and ($JSB < $JSA) and ($JSA < $JLA)) {
				if ((($JSA-1)-($JSB+1)+1) >= (($length_cp-($JLA+1)+1)+(($JLB-1)-1+1))) {
					$IRa=substr($sequence,($JLB-1),($JSB-$JLB+1));
					$LSC=substr($sequence,($JSB+1-1),(($JSA-1)-($JSB+1)+1));
					$IRb=substr($sequence,($JSA-1),($JLA-$JSA+1));
					$SSC=substr($sequence,($JLA+1-1),($length_cp-($JLA+1)+1)).substr($sequence,0,(($JLB-1)-1+1));
					print $output_sequence ">$header\n";
					print $output_sequence $LSC.$IRb.$SSC.$IRa."\n";

					print $output_sequence2 ">$header\_LSC\n";
					print $output_sequence2 $LSC."\n";
					print $output_sequence2 ">$header\_IRb\n";
					print $output_sequence2 $IRb."\n";
					print $output_sequence2 ">$header\_SSC\n";
					print $output_sequence2 $SSC."\n";
					print $output_sequence2 ">$header\_IRa\n";
					print $output_sequence2 $IRa."\n";

					$JLB_adjusted=(($JSA-1)-($JSB+1)+1)+1;
					$JSB_adjusted=(($JSA-1)-($JSB+1)+1)+($JLA-$JSA+1);
					$JSA_adjusted=$length_cp-($JSB-$JLB+1)+1;
					$JLA_adjusted=$length_cp;

					#print $output_coordinate "$filename\tPre-adjustment\t$JSA\t$JLA\t$JLB\t$JSB\n";
					#print $output_coordinate "$filename\tPost-adjustment\t$JLB_adjusted\t$JSB_adjusted\t$JSA_adjusted\t$JLA_adjusted\n";

					my $JLB_temp=$JLB-1;
					my $JSB_temp=$JSB+1;
					my $JSA_temp=$JSA-1;
					my $JLA_temp=$JLA+1;
					my $JLB_adjusted_temp=$JLB_adjusted-1;
					my $JSB_adjusted_temp=$JSB_adjusted+1;
					my $JSA_adjusted_temp=$JSA_adjusted-1;
					my $JLA_adjusted_temp=$JLA_adjusted+1;

					print $output_coordinate "$filename\tPre-adjustment\t$JSB_temp-$JSA_temp\t$JSA-$JLA\t$JLA_temp-$JLB_temp\t$JLB-$JSB\n";
					print $output_coordinate "$filename\tPost-adjustment\t1-$JLB_adjusted_temp\t$JLB_adjusted-$JSB_adjusted\t$JSB_adjusted_temp-$JSA_adjusted_temp\t$JSA_adjusted-$JLA_adjusted\n";
				}elsif ((($JSA-1)-($JSB+1)+1) < (($length_cp-($JLA+1)+1)+(($JLB-1)-1+1))) {
					$IRb=substr($sequence,($JLB-1),($JSB-$JLB+1));
					$SSC=substr($sequence,($JSB+1-1),(($JSA-1)-($JSB+1)+1));
					$IRa=substr($sequence,($JSA-1),($JLA-$JSA+1));
					$LSC=substr($sequence,($JLA+1-1),($length_cp-($JLA+1)+1)).substr($sequence,0,(($JLB-1)-1+1));
					print $output_sequence ">$header\n";
					print $output_sequence $LSC.$IRb.$SSC.$IRa."\n";

					print $output_sequence2 ">$header\_LSC\n";
					print $output_sequence2 $LSC."\n";
					print $output_sequence2 ">$header\_IRb\n";
					print $output_sequence2 $IRb."\n";
					print $output_sequence2 ">$header\_SSC\n";
					print $output_sequence2 $SSC."\n";
					print $output_sequence2 ">$header\_IRa\n";
					print $output_sequence2 $IRa."\n";

					$JLB_adjusted=($length_cp-($JLA+1)+1)+(($JLB-1)-1+1)+1;
					$JSB_adjusted=(($length_cp-($JLA+1)+1)+(($JLB-1)-1+1))+($JSB-$JLB+1);
					$JSA_adjusted=$length_cp-($JLA-$JSA+1)+1;
					$JLA_adjusted=$length_cp;

					#print $output_coordinate "$filename\tPre-adjustment\t$JLB\t$JSB\t$JSA\t$JLA\n";
					#print $output_coordinate "$filename\tPost-adjustment\t$JLB_adjusted\t$JSB_adjusted\t$JSA_adjusted\t$JLA_adjusted\n";

					my $JLB_temp=$JLB-1;
					my $JSB_temp=$JSB+1;
					my $JSA_temp=$JSA-1;
					my $JLA_temp=$JLA+1;
					my $JLB_adjusted_temp=$JLB_adjusted-1;
					my $JSB_adjusted_temp=$JSB_adjusted+1;
					my $JSA_adjusted_temp=$JSA_adjusted-1;
					my $JLA_adjusted_temp=$JLA_adjusted+1;

					print $output_coordinate "$filename\tPre-adjustment\t$JLA_temp-$JLB_temp\t$JLB-$JSB\t$JSB_temp-$JSA_temp\t$JSA-$JLA\n";
					print $output_coordinate "$filename\tPost-adjustment\t1-$JLB_adjusted_temp\t$JLB_adjusted-$JSB_adjusted\t$JSB_adjusted_temp-$JSA_adjusted_temp\t$JSA_adjusted-$JLA_adjusted\n";
				}
			}
		}
		close $output_sequence;
		close $output_sequence2;
		close $output_coordinate;
	}
}




sub argument{
	my @options=("help|h","input|i:s","rc|r:s","same|s:s","length|l:i","output|o:s");
	my %options;
	GetOptions(\%options,@options);
	exec ("pod2usage $0") if ((keys %options)==0 or $options{'h'} or $options{'help'});
	if(!exists $options{'input'}){
		print "***ERROR: Input directory is not assigned!!!\n";
		exec ("pod2usage $0");
	}
	if(!exists $options{'rc'}){
		print "***ERROR: Reverse-complement or not is not assigned!!!\n";
		exec ("pod2usage $0");
	}
	if(!exists $options{'same'}){
		print "***ERROR: Same or not is not assigned!!!\n";
		exec ("pod2usage $0");
	}
	if(!exists $options{'length'}){
		print "***ERROR: Ninimum allowed IR length is not assigned!!!\n";
		exec ("pod2usage $0");
	}
	if(!exists $options{'output'}){
		print "***ERROR: Output directory is not assigned!!!\n";
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

    quadripartite_standardization_v1.pl Standardization of Quadripartite Structure version 1.0

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

    Standardization of Quadripartite Structure version 1.0

=head1 SYNOPSIS

    quadripartite_standardization_v1.pl -i -r -s -l -o
    Standardization of Quadripartite Structure
    Copyright (C) 2025 Xiao-Jian Qu
    Please contact <quxiaojian@sdnu.edu.cn>, if you have any bugs or questions.

    [-h -help]         help information.
    [-i -input]        required: (default: input) input directory name containing FASTA files.
    [-r -rc]           required: (default: N) if N, entire sequence does not require reverse-complement;
                                 if Y, entire sequence requires reverse-complement.
    [-s -same]         required: (default: N) if N, length and nucleotide of IRb and IRa may not be
                                 exactly the same, allowing for a certain degree of sequence difference;
                                 if Y, length and nucleotide of IRb and IRa must be exactly the same.
    [-l -length]       required: (default: 1000) minimum allowed inverted-repeat (IR) length.
    [-o -output]       required: (default: output) output directory name.

=cut
