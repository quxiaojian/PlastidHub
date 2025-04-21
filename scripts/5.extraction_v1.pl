#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Find;
use Data::Dumper;
$|=1;

my $global_options=&argument();
my $input_directory=&default("input","input");
my $pattern=&default("3","pattern");
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

my $maxLenth=14;
my @a = (0..9,'%','a'..'z','A'..'Z','-','+','_');
my $random = "%".join ('',map ($a[int rand @a],0..($maxLenth-1)))."%";


my $pattern1=".gb";
my $pattern2=".gbk";
my $pattern3=".gbf";
my (@filenames1,@filenames2,@filenames3,@filenames4);
find(\&target,$input_directory);
sub target{
	if (/$pattern1/ or /$pattern2/ or /$pattern3/){
		push @filenames1,"$File::Find::name";
	}
	return;
}
@filenames2=@filenames1;
@filenames3=@filenames1;
@filenames4=@filenames1;


#PCGs and regions among PCGs
#genes and regions among genes
#linked CDSs and RNAs and regions among linked CDSs and RNAs
#unlinked CDSs and RNAs and regions among unlinked CDSs and RNAs

#extraction pattern: 1, 2, 3, 4, 5=1+2, 6=3+4, 7=1+2+3+4




####extract PCGs with introns and regions among PCGs
if (($pattern == 1) or ($pattern == 5) or ($pattern == 7)) {
	while (@filenames1) {
		my $filename_gb=shift @filenames1;
		my $target_name=substr($filename_gb,0,rindex($filename_gb,"\."));
		$target_name=~ s/\s/_/g;
		my $filename=substr($target_name,rindex($target_name,"\/")+1);
		#my $target_path=$target_name;
		#$target_path=~ s/$filename//g;
		my $target_name_temp_random1 = $input_directory."/".$random."_".$filename."_temp1";
		my $target_name_temp_random2 = $input_directory."/".$random."_".$filename."_temp2";

		open(my $in_gb1,"<",$filename_gb);
		open(my $out_gb1,">",$target_name_temp_random1);
		while (<$in_gb1>){
			$_=~ s/\r\n/\n/g;
			if ($_=~ /\),\n/){
				$_=~ s/\),\n/\),/g;
			}elsif($_=~ /,\n/){
				$_=~ s/,\n/,/g;
			}
			print $out_gb1 $_;
		}
		close $in_gb1;
		close $out_gb1;

		open(my $in_gb2,"<",$target_name_temp_random1);
		open(my $out_gb2,">",$target_name_temp_random2);
		while (<$in_gb2>){
			$_=~ s/,\s+/,/g;
			print $out_gb2 $_;
		}
		close $in_gb2;
		close $out_gb2;


		#generate_bed_file
		my (@row_array,@start_end,$species_name,$length,$element,@genearray,@output1);
		my $mark=0;
		my $tick=0;
		open (my $in_gb3,"<",$target_name_temp_random2);
		while (<$in_gb3>){
			$_=~ s/\r|\n//g;
			@row_array=split /\s+/,$_;
			my $i=0;
			if (/^LOCUS/i){
				$species_name=$filename;
				$length=$row_array[2];
			}elsif(/ {5}gene {12}/){
				if ($row_array[2]=~ /^\d+..\d+$/){
					$row_array[2]="\+\t$row_array[2]\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~/^complement\((\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
					my @array_start_end=split /\t/,$row_array[2];
					my $first_start_end=$row_array[2];
					my $second_start_end=$row_array[2];
					if ($array_start_end[1] < $array_start_end[2]) {
						$row_array[2]=$row_array[2];
					}elsif ($array_start_end[1] > $array_start_end[2]) {
						$first_start_end=~ s/$array_start_end[2]/$length/g;
						$second_start_end=~ s/$array_start_end[1]/1/g;
						push @start_end,$first_start_end;
						push @start_end,$second_start_end;
						$tick=3;
					}
				}elsif($row_array[2]=~ /^join\((\d+..\d+),(\d+..\d+)\)$/) {
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\((\d+..\d+),(\d+..\d+),(\d+..\d+)\)$/) {
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),(\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\((\d+..\d+),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="+\t$1\t$row_array[1]\t-\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^complement\(join\((\d+..\d+),(\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^complement\(join\((\d+..\d+),(\d+..\d+),(\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\((\d+..\d+),(\d+..\d+)\)$/){
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\((\d+..\d+),(\d+..\d+),(\d+..\d+)\)$/){
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),(\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),(\d+..\d+),(\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),(\d+..\d+),(\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\((\d+..\d+),(\d+..\d+),complement\((\d+..\d+)\)\)$/) {
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^<(\d+..\d+)$/){# positive split no-intron gene
					$row_array[2]="\+\t$1\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
					$tick=2;
				}elsif($row_array[2]=~ /^(\d+..>\d+)$/){# positive split no-intron gene
					$row_array[2]="\+\t$1\t$row_array[1]";
					$row_array[2]=~ s/\..>/\t/g;
					$tick=1;
				}elsif($row_array[2]=~/^complement\(<(\d+..\d+)\)$/){# negative split no-intron gene
					$row_array[2]="-\t$1\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
					$tick=1;
				}elsif($row_array[2]=~/^complement\((\d+..>\d+)\)$/){# negative split no-intron gene
					$row_array[2]="-\t$1\t$row_array[1]";
					$row_array[2]=~ s/\..>/\t/g;
					$tick=2;
				}
				$element=$row_array[2];
				$mark=1;
			}elsif(/^\s+\/gene="(.*)"/ and $mark == 1){
				if ($tick==0) {
					$element=$species_name.":".$1.":".$element;
					push @genearray,$element if ($element!~ /trn/ and $element!~ /rrn/);
					$element=();
					$mark=0;
				}elsif ($tick==1) {
					$element=$species_name.":"."1_".$1.":".$element;#1_trnH-GUG
					push @genearray,$element if ($element!~ /trn/ and $element!~ /rrn/);
					$element=();
					$mark=0;
					$tick=0;
				}elsif ($tick==2) {
					$element=$species_name.":"."2_".$1.":".$element;#2_trnH-GUG
					push @genearray,$element if ($element!~ /trn/ and $element!~ /rrn/);
					$element=();
					$mark=0;
					$tick=0;
				}elsif ($tick==3) {
					push @genearray,$species_name.":"."1_".$1.":".$start_end[1];#1_trnH-GUG
					push @genearray,$species_name.":"."2_".$1.":".$start_end[0];#2_trnH-GUG
					@start_end=();
					$mark=0;
					$tick=0;
				}
			}
		}
		close $in_gb3;
		foreach (@genearray){
			my @array=split /:/,$_;
			push @output1,"$array[0]\t$array[1]\t$array[2]\n";
		}
		@row_array=();
		@genearray=();


		#my (%SP0,%GENE0,%STRAND0,%START0,%END0,%TYPE0,%STRAND_0,%START_0,%END_0,%TYPE_0);
		#my $cnt0=0;
		#foreach (@genearray){
		#	$_=~ s/\r|\n//g;
		#	$cnt0++;
		#	my @array=split /:/,$_;
		#	my $position;
		#	$position="$array[0]\t$array[1]\t$array[2]";
		#	($SP0{$cnt0},$GENE0{$cnt0},$STRAND0{$cnt0},$START0{$cnt0},$END0{$cnt0},$TYPE0{$cnt0},$STRAND_0{$cnt0},$START_0{$cnt0},$END_0{$cnt0},$TYPE_0{$cnt0})=(split /\s+/,$position)[0,1,2,3,4,5,6,7,8,9];
		#	if (defined $STRAND1{$cnt0} eq "") {
		#	}elsif (defined $STRAND1{$cnt0} ne "") {
		#	}
		#}
		#@row_array=();
		#@genearray=();


		#put_fasta_sequence_in_array
		my $flag=0;
		my @sequence;
		my (@fas1,@fas2);
		open(my $in_gb4,"<",$target_name_temp_random2);
		while (<$in_gb4>){
			if ($_=~ /ORIGIN/){
				$flag=1;
			}
			if ($_=~ /\/\//){
				$flag=2;
			}
			if ($flag==1){
				next if ($_=~ /ORIGIN/);
				push @sequence,$_;
			}
		}
		close $in_gb4;
		foreach (@sequence){
			$_=~ s/\r|\n//g;
			$_=~ s/\s*//g;
			$_=~ s/\d+//g;
			push @fas1,$_;
		}
		my $fas1=join "",@fas1;
		my (@fasta1,@fasta2);
		push @fasta1,$species_name,$fas1;
		@fasta2=@fasta1;

		unlink ("$target_name_temp_random1");
		unlink ("$target_name_temp_random2");


		#edit_bed_file
		my (%SP1,%GENE1,%STRAND1,%START1,%END1,%TYPE1,%STRAND2,%START2,%END2,%TYPE2,%STRAND3,%START3,%END3,%TYPE3,@output2);
		my $cnt1=0;
		foreach (@output1) {
			$_=~ s/\r|\n//g;
			$cnt1++;
			($SP1{$cnt1},$GENE1{$cnt1},$STRAND1{$cnt1},$START1{$cnt1},$END1{$cnt1},$TYPE1{$cnt1},$STRAND2{$cnt1},$START2{$cnt1},$END2{$cnt1},$TYPE2{$cnt1})=(split /\s+/,$_)[0,1,2,3,4,5,6,7,8,9];
		}
		foreach (1..$cnt1) {
			if (defined $STRAND2{$_} eq "") {
				push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
			}elsif (defined $STRAND2{$_} ne "") {
				if (($GENE1{$_} ne "rps12")) {
					if (($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START1{$_} < $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START1{$_} > $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START1{$_} < $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START1{$_} > $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}
				}elsif (($GENE1{$_} eq "rps12")) {
					if (($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START1{$_} < $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START1{$_} > $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START1{$_} < $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START1{$_} > $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});

					}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "+") and ($START1{$_} < $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "+") and ($START1{$_} > $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "-") and ($START1{$_} < $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "-") and ($START1{$_} > $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}
				}
			}
		}
		@output1=();


		#Rosa_roxburghii        rps12   -       99457   100247  gene    -       71096   71209   gene
		#Rosa_roxburghii        rps12   -       71096   71209   gene    +       142355  143145  gene

		#Rosa_roxburghii        rps12   +       99457   100247  gene    +       71096   71209   gene
		#Rosa_roxburghii        rps12   +       71096   71209   gene    -       142355  143145  gene

		#Rosa_roxburghii        rps12   -       71096   71209   gene
		#Rosa_roxburghii        rps12   -       99457   100250  gene
		#Rosa_roxburghii        rps12   +       142352  143145  gene

		#Rosa_roxburghii        rps12   +       71096   71209   gene
		#Rosa_roxburghii        rps12   +       99457   100250  gene
		#Rosa_roxburghii        rps12   -       142352  143145  gene


		#sort_bed_file
		my $col=3;
		my (%sort,@output3);
		foreach (@output2){
			my @row=split /\t/,$_;
			$sort{$_}=$row[$col];
		}
		foreach (sort {$sort{$a} <=> $sort{$b}} keys %sort){
			push @output3,"$_\n";
		}
		@output2=();


		#output_bed_file
		open (my $out_bed,">","$output_directory/$filename\_p1_PCG.bed");
		foreach (@output3){
			my @row=split /\t/,$_;
			my $abs=abs($row[3]-$row[4]);
			if ($_!~ /rps12/) {
				print $out_bed $_;
			}
			if (($_=~ /rps12/) and ($abs < 200)) {
				$_=~ s/rps12/rps12-1/;
				print $out_bed $_;
			}
			if (($_=~ /rps12/) and ($abs > 200)) {
				$_=~ s/rps12/rps12-2/;
				print $out_bed $_;
			}
		}
		close $out_bed;


		#extract_gene
		my (%SP2,%GENE2,%STRAND4,%START4,%END4,%TYPE4,$seq1);
		my $cnt2=0;
		open(my $out_coding,">","$output_directory/$filename\_p1_PCG.fasta");
		while (@fasta1){
			my $header=shift @fasta1;
			$seq1=shift @fasta1;
		}
		my ($rps12_header,@rps12_1_sequence,@rps12_2_sequence,@rps12_3_sequence);
		my ($left_header,$left_str,$right_header,$right_str);
		foreach (@output3){
			$_=~ s/\r|\n//g;
			$cnt2++;
			($SP2{$cnt2},$GENE2{$cnt2},$STRAND4{$cnt2},$START4{$cnt2},$END4{$cnt2},$TYPE4{$cnt2})=(split /\s+/,$_)[0,1,2,3,4,5];
			my $str=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
			#if ($STRAND4{$cnt2} eq "-") {
			#	my $rev_com=reverse $str;
			#	$rev_com=~ tr/ACGTacgt/TGCAtgca/;
			#	print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$rev_com."\n";
			#}elsif($STRAND4{$cnt2} eq "+"){
			#	print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$str."\n";
			#}
			if ((($STRAND4{$cnt2} eq "-") and ($GENE2{$cnt2} ne "rps12-1")) and (($STRAND4{$cnt2} eq "-") and ($GENE2{$cnt2} ne "rps12-2"))) {
				my $rev_com=reverse $str;
				$rev_com=~ tr/ACGTacgt/TGCAtgca/;
				if (($GENE2{$cnt2} !~ /1_/) and ($GENE2{$cnt2} !~ /2_/)) {
					print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$rev_com."\n";
				}elsif ($GENE2{$cnt2} =~ /1_(.*)/) {#1_infB
					$left_header=$1;
					$left_str=$rev_com;
				}elsif ($GENE2{$cnt2} =~ /2_(.*)/) {#2_infB
					$right_header=$1;
					$right_str=$rev_com;
				}
			}elsif (($STRAND4{$cnt2} eq "-") and ($GENE2{$cnt2} eq "rps12-1")) {
				my $rev_com=reverse $str;
				$rev_com=~ tr/ACGTacgt/TGCAtgca/;
				push @rps12_1_sequence,$rev_com;
			}elsif (($STRAND4{$cnt2} eq "-") and ($GENE2{$cnt2} eq "rps12-2")) {
				my $rev_com=reverse $str;
				$rev_com=~ tr/ACGTacgt/TGCAtgca/;
				push @rps12_2_sequence,$rev_com;
			}elsif(($STRAND4{$cnt2} eq "+") and ($GENE2{$cnt2} ne "rps12-2")){
				if (($GENE2{$cnt2} !~ /1_/) and ($GENE2{$cnt2} !~ /2_/)) {
					print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$str."\n";
				}elsif ($GENE2{$cnt2} =~ /1_(.*)/) {#1_trnH-GUG
					$left_header=$1;
					$left_str=$str;
				}elsif ($GENE2{$cnt2} =~ /2_(.*)/) {#2_trnH-GUG
					$right_header=$1;
					$right_str=$str;
				}
			}elsif (($STRAND4{$cnt2} eq "+") and ($GENE2{$cnt2} eq "rps12-2")) {
				push @rps12_3_sequence,$str;
			}
			if ((defined $left_header ne "") and (defined $right_header ne "") and ($left_header eq $right_header)) {
				print $out_coding ">".$left_header."_".$SP2{$cnt2}."\n".$left_str.$right_str."\n";
			}
			$rps12_header=">rps12_".$SP2{$cnt2};
		}

		if ((@rps12_1_sequence != 0) and (@rps12_2_sequence != 0) and (@rps12_3_sequence != 0)) {
			my $rps12_1_sequence=shift @rps12_1_sequence;
			my $rps12_2_sequence=shift @rps12_2_sequence;
			my $rps12_3_sequence=shift @rps12_3_sequence;
			print $out_coding $rps12_header."\n".$rps12_1_sequence.$rps12_2_sequence."\n";
			print $out_coding $rps12_header."\n".$rps12_1_sequence.$rps12_3_sequence."\n";
		}elsif ((@rps12_1_sequence != 0) and (@rps12_2_sequence != 0) and (@rps12_3_sequence == 0)) {
			my $rps12_1_sequence=shift @rps12_1_sequence;
			my $rps12_2_sequence=shift @rps12_2_sequence;
			print $out_coding $rps12_header."\n".$rps12_1_sequence.$rps12_2_sequence."\n";
			print $out_coding $rps12_header."\n".$rps12_1_sequence.$rps12_2_sequence."\n";
		}elsif ((@rps12_1_sequence != 0) and (@rps12_2_sequence == 0) and (@rps12_3_sequence != 0)) {
			my $rps12_1_sequence=shift @rps12_1_sequence;
			my $rps12_3_sequence=shift @rps12_3_sequence;
			print $out_coding $rps12_header."\n".$rps12_1_sequence.$rps12_3_sequence."\n";
			print $out_coding $rps12_header."\n".$rps12_1_sequence.$rps12_3_sequence."\n";
		}elsif ((@rps12_1_sequence == 0) and (@rps12_2_sequence != 0) and (@rps12_3_sequence != 0)) {
			my $rps12_2_sequence=shift @rps12_2_sequence;
			my $rps12_3_sequence=shift @rps12_3_sequence;
			print $out_coding $rps12_header."\n".$rps12_2_sequence."\n";
			print $out_coding $rps12_header."\n".$rps12_3_sequence."\n";
		}elsif ((@rps12_1_sequence == 0) and (@rps12_2_sequence == 0) and (@rps12_3_sequence != 0)) {
			my $rps12_3_sequence=shift @rps12_3_sequence;
			print $out_coding $rps12_header."\n".$rps12_3_sequence."\n";
			print $out_coding $rps12_header."\n".$rps12_3_sequence."\n";
		}elsif ((@rps12_1_sequence == 0) and (@rps12_2_sequence != 0) and (@rps12_3_sequence == 0)) {
			my $rps12_2_sequence=shift @rps12_2_sequence;
			print $out_coding $rps12_header."\n".$rps12_2_sequence."\n";
			print $out_coding $rps12_header."\n".$rps12_2_sequence."\n";
		}elsif ((@rps12_1_sequence != 0) and (@rps12_2_sequence == 0) and (@rps12_3_sequence == 0)) {
			my $rps12_1_sequence=shift @rps12_1_sequence;
			print $out_coding $rps12_header."\n".$rps12_1_sequence."\n";
			print $out_coding $rps12_header."\n".$rps12_1_sequence."\n";
		}elsif ((@rps12_1_sequence == 0) and (@rps12_2_sequence == 0) and (@rps12_3_sequence == 0)) {
		}
		close $out_coding;


		#generate_IGS_ranges
		my (%SP3,%GENE3,%STRAND5,%START5,%END5,%TYPE5,$last0,$last1,$last2,@output4);
		my $cnt3=0;
		foreach (@output3){
			$_=~ s/\r|\n//g;
			$cnt3++;
			($SP3{$cnt3},$GENE3{$cnt3},$STRAND5{$cnt3},$START5{$cnt3},$END5{$cnt3},$TYPE5{$cnt3})=(split /\s+/,$_)[0,1,2,3,4,5];
		}
		foreach (keys %SP3){
			if ($_==1 and $START5{$_}!=1){
				unshift @output4,$SP3{$_}."\t"."start".'-'.$GENE3{$_}."\t"."?"."/".$STRAND5{$_}."\t"."1"."\t".($START5{$_}-1)."\t"."?"."/".$TYPE5{$_}."\n";
			}
		}
		foreach (1..($cnt3-1)) {
			$last0=$_-1;
			$last1=$_+1;
			$last2=$_+2;
			next if ((($END5{$_}+1) > ($START5{$last1}-1)) and (($END5{$_}+1) < ($END5{$last1}-1)));
			next if (($_ > 1) and (($END5{$_}+1) < ($END5{$last0}-1)) and (($END5{$_}+1) < ($START5{$last2}-1)));
			if ((($END5{$_}+1) >= ($START5{$last1}-1)) and (($END5{$_}+1) >= ($END5{$last1}-1))){
				push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'.$GENE3{$last2}."\t".$STRAND5{$_}."/".$STRAND5{$last2}."\t".($END5{$_}+1)."\t".($START5{$last2}-1)."\t".$TYPE5{$_}."/".$TYPE5{$last2}."\n";
			}else{
				push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'.$GENE3{$last1}."\t".$STRAND5{$_}."/".$STRAND5{$last1}."\t".($END5{$_}+1)."\t".($START5{$last1}-1)."\t".$TYPE5{$_}."/".$TYPE5{$last1}."\n";
			}
		}
		foreach (keys %SP3){
			if ($_==$cnt3 and $END5{$_}!=$length){
				push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'."end"."\t".$STRAND5{$_}."/"."?"."\t".($END5{$_}+1)."\t".$length."\t".$TYPE5{$_}."/"."?"."\n";
			}
		}
		@output3=();
		#print @output4;


		#extract_IGS
		my (%SP4,%GENE4,%STRAND6,%START6,%END6,%TYPE6,$seq2);
		my $cnt4=0;
		open(my $out_noncoding,">","$output_directory/$filename\_p1_regions_among_PCG.fasta");
		while (@fasta2){
			my $header=shift @fasta2;
			$seq2=shift @fasta2;
		}
		my ($left_header2,$left_str2,$right_header2,$right_str2);
		foreach (@output4){
			$_=~ s/\r|\n//g;
			$cnt4++;
			($SP4{$cnt4},$GENE4{$cnt4},$STRAND6{$cnt4},$START6{$cnt4},$END6{$cnt4},$TYPE6{$cnt4})=(split /\s+/,$_)[0,1,2,3,4,5];
			my $str=substr($seq2,($START6{$cnt4}-1),($END6{$cnt4}-$START6{$cnt4}+1));
			if (($GENE4{$cnt4} !~ /1_/) and ($GENE4{$cnt4} !~ /2_/) and ($GENE4{$cnt4} !~ /start-/) and ($GENE4{$cnt4} !~ /-end/)) {
				print $out_noncoding ">".$GENE4{$cnt4}."_".$SP4{$cnt4}."\n".$str."\n";
			}elsif ($GENE4{$cnt4} =~ /1_/) {#1_infB-psbA
				my $temp_intergenic=$GENE4{$cnt4};
				$temp_intergenic=~ s/1_//g;
				print $out_noncoding ">".$temp_intergenic."_".$SP4{$cnt4}."\n".$str."\n";
			}elsif ($GENE4{$cnt4} =~ /2_/) {#rpl2-2_infB
				my $temp_intergenic=$GENE4{$cnt4};
				$temp_intergenic=~ s/2_//g;
				print $out_noncoding ">".$temp_intergenic."_".$SP4{$cnt4}."\n".$str."\n";
			}elsif ($GENE4{$cnt4} =~ /start-(.*)/) {#start-psbA
				$right_header2=$1;
				$right_str2=$str;
			}elsif ($GENE4{$cnt4} =~ /(.*)-end/) {#rpl2-end
				$left_header2=$1;
				$left_str2=$str;
			}
			if ((defined $left_header2 ne "") and (defined $right_header2 ne "")) {
				print $out_noncoding ">".$left_header2."-".$right_header2."_".$SP4{$cnt4}."\n".$left_str2.$right_str2."\n";
			}
		}
		@output4=();
		close $out_noncoding;
	}
}



####extract genes with introns and intergenic sequences
if (($pattern == 2) or ($pattern == 5) or ($pattern == 7)) {
	while (@filenames2) {
		my $filename_gb=shift @filenames2;
		my $target_name=substr($filename_gb,0,rindex($filename_gb,"\."));
		$target_name=~ s/\s/_/g;
		my $filename=substr($target_name,rindex($target_name,"\/")+1);
		#my $target_path=$target_name;
		#$target_path=~ s/$filename//g;
		my $target_name_temp_random1 = $input_directory."/".$random."_".$filename."_temp1";
		my $target_name_temp_random2 = $input_directory."/".$random."_".$filename."_temp2";

		open(my $in_gb1,"<",$filename_gb);
		open(my $out_gb1,">",$target_name_temp_random1);
		while (<$in_gb1>){
			$_=~ s/\r\n/\n/g;
			if ($_=~ /\),\n/){
				$_=~ s/\),\n/\),/g;
			}elsif($_=~ /,\n/){
				$_=~ s/,\n/,/g;
			}
			print $out_gb1 $_;
		}
		close $in_gb1;
		close $out_gb1;

		open(my $in_gb2,"<",$target_name_temp_random1);
		open(my $out_gb2,">",$target_name_temp_random2);
		while (<$in_gb2>){
			$_=~ s/,\s+/,/g;
			print $out_gb2 $_;
		}
		close $in_gb2;
		close $out_gb2;


		#generate_bed_file
		my (@row_array,@start_end,$species_name,$length,$element,@genearray,@output1);
		my $mark=0;
		my $tick=0;
		open (my $in_gb3,"<",$target_name_temp_random2);
		while (<$in_gb3>){
			$_=~ s/\r|\n//g;
			@row_array=split /\s+/,$_;
			if (/^LOCUS/i){
				$species_name=$filename;
				$length=$row_array[2];
			}elsif(/ {5}gene {12}/){
				if ($row_array[2]=~ /^\d+..\d+$/){
					$row_array[2]="\+\t$row_array[2]\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~/^complement\((\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
					my @array_start_end=split /\t/,$row_array[2];
					my $first_start_end=$row_array[2];
					my $second_start_end=$row_array[2];
					if ($array_start_end[1] < $array_start_end[2]) {
						$row_array[2]=$row_array[2];
					}elsif ($array_start_end[1] > $array_start_end[2]) {
						$first_start_end=~ s/$array_start_end[2]/$length/g;
						$second_start_end=~ s/$array_start_end[1]/1/g;
						push @start_end,$first_start_end;
						push @start_end,$second_start_end;
						$tick=3;
					}
				}elsif($row_array[2]=~ /^join\((\d+..\d+),(\d+..\d+)\)$/) {
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\((\d+..\d+),(\d+..\d+),(\d+..\d+)\)$/) {
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),(\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\((\d+..\d+),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="+\t$1\t$row_array[1]\t-\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^complement\(join\((\d+..\d+),(\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^complement\(join\((\d+..\d+),(\d+..\d+),(\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\((\d+..\d+),(\d+..\d+)\)$/){
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\((\d+..\d+),(\d+..\d+),(\d+..\d+)\)$/){
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),(\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),(\d+..\d+),(\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),(\d+..\d+),(\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\((\d+..\d+),(\d+..\d+),complement\((\d+..\d+)\)\)$/) {
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^<(\d+..\d+)$/){# positive split no-intron gene
					$row_array[2]="\+\t$1\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
					$tick=2;
				}elsif($row_array[2]=~ /^(\d+..>\d+)$/){# positive split no-intron gene
					$row_array[2]="\+\t$1\t$row_array[1]";
					$row_array[2]=~ s/\..>/\t/g;
					$tick=1;
				}elsif($row_array[2]=~/^complement\(<(\d+..\d+)\)$/){# negative split no-intron gene
					$row_array[2]="-\t$1\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
					$tick=1;
				}elsif($row_array[2]=~/^complement\((\d+..>\d+)\)$/){# negative split no-intron gene
					$row_array[2]="-\t$1\t$row_array[1]";
					$row_array[2]=~ s/\..>/\t/g;
					$tick=2;
				}
				$element=$row_array[2];
				$mark=1;
			}elsif(/^\s+\/gene="(.*)"/ and $mark == 1){
				if ($tick==0) {
					$element=$species_name.":".$1.":".$element;
					push @genearray,$element;
					$element=();
					$mark=0;
				}elsif ($tick==1) {
					$element=$species_name.":"."1_".$1.":".$element;#1_trnH-GUG
					push @genearray,$element;
					$element=();
					$mark=0;
					$tick=0;
				}elsif ($tick==2) {
					$element=$species_name.":"."2_".$1.":".$element;#2_trnH-GUG
					push @genearray,$element;
					$element=();
					$mark=0;
					$tick=0;
				}elsif ($tick==3) {
					push @genearray,$species_name.":"."1_".$1.":".$start_end[1];#1_trnH-GUG
					push @genearray,$species_name.":"."2_".$1.":".$start_end[0];#2_trnH-GUG
					@start_end=();
					$mark=0;
					$tick=0;
				}
			}
		}
		close $in_gb3;
		foreach (@genearray){
			my @array=split /:/,$_;
			push @output1,"$array[0]\t$array[1]\t$array[2]\n";
		}
		@row_array=();
		@genearray=();


		#put_fasta_sequence_in_array
		my $flag=0;
		my @sequence;
		my (@fas1,@fas2);
		open(my $in_gb4,"<",$target_name_temp_random2);
		while (<$in_gb4>){
			if ($_=~ /ORIGIN/){
				$flag=1;
			}
			if ($_=~ /\/\//){
				$flag=2;
			}
			if ($flag==1){
				next if ($_=~ /ORIGIN/);
				push @sequence,$_;
			}
		}
		close $in_gb4;
		foreach (@sequence){
			$_=~ s/\r|\n//g;
			$_=~ s/\s*//g;
			$_=~ s/\d+//g;
			push @fas1,$_;
		}
		my $fas1=join "",@fas1;
		my (@fasta1,@fasta2);
		push @fasta1,$species_name,$fas1;
		@fasta2=@fasta1;

		unlink ("$target_name_temp_random1");
		unlink ("$target_name_temp_random2");


		#edit_bed_file
		my (%SP1,%GENE1,%STRAND1,%START1,%END1,%TYPE1,%STRAND2,%START2,%END2,%TYPE2,%STRAND3,%START3,%END3,%TYPE3,@output2);
		my $cnt1=0;
		foreach (@output1) {
			$_=~ s/\r|\n//g;
			$cnt1++;
			($SP1{$cnt1},$GENE1{$cnt1},$STRAND1{$cnt1},$START1{$cnt1},$END1{$cnt1},$TYPE1{$cnt1},$STRAND2{$cnt1},$START2{$cnt1},$END2{$cnt1},$TYPE2{$cnt1})=(split /\s+/,$_)[0,1,2,3,4,5,6,7,8,9];
		}
		foreach (1..$cnt1) {
			if (defined $STRAND2{$_} eq "") {
				push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
			}elsif (defined $STRAND2{$_} ne "") {
				if (($GENE1{$_} ne "rps12")) {
					if (($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START1{$_} < $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START1{$_} > $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START1{$_} < $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START1{$_} > $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}
				}elsif (($GENE1{$_} eq "rps12")) {
					if (($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START1{$_} < $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START1{$_} > $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START1{$_} < $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START1{$_} > $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});

					}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "+") and ($START1{$_} < $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "+") and ($START1{$_} > $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "-") and ($START1{$_} < $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "-") and ($START1{$_} > $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}
				}
			}
		}
		@output1=();


		#Rosa_roxburghii        rps12   -       99457   100247  gene    -       71096   71209   gene
		#Rosa_roxburghii        rps12   -       71096   71209   gene    +       142355  143145  gene

		#Rosa_roxburghii        rps12   +       99457   100247  gene    +       71096   71209   gene
		#Rosa_roxburghii        rps12   +       71096   71209   gene    -       142355  143145  gene

		#Rosa_roxburghii        rps12   -       71096   71209   gene
		#Rosa_roxburghii        rps12   -       99457   100250  gene
		#Rosa_roxburghii        rps12   +       142352  143145  gene

		#Rosa_roxburghii        rps12   +       71096   71209   gene
		#Rosa_roxburghii        rps12   +       99457   100250  gene
		#Rosa_roxburghii        rps12   -       142352  143145  gene


		#sort_bed_file
		my $col=3;
		my (%sort,@output3);
		foreach (@output2){
			my @row=split /\t/,$_;
			$sort{$_}=$row[$col];
		}
		foreach (sort {$sort{$a} <=> $sort{$b}} keys %sort){
			push @output3,"$_\n";
		}
		@output2=();


		##output_bed_file
		#open (my $out_bed,">","$output_directory/$filename\_p2_genes.bed");
		#foreach (@output3){
		#	print $out_bed $_;
		#}
		#close $out_bed;

		#output_bed_file
		open (my $out_bed,">","$output_directory/$filename\_p2_genes.bed");
		foreach (@output3){
			my @row=split /\t/,$_;
			my $abs=abs($row[3]-$row[4]);
			if ($_!~ /rps12/) {
				print $out_bed $_;
			}
			if (($_=~ /rps12/) and ($abs < 200)) {
				$_=~ s/rps12/rps12-1/;
				print $out_bed $_;
			}
			if (($_=~ /rps12/) and ($abs > 200)) {
				$_=~ s/rps12/rps12-2/;
				print $out_bed $_;
			}
		}
		close $out_bed;


		#extract_gene
		my (%SP2,%GENE2,%STRAND4,%START4,%END4,%TYPE4,$seq1);
		my $cnt2=0;
		open(my $out_coding,">","$output_directory/$filename\_p2_genes.fasta");
		while (@fasta1){
			my $header=shift @fasta1;
			$seq1=shift @fasta1;
		}

		my ($rps12_header,@rps12_1_sequence,@rps12_2_sequence,@rps12_3_sequence);
		my ($left_header,$left_str,$right_header,$right_str);
		foreach (@output3){
			$_=~ s/\r|\n//g;
			$cnt2++;
			($SP2{$cnt2},$GENE2{$cnt2},$STRAND4{$cnt2},$START4{$cnt2},$END4{$cnt2},$TYPE4{$cnt2})=(split /\s+/,$_)[0,1,2,3,4,5];
			my $str=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
			#if ($STRAND4{$cnt2} eq "-") {
			#	my $rev_com=reverse $str;
			#	$rev_com=~ tr/ACGTacgt/TGCAtgca/;
			#	print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$rev_com."\n";
			#}elsif($STRAND4{$cnt2} eq "+"){
			#	print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$str."\n";
			#}
			if ((($STRAND4{$cnt2} eq "-") and ($GENE2{$cnt2} ne "rps12-1")) and (($STRAND4{$cnt2} eq "-") and ($GENE2{$cnt2} ne "rps12-2"))) {
				my $rev_com=reverse $str;
				$rev_com=~ tr/ACGTacgt/TGCAtgca/;
				if (($GENE2{$cnt2} !~ /1_/) and ($GENE2{$cnt2} !~ /2_/)) {
					print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$rev_com."\n";
				}elsif ($GENE2{$cnt2} =~ /1_(.*)/) {#1_trnH-GUG
					$left_header=$1;
					$left_str=$rev_com;
				}elsif ($GENE2{$cnt2} =~ /2_(.*)/) {#2_trnH-GUG
					$right_header=$1;
					$right_str=$rev_com;
				}
			}elsif (($STRAND4{$cnt2} eq "-") and ($GENE2{$cnt2} eq "rps12-1")) {
				my $rev_com=reverse $str;
				$rev_com=~ tr/ACGTacgt/TGCAtgca/;
				push @rps12_1_sequence,$rev_com;
			}elsif (($STRAND4{$cnt2} eq "-") and ($GENE2{$cnt2} eq "rps12-2")) {
				my $rev_com=reverse $str;
				$rev_com=~ tr/ACGTacgt/TGCAtgca/;
				push @rps12_2_sequence,$rev_com;
			}elsif(($STRAND4{$cnt2} eq "+") and ($GENE2{$cnt2} ne "rps12-2")){
				if (($GENE2{$cnt2} !~ /1_/) and ($GENE2{$cnt2} !~ /2_/)) {
					print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$str."\n";
				}elsif ($GENE2{$cnt2} =~ /1_(.*)/) {#1_trnH-GUG
					$left_header=$1;
					$left_str=$str;
				}elsif ($GENE2{$cnt2} =~ /2_(.*)/) {#2_trnH-GUG
					$right_header=$1;
					$right_str=$str;
				}
			}elsif (($STRAND4{$cnt2} eq "+") and ($GENE2{$cnt2} eq "rps12-2")) {
				push @rps12_3_sequence,$str;
			}
			if ((defined $left_header ne "") and (defined $right_header ne "") and ($left_header eq $right_header)) {
				print $out_coding ">".$left_header."_".$SP2{$cnt2}."\n".$left_str.$right_str."\n";
			}
			$rps12_header=">rps12_".$SP2{$cnt2};
		}

		if ((@rps12_1_sequence != 0) and (@rps12_2_sequence != 0) and (@rps12_3_sequence != 0)) {
			my $rps12_1_sequence=shift @rps12_1_sequence;
			my $rps12_2_sequence=shift @rps12_2_sequence;
			my $rps12_3_sequence=shift @rps12_3_sequence;
			print $out_coding $rps12_header."\n".$rps12_1_sequence.$rps12_2_sequence."\n";
			print $out_coding $rps12_header."\n".$rps12_1_sequence.$rps12_3_sequence."\n";
		}elsif ((@rps12_1_sequence != 0) and (@rps12_2_sequence != 0) and (@rps12_3_sequence == 0)) {
			my $rps12_1_sequence=shift @rps12_1_sequence;
			my $rps12_2_sequence=shift @rps12_2_sequence;
			print $out_coding $rps12_header."\n".$rps12_1_sequence.$rps12_2_sequence."\n";
			print $out_coding $rps12_header."\n".$rps12_1_sequence.$rps12_2_sequence."\n";
		}elsif ((@rps12_1_sequence != 0) and (@rps12_2_sequence == 0) and (@rps12_3_sequence != 0)) {
			my $rps12_1_sequence=shift @rps12_1_sequence;
			my $rps12_3_sequence=shift @rps12_3_sequence;
			print $out_coding $rps12_header."\n".$rps12_1_sequence.$rps12_3_sequence."\n";
			print $out_coding $rps12_header."\n".$rps12_1_sequence.$rps12_3_sequence."\n";
		}elsif ((@rps12_1_sequence == 0) and (@rps12_2_sequence != 0) and (@rps12_3_sequence != 0)) {
			my $rps12_2_sequence=shift @rps12_2_sequence;
			my $rps12_3_sequence=shift @rps12_3_sequence;
			print $out_coding $rps12_header."\n".$rps12_2_sequence."\n";
			print $out_coding $rps12_header."\n".$rps12_3_sequence."\n";
		}elsif ((@rps12_1_sequence == 0) and (@rps12_2_sequence == 0) and (@rps12_3_sequence != 0)) {
			my $rps12_3_sequence=shift @rps12_3_sequence;
			print $out_coding $rps12_header."\n".$rps12_3_sequence."\n";
			print $out_coding $rps12_header."\n".$rps12_3_sequence."\n";
		}elsif ((@rps12_1_sequence == 0) and (@rps12_2_sequence != 0) and (@rps12_3_sequence == 0)) {
			my $rps12_2_sequence=shift @rps12_2_sequence;
			print $out_coding $rps12_header."\n".$rps12_2_sequence."\n";
			print $out_coding $rps12_header."\n".$rps12_2_sequence."\n";
		}elsif ((@rps12_1_sequence != 0) and (@rps12_2_sequence == 0) and (@rps12_3_sequence == 0)) {
			my $rps12_1_sequence=shift @rps12_1_sequence;
			print $out_coding $rps12_header."\n".$rps12_1_sequence."\n";
			print $out_coding $rps12_header."\n".$rps12_1_sequence."\n";
		}elsif ((@rps12_1_sequence == 0) and (@rps12_2_sequence == 0) and (@rps12_3_sequence == 0)) {
		}

		close $out_coding;


		#generate_IGS_ranges
		my (%SP3,%GENE3,%STRAND5,%START5,%END5,%TYPE5,$last0,$last1,$last2,@output4);
		my $cnt3=0;
		foreach (@output3){
			$_=~ s/\r|\n//g;
			$cnt3++;
			($SP3{$cnt3},$GENE3{$cnt3},$STRAND5{$cnt3},$START5{$cnt3},$END5{$cnt3},$TYPE5{$cnt3})=(split /\s+/,$_)[0,1,2,3,4,5];
		}
		foreach (keys %SP3){
			if ($_==1 and $START5{$_}!=1){
				unshift @output4,$SP3{$_}."\t"."start".'-'.$GENE3{$_}."\t"."?"."/".$STRAND5{$_}."\t"."1"."\t".($START5{$_}-1)."\t"."?"."/".$TYPE5{$_}."\n";
			}
		}
		foreach (1..($cnt3-1)) {
			$last0=$_-1;
			$last1=$_+1;
			$last2=$_+2;
			next if ((($END5{$_}+1) > ($START5{$last1}-1)) and (($END5{$_}+1) < ($END5{$last1}-1)));
			next if (($_ > 1) and (($END5{$_}+1) < ($END5{$last0}-1)) and (($END5{$_}+1) < ($START5{$last2}-1)));
			if ((($END5{$_}+1) >= ($START5{$last1}-1)) and (($END5{$_}+1) >= ($END5{$last1}-1))){
				push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'.$GENE3{$last2}."\t".$STRAND5{$_}."/".$STRAND5{$last2}."\t".($END5{$_}+1)."\t".($START5{$last2}-1)."\t".$TYPE5{$_}."/".$TYPE5{$last2}."\n";
			}else{
				push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'.$GENE3{$last1}."\t".$STRAND5{$_}."/".$STRAND5{$last1}."\t".($END5{$_}+1)."\t".($START5{$last1}-1)."\t".$TYPE5{$_}."/".$TYPE5{$last1}."\n";
			}
		}
		foreach (keys %SP3){
			if ($_==$cnt3 and $END5{$_}!=$length){
				push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'."end"."\t".$STRAND5{$_}."/"."?"."\t".($END5{$_}+1)."\t".$length."\t".$TYPE5{$_}."/"."?"."\n";
			}
		}
		@output3=();
		#print @output4;


		#extract_IGS
		my (%SP4,%GENE4,%STRAND6,%START6,%END6,%TYPE6,$seq2);
		my $cnt4=0;
		open(my $out_noncoding,">","$output_directory/$filename\_p2_regions_among_genes.fasta");
		while (@fasta2){
			my $header=shift @fasta2;
			$seq2=shift @fasta2;
		}
		my ($left_header2,$left_str2,$right_header2,$right_str2);
		foreach (@output4){
			$_=~ s/\r|\n//g;
			$cnt4++;
			($SP4{$cnt4},$GENE4{$cnt4},$STRAND6{$cnt4},$START6{$cnt4},$END6{$cnt4},$TYPE6{$cnt4})=(split /\s+/,$_)[0,1,2,3,4,5];
			my $str=substr($seq2,($START6{$cnt4}-1),($END6{$cnt4}-$START6{$cnt4}+1));
			if (($GENE4{$cnt4} !~ /1_/) and ($GENE4{$cnt4} !~ /2_/) and ($GENE4{$cnt4} !~ /start-/) and ($GENE4{$cnt4} !~ /-end/)) {
				print $out_noncoding ">".$GENE4{$cnt4}."_".$SP4{$cnt4}."\n".$str."\n";
			}elsif ($GENE4{$cnt4} =~ /1_/) {#1_trnH-GUG-psbA
				my $temp_intergenic=$GENE4{$cnt4};
				$temp_intergenic=~ s/1_//g;
				print $out_noncoding ">".$temp_intergenic."_".$SP4{$cnt4}."\n".$str."\n";
			}elsif ($GENE4{$cnt4} =~ /2_/) {#rpl2-2_trnH-GUG
				my $temp_intergenic=$GENE4{$cnt4};
				$temp_intergenic=~ s/2_//g;
				print $out_noncoding ">".$temp_intergenic."_".$SP4{$cnt4}."\n".$str."\n";
			}elsif ($GENE4{$cnt4} =~ /start-(.*)/) {#start-trnH-GUG
				$right_header2=$1;
				$right_str2=$str;
			}elsif ($GENE4{$cnt4} =~ /(.*)-end/) {#rpl2-end
				$left_header2=$1;
				$left_str2=$str;
			}
			if ((defined $left_header2 ne "") and (defined $right_header2 ne "")) {
				print $out_noncoding ">".$left_header2."-".$right_header2."_".$SP4{$cnt4}."\n".$left_str2.$right_str2."\n";
			}
		}
		@output4=();
		close $out_noncoding;
	}
}




####extract coding and noncoding, including linked CDSs and tRNAs
if (($pattern == 3) or ($pattern == 6) or ($pattern == 7)) {
	while (@filenames3) {
		my $filename_gb=shift @filenames3;
		my $target_name=substr($filename_gb,0,rindex($filename_gb,"\."));
		$target_name=~ s/\s/_/g;
		my $filename=substr($target_name,rindex($target_name,"\/")+1);
		#my $target_path=$target_name;
		#$target_path=~ s/$filename//g;
		my $target_name_temp_random1 = $input_directory."/".$random."_".$filename."_temp1";
		my $target_name_temp_random2 = $input_directory."/".$random."_".$filename."_temp2";

		open(my $in_gb1,"<",$filename_gb);
		open(my $out_gb1,">",$target_name_temp_random1);
		while (<$in_gb1>){
			$_=~ s/\r\n/\n/g;
			if ($_=~ /\),\n/){
				$_=~ s/\),\n/\),/g;
			}elsif($_=~ /,\n/){
				$_=~ s/,\n/,/g;
			}
			print $out_gb1 $_;
		}
		close $in_gb1;
		close $out_gb1;

		open(my $in_gb2,"<",$target_name_temp_random1);
		open(my $out_gb2,">",$target_name_temp_random2);
		while (<$in_gb2>){
			$_=~ s/,\s+/,/g;
			print $out_gb2 $_;
		}
		close $in_gb2;
		close $out_gb2;


		#generate_bed_file
		my (@row_array,@start_end,$species_name,$length,$element,@genearray,@output1);
		my $mark=0;
		my $tick=0;
		open (my $in_gb3,"<",$target_name_temp_random2);
		while (<$in_gb3>){
			$_=~ s/\r|\n//g;
			@row_array=split /\s+/,$_;
			if (/^LOCUS/i){
				$species_name=$filename;
				$length=$row_array[2];
			}elsif(/ {5}CDS {13}/ or / {5}tRNA {12}/ or / {5}rRNA {12}/){
				if ($row_array[2]=~ /^\d+..\d+$/){
					$row_array[2]="\+\t$row_array[2]\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~/^complement\((\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
					my @array_start_end=split /\t/,$row_array[2];
					my $first_start_end=$row_array[2];
					my $second_start_end=$row_array[2];
					if ($array_start_end[1] < $array_start_end[2]) {
						$row_array[2]=$row_array[2];
					}elsif ($array_start_end[1] > $array_start_end[2]) {
						$first_start_end=~ s/$array_start_end[2]/$length/g;
						$second_start_end=~ s/$array_start_end[1]/1/g;
						push @start_end,$first_start_end;
						push @start_end,$second_start_end;
						$tick=3;
					}
				}elsif($row_array[2]=~ /^join\((\d+..\d+),(\d+..\d+)\)$/) {
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\((\d+..\d+),(\d+..\d+),(\d+..\d+)\)$/) {
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),(\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\((\d+..\d+),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="+\t$1\t$row_array[1]\t-\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^complement\(join\((\d+..\d+),(\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^complement\(join\((\d+..\d+),(\d+..\d+),(\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\((\d+..\d+),(\d+..\d+)\)$/){
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\((\d+..\d+),(\d+..\d+),(\d+..\d+)\)$/){
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),(\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),(\d+..\d+),(\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),(\d+..\d+),(\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\((\d+..\d+),(\d+..\d+),complement\((\d+..\d+)\)\)$/) {
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^<(\d+..\d+)$/){# positive split no-intron gene
					$row_array[2]="\+\t$1\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
					$tick=2;
				}elsif($row_array[2]=~ /^(\d+..>\d+)$/){# positive split no-intron gene
					$row_array[2]="\+\t$1\t$row_array[1]";
					$row_array[2]=~ s/\..>/\t/g;
					$tick=1;
				}elsif($row_array[2]=~/^complement\(<(\d+..\d+)\)$/){# negative split no-intron gene
					$row_array[2]="-\t$1\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
					$tick=1;
				}elsif($row_array[2]=~/^complement\((\d+..>\d+)\)$/){# negative split no-intron gene
					$row_array[2]="-\t$1\t$row_array[1]";
					$row_array[2]=~ s/\..>/\t/g;
					$tick=2;
				}
				$element=$row_array[2];
				$mark=1;
			}elsif(/^\s+\/gene="(.*)"/ and $mark == 1){
				if ($tick==0) {
					$element=$species_name.":".$1.":".$element;
					push @genearray,$element;
					$element=();
					$mark=0;
				}elsif ($tick==1) {
					$element=$species_name.":"."1_".$1.":".$element;#1_trnH-GUG
					push @genearray,$element;
					$element=();
					$mark=0;
					$tick=0;
				}elsif ($tick==2) {
					$element=$species_name.":"."2_".$1.":".$element;#2_trnH-GUG
					push @genearray,$element;
					$element=();
					$mark=0;
					$tick=0;
				}elsif ($tick==3) {
					push @genearray,$species_name.":"."1_".$1.":".$start_end[1];#1_trnH-GUG
					push @genearray,$species_name.":"."2_".$1.":".$start_end[0];#2_trnH-GUG
					@start_end=();
					$mark=0;
					$tick=0;
				}
			}
		}
		close $in_gb3;
		foreach (@genearray){
			my @array=split /:/,$_;
			push @output1,"$array[0]\t$array[1]\t$array[2]\n";
		}
		@row_array=();
		@genearray=();


		#put_fasta_sequence_in_array
		my $flag=0;
		my @sequence;
		my (@fas1,@fas2);
		open(my $in_gb4,"<",$target_name_temp_random2);
		while (<$in_gb4>){
			if ($_=~ /ORIGIN/){
				$flag=1;
			}
			if ($_=~ /\/\//){
				$flag=2;
			}
			if ($flag==1){
				next if ($_=~ /ORIGIN/);
				push @sequence,$_;
			}
		}
		close $in_gb4;
		foreach (@sequence){
			$_=~ s/\r|\n//g;
			$_=~ s/\s*//g;
			$_=~ s/\d+//g;
			push @fas1,$_;
		}
		my $fas1=join "",@fas1;
		my (@fasta1,@fasta2);
		push @fasta1,$species_name,$fas1;
		@fasta2=@fasta1;

		unlink ("$target_name_temp_random1");
		unlink ("$target_name_temp_random2");


		#edit_bed_file
		my (%SP1,%GENE1,%STRAND1,%START1,%END1,%TYPE1,%STRAND2,%START2,%END2,%TYPE2,%STRAND3,%START3,%END3,%TYPE3,@output2);
		my $cnt1=0;
		foreach (@output1) {
			$_=~ s/\r|\n//g;
			$cnt1++;
			($SP1{$cnt1},$GENE1{$cnt1},$STRAND1{$cnt1},$START1{$cnt1},$END1{$cnt1},$TYPE1{$cnt1},$STRAND2{$cnt1},$START2{$cnt1},$END2{$cnt1},$TYPE2{$cnt1},$STRAND3{$cnt1},$START3{$cnt1},$END3{$cnt1},$TYPE3{$cnt1})=(split /\s+/,$_)[0,1,2,3,4,5,6,7,8,9,10,11,12,13];
		}
		foreach (1..$cnt1) {
			if (defined $STRAND2{$_} eq "") {
				push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
			}elsif ((defined $STRAND2{$_} ne "") and (defined $STRAND3{$_} eq "")) {
				if ($GENE1{$_} ne "rps12") {
					if (($STRAND1{$_} eq "-") and ($START1{$_} < $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}elsif(($STRAND1{$_} eq "-") and ($START1{$_} > $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}elsif(($STRAND1{$_} eq "+") and ($START1{$_} < $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}elsif(($STRAND1{$_} eq "+") and ($START1{$_} > $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}
				}elsif ($GENE1{$_} eq "rps12") {
					if (abs($START1{$_}-$START2{$_}) < 5000) {
						if (($STRAND1{$_} eq "-") and ($START1{$_} < $START2{$_})){
							push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
							push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						}elsif(($STRAND1{$_} eq "-") and ($START1{$_} > $START2{$_})){
							push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
							push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						}elsif(($STRAND1{$_} eq "+") and ($START1{$_} < $START2{$_})){
							push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
							push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						}elsif(($STRAND1{$_} eq "+") and ($START1{$_} > $START2{$_})){
							push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
							push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						}
					}elsif (abs($START1{$_}-$START2{$_}) > 5000) {
						if (($STRAND1{$_} eq "-") and ($START1{$_} < $START2{$_})){
							push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
							push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						}elsif(($STRAND1{$_} eq "-") and ($START1{$_} > $START2{$_})){
							push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
							push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						}elsif(($STRAND1{$_} eq "+") and ($START1{$_} < $START2{$_})){
							push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
							push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						}elsif(($STRAND1{$_} eq "+") and ($START1{$_} > $START2{$_})){
							push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
							push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						}
					}
				}
			}elsif ((defined $STRAND2{$_} ne "") and (defined $STRAND3{$_} ne "")) {
				if ($GENE1{$_} ne "rps12") {
					if (($STRAND1{$_} eq "-") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-3"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "-") and ($START1{$_} > $START2{$_}) and ($START2{$_} > $START3{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-3"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "-") and ($START1{$_} > $START3{$_}) and ($START3{$_} > $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-3"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "+") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-3"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "+") and ($START1{$_} > $START2{$_}) and ($START2{$_} > $START3{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-3"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "+") and ($START1{$_} > $START3{$_}) and ($START3{$_} > $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-3"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}
				}elsif ($GENE1{$_} eq "rps12") {
					if (($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_})){#1<2<3
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START2{$_} > $START3{$_}) and ($START3{$_} > $START1{$_})){#1<3<2
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START1{$_} > $START2{$_}) and ($START3{$_} > $START1{$_})){#2<1<3
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START1{$_} > $START3{$_}) and ($START3{$_} > $START2{$_})){#2<3<1
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START2{$_} > $START1{$_}) and ($START1{$_} > $START3{$_})){#3<1<2
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START1{$_} > $START2{$_}) and ($START2{$_} > $START3{$_})){#3<2<1
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});

					}elsif (($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "+") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_})){#1<2<3
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "+") and ($START2{$_} > $START3{$_}) and ($START3{$_} > $START1{$_})){#1<3<2
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "+") and ($START1{$_} > $START2{$_}) and ($START3{$_} > $START1{$_})){#2<1<3
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "+") and ($START1{$_} > $START3{$_}) and ($START3{$_} > $START2{$_})){#2<3<1
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "+") and ($START2{$_} > $START1{$_}) and ($START1{$_} > $START3{$_})){#3<1<2
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "+") and ($START1{$_} > $START2{$_}) and ($START2{$_} > $START3{$_})){#3<2<1
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});

					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_})){#1<2<3
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START1{$_} < $START3{$_}) and ($START3{$_} < $START2{$_})){#1<3<2
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START2{$_} < $START1{$_}) and ($START1{$_} < $START3{$_})){#2<1<3
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START1{$_} > $START3{$_}) and ($START3{$_} > $START2{$_})){#2<3<1
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START2{$_} > $START1{$_}) and ($START1{$_} > $START3{$_})){#3<1<2
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START1{$_} > $START2{$_}) and ($START2{$_} > $START3{$_})){#3<2<1
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});

					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "-") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_})){#1<2<3
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "-") and ($START1{$_} < $START3{$_}) and ($START3{$_} < $START2{$_})){#1<3<2
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "-") and ($START2{$_} < $START1{$_}) and ($START1{$_} < $START3{$_})){#2<1<3
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "-") and ($START1{$_} > $START3{$_}) and ($START3{$_} > $START2{$_})){#2<3<1
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "-") and ($START2{$_} > $START1{$_}) and ($START1{$_} > $START3{$_})){#3<1<2
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "-") and ($START1{$_} > $START2{$_}) and ($START2{$_} > $START3{$_})){#3<2<1
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}
				}
			}
		}


		#sort_bed_file
		my $col=3;
		my (%sort,@output3);
		foreach (@output2){
			my @row=split /\t/,$_;
			$sort{$_}=$row[$col];
		}
		foreach (sort {$sort{$a} <=> $sort{$b}} keys %sort){
			push @output3,"$_\n";
		}
		@output2=();


		#output_bed_file
		#open (my $out_bed,">","$output_directory/$filename\_p3_linked_CDS_RNA.bed");
		#foreach (@output3){
		#	print $out_bed $_;
		#}
		#close $out_bed;


		#output_bed_file
		open (my $out_bed,">","$output_directory/$filename\_p3_linked_CDS_RNA.bed");
		foreach (@output3){
			my @row=split /\t/,$_;
			my $abs=abs($row[3]-$row[4]);

			if (($row[1] ne "rps12") and ($row[1] ne "rps12-1") and ($row[1] ne "rps12-2")) {
				print $out_bed $_;
			}
			if (($row[1] eq "rps12") and ($abs < 200)) {
				$_=~ s/rps12/rps12-1/;
				print $out_bed $_;
			}
			if (($row[1] eq "rps12") and ($abs > 200)) {
				$_=~ s/rps12/rps12-2/;
				print $out_bed $_;
			}
			if ($row[1] eq "rps12-1") {
				$_=~ s/rps12-1/rps12-2/;
				print $out_bed $_;
			}
			if ($row[1] eq "rps12-2") {
				$_=~ s/rps12-2/rps12-3/;
				print $out_bed $_;
			}
		}
		close $out_bed;


		####for rps12 in @output1 to @output2
		####unlink of rps12, with an intron between exon2 and exon3
		#Rosa_roxburghii2        rps12   -       71096   71209   CDS
		#Rosa_roxburghii2        rps12   -       99457   99482   CDS     -       100019  100250  CDS
		#Rosa_roxburghii2        rps12   +       142352  142583  CDS     +       143120  143145  CDS

		#Rosa_roxburghii2        rps12   -       99457   99482   CDS     -       100019  100250  CDS
		#Rosa_roxburghii2        rps12   +       142352  142583  CDS     +       143120  143145  CDS

		#Rosa_roxburghii2        rps12   -       71096   71209   CDS
		#Rosa_roxburghii2        rps12   +       142352  142583  CDS     +       143120  143145  CDS

		#Rosa_roxburghii2        rps12   -       71096   71209   CDS
		#Rosa_roxburghii2        rps12   -       99457   99482   CDS     -       100019  100250  CDS

		#Rosa_roxburghii2        rps12   +       142352  142583  CDS     +       143120  143145  CDS

		#Rosa_roxburghii2        rps12   -       99457   99482   CDS     -       100019  100250  CDS

		#Rosa_roxburghii2        rps12   -       71096   71209   CDS


		####unlink of rps12, without an intron between exon2 and exon3
		#Rosa_roxburghii2        rps12   -       71096   71209   CDS
		#Rosa_roxburghii2        rps12   -       99457   100250  CDS
		#Rosa_roxburghii2        rps12   +       142352  143145  CDS

		#Rosa_roxburghii2        rps12   -       99457   100250  CDS
		#Rosa_roxburghii2        rps12   +       142352  143145  CDS

		#Rosa_roxburghii2        rps12   -       71096   71209   CDS
		#Rosa_roxburghii2        rps12   +       142352  143145  CDS

		#Rosa_roxburghii2        rps12   -       71096   71209   CDS
		#Rosa_roxburghii2        rps12   -       99457   100250  CDS

		#Rosa_roxburghii2        rps12   +       142352  143145  CDS

		#Rosa_roxburghii2        rps12   -       99457   100250  CDS

		#Rosa_roxburghii2        rps12   -       71096   71209   CDS


		####link of rps12, with an intron between exon2 and exon3
		#Rosa_roxburghii2        rps12   -       71096   71209   CDS	 -       99457   99482   CDS     -       100019  100250  CDS
		#Rosa_roxburghii2        rps12   -       71096   71209   CDS     +       142352  142583  CDS     +       143120  143145  CDS

		#Rosa_roxburghii2        rps12   -       71096   71209   CDS     +       142352  142583  CDS     +       143120  143145  CDS

		#Rosa_roxburghii2        rps12   -       71096   71209   CDS	 -       99457   99482   CDS     -       100019  100250  CDS


		####link of rps12, without an intron between exon2 and exon3
		#Rosa_roxburghii2        rps12   -       71096   71209   CDS	 -       99457   100250   CDS
		#Rosa_roxburghii2        rps12   -       71096   71209   CDS     +       142352  143145   CDS

		#Rosa_roxburghii2        rps12   -       71096   71209   CDS     +       142352  143145   CDS

		#Rosa_roxburghii2        rps12   -       71096   71209   CDS	 -       99457   100250   CDS


		#extract_gene
		my (%SP2,%GENE2,%STRAND4,%START4,%END4,%TYPE4,%STRAND5,%START5,%END5,%TYPE5,%STRAND6,%START6,%END6,%TYPE6,$seq1);
		my $cnt2=0;
		open(my $out_coding,">","$output_directory/$filename\_p3_linked_CDS_RNA.fasta");
		while (@fasta1){
			my $header=shift @fasta1;
			$seq1=shift @fasta1;
		}

		#sort_bed_file
		my $col2=3;
		my (%sort2,@output5);
		foreach (@output1){
			my @row=split /\t/,$_;
			$sort2{$_}=$row[$col2];
		}
		foreach (sort {$sort2{$a} <=> $sort2{$b}} keys %sort2){
			push @output5,"$_\n";
		}

		####for rps12 in @output1
		####unlink of rps12, with an intron between exon2 and exon3
		####unlink of rps12, without an intron between exon2 and exon3
		####link of rps12, with an intron between exon2 and exon3
		####link of rps12, without an intron between exon2 and exon3

		my ($rps12_header,@rps12_1_sequence,@rps12_2_sequence,@rps12_3_sequence);
		my ($left_header,$left_str,$right_header,$right_str);
		foreach (@output5){
			$_=~ s/\r|\n//g;
			$cnt2++;
			($SP2{$cnt2},$GENE2{$cnt2},$STRAND4{$cnt2},$START4{$cnt2},$END4{$cnt2},$TYPE4{$cnt2},$STRAND5{$cnt2},$START5{$cnt2},$END5{$cnt2},$TYPE5{$cnt2},$STRAND6{$cnt2},$START6{$cnt2},$END6{$cnt2},$TYPE6{$cnt2})=(split /\s+/,$_)[0,1,2,3,4,5,6,7,8,9,10,11,12,13];
			if (defined $STRAND5{$cnt2} eq "") {
				my $str1=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
				my $rev_com1=reverse $str1;
				$rev_com1=~ tr/ACGTacgt/TGCAtgca/;
				if ($GENE2{$cnt2} ne "rps12") {
					if ($STRAND4{$cnt2} eq "-") {
						if (($GENE2{$cnt2} !~ /1_/) and ($GENE2{$cnt2} !~ /2_/)) {
							print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$rev_com1."\n";
						}elsif ($GENE2{$cnt2} =~ /1_(.*)/) {#1_trnH-GUG
							$left_header=$1;
							$left_str=$rev_com1;
						}elsif ($GENE2{$cnt2} =~ /2_(.*)/) {#2_trnH-GUG
							$right_header=$1;
							$right_str=$rev_com1;
						}
					}elsif($STRAND4{$cnt2} eq "+"){
						if (($GENE2{$cnt2} !~ /1_/) and ($GENE2{$cnt2} !~ /2_/)) {
							print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$str1."\n";
						}elsif ($GENE2{$cnt2} =~ /1_(.*)/) {#1_trnH-GUG
							$left_header=$1;
							$left_str=$str1;
						}elsif ($GENE2{$cnt2} =~ /2_(.*)/) {#2_trnH-GUG
							$right_header=$1;
							$right_str=$str1;
						}
					}
				}elsif ($GENE2{$cnt2} eq "rps12") {
					my $length_rps12=abs($END4{$cnt2}-$START4{$cnt2}+1);
					if (($STRAND4{$cnt2} eq "-")) {
						if ($length_rps12 < 200) {
							push @rps12_1_sequence,$rev_com1;
						}elsif ($length_rps12 > 200) {
							push @rps12_2_sequence,$rev_com1;
						}
					}elsif($STRAND4{$cnt2} eq "+"){
						if ($length_rps12 < 200) {
							push @rps12_1_sequence,$str1;
						}elsif ($length_rps12 > 200) {
							push @rps12_2_sequence,$str1;
						}
					}
				}

			##unlink, rps12-1, rps12-2 with intron
			#Rosa_roxburghii1        rps12   -       71096   71209   CDS	<200
			#Rosa_roxburghii1        rps12   +       71096   71209   CDS	<200

			##unlink, rps12-1, rps12-2 without intron
			#Rosa_roxburghii1        rps12   -       71096   71209   CDS	<200
			#Rosa_roxburghii1        rps12   -       99457   100250  CDS	>200
			#Rosa_roxburghii1        rps12   +       142352  143145  CDS	>200
			#Rosa_roxburghii1        rps12   +       99457   100250  CDS	>200
			#Rosa_roxburghii1        rps12   -       142352  143145  CDS	>200

			#Rosa_roxburghii1        rps12   +       71096   71209   CDS	<200
			#Rosa_roxburghii1        rps12   +       99457   100250  CDS	>200
			#Rosa_roxburghii1        rps12   -       142352  143145  CDS	>200
			#Rosa_roxburghii1        rps12   -       99457   100250  CDS	>200
			#Rosa_roxburghii1        rps12   +       142352  143145  CDS	>200

			}elsif((defined $STRAND5{$cnt2} ne "") and (defined $STRAND6{$cnt2} eq "")) {
				if ($GENE2{$cnt2} ne "rps12") {
					if (($STRAND4{$cnt2} eq "-") and ($START4{$cnt2} < $START5{$cnt2})) {
						my $str2=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1)).substr($seq1,($START5{$cnt2}-1),($END5{$cnt2}-$START5{$cnt2}+1));
						my $rev_com2=reverse $str2;
						$rev_com2=~ tr/ACGTacgt/TGCAtgca/;
						print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$rev_com2."\n";
					}elsif(($STRAND4{$cnt2} eq "-") and ($START4{$cnt2} > $START5{$cnt2})) {
						my $str3=substr($seq1,($START5{$cnt2}-1),($END5{$cnt2}-$START5{$cnt2}+1)).substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
						my $rev_com3=reverse $str3;
						$rev_com3=~ tr/ACGTacgt/TGCAtgca/;
						print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$rev_com3."\n";
					}elsif(($STRAND4{$cnt2} eq "+") and ($START4{$cnt2} < $START5{$cnt2})){
						my $str4=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1)).substr($seq1,($START5{$cnt2}-1),($END5{$cnt2}-$START5{$cnt2}+1));
						print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$str4."\n";
					}elsif(($STRAND4{$cnt2} eq "+") and ($START4{$cnt2} > $START5{$cnt2})){
						my $str5=substr($seq1,($START5{$cnt2}-1),($END5{$cnt2}-$START5{$cnt2}+1)).substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
						print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$str5."\n";
					}

				##link, rps12-1, rps12-2 with intron
				#166675730931123         rps12   -       71096   71209   CDS     +       142355  142583  CDS     +       143090  143145  CDS
				#166675730931123         rps12   -       99457   99512   CDS     -       100019  100247  CDS     -       71096   71209   CDS


				##link, rps12-1, rps12-2 without intron
				#166675730931123         rps12   -       71096   71209   CDS     -       99457   100247  CDS	<200	>200
				#166675730931123         rps12   -       99457   100247  CDS     -       71096   71209   CDS	>200	<200

				#166675730931123         rps12   +       71096   71209   CDS     +       99457   100247  CDS	<200	>200
				#166675730931123         rps12   +       99457   100247  CDS     +       71096   71209   CDS	>200	<200

				#166675730931123         rps12   -       71096   71209   CDS     +       142355  143145  CDS	<200	>200
				#166675730931123         rps12   -       142355  143145  CDS     +       71096   71209   CDS	>200	<200
				#166675730931123         rps12   +       71096   71209   CDS     -       142355  143145  CDS	<200	>200
				#166675730931123         rps12   +       142355  143145  CDS     -       71096   71209   CDS	>200	<200


				##unlink, rps12-1, rps12-2 with intron
				#Rosa_roxburghii1        rps12   -       71096   71209   CDS	<200
				#Rosa_roxburghii1        rps12   -       99457   99482   CDS     -       100019  100250  CDS	<80		>200
				#Rosa_roxburghii1        rps12   -       100019  100250  CDS     -       99457   99482   CDS	>200	<80
				#Rosa_roxburghii1        rps12   +       142352  142583  CDS     +       143120  143145  CDS	>200	<80
				#Rosa_roxburghii1        rps12   +       143120  143145  CDS     +       142352  142583  CDS	<80		>200


				##unlink, rps12-1, rps12-2 without intron
				#Rosa_roxburghii1        rps12   -       71096   71209   CDS	<200
				#Rosa_roxburghii1        rps12   -       99457   100250  CDS	>200
				#Rosa_roxburghii1        rps12   +       142352  143145  CDS	>200

				#Rosa_roxburghii1        rps12   +       71096   71209   CDS	<200
				#Rosa_roxburghii1        rps12   +       99457   100250  CDS	>200
				#Rosa_roxburghii1        rps12   -       142352  143145  CDS	>200

				}elsif ($GENE2{$cnt2} eq "rps12") {
					my $length_rps12_1=abs($END4{$cnt2}-$START4{$cnt2}+1);
					my $length_rps12_2=abs($END5{$cnt2}-$START5{$cnt2}+1);
					my $str_rps12_1=substr($seq1,($START4{$cnt2}-1),$length_rps12_1);
					my $str_rps12_2=substr($seq1,($START5{$cnt2}-1),$length_rps12_2);
					my $str2=substr($seq1,($START4{$cnt2}-1),$length_rps12_1).substr($seq1,($START5{$cnt2}-1),$length_rps12_2);
					my $str3=substr($seq1,($START5{$cnt2}-1),$length_rps12_2).substr($seq1,($START4{$cnt2}-1),$length_rps12_1);
					my $rev_com1=reverse $str_rps12_1;
					my $rev_com2=reverse $str_rps12_2;
					my $rev_com3=reverse $str2;
					my $rev_com4=reverse $str3;
					$rev_com1=~ tr/ACGTacgt/TGCAtgca/;
					$rev_com2=~ tr/ACGTacgt/TGCAtgca/;
					$rev_com3=~ tr/ACGTacgt/TGCAtgca/;
					$rev_com4=~ tr/ACGTacgt/TGCAtgca/;
					if (($STRAND4{$cnt2} eq "-") and ($STRAND5{$cnt2} eq "-")) {
						if (($START4{$cnt2} < $START5{$cnt2}) and ($length_rps12_1 > 80) and ($length_rps12_1 < 200) and ($length_rps12_2 > 200)) {
							print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$rev_com3."\n";
						}elsif (($START4{$cnt2} > $START5{$cnt2}) and ($length_rps12_1 > 200) and ($length_rps12_2 > 80) and ($length_rps12_2 < 200)) {
							print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$rev_com4."\n";
						}elsif (($START4{$cnt2} < $START5{$cnt2}) and ($length_rps12_1 < 80) and ($length_rps12_2 > 200)) {
							push @rps12_2_sequence,$rev_com3;
						}elsif (($START4{$cnt2} > $START5{$cnt2}) and ($length_rps12_1 > 200) and ($length_rps12_2 < 80)) {
							push @rps12_2_sequence,$rev_com4;
						}
					}elsif(($STRAND4{$cnt2} eq "+") and ($STRAND5{$cnt2} eq "+")){
						if (($START4{$cnt2} < $START5{$cnt2}) and ($length_rps12_1 > 80) and ($length_rps12_1 < 200) and ($length_rps12_2 > 200)) {
							print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$str2."\n";
						}elsif (($START4{$cnt2} > $START5{$cnt2}) and ($length_rps12_1 > 200) and ($length_rps12_2 > 80) and ($length_rps12_2 < 200)) {
							print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$str3."\n";
						}elsif (($START4{$cnt2} < $START5{$cnt2}) and ($length_rps12_1 > 200) and ($length_rps12_2 < 80)) {
							push @rps12_2_sequence,$str2;
						}elsif (($START4{$cnt2} > $START5{$cnt2}) and ($length_rps12_1 < 80) and ($length_rps12_2 >200)) {
							push @rps12_2_sequence,$str3;
						}
					}elsif(($STRAND4{$cnt2} eq "-") and ($STRAND5{$cnt2} eq "+")){
						if (($START4{$cnt2} < $START5{$cnt2}) and ($length_rps12_1 > 80) and ($length_rps12_1 < 200) and ($length_rps12_2 > 200)) {
							print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$rev_com1.$str_rps12_2."\n";
						}elsif (($START4{$cnt2} > $START5{$cnt2}) and ($length_rps12_1 > 200) and ($length_rps12_2 > 80) and ($length_rps12_2 < 200)) {
							print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$str_rps12_2.$rev_com1."\n";
						}
					}elsif (($STRAND4{$cnt2} eq "+") and ($STRAND5{$cnt2} eq "-")) {
						if (($START4{$cnt2} < $START5{$cnt2}) and ($length_rps12_1 > 80) and ($length_rps12_1 < 200) and ($length_rps12_2 > 200)) {
							print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$str_rps12_1.$rev_com2."\n";
						}elsif (($START4{$cnt2} > $START5{$cnt2}) and ($length_rps12_1 > 200) and ($length_rps12_2 > 80) and ($length_rps12_2 < 200)) {
							print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$rev_com2.$str_rps12_1."\n";
						}
					}
				}
			}elsif ((defined $STRAND5{$cnt2} ne "") and (defined $STRAND6{$cnt2} ne "")) {
				if ($GENE2{$cnt2} ne "rps12") {
					if (($STRAND4{$cnt2} eq "-") and ($START4{$cnt2} < $START5{$cnt2}) and ($START5{$cnt2} < $START6{$cnt2})) {
						my $str6=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1)).substr($seq1,($START5{$cnt2}-1),($END5{$cnt2}-$START5{$cnt2}+1)).substr($seq1,($START6{$cnt2}-1),($END6{$cnt2}-$START6{$cnt2}+1));
						my $rev_com4=reverse $str6;
						$rev_com4=~ tr/ACGTacgt/TGCAtgca/;
						print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$rev_com4."\n";
					}elsif(($STRAND4{$cnt2} eq "-") and ($START4{$cnt2} > $START5{$cnt2}) and ($START5{$cnt2} > $START6{$cnt2})) {
						my $str7=substr($seq1,($START6{$cnt2}-1),($END6{$cnt2}-$START6{$cnt2}+1)).substr($seq1,($START5{$cnt2}-1),($END5{$cnt2}-$START5{$cnt2}+1)).substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
						my $rev_com5=reverse $str7;
						$rev_com5=~ tr/ACGTacgt/TGCAtgca/;
						print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$rev_com5."\n";
					#}elsif(($STRAND4{$cnt2} eq "-") and ($START4{$cnt2} > $START6{$cnt2}) and ($START6{$cnt2} > $START5{$cnt2})) {
					#	my $str8=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
					#	my $str9=substr($seq1,($START5{$cnt2}-1),($END5{$cnt2}-$START5{$cnt2}+1)).substr($seq1,($START6{$cnt2}-1),($END6{$cnt2}-$START6{$cnt2}+1));
					#	my $rev_com6=reverse $str8;
					#	$rev_com6=~ tr/ACGTacgt/TGCAtgca/;
					#	print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$rev_com6.$str9."\n";
					}elsif(($STRAND4{$cnt2} eq "+") and ($START4{$cnt2} < $START5{$cnt2}) and ($START5{$cnt2} < $START6{$cnt2})){
						my $str10=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1)).substr($seq1,($START5{$cnt2}-1),($END5{$cnt2}-$START5{$cnt2}+1)).substr($seq1,($START6{$cnt2}-1),($END6{$cnt2}-$START6{$cnt2}+1));
						print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$str10."\n";
					}elsif(($STRAND4{$cnt2} eq "+") and ($START4{$cnt2} > $START5{$cnt2}) and ($START5{$cnt2} > $START6{$cnt2})){
						my $str11=substr($seq1,($START6{$cnt2}-1),($END6{$cnt2}-$START6{$cnt2}+1)).substr($seq1,($START5{$cnt2}-1),($END5{$cnt2}-$START5{$cnt2}+1)).substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
						print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$str11."\n";
					#}elsif(($STRAND4{$cnt2} eq "+") and ($START4{$cnt2} > $START6{$cnt2}) and ($START6{$cnt2} > $START5{$cnt2})) {
					#	my $str12=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1)).substr($seq1,($START5{$cnt2}-1),($END5{$cnt2}-$START5{$cnt2}+1)).substr($seq1,($START6{$cnt2}-1),($END6{$cnt2}-$START6{$cnt2}+1));
					#	print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$str12."\n";
					}

				##link, rps12-1, rps12-2 with intron
				#166675730931123         rps12   -       71096   71209   CDS     -       99457   99512   CDS     -       100019  100247  CDS
				#166675730931123         rps12   -       71096   71209   CDS     +       142355  142583  CDS     +       143090  143145  CDS

				#Rosa_roxburghii1        rps12   -       99457   99512   CDS     -       100019  100247  CDS     -       71096   71209   CDS
				#Rosa_roxburghii1        rps12   -       71096   71209   CDS     +       142355  142583  CDS     +       143090  143145  CDS

				#166675730931123         rps12   -       71096   71209   CDS     -       99457   99512   CDS     -       100019  100247  CDS
				#166675730931123         rps12   +       142355  142583  CDS     +       143090  143145  CDS     -       71096   71209   CDS

				#166675730931123         rps12   -       99457   99512   CDS     -       100019  100247  CDS     -       71096   71209   CDS
				#166675730931123         rps12   +       142355  142583  CDS     +       143090  143145  CDS     -       71096   71209   CDS


				#166675730931123         rps12   +       71096   71209   CDS     +       99457   99512   CDS     +       100019  100247  CDS
				#166675730931123         rps12   +       71096   71209   CDS     -       142355  142583  CDS     -       143090  143145  CDS

				#Rosa_roxburghii1        rps12   +       99457   99512   CDS     +       100019  100247  CDS     +       71096   71209   CDS
				#Rosa_roxburghii1        rps12   +       71096   71209   CDS     -       142355  142583  CDS     -       143090  143145  CDS

				#166675730931123         rps12   +       71096   71209   CDS     +       99457   99512   CDS     +       100019  100247  CDS
				#166675730931123         rps12   -       142355  142583  CDS     -       143090  143145  CDS     +       71096   71209   CDS

				#166675730931123         rps12   +       99457   99512   CDS     +       100019  100247  CDS     +       71096   71209   CDS
				#166675730931123         rps12   -       142355  142583  CDS     -       143090  143145  CDS     +       71096   71209   CDS


				##link, rps12-1, rps12-2 with intron
				#Rosa_roxburghii1         rps12   -       71096   71209   CDS     -       99457   99512   CDS     -       100019  100247  CDS
				#Rosa_roxburghii1         rps12   -       99457   99512   CDS     -       100019  100247  CDS     -       71096   71209   CDS
				#Rosa_roxburghii1         rps12   -       142355  142583  CDS     -       143090  143145  CDS     +       71096   71209   CDS
				#Rosa_roxburghii1         rps12   +       71096   71209   CDS     -       142355  142583  CDS     -       143090  143145  CDS
				#Rosa_roxburghii1         rps12   +       142355  142583  CDS     +       143090  143145  CDS     -       71096   71209   CDS
				#Rosa_roxburghii1         rps12   -       71096   71209   CDS     +       142355  142583  CDS     +       143090  143145  CDS
				#Rosa_roxburghii1         rps12   +       71096   71209   CDS     +       99457   99512   CDS     +       100019  100247  CDS
				#Rosa_roxburghii1         rps12   +       99457   99512   CDS     +       100019  100247  CDS     +       71096   71209   CDS

				}elsif ($GENE2{$cnt2} eq "rps12") {
					my $str6=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
					my $str7=substr($seq1,($START5{$cnt2}-1),($END5{$cnt2}-$START5{$cnt2}+1));
					my $str8=substr($seq1,($START6{$cnt2}-1),($END6{$cnt2}-$START6{$cnt2}+1));
					my $rev_com6=reverse $str6;
					my $rev_com7=reverse $str7;
					my $rev_com8=reverse $str8;
					$rev_com6=~ tr/ACGTacgt/TGCAtgca/;
					$rev_com7=~ tr/ACGTacgt/TGCAtgca/;
					$rev_com8=~ tr/ACGTacgt/TGCAtgca/;
					if (($STRAND4{$cnt2} eq "-") and ($STRAND5{$cnt2} eq "-") and ($STRAND6{$cnt2} eq "-")) {
						if (($START4{$cnt2} < $START5{$cnt2}) and ($START5{$cnt2} < $START6{$cnt2})) {
							print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$rev_com6.$rev_com8.$rev_com7."\n";
						}elsif (($START4{$cnt2} < $START5{$cnt2}) and ($START6{$cnt2} < $START4{$cnt2})) {
							print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$rev_com8.$rev_com7.$rev_com6."\n";
						}
					}elsif(($STRAND4{$cnt2} eq "-") and ($STRAND5{$cnt2} eq "-") and ($STRAND6{$cnt2} eq "+")) {
						if (($START4{$cnt2} < $START5{$cnt2}) and ($START6{$cnt2} < $START4{$cnt2})) {
							print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$rev_com8.$rev_com7.$rev_com6."\n";
						}
					}elsif(($STRAND4{$cnt2} eq "+") and ($STRAND5{$cnt2} eq "-") and ($STRAND6{$cnt2} eq "-")) {
						if (($START4{$cnt2} < $START5{$cnt2}) and ($START5{$cnt2} < $START6{$cnt2})) {
							print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$str6.$rev_com8.$rev_com7."\n";
						}
					}elsif(($STRAND4{$cnt2} eq "+") and ($STRAND5{$cnt2} eq "+") and ($STRAND6{$cnt2} eq "-")) {
						if (($START4{$cnt2} < $START5{$cnt2}) and ($START6{$cnt2} < $START4{$cnt2})) {
							print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$rev_com8.$str6.$str7."\n";
						}
					}elsif(($STRAND4{$cnt2} eq "-") and ($STRAND5{$cnt2} eq "+") and ($STRAND6{$cnt2} eq "+")) {
						if (($START4{$cnt2} < $START5{$cnt2}) and ($START5{$cnt2} < $START6{$cnt2})) {
							print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$rev_com6.$str7.$str8."\n";
						}
					}elsif(($STRAND4{$cnt2} eq "+") and ($STRAND5{$cnt2} eq "+") and ($STRAND6{$cnt2} eq "+")) {
						if (($START4{$cnt2} < $START5{$cnt2}) and ($START5{$cnt2} < $START6{$cnt2})) {
							print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$str6.$str7.$str8."\n";
						}elsif (($START4{$cnt2} < $START5{$cnt2}) and ($START6{$cnt2} < $START4{$cnt2})) {
							print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$str8.$str6.$str7."\n";
						}
					}
				}
			}
			if ((defined $left_header ne "") and (defined $right_header ne "") and ($left_header eq $right_header)) {
				print $out_coding ">".$left_header."_".$SP2{$cnt2}."\n".$left_str.$right_str."\n";
			}
			$rps12_header=">rps12_".$SP2{$cnt2};
		}


		if ((@rps12_1_sequence != 0) and (@rps12_2_sequence != 0)) {
			if (@rps12_2_sequence == 1) {
				my $rps12_1_sequence=shift @rps12_1_sequence;
				my $rps12_2_sequence=shift @rps12_2_sequence;
				print $out_coding $rps12_header."\n".$rps12_1_sequence.$rps12_2_sequence."\n";
			}elsif (@rps12_2_sequence == 2) {
				my $rps12_1_sequence=shift @rps12_1_sequence;
				my $rps12_2_sequence=shift @rps12_2_sequence;
				my $rps12_3_sequence=shift @rps12_2_sequence;
				print $out_coding $rps12_header."\n".$rps12_1_sequence.$rps12_2_sequence."\n";
				print $out_coding $rps12_header."\n".$rps12_1_sequence.$rps12_3_sequence."\n";
			}
		}elsif ((@rps12_1_sequence != 0) and (@rps12_2_sequence == 0)) {
			my $rps12_1_sequence=shift @rps12_1_sequence;
			print $out_coding $rps12_header."\n".$rps12_1_sequence."\n";
			print $out_coding $rps12_header."\n".$rps12_1_sequence."\n";
		}elsif ((@rps12_1_sequence == 0) and (@rps12_2_sequence != 0)) {
			if (@rps12_2_sequence == 1) {
				my $rps12_2_sequence=shift @rps12_2_sequence;
				print $out_coding $rps12_header."\n".$rps12_2_sequence."\n";
			}elsif (@rps12_2_sequence == 2) {
				my $rps12_2_sequence=shift @rps12_2_sequence;
				my $rps12_3_sequence=shift @rps12_2_sequence;
				print $out_coding $rps12_header."\n".$rps12_2_sequence."\n";
				print $out_coding $rps12_header."\n".$rps12_3_sequence."\n";
			}
		}elsif ((@rps12_1_sequence == 0) and (@rps12_2_sequence == 0)) {
		}

		close $out_coding;
		@output1=();




		#generate_IGS_ranges
		my (%SP3,%GENE3,%STRAND7,%START7,%END7,%TYPE7,$last0,$last1,$last2,@output4);
		my $cnt3=0;
		foreach (@output3){
			$_=~ s/\r|\n//g;
			$cnt3++;
			($SP3{$cnt3},$GENE3{$cnt3},$STRAND7{$cnt3},$START7{$cnt3},$END7{$cnt3},$TYPE7{$cnt3})=(split /\s+/,$_)[0,1,2,3,4,5];
		}
		foreach (keys %SP3){
			if ($_==1 and $START7{$_}!=1){
				unshift @output4,$SP3{$_}."\t"."start".'-'.$GENE3{$_}."\t"."?"."/".$STRAND7{$_}."\t"."1"."\t".($START7{$_}-1)."\t"."?"."/".$TYPE7{$_}."\n";
			}
		}
		foreach (1..($cnt3-1)) {
			$last0=$_-1;
			$last1=$_+1;
			$last2=$_+2;
			next if ((($END7{$_}+1) > ($START7{$last1}-1)) and (($END7{$_}+1) < ($END7{$last1}-1)));
			next if (($_ > 1) and (($END7{$_}+1) < ($END7{$last0}-1)) and (($END7{$_}+1) < ($START7{$last2}-1)));
			if ((($END7{$_}+1) >= ($START7{$last1}-1)) and (($END7{$_}+1) >= ($END7{$last1}-1))){
				push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'.$GENE3{$last2}."\t".$STRAND7{$_}."/".$STRAND7{$last2}."\t".($END7{$_}+1)."\t".($START7{$last2}-1)."\t".$TYPE7{$_}."/".$TYPE7{$last2}."\n";
			}else{
				push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'.$GENE3{$last1}."\t".$STRAND7{$_}."/".$STRAND7{$last1}."\t".($END7{$_}+1)."\t".($START7{$last1}-1)."\t".$TYPE7{$_}."/".$TYPE7{$last1}."\n";
			}
		}
		foreach (keys %SP3){
			if ($_==$cnt3 and $END7{$_}!=$length){
				push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'."end"."\t".$STRAND7{$_}."/"."?"."\t".($END7{$_}+1)."\t".$length."\t".$TYPE7{$_}."/"."?"."\n";
			}
		}
		@output3=();
		#print @output4;


		#extract_IGS
		my (%SP4,%GENE4,%STRAND8,%START8,%END8,%TYPE8,$seq2);
		my $cnt4=0;
		open(my $out_noncoding,">","$output_directory/$filename\_p3_regions_among_linked_CDS_RNA.fasta");
		while (@fasta2){
			my $header=shift @fasta2;
			$seq2=shift @fasta2;
		}
		my ($left_header2,$left_str2,$right_header2,$right_str2);
		foreach (@output4){
			$_=~ s/\r|\n//g;
			$cnt4++;
			($SP4{$cnt4},$GENE4{$cnt4},$STRAND8{$cnt4},$START8{$cnt4},$END8{$cnt4},$TYPE8{$cnt4})=(split /\s+/,$_)[0,1,2,3,4,5];
			my $str=substr($seq2,($START8{$cnt4}-1),($END8{$cnt4}-$START8{$cnt4}+1));
			if (($GENE4{$cnt4} !~ /1_/) and ($GENE4{$cnt4} !~ /2_/) and ($GENE4{$cnt4} !~ /start-/) and ($GENE4{$cnt4} !~ /-end/)) {
				print $out_noncoding ">".$GENE4{$cnt4}."_".$SP4{$cnt4}."\n".$str."\n";
			}elsif ($GENE4{$cnt4} =~ /1_/) {#1_trnH-GUG-psbA
				my $temp_intergenic=$GENE4{$cnt4};
				$temp_intergenic=~ s/1_//g;
				print $out_noncoding ">".$temp_intergenic."_".$SP4{$cnt4}."\n".$str."\n";
			}elsif ($GENE4{$cnt4} =~ /2_/) {#rpl2-2-2_trnH-GUG
				my $temp_intergenic=$GENE4{$cnt4};
				$temp_intergenic=~ s/2_//g;
				print $out_noncoding ">".$temp_intergenic."_".$SP4{$cnt4}."\n".$str."\n";
			}elsif ($GENE4{$cnt4} =~ /start-(.*)/) {#start-trnH-GUG
				$right_header2=$1;
				$right_str2=$str;
			}elsif ($GENE4{$cnt4} =~ /(.*)-end/) {#rpl2-2-end
				$left_header2=$1;
				$left_str2=$str;
			}
			if ((defined $left_header2 ne "") and (defined $right_header2 ne "")) {
				print $out_noncoding ">".$left_header2."-".$right_header2."_".$SP4{$cnt4}."\n".$left_str2.$right_str2."\n";
			}
		}
		@output4=();
		close $out_noncoding;
	}
}




####extract coding and noncoding, including unlinked CDSs and tRNAs
if (($pattern == 4) or ($pattern == 6) or ($pattern == 7)) {
	while (@filenames4) {
		my $filename_gb=shift @filenames4;
		my $target_name=substr($filename_gb,0,rindex($filename_gb,"\."));
		$target_name=~ s/\s/_/g;
		my $filename=substr($target_name,rindex($target_name,"\/")+1);
		#my $target_path=$target_name;
		#$target_path=~ s/$filename//g;
		my $target_name_temp_random1 = $input_directory."/".$random."_".$filename."_temp1";
		my $target_name_temp_random2 = $input_directory."/".$random."_".$filename."_temp2";

		open(my $in_gb1,"<",$filename_gb);
		open(my $out_gb1,">",$target_name_temp_random1);
		while (<$in_gb1>){
			$_=~ s/\r\n/\n/g;
			if ($_=~ /\),\n/){
				$_=~ s/\),\n/\),/g;
			}elsif($_=~ /,\n/){
				$_=~ s/,\n/,/g;
			}
			print $out_gb1 $_;
		}
		close $in_gb1;
		close $out_gb1;

		open(my $in_gb2,"<",$target_name_temp_random1);
		open(my $out_gb2,">",$target_name_temp_random2);
		while (<$in_gb2>){
			$_=~ s/,\s+/,/g;
			print $out_gb2 $_;
		}
		close $in_gb2;
		close $out_gb2;


		#generate_bed_file
		my (@row_array,@start_end,$species_name,$length,$element,@genearray,@output1);
		my $mark=0;
		my $tick=0;
		open (my $in_gb3,"<",$target_name_temp_random2);
		while (<$in_gb3>){
			$_=~ s/\r|\n//g;
			@row_array=split /\s+/,$_;
			if (/^LOCUS/i){
				$species_name=$filename;
				$length=$row_array[2];
			}elsif(/ {5}CDS {13}/ or / {5}tRNA {12}/ or / {5}rRNA {12}/){
				if ($row_array[2]=~ /^\d+..\d+$/){
					$row_array[2]="\+\t$row_array[2]\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~/^complement\((\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
					my @array_start_end=split /\t/,$row_array[2];
					my $first_start_end=$row_array[2];
					my $second_start_end=$row_array[2];
					if ($array_start_end[1] < $array_start_end[2]) {
						$row_array[2]=$row_array[2];
					}elsif ($array_start_end[1] > $array_start_end[2]) {
						$first_start_end=~ s/$array_start_end[2]/$length/g;
						$second_start_end=~ s/$array_start_end[1]/1/g;
						push @start_end,$first_start_end;
						push @start_end,$second_start_end;
						$tick=3;
					}
				}elsif($row_array[2]=~ /^join\((\d+..\d+),(\d+..\d+)\)$/) {
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\((\d+..\d+),(\d+..\d+),(\d+..\d+)\)$/) {
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),(\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\((\d+..\d+),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="+\t$1\t$row_array[1]\t-\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^complement\(join\((\d+..\d+),(\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^complement\(join\((\d+..\d+),(\d+..\d+),(\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\((\d+..\d+),(\d+..\d+)\)$/){
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\((\d+..\d+),(\d+..\d+),(\d+..\d+)\)$/){
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),(\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),(\d+..\d+),(\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),(\d+..\d+),(\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\((\d+..\d+),(\d+..\d+),complement\((\d+..\d+)\)\)$/) {
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^<(\d+..\d+)$/){# positive split no-intron gene
					$row_array[2]="\+\t$1\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
					$tick=2;
				}elsif($row_array[2]=~ /^(\d+..>\d+)$/){# positive split no-intron gene
					$row_array[2]="\+\t$1\t$row_array[1]";
					$row_array[2]=~ s/\..>/\t/g;
					$tick=1;
				}elsif($row_array[2]=~/^complement\(<(\d+..\d+)\)$/){# negative split no-intron gene
					$row_array[2]="-\t$1\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
					$tick=1;
				}elsif($row_array[2]=~/^complement\((\d+..>\d+)\)$/){# negative split no-intron gene
					$row_array[2]="-\t$1\t$row_array[1]";
					$row_array[2]=~ s/\..>/\t/g;
					$tick=2;
				}
				$element=$row_array[2];
				$mark=1;
			}elsif(/^\s+\/gene="(.*)"/ and $mark == 1){
				if ($tick==0) {
					$element=$species_name.":".$1.":".$element;
					push @genearray,$element;
					$element=();
					$mark=0;
				}elsif ($tick==1) {
					$element=$species_name.":"."1_".$1.":".$element;#1_trnH-GUG
					push @genearray,$element;
					$element=();
					$mark=0;
					$tick=0;
				}elsif ($tick==2) {
					$element=$species_name.":"."2_".$1.":".$element;#2_trnH-GUG
					push @genearray,$element;
					$element=();
					$mark=0;
					$tick=0;
				}elsif ($tick==3) {
					push @genearray,$species_name.":"."1_".$1.":".$start_end[1];#1_trnH-GUG
					push @genearray,$species_name.":"."2_".$1.":".$start_end[0];#2_trnH-GUG
					@start_end=();
					$mark=0;
					$tick=0;
				}
			}
		}
		close $in_gb3;
		foreach (@genearray){
			my @array=split /:/,$_;
			push @output1,"$array[0]\t$array[1]\t$array[2]\n";
		}
		@row_array=();
		@genearray=();


		#put_fasta_sequence_in_array
		my $flag=0;
		my @sequence;
		my (@fas1,@fas2);
		open(my $in_gb4,"<",$target_name_temp_random2);
		while (<$in_gb4>){
			if ($_=~ /ORIGIN/){
				$flag=1;
			}
			if ($_=~ /\/\//){
				$flag=2;
			}
			if ($flag==1){
				next if ($_=~ /ORIGIN/);
				push @sequence,$_;
			}
		}
		close $in_gb4;
		foreach (@sequence){
			$_=~ s/\r|\n//g;
			$_=~ s/\s*//g;
			$_=~ s/\d+//g;
			push @fas1,$_;
		}
		my $fas1=join "",@fas1;
		my (@fasta1,@fasta2);
		push @fasta1,$species_name,$fas1;
		@fasta2=@fasta1;

		unlink ("$target_name_temp_random1");
		unlink ("$target_name_temp_random2");


		#edit_bed_file
		my (%SP1,%GENE1,%STRAND1,%START1,%END1,%TYPE1,%STRAND2,%START2,%END2,%TYPE2,%STRAND3,%START3,%END3,%TYPE3,@output2);
		my $cnt1=0;
		foreach (@output1) {
			$_=~ s/\r|\n//g;
			$cnt1++;
			($SP1{$cnt1},$GENE1{$cnt1},$STRAND1{$cnt1},$START1{$cnt1},$END1{$cnt1},$TYPE1{$cnt1},$STRAND2{$cnt1},$START2{$cnt1},$END2{$cnt1},$TYPE2{$cnt1},$STRAND3{$cnt1},$START3{$cnt1},$END3{$cnt1},$TYPE3{$cnt1})=(split /\s+/,$_)[0,1,2,3,4,5,6,7,8,9,10,11,12,13];
		}

		#foreach (1..$cnt1) {
		#	if (defined $STRAND2{$_} eq "") {
		#		push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
		#	}elsif ((defined $STRAND2{$_} ne "") and (defined $STRAND3{$_} eq "")) {
		#		if (($STRAND1{$_} eq "-") and ($START1{$_} < $START2{$_})){
		#			push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
		#			push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
		#		}elsif(($STRAND1{$_} eq "-") and ($START1{$_} > $START2{$_})){
		#			push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
		#			push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
		#		}elsif(($STRAND1{$_} eq "+") and ($START1{$_} < $START2{$_})){
		#			push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
		#			push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
		#		}elsif(($STRAND1{$_} eq "+") and ($START1{$_} > $START2{$_})){
		#			push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
		#			push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
		#		}
		#	}elsif ((defined $STRAND2{$_} ne "") and (defined $STRAND3{$_} ne "")) {
		#		if (($STRAND1{$_} eq "-") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_})){
		#			push @output2,($SP1{$_}."\t".$GENE1{$_}."-3"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
		#			push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
		#			push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
		#		}elsif(($STRAND1{$_} eq "-") and ($START1{$_} > $START2{$_}) and ($START2{$_} > $START3{$_})){
		#			push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
		#			push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
		#			push @output2,($SP1{$_}."\t".$GENE1{$_}."-3"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
		#		}elsif(($STRAND1{$_} eq "-") and ($START1{$_} > $START3{$_}) and ($START3{$_} > $START2{$_})){
		#			push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
		#			push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
		#			push @output2,($SP1{$_}."\t".$GENE1{$_}."-3"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
		#		}elsif(($STRAND1{$_} eq "+") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_})){
		#			push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
		#			push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
		#			push @output2,($SP1{$_}."\t".$GENE1{$_}."-3"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
		#		}elsif(($STRAND1{$_} eq "+") and ($START1{$_} > $START2{$_}) and ($START2{$_} > $START3{$_})){
		#			push @output2,($SP1{$_}."\t".$GENE1{$_}."-3"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
		#			push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
		#			push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
		#		}elsif(($STRAND1{$_} eq "+") and ($START1{$_} > $START3{$_}) and ($START3{$_} > $START2{$_})){
		#			push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
		#			push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
		#			push @output2,($SP1{$_}."\t".$GENE1{$_}."-3"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
		#		}
		#	}
		#}

		foreach (1..$cnt1) {
			if (defined $STRAND2{$_} eq "") {
				push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
			}elsif ((defined $STRAND2{$_} ne "") and (defined $STRAND3{$_} eq "")) {
				if ($GENE1{$_} ne "rps12") {
					if (($STRAND1{$_} eq "-") and ($START1{$_} < $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}elsif(($STRAND1{$_} eq "-") and ($START1{$_} > $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}elsif(($STRAND1{$_} eq "+") and ($START1{$_} < $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}elsif(($STRAND1{$_} eq "+") and ($START1{$_} > $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					}
				}elsif ($GENE1{$_} eq "rps12") {
					if (abs($START1{$_}-$START2{$_}) < 5000) {
						if (($STRAND1{$_} eq "-") and ($START1{$_} < $START2{$_})){
							push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
							push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						}elsif(($STRAND1{$_} eq "-") and ($START1{$_} > $START2{$_})){
							push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
							push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						}elsif(($STRAND1{$_} eq "+") and ($START1{$_} < $START2{$_})){
							push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
							push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						}elsif(($STRAND1{$_} eq "+") and ($START1{$_} > $START2{$_})){
							push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
							push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						}
					}elsif (abs($START1{$_}-$START2{$_}) > 5000) {
						if (($STRAND1{$_} eq "-") and ($START1{$_} < $START2{$_})){
							push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
							push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						}elsif(($STRAND1{$_} eq "-") and ($START1{$_} > $START2{$_})){
							push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
							push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						}elsif(($STRAND1{$_} eq "+") and ($START1{$_} < $START2{$_})){
							push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
							push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						}elsif(($STRAND1{$_} eq "+") and ($START1{$_} > $START2{$_})){
							push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
							push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						}
					}
				}
			}elsif ((defined $STRAND2{$_} ne "") and (defined $STRAND3{$_} ne "")) {
				if ($GENE1{$_} ne "rps12") {
					if (($STRAND1{$_} eq "-") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-3"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "-") and ($START1{$_} > $START2{$_}) and ($START2{$_} > $START3{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-3"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "-") and ($START1{$_} > $START3{$_}) and ($START3{$_} > $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-3"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "+") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-3"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "+") and ($START1{$_} > $START2{$_}) and ($START2{$_} > $START3{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-3"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "+") and ($START1{$_} > $START3{$_}) and ($START3{$_} > $START2{$_})){
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-3"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}
				}elsif ($GENE1{$_} eq "rps12") {
					if (($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_})){#1<2<3
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START2{$_} > $START3{$_}) and ($START3{$_} > $START1{$_})){#1<3<2
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START1{$_} > $START2{$_}) and ($START3{$_} > $START1{$_})){#2<1<3
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START1{$_} > $START3{$_}) and ($START3{$_} > $START2{$_})){#2<3<1
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START2{$_} > $START1{$_}) and ($START1{$_} > $START3{$_})){#3<1<2
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START1{$_} > $START2{$_}) and ($START2{$_} > $START3{$_})){#3<2<1
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});

					}elsif (($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "+") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_})){#1<2<3
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "+") and ($START2{$_} > $START3{$_}) and ($START3{$_} > $START1{$_})){#1<3<2
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "+") and ($START1{$_} > $START2{$_}) and ($START3{$_} > $START1{$_})){#2<1<3
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "+") and ($START1{$_} > $START3{$_}) and ($START3{$_} > $START2{$_})){#2<3<1
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "+") and ($START2{$_} > $START1{$_}) and ($START1{$_} > $START3{$_})){#3<1<2
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "+") and ($START1{$_} > $START2{$_}) and ($START2{$_} > $START3{$_})){#3<2<1
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});

					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_})){#1<2<3
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START1{$_} < $START3{$_}) and ($START3{$_} < $START2{$_})){#1<3<2
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START2{$_} < $START1{$_}) and ($START1{$_} < $START3{$_})){#2<1<3
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START1{$_} > $START3{$_}) and ($START3{$_} > $START2{$_})){#2<3<1
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START2{$_} > $START1{$_}) and ($START1{$_} > $START3{$_})){#3<1<2
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START1{$_} > $START2{$_}) and ($START2{$_} > $START3{$_})){#3<2<1
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});

					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "-") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_})){#1<2<3
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "-") and ($START1{$_} < $START3{$_}) and ($START3{$_} < $START2{$_})){#1<3<2
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "-") and ($START2{$_} < $START1{$_}) and ($START1{$_} < $START3{$_})){#2<1<3
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "-") and ($START1{$_} > $START3{$_}) and ($START3{$_} > $START2{$_})){#2<3<1
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "-") and ($START2{$_} > $START1{$_}) and ($START1{$_} > $START3{$_})){#3<1<2
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "-") and ($START1{$_} > $START2{$_}) and ($START2{$_} > $START3{$_})){#3<2<1
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
						push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
					}
				}
			}
		}
		@output1=();


		#sort_bed_file
		my $col=3;
		my (%sort,@output3);
		foreach (@output2){
			my @row=split /\t/,$_;
			$sort{$_}=$row[$col];
		}
		foreach (sort {$sort{$a} <=> $sort{$b}} keys %sort){
			push @output3,"$_\n";
		}
		@output2=();


		#output_bed_file
		#open (my $out_bed,">","$output_directory/$filename\_p4_unlinked_CDS_RNA.bed");
		#foreach (@output3){
		#	print $out_bed $_;
		#}
		#close $out_bed;


		#output_bed_file
		open (my $out_bed,">","$output_directory/$filename\_p4_unlinked_CDS_RNA.bed");
		foreach (@output3){
			my @row=split /\t/,$_;
			my $abs=abs($row[3]-$row[4]);

			if (($row[1] ne "rps12") and ($row[1] ne "rps12-1") and ($row[1] ne "rps12-2")) {
				print $out_bed $_;
			}
			if (($row[1] eq "rps12") and ($abs < 200)) {
				$_=~ s/rps12/rps12-1/;
				print $out_bed $_;
			}
			if (($row[1] eq "rps12") and ($abs > 200)) {
				$_=~ s/rps12/rps12-2/;
				print $out_bed $_;
			}
			if ($row[1] eq "rps12-1") {
				$_=~ s/rps12-1/rps12-2/;
				print $out_bed $_;
			}
			if ($row[1] eq "rps12-2") {
				$_=~ s/rps12-2/rps12-3/;
				print $out_bed $_;
			}
		}
		close $out_bed;


		#extract_coding_regions
		my (%SP2,%GENE2,%STRAND4,%START4,%END4,%TYPE4,$seq1);
		my $cnt2=0;
		open(my $out_coding,">","$output_directory/$filename\_p4_unlinked_CDS_RNA.fasta");
		while (@fasta1){
			my $header=shift @fasta1;
			$seq1=shift @fasta1;
		}
		my ($left_header,$left_str,$right_header,$right_str);
		foreach (@output3){
			$_=~ s/\r|\n//g;
			$cnt2++;
			($SP2{$cnt2},$GENE2{$cnt2},$STRAND4{$cnt2},$START4{$cnt2},$END4{$cnt2},$TYPE4{$cnt2})=(split /\s+/,$_)[0,1,2,3,4,5];
			my $str=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
			if ($STRAND4{$cnt2} eq "-") {
				my $rev_com=reverse $str;
				$rev_com=~ tr/ACGTacgt/TGCAtgca/;
				if (($GENE2{$cnt2} !~ /1_/) and ($GENE2{$cnt2} !~ /2_/)) {
					print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$rev_com."\n";
				}elsif ($GENE2{$cnt2} =~ /1_(.*)/) {#1_trnH-GUG
					$left_header=$1;
					$left_str=$rev_com;
				}elsif ($GENE2{$cnt2} =~ /2_(.*)/) {#2_trnH-GUG
					$right_header=$1;
					$right_str=$rev_com;
				}
			}elsif($STRAND4{$cnt2} eq "+"){
				if (($GENE2{$cnt2} !~ /1_/) and ($GENE2{$cnt2} !~ /2_/)) {
					print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$str."\n";
				}elsif ($GENE2{$cnt2} =~ /1_(.*)/) {#1_trnH-GUG
					$left_header=$1;
					$left_str=$str;
				}elsif ($GENE2{$cnt2} =~ /2_(.*)/) {#2_trnH-GUG
					$right_header=$1;
					$right_str=$str;
				}
			}
			if ((defined $left_header ne "") and (defined $right_header ne "") and ($left_header eq $right_header)) {
				print $out_coding ">".$left_header."_".$SP2{$cnt2}."\n".$left_str.$right_str."\n";
			}
		}
		close $out_coding;


		#generate_noncoding_ranges
		my (%SP3,%GENE3,%STRAND5,%START5,%END5,%TYPE5,$last,@output4);
		my $cnt3=0;
		foreach (@output3){
			$_=~ s/\r|\n//g;
			$cnt3++;
			($SP3{$cnt3},$GENE3{$cnt3},$STRAND5{$cnt3},$START5{$cnt3},$END5{$cnt3},$TYPE5{$cnt3})=(split /\s+/,$_)[0,1,2,3,4,5];
		}
		foreach (keys %SP3){
			if ($_==1 and $START5{$_}!=1){
				unshift @output4,$SP3{$_}."\t"."start".'-'.$GENE3{$_}."\t"."?"."/".$STRAND5{$_}."\t"."1"."\t".($START5{$_}-1)."\t"."?"."/".$TYPE5{$_}."\n";
			}
		}
		foreach (1..($cnt3-1)) {
			$last=$_+1;
			next if ($END5{$_}+1) > ($START5{$last}-1);
			push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'.$GENE3{$last}."\t".$STRAND5{$_}."/".$STRAND5{$last}."\t".($END5{$_}+1)."\t".($START5{$last}-1)."\t".$TYPE5{$_}."/".$TYPE5{$last}."\n";
		}
		foreach (keys %SP3){
			if ($_==$cnt3 and $END5{$_}!=$length){
				push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'."end"."\t".$STRAND5{$_}."/"."?"."\t".($END5{$_}+1)."\t".$length."\t".$TYPE5{$_}."/"."?"."\n";
			}
		}
		@output3=();
		#print @output4;


		#extract_noncoding_regions
		my (%SP4,%GENE4,%STRAND6,%START6,%END6,%TYPE6,$seq2);
		my $cnt4=0;
		open(my $out_noncoding,">","$output_directory/$filename\_p4_regions_among_unlinked_CDS_RNA.fasta");
		while (@fasta2){
			my $header=shift @fasta2;
			$seq2=shift @fasta2;
		}
		my ($left_header2,$left_str2,$right_header2,$right_str2);
		foreach (@output4){
			$_=~ s/\r|\n//g;
			$cnt4++;
			($SP4{$cnt4},$GENE4{$cnt4},$STRAND6{$cnt4},$START6{$cnt4},$END6{$cnt4},$TYPE6{$cnt4})=(split /\s+/,$_)[0,1,2,3,4,5];
			my $str=substr($seq2,($START6{$cnt4}-1),($END6{$cnt4}-$START6{$cnt4}+1));
			if (($GENE4{$cnt4} !~ /1_/) and ($GENE4{$cnt4} !~ /2_/) and ($GENE4{$cnt4} !~ /start-/) and ($GENE4{$cnt4} !~ /-end/)) {
				print $out_noncoding ">".$GENE4{$cnt4}."_".$SP4{$cnt4}."\n".$str."\n";
			}elsif ($GENE4{$cnt4} =~ /1_/) {#1_trnH-GUG-psbA
				my $temp_intergenic=$GENE4{$cnt4};
				$temp_intergenic=~ s/1_//g;
				print $out_noncoding ">".$temp_intergenic."_".$SP4{$cnt4}."\n".$str."\n";
			}elsif ($GENE4{$cnt4} =~ /2_/) {#rpl2-2_trnH-GUG
				my $temp_intergenic=$GENE4{$cnt4};
				$temp_intergenic=~ s/2_//g;
				print $out_noncoding ">".$temp_intergenic."_".$SP4{$cnt4}."\n".$str."\n";
			}elsif ($GENE4{$cnt4} =~ /start-(.*)/) {#start-trnH-GUG
				$right_header2=$1;
				$right_str2=$str;
			}elsif ($GENE4{$cnt4} =~ /(.*)-end/) {#rpl2-2-end
				$left_header2=$1;
				$left_str2=$str;
			}
			if ((defined $left_header2 ne "") and (defined $right_header2 ne "")) {
				print $out_noncoding ">".$left_header2."-".$right_header2."_".$SP4{$cnt4}."\n".$left_str2.$right_str2."\n";
			}
		}
		@output4=();
		close $out_noncoding;
	}
}




####function
sub argument{
	my @options=("help|h","input|i:s","pattern|p:i","output|o:s");
	my %options;
	GetOptions(\%options,@options);
	exec ("pod2usage $0") if ((keys %options)==0 or $options{'h'} or $options{'help'});
	if(!exists $options{'input'}){
		print "***ERROR: No input directory is assigned!!!\n";
		exec ("pod2usage $0");
	}elsif(!exists $options{'pattern'}){
		print "***ERROR: No extraction pattern is assigned!!!\n";
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

    extraction_v1.pl version 1.0

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

    extract coding and noncoding sequences from .gb/.gbf/.gbk file(s)

=head1 SYNOPSIS

    extraction_v1.pl -i -p -o
    Copyright (C) 2025 Xiao-Jian Qu
    Please contact <quxiaojian@sdnu.edu.cn>, if you have any bugs or questions.

    [-h -help]    help information.
    [-i -input]   required: (default: input) input directory containing GenBank-formatted files.
    [-p -pattern] required: (default: 3) 1/2/3/4/5/6/7 extraction pattern, 1=PCGs with introns,
                  and sequences among PCGs, 2=genes with introns, and intergenic sequences,
                  3=coding sequences with linked CDSs and tRNAs, and noncoding sequences,
                  4=coding sequences with unlinked CDSs and tRNAs, and noncoding sequences,
                  5=1+2, 6=3+4, 7=1+2+3+4.
    [-o -output]  required: (default: output) output directory.

=cut
