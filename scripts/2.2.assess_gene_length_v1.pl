#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Find;
use Data::Dumper;
$|=1;

my $global_options=&argument();
my $input_reference=&default("reference","reference");
my $input_target=&default("target","target");
my $output_directory=&default("outptu","output");

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


#PCGs and regions among PCGs
#genes and regions among genes
#linked CDSs and RNAs and regions among linked CDSs and RNAs
#unlinked CDSs and RNAs and regions among unlinked CDSs and RNAs




####extract coding and noncoding, including linked CDSs and tRNAs
my $pattern1=".gb";
my $pattern2=".gbf";
my $pattern3=".gbk";
my (@filenames1,@filenames2);
find(\&target1,$input_reference);
sub target1{
	if (/$pattern1/ or /$pattern2/ or /$pattern3/){
        push @filenames1,"$File::Find::name";
    }
    return;
}
@filenames2=@filenames1;

my $filename_reference;
while (@filenames1) {
	my $filename_gb=shift @filenames1;
	my $target_name=$filename_gb;#target/Amborella_trichopoda.gb
	$target_name=substr($target_name,0,rindex($target_name,"\."));#target/Amborella_trichopoda
	$target_name=~ s/\s/_/g;
	my $filename=substr($target_name,rindex($target_name,"\/")+1);#Amborella_trichopoda
	$filename_reference=$filename;
	my $target_path=$target_name;
	#$target_path=~ s/$filename//g;
	$target_path=substr($target_name,0,rindex($target_name,"\/"));#target

#{
	#my $target_name=substr($input_reference,0,rindex($input_reference,"\."));
	#$target_name=~ s/\s/_/g;
	#my $filename=substr($target_name,rindex($target_name,"\/")+1);
	#$filename_reference=$filename;
	#my $target_path=$target_name;
	##$target_path=~ s/$filename//g;
	#$target_path=substr($target_name,0,rindex($target_name,"\/"));

	open(my $in_gb1,"<",$filename_gb);
	open(my $out_gb1,">","$target_path\_temp1");
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

	open(my $in_gb2,"<","$target_path\_temp1");
	open(my $out_gb2,">","$target_path\_temp2");
	while (<$in_gb2>){
		$_=~ s/,\s+/,/g;
		print $out_gb2 $_;
	}
	close $in_gb2;
	close $out_gb2;


	#generate_bed_file
	my (@row_array,$species_name,$length,$element,@genearray,@output1);
	my $mark=0;
	my $tick=0;
	open (my $in_gb3,"<","$target_path\_temp2");
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
	open(my $in_gb4,"<","$target_path\_temp2");
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

	unlink "$target_path\_temp1";
	unlink "$target_path\_temp2";


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
	#open (my $out_bed,">","$filename.bed");
	#foreach (@output3){
	#	print $out_bed $_;
	#}
	#close $out_bed;


	#output_bed_file
	open (my $out_bed,">","$target_path\\$filename.bed");
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
	unlink("$target_path\\$filename.bed");


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
	open(my $out_coding,">","$target_path\\$filename.fasta");
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




	##generate_IGS_ranges
	#my (%SP3,%GENE3,%STRAND7,%START7,%END7,%TYPE7,$last0,$last1,$last2,@output4);
	#my $cnt3=0;
	#foreach (@output3){
	#	$_=~ s/\r|\n//g;
	#	$cnt3++;
	#	($SP3{$cnt3},$GENE3{$cnt3},$STRAND7{$cnt3},$START7{$cnt3},$END7{$cnt3},$TYPE7{$cnt3})=(split /\s+/,$_)[0,1,2,3,4,5];
	#}
	#foreach (keys %SP3){
	#	if ($_==1 and $START7{$_}!=1){
	#		unshift @output4,$SP3{$_}."\t"."start".'-'.$GENE3{$_}."\t"."?"."/".$STRAND7{$_}."\t"."1"."\t".($START7{$_}-1)."\t"."?"."/".$TYPE7{$_}."\n";
	#	}
	#}
	#foreach (1..($cnt3-1)) {
	#	$last0=$_-1;
	#	$last1=$_+1;
	#	$last2=$_+2;
	#	next if ((($END7{$_}+1) >= ($START7{$last1}-1)) and (($END7{$_}+1) < ($END7{$last1}-1)));
	#	next if (($_ > 1) and (($END7{$_}+1) < ($END7{$last0}-1)) and (($END7{$_}+1) < ($START7{$last2}-1)));
	#	if ((($END7{$_}+1) >= ($START7{$last1}-1)) and (($END7{$_}+1) >= ($END7{$last1}-1))){
	#		push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'.$GENE3{$last2}."\t".$STRAND7{$_}."/".$STRAND7{$last2}."\t".($END7{$_}+1)."\t".($START7{$last2}-1)."\t".$TYPE7{$_}."/".$TYPE7{$last2}."\n";
	#	}else{
	#		push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'.$GENE3{$last1}."\t".$STRAND7{$_}."/".$STRAND7{$last1}."\t".($END7{$_}+1)."\t".($START7{$last1}-1)."\t".$TYPE7{$_}."/".$TYPE7{$last1}."\n";
	#	}
	#}
	#foreach (keys %SP3){
	#	if ($_==$cnt3 and $END7{$_}!=$length){
	#		push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'."end"."\t".$STRAND7{$_}."/"."?"."\t".($END7{$_}+1)."\t".$length."\t".$TYPE7{$_}."/"."?"."\n";
	#	}
	#}
	#@output3=();
	##print @output4;
	#
	#
	##extract_IGS
	#my (%SP4,%GENE4,%STRAND8,%START8,%END8,%TYPE8,$seq2);
	#my $cnt4=0;
	#open(my $out_noncoding,">","$filename\_regions_among_linked_CDS_RNA.fasta");
	#while (@fasta2){
	#	my $header=shift @fasta2;
	#	$seq2=shift @fasta2;
	#}
	#my ($left_header2,$left_str2,$right_header2,$right_str2);
	#foreach (@output4){
	#	$_=~ s/\r|\n//g;
	#	$cnt4++;
	#	($SP4{$cnt4},$GENE4{$cnt4},$STRAND8{$cnt4},$START8{$cnt4},$END8{$cnt4},$TYPE8{$cnt4})=(split /\s+/,$_)[0,1,2,3,4,5];
	#	my $str=substr($seq2,($START8{$cnt4}-1),($END8{$cnt4}-$START8{$cnt4}+1));
	#	if (($GENE4{$cnt4} !~ /1_/) and ($GENE4{$cnt4} !~ /2_/) and ($GENE4{$cnt4} !~ /start-/) and ($GENE4{$cnt4} !~ /-end/)) {
	#		print $out_noncoding ">".$GENE4{$cnt4}."_".$SP4{$cnt4}."\n".$str."\n";
	#	}elsif ($GENE4{$cnt4} =~ /1_/) {#1_trnH-GUG-psbA
	#		my $temp_intergenic=$GENE4{$cnt4};
	#		$temp_intergenic=~ s/1_//g;
	#		print $out_noncoding ">".$temp_intergenic."_".$SP4{$cnt4}."\n".$str."\n";
	#	}elsif ($GENE4{$cnt4} =~ /2_/) {#rpl2-2-2_trnH-GUG
	#		my $temp_intergenic=$GENE4{$cnt4};
	#		$temp_intergenic=~ s/2_//g;
	#		print $out_noncoding ">".$temp_intergenic."_".$SP4{$cnt4}."\n".$str."\n";
	#	}elsif ($GENE4{$cnt4} =~ /start-(.*)/) {#start-trnH-GUG
	#		$right_header2=$1;
	#		$right_str2=$str;
	#	}elsif ($GENE4{$cnt4} =~ /(.*)-end/) {#rpl2-2-end
	#		$left_header2=$1;
	#		$left_str2=$str;
	#	}
	#	if ((defined $left_header2 ne "") and (defined $right_header2 ne "")) {
	#		print $out_noncoding ">".$left_header2."-".$right_header2."_".$SP4{$cnt4}."\n".$left_str2.$right_str2."\n";
	#	}
	#}
	#@output4=();
	#close $out_noncoding;
}




my $pattern4=".gb";
my $pattern5=".gbf";
my $pattern6=".gbk";
my (@filenames3,@filenames4);
find(\&target2,$input_target);
sub target2{
	if (/$pattern4/ or /$pattern5/ or /$pattern6/){
        push @filenames3,"$File::Find::name";
    }
    return;
}
@filenames4=@filenames3;

while (@filenames3) {
	my $filename_gb=shift @filenames3;
	my $target_name=$filename_gb;#target/Amborella_trichopoda2.gb #target/Rosa_roxburghii.gb
	$target_name=substr($target_name,0,rindex($target_name,"\."));#target/Amborella_trichopoda2 #target/Rosa_roxburghii
	$target_name=~ s/\s/_/g;
	my $filename=substr($target_name,rindex($target_name,"\/")+1);#Amborella_trichopoda2 #Rosa_roxburghii
	my $target_path=$target_name;
	#$target_path=~ s/$filename//g;
	$target_path=substr($target_name,0,rindex($target_name,"\/"));#target

#{
	#my $target_name=substr($input_target,0,rindex($input_target,"\."));
	#$target_name=~ s/\s/_/g;
	#my $filename=substr($target_name,rindex($target_name,"\\")+1);
	#my $target_path=$target_name;
	##$target_path=~ s/$filename//g;
	#$target_path=substr($target_name,0,rindex($target_name,"\\"));

	open(my $in_gb1,"<",$filename_gb);
	open(my $out_gb1,">","$target_path\_temp1");
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

	open(my $in_gb2,"<","$target_path\_temp1");
	open(my $out_gb2,">","$target_path\_temp2");
	while (<$in_gb2>){
		$_=~ s/,\s+/,/g;
		print $out_gb2 $_;
	}
	close $in_gb2;
	close $out_gb2;


	#generate_bed_file
	my (@row_array,$species_name,$length,$element,@genearray,@output1);
	my $mark=0;
	my $tick=0;
	open (my $in_gb3,"<","$target_path\_temp2");
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
	open(my $in_gb4,"<","$target_path\_temp2");
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

	unlink "$target_path\_temp1";
	unlink "$target_path\_temp2";


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
	#open (my $out_bed,">","$filename.bed");
	#foreach (@output3){
	#	print $out_bed $_;
	#}
	#close $out_bed;


	#output_bed_file
	open (my $out_bed,">","$target_path\\$filename.bed");
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
	unlink("$target_path\\$filename.bed");


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
	open(my $out_coding,">","$target_path\\$filename.fasta");
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




	##generate_IGS_ranges
	#my (%SP3,%GENE3,%STRAND7,%START7,%END7,%TYPE7,$last0,$last1,$last2,@output4);
	#my $cnt3=0;
	#foreach (@output3){
	#	$_=~ s/\r|\n//g;
	#	$cnt3++;
	#	($SP3{$cnt3},$GENE3{$cnt3},$STRAND7{$cnt3},$START7{$cnt3},$END7{$cnt3},$TYPE7{$cnt3})=(split /\s+/,$_)[0,1,2,3,4,5];
	#}
	#foreach (keys %SP3){
	#	if ($_==1 and $START7{$_}!=1){
	#		unshift @output4,$SP3{$_}."\t"."start".'-'.$GENE3{$_}."\t"."?"."/".$STRAND7{$_}."\t"."1"."\t".($START7{$_}-1)."\t"."?"."/".$TYPE7{$_}."\n";
	#	}
	#}
	#foreach (1..($cnt3-1)) {
	#	$last0=$_-1;
	#	$last1=$_+1;
	#	$last2=$_+2;
	#	next if ((($END7{$_}+1) >= ($START7{$last1}-1)) and (($END7{$_}+1) < ($END7{$last1}-1)));
	#	next if (($_ > 1) and (($END7{$_}+1) < ($END7{$last0}-1)) and (($END7{$_}+1) < ($START7{$last2}-1)));
	#	if ((($END7{$_}+1) >= ($START7{$last1}-1)) and (($END7{$_}+1) >= ($END7{$last1}-1))){
	#		push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'.$GENE3{$last2}."\t".$STRAND7{$_}."/".$STRAND7{$last2}."\t".($END7{$_}+1)."\t".($START7{$last2}-1)."\t".$TYPE7{$_}."/".$TYPE7{$last2}."\n";
	#	}else{
	#		push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'.$GENE3{$last1}."\t".$STRAND7{$_}."/".$STRAND7{$last1}."\t".($END7{$_}+1)."\t".($START7{$last1}-1)."\t".$TYPE7{$_}."/".$TYPE7{$last1}."\n";
	#	}
	#}
	#foreach (keys %SP3){
	#	if ($_==$cnt3 and $END7{$_}!=$length){
	#		push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'."end"."\t".$STRAND7{$_}."/"."?"."\t".($END7{$_}+1)."\t".$length."\t".$TYPE7{$_}."/"."?"."\n";
	#	}
	#}
	#@output3=();
	##print @output4;
	#
	#
	##extract_IGS
	#my (%SP4,%GENE4,%STRAND8,%START8,%END8,%TYPE8,$seq2);
	#my $cnt4=0;
	#open(my $out_noncoding,">","$filename\_regions_among_linked_CDS_RNA.fasta");
	#while (@fasta2){
	#	my $header=shift @fasta2;
	#	$seq2=shift @fasta2;
	#}
	#my ($left_header2,$left_str2,$right_header2,$right_str2);
	#foreach (@output4){
	#	$_=~ s/\r|\n//g;
	#	$cnt4++;
	#	($SP4{$cnt4},$GENE4{$cnt4},$STRAND8{$cnt4},$START8{$cnt4},$END8{$cnt4},$TYPE8{$cnt4})=(split /\s+/,$_)[0,1,2,3,4,5];
	#	my $str=substr($seq2,($START8{$cnt4}-1),($END8{$cnt4}-$START8{$cnt4}+1));
	#	if (($GENE4{$cnt4} !~ /1_/) and ($GENE4{$cnt4} !~ /2_/) and ($GENE4{$cnt4} !~ /start-/) and ($GENE4{$cnt4} !~ /-end/)) {
	#		print $out_noncoding ">".$GENE4{$cnt4}."_".$SP4{$cnt4}."\n".$str."\n";
	#	}elsif ($GENE4{$cnt4} =~ /1_/) {#1_trnH-GUG-psbA
	#		my $temp_intergenic=$GENE4{$cnt4};
	#		$temp_intergenic=~ s/1_//g;
	#		print $out_noncoding ">".$temp_intergenic."_".$SP4{$cnt4}."\n".$str."\n";
	#	}elsif ($GENE4{$cnt4} =~ /2_/) {#rpl2-2-2_trnH-GUG
	#		my $temp_intergenic=$GENE4{$cnt4};
	#		$temp_intergenic=~ s/2_//g;
	#		print $out_noncoding ">".$temp_intergenic."_".$SP4{$cnt4}."\n".$str."\n";
	#	}elsif ($GENE4{$cnt4} =~ /start-(.*)/) {#start-trnH-GUG
	#		$right_header2=$1;
	#		$right_str2=$str;
	#	}elsif ($GENE4{$cnt4} =~ /(.*)-end/) {#rpl2-2-end
	#		$left_header2=$1;
	#		$left_str2=$str;
	#	}
	#	if ((defined $left_header2 ne "") and (defined $right_header2 ne "")) {
	#		print $out_noncoding ">".$left_header2."-".$right_header2."_".$SP4{$cnt4}."\n".$left_str2.$right_str2."\n";
	#	}
	#}
	#@output4=();
	#close $out_noncoding;
}




####
my (%rrn_reference,%trn_reference,%pcg_reference);
while (@filenames2) {
	my $filename_gb=shift @filenames2;
	my $reference_name=$filename_gb;
	$reference_name=substr($reference_name,0,rindex($reference_name,"\."));
	$reference_name=~ s/\s/_/g;
	#my $filename=substr($reference_name,rindex($reference_name,"\/")+1);

	#my $reference_name=substr($input_reference,0,rindex($input_reference,"\."));
	#$reference_name=~ s/\s/_/g;
	#my (%rrn_reference,%trn_reference,%pcg_reference);
	open(my $in_fasta_reference,"<","$reference_name.fasta");
	my ($header,$sequence);
	while (defined ($header=<$in_fasta_reference>) and defined ($sequence=<$in_fasta_reference>)) {
		$header=~ s/\r|\n//g;
		$sequence=~ s/\r|\n//g;
		my $length=length $sequence;
		if ($header=~ /^>(rrn.*?)(_.*)/){
			my $genes=$1;
			$rrn_reference{$genes}=$length;
		}elsif ($header=~ /^>(trn.*?)(_.*)/){
			my $genes=$1;
			$trn_reference{$genes}=$length;
		}elsif ($header=~ /^>(?!trn)(.*?)(_.*)/ and $header=~ /^>(?!rrn)(.*?)(_.*)/){
			my $genes=$1;
			$genes=~ s/rps12\-1/rps12/;
			$genes=~ s/rps12\-2/rps12/;
			$pcg_reference{$genes}=$length;
		}
	}
	close $in_fasta_reference;
	unlink("$reference_name.fasta");
}




while (@filenames4) {
	my $filename_gb=shift @filenames4;
	my $target_name=$filename_gb;
	$target_name=substr($target_name,0,rindex($target_name,"\."));
	$target_name=~ s/\s/_/g;
	my $filename=substr($target_name,rindex($target_name,"\/")+1);

	#my $target_name=substr($input_target,0,rindex($input_target,"\."));
	#$target_name=~ s/\s/_/g;
	my (%rrn_target,%trn_target,%pcg_target);
	open(my $in_fasta_target,"<","$target_name.fasta");
	my ($header2,$sequence2);
	while (defined ($header2=<$in_fasta_target>) and defined ($sequence2=<$in_fasta_target>)) {
		$header2=~ s/\r|\n//g;
		$sequence2=~ s/\r|\n//g;
		my $length=length $sequence2;
		if ($header2=~ /^>(rrn.*?)(_.*)/){
			my $genes=$1;
			$rrn_target{$genes}=$length;
		}elsif ($header2=~ /^>(trn.*?)(_.*)/){
			my $genes=$1;
			$trn_target{$genes}=$length;
		}elsif ($header2=~ /^>(?!trn)(.*?)(_.*)/ and $header2=~ /^>(?!rrn)(.*?)(_.*)/){
			my $genes=$1;
			$genes=~ s/rps12\-1/rps12/;
			$genes=~ s/rps12\-2/rps12/;
			$pcg_target{$genes}=$length;
		}
	}
	close $in_fasta_target;
	unlink("$target_name.fasta");




	####output the statistic results
	open(my $assessment,">","$output_directory/$filename_reference\_$filename.txt");
	print $assessment "type"."\t"."gene_name"."\t"."reference_length"."\t"."target_length"."\t"."length_difference"."\n";
	my @array_PCGs;
	foreach my $key (keys %pcg_reference) {
		if (exists $pcg_target{$key}) {
			my $difference=$pcg_target{$key}-$pcg_reference{$key};
			push @array_PCGs,(["PCGs",$key,$pcg_reference{$key},$pcg_target{$key},$difference]);
			#print $assessment "PCGs"."\t".$key."\t\t".$pcg_reference{$key}."\t".$pcg_target{$key}."\t".$difference."\n";
		}elsif (!exists $pcg_target{$key}) {
			my $difference=0-$pcg_reference{$key};
			push @array_PCGs,(["PCGs",$key,$pcg_reference{$key},"NA",$difference]);
			#print $assessment "PCGs"."\t".$key."\t\t".$pcg_reference{$key}."\t"."NA"."\t".$difference."\n";
		}
	}
	foreach my $key (keys %pcg_target) {
		if (!exists $pcg_reference{$key}) {
			my $difference=$pcg_target{$key}-0;
			push @array_PCGs,(["PCGs",$key,"NA",$pcg_target{$key},$difference]);
			#print $assessment "PCGs"."\t".$key."\t\t"."NA"."\t".$pcg_target{$key}."\t".$difference."\n";
		}
	}
	foreach my $item (sort {abs($b->[4]) <=> abs($a->[4]) or $a->[1] cmp $b->[1]} @array_PCGs){
		print $assessment "$item->[0]\t$item->[1]\t$item->[2]\t$item->[3]\t$item->[4]\n";
	}

	my @array_tRNAs;
	foreach my $key (keys %trn_reference) {
		if (exists $trn_target{$key}) {
			my $difference=$trn_target{$key}-$trn_reference{$key};
			push @array_tRNAs,(["tRNAs",$key,$trn_reference{$key},$trn_target{$key},$difference]);
			#print $assessment "tRNAs"."\t".$key."\t".$trn_reference{$key}."\t".$trn_target{$key}."\t".$difference."\n";
		}elsif (!exists $trn_target{$key}) {
			my $difference=0-$trn_reference{$key};
			push @array_tRNAs,(["tRNAs",$key,$trn_reference{$key},"NA",$difference]);
			#print $assessment "tRNAs"."\t".$key."\t".$trn_reference{$key}."\t"."NA"."\t".$difference."\n";
		}
	}
	foreach my $key (keys %trn_target) {
		if (!exists $trn_reference{$key}) {
			my $difference=$trn_target{$key}-0;
			push @array_tRNAs,(["tRNAs",$key,"NA",$trn_target{$key},$difference]);
			#print $assessment "tRNAs"."\t".$key."\t"."NA"."\t".$trn_target{$key}."\t".$difference."\n";
		}
	}
	foreach my $item (sort {abs($b->[4]) <=> abs($a->[4]) or $a->[1] cmp $b->[1]} @array_tRNAs){
		print $assessment "$item->[0]\t$item->[1]\t$item->[2]\t$item->[3]\t$item->[4]\n";
	}

	my @array_rRNAs;
	foreach my $key (keys %rrn_reference) {
		if (exists $rrn_target{$key}) {
			my $difference=$rrn_target{$key}-$rrn_reference{$key};
			push @array_rRNAs,(["rRNAs",$key,$rrn_reference{$key},$rrn_target{$key},$difference]);
			#print $assessment "rRNAs"."\t".$key."\t".$rrn_reference{$key}."\t".$rrn_target{$key}."\t".$difference."\n";
		}elsif (!exists $rrn_target{$key}) {
			my $difference=0-$rrn_reference{$key};
			push @array_rRNAs,(["rRNAs",$key,$rrn_reference{$key},"NA",$difference]);
			#print $assessment "rRNAs"."\t".$key."\t".$rrn_reference{$key}."\t"."NA"."\t".$difference."\n";
		}
	}
	foreach my $key (keys %rrn_target) {
		if (!exists $rrn_reference{$key}) {
			my $difference=$rrn_target{$key}-0;
			push @array_rRNAs,(["rRNAs",$key,"NA",$rrn_target{$key},$difference]);
			#print $assessment "rRNAs"."\t".$key."\t"."NA"."\t".$rrn_target{$key}."\t".$difference."\n";
		}
	}
	foreach my $item (sort {abs($b->[4]) <=> abs($a->[4]) or $a->[1] cmp $b->[1]} @array_rRNAs){
		print $assessment "$item->[0]\t$item->[1]\t$item->[2]\t$item->[3]\t$item->[4]\n";
	}

	close $assessment;
}








####function
sub argument{
	my @options=("help|h","reference|r:s","target|t:s","output|o:s");
	my %options;
	GetOptions(\%options,@options);
	exec ("pod2usage $0") if ((keys %options)==0 or $options{'h'} or $options{'help'});
	if(!exists $options{'reference'}){
		print "***ERROR: No reference directory is assigned!!!\n";
		exec ("pod2usage $0");
	}elsif(!exists $options{'target'}){
		print "***ERROR: No target directory is assigned!!!\n";
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

    assess_gene_length_v1.pl version 1.0

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

    assess gene length between reference and target files

=head1 SYNOPSIS

    assess_gene_length_v1.pl -r -t -o
    Copyright (C) 2025 Xiao-Jian Qu
    Please contact <quxiaojian@sdnu.edu.cn>, if you have any bugs or questions.

    [-h -help]        help information.
    [-r -reference]   required: (default: reference) input reference directory containing one GenBank-formatted file.
    [-t -target]      required: (default: target) input target directory containing GenBank-formatted file(s).
    [-o -output]      required: (default: output) output directory containing TXT-format file(s) with statistic result.

=cut
