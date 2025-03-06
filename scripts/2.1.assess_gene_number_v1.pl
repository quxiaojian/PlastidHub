#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Find;
use Data::Dumper;
$|=1;

my $global_options=&argument();
my $input_reference=&default("reference","reference");
my $input_target=&default("target","target");
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


#PCGs and regions among PCGs
#genes and regions among genes
#linked CDSs and RNAs and regions among linked CDSs and RNAs
#unlinked CDSs and RNAs and regions among unlinked CDSs and RNAs




####extract genes with introns and ingergenic sequences
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
		}elsif(/ {5}gene {12}/){
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
	#open (my $out_bed,">","$filename\.bed");
	#foreach (@output3){
	#	print $out_bed $_;
	#}
	#close $out_bed;

	#output_bed_file
	open (my $out_bed,">","$target_path\\$filename\.bed");
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


	##extract_gene
	#my (%SP2,%GENE2,%STRAND4,%START4,%END4,%TYPE4,$seq1);
	#my $cnt2=0;
	#open(my $out_coding,">","$filename\_genes.fasta");
	#while (@fasta1){
	#	my $header=shift @fasta1;
	#	$seq1=shift @fasta1;
	#}
	#
	#my ($rps12_header,@rps12_1_sequence,@rps12_2_sequence,@rps12_3_sequence);
	#my ($left_header,$left_str,$right_header,$right_str);
	#foreach (@output3){
	#	$_=~ s/\r|\n//g;
	#	$cnt2++;
	#	($SP2{$cnt2},$GENE2{$cnt2},$STRAND4{$cnt2},$START4{$cnt2},$END4{$cnt2},$TYPE4{$cnt2})=(split /\s+/,$_)[0,1,2,3,4,5];
	#	my $str=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
	#	#if ($STRAND4{$cnt2} eq "-") {
	#	#	my $rev_com=reverse $str;
	#	#	$rev_com=~ tr/ACGTacgt/TGCAtgca/;
	#	#	print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$rev_com."\n";
	#	#}elsif($STRAND4{$cnt2} eq "+"){
	#	#	print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$str."\n";
	#	#}
	#	if ((($STRAND4{$cnt2} eq "-") and ($GENE2{$cnt2} ne "rps12-1")) and (($STRAND4{$cnt2} eq "-") and ($GENE2{$cnt2} ne "rps12-2"))) {
	#		my $rev_com=reverse $str;
	#		$rev_com=~ tr/ACGTacgt/TGCAtgca/;
	#		if (($GENE2{$cnt2} !~ /1_/) and ($GENE2{$cnt2} !~ /2_/)) {
	#			print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$rev_com."\n";
	#		}elsif ($GENE2{$cnt2} =~ /1_(.*)/) {#1_trnH-GUG
	#			$left_header=$1;
	#			$left_str=$rev_com;
	#		}elsif ($GENE2{$cnt2} =~ /2_(.*)/) {#2_trnH-GUG
	#			$right_header=$1;
	#			$right_str=$rev_com;
	#		}
	#	}elsif (($STRAND4{$cnt2} eq "-") and ($GENE2{$cnt2} eq "rps12-1")) {
	#		my $rev_com=reverse $str;
	#		$rev_com=~ tr/ACGTacgt/TGCAtgca/;
	#		push @rps12_1_sequence,$rev_com;
	#	}elsif (($STRAND4{$cnt2} eq "-") and ($GENE2{$cnt2} eq "rps12-2")) {
	#		my $rev_com=reverse $str;
	#		$rev_com=~ tr/ACGTacgt/TGCAtgca/;
	#		push @rps12_2_sequence,$rev_com;
	#	}elsif(($STRAND4{$cnt2} eq "+") and ($GENE2{$cnt2} ne "rps12-2")){
	#		if (($GENE2{$cnt2} !~ /1_/) and ($GENE2{$cnt2} !~ /2_/)) {
	#			print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$str."\n";
	#		}elsif ($GENE2{$cnt2} =~ /1_(.*)/) {#1_trnH-GUG
	#			$left_header=$1;
	#			$left_str=$str;
	#		}elsif ($GENE2{$cnt2} =~ /2_(.*)/) {#2_trnH-GUG
	#			$right_header=$1;
	#			$right_str=$str;
	#		}
	#	}elsif (($STRAND4{$cnt2} eq "+") and ($GENE2{$cnt2} eq "rps12-2")) {
	#		push @rps12_3_sequence,$str;
	#	}
	#	if ((defined $left_header ne "") and (defined $right_header ne "") and ($left_header eq $right_header)) {
	#		print $out_coding ">".$left_header."_".$SP2{$cnt2}."\n".$left_str.$right_str."\n";
	#	}
	#	$rps12_header=">rps12_".$SP2{$cnt2};
	#}
	#
	#if ((@rps12_1_sequence != 0) and (@rps12_2_sequence != 0) and (@rps12_3_sequence != 0)) {
	#	my $rps12_1_sequence=shift @rps12_1_sequence;
	#	my $rps12_2_sequence=shift @rps12_2_sequence;
	#	my $rps12_3_sequence=shift @rps12_3_sequence;
	#	print $out_coding $rps12_header."\n".$rps12_1_sequence.$rps12_2_sequence."\n";
	#	print $out_coding $rps12_header."\n".$rps12_1_sequence.$rps12_3_sequence."\n";
	#}elsif ((@rps12_1_sequence != 0) and (@rps12_2_sequence != 0) and (@rps12_3_sequence == 0)) {
	#	my $rps12_1_sequence=shift @rps12_1_sequence;
	#	my $rps12_2_sequence=shift @rps12_2_sequence;
	#	print $out_coding $rps12_header."\n".$rps12_1_sequence.$rps12_2_sequence."\n";
	#	print $out_coding $rps12_header."\n".$rps12_1_sequence.$rps12_2_sequence."\n";
	#}elsif ((@rps12_1_sequence != 0) and (@rps12_2_sequence == 0) and (@rps12_3_sequence != 0)) {
	#	my $rps12_1_sequence=shift @rps12_1_sequence;
	#	my $rps12_3_sequence=shift @rps12_3_sequence;
	#	print $out_coding $rps12_header."\n".$rps12_1_sequence.$rps12_3_sequence."\n";
	#	print $out_coding $rps12_header."\n".$rps12_1_sequence.$rps12_3_sequence."\n";
	#}elsif ((@rps12_1_sequence == 0) and (@rps12_2_sequence != 0) and (@rps12_3_sequence != 0)) {
	#	my $rps12_2_sequence=shift @rps12_2_sequence;
	#	my $rps12_3_sequence=shift @rps12_3_sequence;
	#	print $out_coding $rps12_header."\n".$rps12_2_sequence."\n";
	#	print $out_coding $rps12_header."\n".$rps12_3_sequence."\n";
	#}elsif ((@rps12_1_sequence == 0) and (@rps12_2_sequence == 0) and (@rps12_3_sequence != 0)) {
	#	my $rps12_3_sequence=shift @rps12_3_sequence;
	#	print $out_coding $rps12_header."\n".$rps12_3_sequence."\n";
	#	print $out_coding $rps12_header."\n".$rps12_3_sequence."\n";
	#}elsif ((@rps12_1_sequence == 0) and (@rps12_2_sequence != 0) and (@rps12_3_sequence == 0)) {
	#	my $rps12_2_sequence=shift @rps12_2_sequence;
	#	print $out_coding $rps12_header."\n".$rps12_2_sequence."\n";
	#	print $out_coding $rps12_header."\n".$rps12_2_sequence."\n";
	#}elsif ((@rps12_1_sequence != 0) and (@rps12_2_sequence == 0) and (@rps12_3_sequence == 0)) {
	#	my $rps12_1_sequence=shift @rps12_1_sequence;
	#	print $out_coding $rps12_header."\n".$rps12_1_sequence."\n";
	#	print $out_coding $rps12_header."\n".$rps12_1_sequence."\n";
	#}elsif ((@rps12_1_sequence == 0) and (@rps12_2_sequence == 0) and (@rps12_3_sequence == 0)) {
	#}
	#
	#close $out_coding;
	#
	#
	##generate_IGS_ranges
	#my (%SP3,%GENE3,%STRAND5,%START5,%END5,%TYPE5,$last0,$last1,$last2,@output4);
	#my $cnt3=0;
	#foreach (@output3){
	#	$_=~ s/\r|\n//g;
	#	$cnt3++;
	#	($SP3{$cnt3},$GENE3{$cnt3},$STRAND5{$cnt3},$START5{$cnt3},$END5{$cnt3},$TYPE5{$cnt3})=(split /\s+/,$_)[0,1,2,3,4,5];
	#}
	#foreach (keys %SP3){
	#	if ($_==1 and $START5{$_}!=1){
	#		unshift @output4,$SP3{$_}."\t"."start".'-'.$GENE3{$_}."\t"."?"."/".$STRAND5{$_}."\t"."1"."\t".($START5{$_}-1)."\t"."?"."/".$TYPE5{$_}."\n";
	#	}
	#}
	#foreach (1..($cnt3-1)) {
	#	$last0=$_-1;
	#	$last1=$_+1;
	#	$last2=$_+2;
	#	next if ((($END5{$_}+1) >= ($START5{$last1}-1)) and (($END5{$_}+1) < ($END5{$last1}-1)));
	#	next if (($_ > 1) and (($END5{$_}+1) < ($END5{$last0}-1)) and (($END5{$_}+1) < ($START5{$last2}-1)));
	#	if ((($END5{$_}+1) >= ($START5{$last1}-1)) and (($END5{$_}+1) >= ($END5{$last1}-1))){
	#		push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'.$GENE3{$last2}."\t".$STRAND5{$_}."/".$STRAND5{$last2}."\t".($END5{$_}+1)."\t".($START5{$last2}-1)."\t".$TYPE5{$_}."/".$TYPE5{$last2}."\n";
	#	}else{
	#		push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'.$GENE3{$last1}."\t".$STRAND5{$_}."/".$STRAND5{$last1}."\t".($END5{$_}+1)."\t".($START5{$last1}-1)."\t".$TYPE5{$_}."/".$TYPE5{$last1}."\n";
	#	}
	#}
	#foreach (keys %SP3){
	#	if ($_==$cnt3 and $END5{$_}!=$length){
	#		push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'."end"."\t".$STRAND5{$_}."/"."?"."\t".($END5{$_}+1)."\t".$length."\t".$TYPE5{$_}."/"."?"."\n";
	#	}
	#}
	#@output3=();
	##print @output4;
	#
	#
	##extract_IGS
	#my (%SP4,%GENE4,%STRAND6,%START6,%END6,%TYPE6,$seq2);
	#my $cnt4=0;
	#open(my $out_noncoding,">","$filename\_regions_among_genes.fasta");
	#while (@fasta2){
	#	my $header=shift @fasta2;
	#	$seq2=shift @fasta2;
	#}
	#my ($left_header2,$left_str2,$right_header2,$right_str2);
	#foreach (@output4){
	#	$_=~ s/\r|\n//g;
	#	$cnt4++;
	#	($SP4{$cnt4},$GENE4{$cnt4},$STRAND6{$cnt4},$START6{$cnt4},$END6{$cnt4},$TYPE6{$cnt4})=(split /\s+/,$_)[0,1,2,3,4,5];
	#	my $str=substr($seq2,($START6{$cnt4}-1),($END6{$cnt4}-$START6{$cnt4}+1));
	#	if (($GENE4{$cnt4} !~ /1_/) and ($GENE4{$cnt4} !~ /2_/) and ($GENE4{$cnt4} !~ /start-/) and ($GENE4{$cnt4} !~ /-end/)) {
	#		print $out_noncoding ">".$GENE4{$cnt4}."_".$SP4{$cnt4}."\n".$str."\n";
	#	}elsif ($GENE4{$cnt4} =~ /1_/) {#1_trnH-GUG-psbA
	#		my $temp_intergenic=$GENE4{$cnt4};
	#		$temp_intergenic=~ s/1_//g;
	#		print $out_noncoding ">".$temp_intergenic."_".$SP4{$cnt4}."\n".$str."\n";
	#	}elsif ($GENE4{$cnt4} =~ /2_/) {#rpl2-2_trnH-GUG
	#		my $temp_intergenic=$GENE4{$cnt4};
	#		$temp_intergenic=~ s/2_//g;
	#		print $out_noncoding ">".$temp_intergenic."_".$SP4{$cnt4}."\n".$str."\n";
	#	}elsif ($GENE4{$cnt4} =~ /start-(.*)/) {#start-trnH-GUG
	#		$right_header2=$1;
	#		$right_str2=$str;
	#	}elsif ($GENE4{$cnt4} =~ /(.*)-end/) {#rpl2-end
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
		}elsif(/ {5}gene {12}/){
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
	#open (my $out_bed,">","$filename\_genes.bed");
	#foreach (@output3){
	#	print $out_bed $_;
	#}
	#close $out_bed;

	#output_bed_file
	open (my $out_bed,">","$target_path\\$filename\.bed");

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


	##extract_gene
	#my (%SP2,%GENE2,%STRAND4,%START4,%END4,%TYPE4,$seq1);
	#my $cnt2=0;
	#open(my $out_coding,">","$filename\_genes.fasta");
	#while (@fasta1){
	#	my $header=shift @fasta1;
	#	$seq1=shift @fasta1;
	#}
	#
	#my ($rps12_header,@rps12_1_sequence,@rps12_2_sequence,@rps12_3_sequence);
	#my ($left_header,$left_str,$right_header,$right_str);
	#foreach (@output3){
	#	$_=~ s/\r|\n//g;
	#	$cnt2++;
	#	($SP2{$cnt2},$GENE2{$cnt2},$STRAND4{$cnt2},$START4{$cnt2},$END4{$cnt2},$TYPE4{$cnt2})=(split /\s+/,$_)[0,1,2,3,4,5];
	#	my $str=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
	#	#if ($STRAND4{$cnt2} eq "-") {
	#	#	my $rev_com=reverse $str;
	#	#	$rev_com=~ tr/ACGTacgt/TGCAtgca/;
	#	#	print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$rev_com."\n";
	#	#}elsif($STRAND4{$cnt2} eq "+"){
	#	#	print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$str."\n";
	#	#}
	#	if ((($STRAND4{$cnt2} eq "-") and ($GENE2{$cnt2} ne "rps12-1")) and (($STRAND4{$cnt2} eq "-") and ($GENE2{$cnt2} ne "rps12-2"))) {
	#		my $rev_com=reverse $str;
	#		$rev_com=~ tr/ACGTacgt/TGCAtgca/;
	#		if (($GENE2{$cnt2} !~ /1_/) and ($GENE2{$cnt2} !~ /2_/)) {
	#			print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$rev_com."\n";
	#		}elsif ($GENE2{$cnt2} =~ /1_(.*)/) {#1_trnH-GUG
	#			$left_header=$1;
	#			$left_str=$rev_com;
	#		}elsif ($GENE2{$cnt2} =~ /2_(.*)/) {#2_trnH-GUG
	#			$right_header=$1;
	#			$right_str=$rev_com;
	#		}
	#	}elsif (($STRAND4{$cnt2} eq "-") and ($GENE2{$cnt2} eq "rps12-1")) {
	#		my $rev_com=reverse $str;
	#		$rev_com=~ tr/ACGTacgt/TGCAtgca/;
	#		push @rps12_1_sequence,$rev_com;
	#	}elsif (($STRAND4{$cnt2} eq "-") and ($GENE2{$cnt2} eq "rps12-2")) {
	#		my $rev_com=reverse $str;
	#		$rev_com=~ tr/ACGTacgt/TGCAtgca/;
	#		push @rps12_2_sequence,$rev_com;
	#	}elsif(($STRAND4{$cnt2} eq "+") and ($GENE2{$cnt2} ne "rps12-2")){
	#		if (($GENE2{$cnt2} !~ /1_/) and ($GENE2{$cnt2} !~ /2_/)) {
	#			print $out_coding ">".$GENE2{$cnt2}."_".$SP2{$cnt2}."\n".$str."\n";
	#		}elsif ($GENE2{$cnt2} =~ /1_(.*)/) {#1_trnH-GUG
	#			$left_header=$1;
	#			$left_str=$str;
	#		}elsif ($GENE2{$cnt2} =~ /2_(.*)/) {#2_trnH-GUG
	#			$right_header=$1;
	#			$right_str=$str;
	#		}
	#	}elsif (($STRAND4{$cnt2} eq "+") and ($GENE2{$cnt2} eq "rps12-2")) {
	#		push @rps12_3_sequence,$str;
	#	}
	#	if ((defined $left_header ne "") and (defined $right_header ne "") and ($left_header eq $right_header)) {
	#		print $out_coding ">".$left_header."_".$SP2{$cnt2}."\n".$left_str.$right_str."\n";
	#	}
	#	$rps12_header=">rps12_".$SP2{$cnt2};
	#}
	#
	#if ((@rps12_1_sequence != 0) and (@rps12_2_sequence != 0) and (@rps12_3_sequence != 0)) {
	#	my $rps12_1_sequence=shift @rps12_1_sequence;
	#	my $rps12_2_sequence=shift @rps12_2_sequence;
	#	my $rps12_3_sequence=shift @rps12_3_sequence;
	#	print $out_coding $rps12_header."\n".$rps12_1_sequence.$rps12_2_sequence."\n";
	#	print $out_coding $rps12_header."\n".$rps12_1_sequence.$rps12_3_sequence."\n";
	#}elsif ((@rps12_1_sequence != 0) and (@rps12_2_sequence != 0) and (@rps12_3_sequence == 0)) {
	#	my $rps12_1_sequence=shift @rps12_1_sequence;
	#	my $rps12_2_sequence=shift @rps12_2_sequence;
	#	print $out_coding $rps12_header."\n".$rps12_1_sequence.$rps12_2_sequence."\n";
	#	print $out_coding $rps12_header."\n".$rps12_1_sequence.$rps12_2_sequence."\n";
	#}elsif ((@rps12_1_sequence != 0) and (@rps12_2_sequence == 0) and (@rps12_3_sequence != 0)) {
	#	my $rps12_1_sequence=shift @rps12_1_sequence;
	#	my $rps12_3_sequence=shift @rps12_3_sequence;
	#	print $out_coding $rps12_header."\n".$rps12_1_sequence.$rps12_3_sequence."\n";
	#	print $out_coding $rps12_header."\n".$rps12_1_sequence.$rps12_3_sequence."\n";
	#}elsif ((@rps12_1_sequence == 0) and (@rps12_2_sequence != 0) and (@rps12_3_sequence != 0)) {
	#	my $rps12_2_sequence=shift @rps12_2_sequence;
	#	my $rps12_3_sequence=shift @rps12_3_sequence;
	#	print $out_coding $rps12_header."\n".$rps12_2_sequence."\n";
	#	print $out_coding $rps12_header."\n".$rps12_3_sequence."\n";
	#}elsif ((@rps12_1_sequence == 0) and (@rps12_2_sequence == 0) and (@rps12_3_sequence != 0)) {
	#	my $rps12_3_sequence=shift @rps12_3_sequence;
	#	print $out_coding $rps12_header."\n".$rps12_3_sequence."\n";
	#	print $out_coding $rps12_header."\n".$rps12_3_sequence."\n";
	#}elsif ((@rps12_1_sequence == 0) and (@rps12_2_sequence != 0) and (@rps12_3_sequence == 0)) {
	#	my $rps12_2_sequence=shift @rps12_2_sequence;
	#	print $out_coding $rps12_header."\n".$rps12_2_sequence."\n";
	#	print $out_coding $rps12_header."\n".$rps12_2_sequence."\n";
	#}elsif ((@rps12_1_sequence != 0) and (@rps12_2_sequence == 0) and (@rps12_3_sequence == 0)) {
	#	my $rps12_1_sequence=shift @rps12_1_sequence;
	#	print $out_coding $rps12_header."\n".$rps12_1_sequence."\n";
	#	print $out_coding $rps12_header."\n".$rps12_1_sequence."\n";
	#}elsif ((@rps12_1_sequence == 0) and (@rps12_2_sequence == 0) and (@rps12_3_sequence == 0)) {
	#}
	#
	#close $out_coding;
	#
	#
	##generate_IGS_ranges
	#my (%SP3,%GENE3,%STRAND5,%START5,%END5,%TYPE5,$last0,$last1,$last2,@output4);
	#my $cnt3=0;
	#foreach (@output3){
	#	$_=~ s/\r|\n//g;
	#	$cnt3++;
	#	($SP3{$cnt3},$GENE3{$cnt3},$STRAND5{$cnt3},$START5{$cnt3},$END5{$cnt3},$TYPE5{$cnt3})=(split /\s+/,$_)[0,1,2,3,4,5];
	#}
	#foreach (keys %SP3){
	#	if ($_==1 and $START5{$_}!=1){
	#		unshift @output4,$SP3{$_}."\t"."start".'-'.$GENE3{$_}."\t"."?"."/".$STRAND5{$_}."\t"."1"."\t".($START5{$_}-1)."\t"."?"."/".$TYPE5{$_}."\n";
	#	}
	#}
	#foreach (1..($cnt3-1)) {
	#	$last0=$_-1;
	#	$last1=$_+1;
	#	$last2=$_+2;
	#	next if ((($END5{$_}+1) >= ($START5{$last1}-1)) and (($END5{$_}+1) < ($END5{$last1}-1)));
	#	next if (($_ > 1) and (($END5{$_}+1) < ($END5{$last0}-1)) and (($END5{$_}+1) < ($START5{$last2}-1)));
	#	if ((($END5{$_}+1) >= ($START5{$last1}-1)) and (($END5{$_}+1) >= ($END5{$last1}-1))){
	#		push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'.$GENE3{$last2}."\t".$STRAND5{$_}."/".$STRAND5{$last2}."\t".($END5{$_}+1)."\t".($START5{$last2}-1)."\t".$TYPE5{$_}."/".$TYPE5{$last2}."\n";
	#	}else{
	#		push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'.$GENE3{$last1}."\t".$STRAND5{$_}."/".$STRAND5{$last1}."\t".($END5{$_}+1)."\t".($START5{$last1}-1)."\t".$TYPE5{$_}."/".$TYPE5{$last1}."\n";
	#	}
	#}
	#foreach (keys %SP3){
	#	if ($_==$cnt3 and $END5{$_}!=$length){
	#		push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'."end"."\t".$STRAND5{$_}."/"."?"."\t".($END5{$_}+1)."\t".$length."\t".$TYPE5{$_}."/"."?"."\n";
	#	}
	#}
	#@output3=();
	##print @output4;
	#
	#
	##extract_IGS
	#my (%SP4,%GENE4,%STRAND6,%START6,%END6,%TYPE6,$seq2);
	#my $cnt4=0;
	#open(my $out_noncoding,">","$filename\_regions_among_genes.fasta");
	#while (@fasta2){
	#	my $header=shift @fasta2;
	#	$seq2=shift @fasta2;
	#}
	#my ($left_header2,$left_str2,$right_header2,$right_str2);
	#foreach (@output4){
	#	$_=~ s/\r|\n//g;
	#	$cnt4++;
	#	($SP4{$cnt4},$GENE4{$cnt4},$STRAND6{$cnt4},$START6{$cnt4},$END6{$cnt4},$TYPE6{$cnt4})=(split /\s+/,$_)[0,1,2,3,4,5];
	#	my $str=substr($seq2,($START6{$cnt4}-1),($END6{$cnt4}-$START6{$cnt4}+1));
	#	if (($GENE4{$cnt4} !~ /1_/) and ($GENE4{$cnt4} !~ /2_/) and ($GENE4{$cnt4} !~ /start-/) and ($GENE4{$cnt4} !~ /-end/)) {
	#		print $out_noncoding ">".$GENE4{$cnt4}."_".$SP4{$cnt4}."\n".$str."\n";
	#	}elsif ($GENE4{$cnt4} =~ /1_/) {#1_trnH-GUG-psbA
	#		my $temp_intergenic=$GENE4{$cnt4};
	#		$temp_intergenic=~ s/1_//g;
	#		print $out_noncoding ">".$temp_intergenic."_".$SP4{$cnt4}."\n".$str."\n";
	#	}elsif ($GENE4{$cnt4} =~ /2_/) {#rpl2-2_trnH-GUG
	#		my $temp_intergenic=$GENE4{$cnt4};
	#		$temp_intergenic=~ s/2_//g;
	#		print $out_noncoding ">".$temp_intergenic."_".$SP4{$cnt4}."\n".$str."\n";
	#	}elsif ($GENE4{$cnt4} =~ /start-(.*)/) {#start-trnH-GUG
	#		$right_header2=$1;
	#		$right_str2=$str;
	#	}elsif ($GENE4{$cnt4} =~ /(.*)-end/) {#rpl2-end
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
	open(my $in_bed_reference,"<","$reference_name\.bed");
	while (<$in_bed_reference>) {
		my @array=split /\t/,$_;
		if ($array[1]=~ /((rrn)(.+))/){
			my $genes=$1;
			$rrn_reference{$genes}++;
		}elsif ($array[1]=~ /((trn)(.+))/){
			my $genes=$1;
			$trn_reference{$genes}++;
		}elsif ($array[1]=~ /((?!trn)(.+))/ and $array[1]=~ /((?!rrn)(.+))/){
			my $genes=$1;
			$genes=~ s/rps12\-1/rps12/;
			$genes=~ s/rps12\-2/rps12/;
			$pcg_reference{$genes}++;
		}
	}
	close $in_bed_reference;
	unlink("$reference_name\.bed");
}

my $rrn_number_reference=keys %rrn_reference;
my $trn_number_reference=keys %trn_reference;
my $pcg_number_reference=keys %pcg_reference;
#print $rrn_number_reference."\n".$trn_number_reference."\n".$pcg_number_reference."\n";




while (@filenames4) {
	my $filename_gb=shift @filenames4;
	my $target_name=$filename_gb;
	$target_name=substr($target_name,0,rindex($target_name,"\."));
	$target_name=~ s/\s/_/g;
	my $filename=substr($target_name,rindex($target_name,"\/")+1);

	#my $target_name=substr($input_target,0,rindex($input_target,"\."));
	#$target_name=~ s/\s/_/g;
	my (%rrn_target,%trn_target,%pcg_target);
	open(my $in_bed_target,"<","$target_name\.bed");
	while (<$in_bed_target>) {
		my @array=split /\t/,$_;
		if ($array[1]=~ /((rrn)(.+))/){
			my $genes=$1;
			$rrn_target{$genes}++;
		}elsif ($array[1]=~ /((trn)(.+))/){
			my $genes=$1;
			$trn_target{$genes}++;
		}elsif ($array[1]=~ /((?!trn)(.+))/ and $array[1]=~ /((?!rrn)(.+))/){
			my $genes=$1;
			$genes=~ s/rps12\-1/rps12/;
			$genes=~ s/rps12\-2/rps12/;
			$pcg_target{$genes}++;
		}
	}
	close $in_bed_target;
	unlink("$target_name\.bed");

	my $rrn_number_target=keys %rrn_target;
	my $trn_number_target=keys %trn_target;
	my $pcg_number_target=keys %pcg_target;
	#print $rrn_number_target."\n".$trn_number_target."\n".$pcg_number_target."\n";






	####total reference, total target, total hitting, total missing, total redundant
	my $total_reference=$rrn_number_reference+$trn_number_reference+$pcg_number_reference;
	my $total_target=$rrn_number_target+$trn_number_target+$pcg_number_target;

	my @total_hitting;
	my @total_missing;
	my @total_redundant;
	my %total_reference=(%rrn_reference,%trn_reference,%pcg_reference);
	my %total_target=(%rrn_target,%trn_target,%pcg_target);

	foreach my $key (keys %total_reference) {
		if (exists $total_target{$key}) {
			push @total_hitting,$key;
		}
		if (!exists $total_target{$key}) {
			push @total_missing,$key;
		}
	}
	my $total_hitting=@total_hitting;
	my $total_missing=@total_missing;
	foreach my $key (keys %total_target) {
		if (!exists $total_reference{$key}) {
			push @total_redundant,$key;
		}
	}
	my $total_redundant=@total_redundant;

	my $percentage_hitting;
	my $percentage_hitting_s;
	if ($total_reference != 0) {
		$percentage_hitting=$total_hitting/$total_reference;
		$percentage_hitting_s=sprintf("%.2f",$percentage_hitting)*100;
	}elsif ($total_reference == 0) {
		$percentage_hitting=$total_hitting."/".$total_reference;
		$percentage_hitting_s=$percentage_hitting."*".100;
	}

	my $percentage_missing;
	my $percentage_missing_s;
	if ($total_reference != 0) {
		$percentage_missing=$total_missing/$total_reference;
		$percentage_missing_s=sprintf("%.2f",$percentage_missing)*100;
	}elsif ($total_reference == 0) {
		$percentage_missing=$total_missing."/".$total_reference;
		$percentage_missing_s=$percentage_missing."*".100;
	}

	my $percentage_redundant;
	my $percentage_redundant_s;
	if ($total_reference != 0) {
		$percentage_redundant=$total_redundant/$total_reference;
		$percentage_redundant_s=sprintf("%.2f",$percentage_redundant)*100;
	}elsif ($total_reference == 0) {
		$percentage_redundant=$total_redundant."/".$total_reference;
		$percentage_redundant_s=$percentage_redundant."*".100;
	}



	####rrn/trn/pcg reference, rrn/trn/pcg target, rrn/trn/pcg hitting, rrn/trn/pcg missing, rrn/trn/pcg redundant
	my @rrn_hitting;
	my @rrn_missing;
	my @rrn_redundant;
	my @trn_hitting;
	my @trn_missing;
	my @trn_redundant;
	my @pcg_hitting;
	my @pcg_missing;
	my @pcg_redundant;

	foreach my $key (keys %rrn_reference) {
		if (exists $rrn_target{$key}) {
			push @rrn_hitting,$key;
		}
		if (!exists $rrn_target{$key}) {
			push @rrn_missing,$key;
		}
	}
	my $rrn_hitting=@rrn_hitting;
	my $rrn_missing=@rrn_missing;
	foreach my $key (keys %rrn_target) {
		if (!exists $rrn_reference{$key}) {
			push @rrn_redundant,$key;
		}
	}
	my $rrn_redundant=@rrn_redundant;

	foreach my $key (keys %trn_reference) {
		if (exists $trn_target{$key}) {
			push @trn_hitting,$key;
		}
		if (!exists $trn_target{$key}) {
			push @trn_missing,$key;
		}
	}
	my $trn_hitting=@trn_hitting;
	my $trn_missing=@trn_missing;
	foreach my $key (keys %trn_target) {
		if (!exists $trn_reference{$key}) {
			push @trn_redundant,$key;
		}
	}
	my $trn_redundant=@trn_redundant;

	foreach my $key (keys %pcg_reference) {
		if (exists $pcg_target{$key}) {
			push @pcg_hitting,$key;
		}
		if (!exists $pcg_target{$key}) {
			push @pcg_missing,$key;
		}
	}
	my $pcg_hitting=@pcg_hitting;
	my $pcg_missing=@pcg_missing;
	foreach my $key (keys %pcg_target) {
		if (!exists $pcg_reference{$key}) {
			push @pcg_redundant,$key;
		}
	}
	my $pcg_redundant=@pcg_redundant;




	####output the statistic results
	open(my $assessment,">","$output_directory/$filename_reference\_$filename.txt");
	print $assessment "H:hitting, M:missing, R:redundant, n:number, p:percentage, P:PCGs, T:tRNAs, R:rRNAs"."\n\n";
	print $assessment "1. Number of hitting, missing, and redundant genes:"."\n";
	print $assessment "|$total_reference total gene numbers in reference|"."\n";
	print $assessment "|$total_target total gene numbers in target|"."\n\n";

	print $assessment "|$total_hitting number of hitting genes (nH)|"."\n";
	print $assessment "|$total_missing number of missing genes (nM)|"."\n";
	print $assessment "|$total_redundant number of redundant genes (nR)|"."\n\n";

	print $assessment "|$percentage_hitting_s% percentage of hitting gene (pH)|"."\n";
	print $assessment "|$percentage_missing_s% percentage of missing genes (pM)|"."\n";
	print $assessment "|$percentage_redundant_s% percentage of redundant genes (pR)|"."\n\n";

	print $assessment "2. Number of hitting, missing, and redundant PCGs, tRNAs and rRNAs:"."\n";
	print $assessment "|$pcg_hitting number of hitting PCGs (nHP)|"."\n";
	print $assessment "|$pcg_missing number of missing PCGs (nMP)|"."\n";
	print $assessment "|$pcg_redundant number of redundant PCGs (nRP)|"."\n\n";

	print $assessment "|$trn_hitting number of hitting tRNAs (nHT)|"."\n";
	print $assessment "|$trn_missing number of missing tRNAs (nMT)|"."\n";
	print $assessment "|$trn_redundant number of redundant tRNAs (nRT)|"."\n\n";

	print $assessment "|$rrn_hitting number of hitting rRNAs (nHR)|"."\n";
	print $assessment "|$rrn_missing number of missing rRNAs (nMR)|"."\n";
	print $assessment "|$rrn_redundant number of redundant rRNAs (nRR)|"."\n\n";

	print $assessment "3. Name of hitting, missing, and redundant PCGs, tRNAs and rRNAs:"."\n";
	print $assessment "|name of hitting PCGs|"."\n";
	if ($pcg_hitting!=0) {
		my @pcg_hitting_sorted=sort @pcg_hitting;
		print $assessment "@pcg_hitting_sorted"."\n";
	}elsif ($pcg_hitting==0) {
		print $assessment "None"."\n";
	}
	print $assessment "|name of missing PCGs|"."\n";
	if ($pcg_missing!=0) {
		my @pcg_missing_sorted=sort @pcg_missing;
		print $assessment "@pcg_missing_sorted"."\n";
	}elsif ($pcg_missing==0) {
		print $assessment "None"."\n";
	}
	print $assessment "|name of redundant PCGs|"."\n";
	if ($pcg_redundant!=0) {
		my @pcg_redundant_sorted=sort @pcg_redundant;
		print $assessment "@pcg_redundant_sorted"."\n\n";
	}elsif ($pcg_redundant==0) {
		print $assessment "None"."\n\n";
	}

	print $assessment "|name of hitting tRNAs|"."\n";
	if ($trn_hitting!=0) {
		my @trn_hitting_sorted=sort @trn_hitting;
		print $assessment "@trn_hitting_sorted"."\n";
	}elsif ($trn_hitting==0) {
		print $assessment "None"."\n";
	}
	print $assessment "|name of missing tRNAs|"."\n";
	if ($trn_missing!=0) {
		my @trn_missing_sorted=sort @trn_missing;
		print $assessment "@trn_missing_sorted"."\n";
	}elsif ($trn_missing==0) {
		print $assessment "None"."\n";
	}
	print $assessment "|name of redundant tRNAs|"."\n";
	if ($trn_redundant!=0) {
		my @trn_redundant_sorted=sort @trn_redundant;
		print $assessment "@trn_redundant_sorted"."\n\n";
	}elsif ($trn_redundant==0) {
		print $assessment "None"."\n\n";
	}

	print $assessment "|name of hitting rRNAs|"."\n";
	if ($rrn_hitting!=0) {
		my @rrn_hitting_sorted=sort @rrn_hitting;
		print $assessment "@rrn_hitting_sorted"."\n";
	}elsif ($rrn_hitting==0) {
		print $assessment "None"."\n";
	}
	print $assessment "|name of missing rRNAs|"."\n";
	if ($rrn_missing!=0) {
		my @rrn_missing_sorted=sort @rrn_missing;
		print $assessment "@rrn_missing_sorted"."\n";
	}elsif ($rrn_missing==0) {
		print $assessment "None"."\n";
	}
	print $assessment "|name of redundant rRNAs|"."\n";
	if ($rrn_redundant!=0) {
		my @rrn_redundant_sorted=sort @rrn_redundant;
		print $assessment "@rrn_redundant_sorted"."\n";
	}elsif ($rrn_redundant==0) {
		print $assessment "None"."\n";
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

    assess_gene_number_v1.pl version 1.0

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

    assess gene number and name between reference and target files

=head1 SYNOPSIS

    assess_gene_number_v1.pl -r -t -o
    Copyright (C) 2025 Xiao-Jian Qu
    Please contact <quxiaojian@sdnu.edu.cn>, if you have any bugs or questions.

    [-h -help]        help information.
    [-r -reference]   required: (default: reference) input reference directory containing one GenBank-formatted file.
    [-t -target]      required: (default: target) input target directory containing GenBank-formatted file(s).
    [-o -output]      required: (default: output) output directory containing TXT-format file(s) with statistic result.

=cut
