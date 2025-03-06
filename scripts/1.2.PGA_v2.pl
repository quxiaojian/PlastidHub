#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Find;
no warnings "uninitialized";
use Data::Dumper;
$|=1;

my $global_options=&argument();
my $reference_directory=&default("reference","reference");
my $sequence_directory=&default("target","target");
my $inverted_repeat=&default("1000","ir");
my $similarity_value=&default("40","percent");
my $link=&default("Y","link");
my $redundancy=&default("N","redundancy");
my $qcoverage=&default("0.5,2.0","qcoverage");
my ($short,$long)=split /,/,$qcoverage;
$short=eval($short);
$long=eval($long);
my $output_directory=&default("gb","out");
my $type=&default("circular","form");
my $warning=&default("warning","warning");


print "\nPGA_v2.pl Plastid Genome Annotator version 2.0
Copyright (C) 2025 Xiao-Jian Qu
Email: quxiaojian\@sdnu.edu.cn\n\n";

my $now1=time;
my $now2=&gettime;
print "$now2 || Begin extracting annotations from reference!";
print "\n";

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
my @a = (0..9,'$','%','a'..'z','A'..'Z','-','+','_');
my $random = "%".join ('',map ($a[int rand @a],0..($maxLenth-1)))."%";

############################################################
## extract_gene_coding_CDS
############################################################
my $pattern1=".gb";
my $pattern2=".gbf";
my $pattern3=".gbk";
my @filenames;
find(\&target1,$reference_directory);
sub target1{
	if (/$pattern1/ or /$pattern2/ or /$pattern3/){
        push @filenames,"$File::Find::name";
    }
    return;
}

my $reference_path;
my $filename_base;
my @reference_filename;
while (@filenames) {
	my $filename_gb=shift @filenames;
	$filename_base=$filename_gb;
	#$filename_base=~ s/(.*).gb/$1/g;
	$filename_base=~ s/\s+/_/gg;
	$filename_base=substr($filename_base,0,rindex($filename_base,"\."));
	push @reference_filename,$filename_base;
	my $latin_name=substr ($filename_base,rindex($filename_base,"\/")+1);
	$reference_path=$filename_base;
	$reference_path=~ s/$latin_name//g;
	open(my $in_gb,"<",$filename_gb);
	open(my $out_gb,">","$filename_base\_temp1");
	while (<$in_gb>){
		$_=~ s/\r\n/\n/g;
		if ($_=~ /\),\n/){
			$_=~ s/\),\n/\),/g;
		}elsif($_=~ /,\n/){
			$_=~ s/,\n/,/g;
		}
		print $out_gb $_;
	}
	close $in_gb;
	close $out_gb;

	open(my $in_gbk,"<","$filename_base\_temp1");
	open(my $out_gbk,">","$filename_base\_temp2");
	while (<$in_gbk>){
		$_=~ s/,\s+/,/g;
		print $out_gbk $_;
	}
	close $in_gbk;
	close $out_gbk;

	############################################################
	## extract_bed_gene
	############################################################
#	if (!-e "$filename_base\_gene.bed"){
	{
		#generate_bed_file
		my (@row_array,$species_name,$length,$element,@genearray,@output1);
		my $mark=0;
		my $tick=0;
		open (my $in_genbank,"<","$filename_base\_temp2");
		while (<$in_genbank>){
			$_=~ s/\r|\n//g;
			@row_array=split /\s+/,$_;
			if (/^LOCUS/i){
				$species_name=$latin_name;
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
					$element=$species_name.":"."1_".$1.":".$element;
					push @genearray,$element;
					$element=();
					$mark=0;
					$tick=0;
				}elsif ($tick==2) {
					$element=$species_name.":"."2_".$1.":".$element;
					push @genearray,$element;
					$element=();
					$mark=0;
					$tick=0;
				}
			}
		}
		close $in_genbank;
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
		open(my $in_genebank,"<","$filename_base\_temp2");
	    while (<$in_genebank>){
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
		close $in_genebank;
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
			}elsif (defined $STRAND2{$_} ne "") {# rps12_gene,including rps12+1_gene and rps12+2_gene
				push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
				push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
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
		open (my $out_bed,">","$filename_base\_gene.bed");
		foreach (@output3){
			my @row=split /\t/,$_;
			my $abs=abs($row[3]-$row[4]);
			if ($_!~ /rps12/) {
				print $out_bed $_;
			}
			if (($_=~ /rps12/) and ($abs < 200)) {
				$_=~ s/rps12/rps12+1/;
				print $out_bed $_;
			}
			if (($_=~ /rps12/) and ($abs > 200)) {
				$_=~ s/rps12/rps12+2/;
				print $out_bed $_;
			}
		}
		close $out_bed;

		#extract_gene
		my (%SP2,%GENE2,%STRAND4,%START4,%END4,%TYPE4,$seq1);
		my $cnt2=0;
		open(my $out_coding,">","$filename_base\_gene.fasta");
		while (@fasta1){
			my $header=shift @fasta1;
			$seq1=shift @fasta1;
		}
		my %hash_join;
		foreach (@output3){
			$_=~ s/\r|\n//g;
			$cnt2++;
			($SP2{$cnt2},$GENE2{$cnt2},$STRAND4{$cnt2},$START4{$cnt2},$END4{$cnt2},$TYPE4{$cnt2})=(split /\s+/,$_)[0,1,2,3,4,5];
			my $str=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
			my $rev_com=reverse $str;
			$rev_com=~ tr/ACGTacgt/TGCAtgca/;
			if (($STRAND4{$cnt2} eq "-") and (($GENE2{$cnt2}!~ /1_/) and ($GENE2{$cnt2}!~ /2_/))) {
				print $out_coding ">".$GENE2{$cnt2}."_gene_".$SP2{$cnt2}."\n".$rev_com."\n";
			}elsif (($STRAND4{$cnt2} eq "-") and ($GENE2{$cnt2}=~ /(^1_(.+))/)) {
				my $header1=$1."_gene_".$SP2{$cnt2};
				$hash_join{$header1}=$rev_com;
			}elsif (($STRAND4{$cnt2} eq "-") and ($GENE2{$cnt2}=~ /(^2_(.+))/)) {
				my $header2=$1."_gene_".$SP2{$cnt2};
				$hash_join{$header2}=$rev_com;
			}elsif(($STRAND4{$cnt2} eq "+") and (($GENE2{$cnt2}!~ /1_/) and ($GENE2{$cnt2}!~ /2_/))){
				print $out_coding ">".$GENE2{$cnt2}."_gene_".$SP2{$cnt2}."\n".$str."\n";
			}elsif (($STRAND4{$cnt2} eq "+") and ($GENE2{$cnt2}=~ /(^1_(.+))/)) {
				my $header1=$1."_gene_".$SP2{$cnt2};
				$hash_join{$header1}=$str;
			}elsif (($STRAND4{$cnt2} eq "+") and ($GENE2{$cnt2}=~ /(^2_(.+))/)) {
				my $header2=$1."_gene_".$SP2{$cnt2};
				$hash_join{$header2}=$str;
			}
		}
		my ($join1,$join2,$header_join);
		foreach my $key (keys %hash_join) {
			$key=~ /^(\d)_((.+)_gene_(.+))/;
			$header_join=$2;
			if ($1==1) {
				$join1=$hash_join{$key};
			}elsif($1==2) {
				$join2=$hash_join{$key};
			}
		}
		print $out_coding ">".$header_join."\n".$join1.$join2."\n";
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
			next if ((($END5{$_}+1) >= ($START5{$last1}-1)) and (($END5{$_}+1) < ($END5{$last1}-1)));
			next if (($_ > 1) and (($END5{$_}+1) < ($END5{$last0}-1)) and (($END5{$_}+1) < ($START5{$last2}-1)));
			if ((($END5{$_}+1) >= ($START5{$last1}-1)) and (($END5{$_}+1) >= ($END5{$last1}-1))){
	    		push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'.$GENE3{$last2}."\t".$STRAND5{$_}."/".$STRAND5{$last2}."\t".($END5{$_}+1)."\t".($START5{$last2}-1)."\t".$TYPE5{$_}."/".$TYPE5{$last2}."\n";
	        }else{
	    		push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'.$GENE3{$last1}."\t".$STRAND5{$_}."/".$STRAND5{$last1}."\t".($END5{$_}+1)."\t".($START5{$last1}-1)."\t".$TYPE5{$_}."/".$TYPE5{$last1}."\n";
	        }
		}
		foreach (keys %SP3){
			if ($_==$cnt3){
				push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'."end"."\t".$STRAND5{$_}."/"."?"."\t".($END5{$_}+1)."\t".$length."\t".$TYPE5{$_}."/"."?"."\n";
			}
		}
		@output3=();

		#extract_IGS
		my (%SP4,%GENE4,%STRAND6,%START6,%END6,%TYPE6,$seq2);
		my $cnt4=0;
		open(my $out_noncoding,">","$filename_base\_IGS.fasta");
		while (@fasta2){
			my $header=shift @fasta2;
			$seq2=shift @fasta2;
		}
		foreach (@output4){
			$_=~ s/\r|\n//g;
			$cnt4++;
			($SP4{$cnt4},$GENE4{$cnt4},$STRAND6{$cnt4},$START6{$cnt4},$END6{$cnt4},$TYPE6{$cnt4})=(split /\s+/,$_)[0,1,2,3,4,5];
			my $str=substr($seq2,($START6{$cnt4}-1),($END6{$cnt4}-$START6{$cnt4}+1));
			print $out_noncoding ">".$GENE4{$cnt4}."_gene_".$SP4{$cnt4}."\n".$str."\n";
		}
		@output4=();
		close $out_noncoding;
	    unlink "$filename_base\_IGS.fasta";
    }


	############################################################
	## extract_bed_coding
	############################################################
#	if (!-e "$filename_base\_coding.bed"){
	{
		#generate_bed_file
		my (@row_array,$species_name,$length,$element,@genearray,@output1);
		my $mark=0;
		my $tick=0;
		open (my $in_genbank,"<","$filename_base\_temp2");
		while (<$in_genbank>){
			$_=~ s/\r|\n//g;
			@row_array=split /\s+/,$_;
			if (/^LOCUS/i){
				$species_name=$latin_name;
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
					$element=$species_name.":"."1_".$1.":".$element;
					push @genearray,$element;
					$element=();
					$mark=0;
					$tick=0;
				}elsif ($tick==2) {
					$element=$species_name.":"."2_".$1.":".$element;
					push @genearray,$element;
					$element=();
					$mark=0;
					$tick=0;
				}
			}
		}
		close $in_genbank;
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
		open(my $in_genebank,"<","$filename_base\_temp2");
	    while (<$in_genebank>){
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
		close $in_genebank;
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
			}elsif ((defined $STRAND2{$_} ne "") and (defined $STRAND3{$_} eq "") and ($GENE1{$_} ne "rps12")) {# non-rps12
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
			}elsif ((defined $STRAND2{$_} ne "") and (defined $STRAND3{$_} eq "") and ($GENE1{$_} eq "rps12")) {# rps12 in IRb and IRa
				push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
				push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
			}elsif ((defined $STRAND2{$_} ne "") and (defined $STRAND3{$_} ne "") and ($GENE1{$_} ne "rps12")) {# ycf3 or clpP
				if (($STRAND1{$_} eq "-") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_})){
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-3"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
				}elsif(($STRAND1{$_} eq "-") and ($START1{$_} > $START2{$_}) and ($START2{$_} > $START3{$_})){
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
				}
			}elsif ((defined $STRAND2{$_} ne "") and (defined $STRAND3{$_} ne "") and ($GENE1{$_} eq "rps12")) {# rps12
					push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
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
		open (my $out_bed,">","$filename_base\_coding.bed");
		foreach (@output3){
			my @row=split /\t/,$_;
			my $abs=abs($row[3]-$row[4]);
			if ($_!~ /rps12/) {
				print $out_bed $_;
			}
			if (($_=~ /rps12/) and ($abs < 50)) {
				$_=~ s/rps12/rps12+2-2/;
				print $out_bed $_;
			}
			if (($_=~ /rps12/) and ($abs > 50) and ($abs < 200)) {
				$_=~ s/rps12/rps12+1/;
				print $out_bed $_;
			}
			if (($_=~ /rps12/) and ($abs > 200) and ($abs < 245)) {
				$_=~ s/rps12/rps12+2-1/;
				print $out_bed $_;
			}
			if (($_=~ /rps12/) and ($abs > 245)) {
				$_=~ s/rps12/rps12+2/;
				print $out_bed $_;
			}
		}
		close $out_bed;

		#extract_coding_regions
		my (%SP2,%GENE2,%STRAND4,%START4,%END4,%TYPE4,$seq1);
		my $cnt2=0;
		open(my $out_coding,">","$filename_base\_coding.fasta");
		while (@fasta1){
			my $header=shift @fasta1;
			$seq1=shift @fasta1;
		}
		my %hash_join;
		foreach (@output3){
			$_=~ s/\r|\n//g;
			$cnt2++;
			($SP2{$cnt2},$GENE2{$cnt2},$STRAND4{$cnt2},$START4{$cnt2},$END4{$cnt2},$TYPE4{$cnt2})=(split /\s+/,$_)[0,1,2,3,4,5];
			my $str=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
			my $rev_com=reverse $str;
			$rev_com=~ tr/ACGTacgt/TGCAtgca/;
			if (($STRAND4{$cnt2} eq "-") and (($GENE2{$cnt2}!~ /1_/) and ($GENE2{$cnt2}!~ /2_/))) {
				print $out_coding ">".$GENE2{$cnt2}."_coding_".$SP2{$cnt2}."\n".$rev_com."\n";
			}elsif (($STRAND4{$cnt2} eq "-") and ($GENE2{$cnt2}=~ /(^1_(.+))/)) {
				my $header1=$1."_coding_".$SP2{$cnt2};
				$hash_join{$header1}=$rev_com;
			}elsif (($STRAND4{$cnt2} eq "-") and ($GENE2{$cnt2}=~ /(^2_(.+))/)) {
				my $header2=$1."_coding_".$SP2{$cnt2};
				$hash_join{$header2}=$rev_com;
			}elsif(($STRAND4{$cnt2} eq "+") and (($GENE2{$cnt2}!~ /1_/) and ($GENE2{$cnt2}!~ /2_/))){
				print $out_coding ">".$GENE2{$cnt2}."_coding_".$SP2{$cnt2}."\n".$str."\n";
			}elsif (($STRAND4{$cnt2} eq "+") and ($GENE2{$cnt2}=~ /(^1_(.+))/)) {
				my $header1=$1."_coding_".$SP2{$cnt2};
				$hash_join{$header1}=$str;
			}elsif (($STRAND4{$cnt2} eq "+") and ($GENE2{$cnt2}=~ /(^2_(.+))/)) {
				my $header2=$1."_coding_".$SP2{$cnt2};
				$hash_join{$header2}=$str;
			}
		}
		my ($join1,$join2,$header_join);
		foreach my $key (keys %hash_join) {
			$key=~ /^(\d)_((.+)_coding_(.+))/;
			$header_join=$2;
			if ($1==1) {
				$join1=$hash_join{$key};
			}elsif($1==2) {
				$join2=$hash_join{$key};
			}
		}
		print $out_coding ">".$header_join."\n".$join1.$join2."\n";
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
			next if ($END5{$_}+1) >= ($START5{$last}-1);
			push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'.$GENE3{$last}."\t".$STRAND5{$_}."/".$STRAND5{$last}."\t".($END5{$_}+1)."\t".($START5{$last}-1)."\t".$TYPE5{$_}."/".$TYPE5{$last}."\n";
		}
		foreach (keys %SP3){
			if ($_==$cnt3){
				push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'."end"."\t".$STRAND5{$_}."/"."?"."\t".($END5{$_}+1)."\t".$length."\t".$TYPE5{$_}."/"."?"."\n";
			}
		}
		@output3=();

		#extract_noncoding_regions
		my (%SP4,%GENE4,%STRAND6,%START6,%END6,%TYPE6,$seq2);
		my $cnt4=0;
		open(my $out_noncoding,">","$filename_base\_noncoding.fasta");
		while (@fasta2){
			my $header=shift @fasta2;
			$seq2=shift @fasta2;
		}
		foreach (@output4){
			$_=~ s/\r|\n//g;
			$cnt4++;
			($SP4{$cnt4},$GENE4{$cnt4},$STRAND6{$cnt4},$START6{$cnt4},$END6{$cnt4},$TYPE6{$cnt4})=(split /\s+/,$_)[0,1,2,3,4,5];
			my $str=substr($seq2,($START6{$cnt4}-1),($END6{$cnt4}-$START6{$cnt4}+1));
			print $out_noncoding ">".$GENE4{$cnt4}."_coding_".$SP4{$cnt4}."\n".$str."\n";
		}
		@output4=();
		close $out_noncoding;
	    unlink "$filename_base\_noncoding.fasta";
    }


	############################################################
	## extract_bed_CDS
	############################################################
#	if (!-e "$filename_base\_CDS.bed"){
	{
		#generate_bed_file
		my (@row_array,$species_name,$length,$element,@genearray,@output1);
		my $mark=0;
		my $tick=0;
		open (my $in_genbank,"<","$filename_base\_temp2");
		while (<$in_genbank>){
			$_=~ s/\r|\n//g;
			@row_array=split /\s+/,$_;
			if (/^LOCUS/i){
				$species_name=$latin_name;
				$length=$row_array[2];
			}elsif(/ {5}CDS {13}/){
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
					$element=$species_name.":"."1_".$1.":".$element;
					push @genearray,$element;
					$element=();
					$mark=0;
					$tick=0;
				}elsif ($tick==2) {
					$element=$species_name.":"."2_".$1.":".$element;
					push @genearray,$element;
					$element=();
					$mark=0;
					$tick=0;
				}
			}
		}
		close $in_genbank;
		foreach (@genearray){
			my @array=split /:/,$_;
			push @output1,"$array[0]\t$array[1]\t$array[2]\n";
		}
		@row_array=();
		@genearray=();
		my @output5=@output1;

		#put_fasta_sequence_in_array
	    my $flag=0;
	    my @sequence;
		my (@fas1,@fas2);
		open(my $in_genebank,"<","$filename_base\_temp2");
	    while (<$in_genebank>){
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
		close $in_genebank;
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
			}elsif ((defined $STRAND2{$_} ne "") and (defined $STRAND3{$_} eq "") and ($GENE1{$_} ne "rps12")) {# non-rps12
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
			}elsif ((defined $STRAND2{$_} ne "") and (defined $STRAND3{$_} eq "") and ($GENE1{$_} eq "rps12")) {# rps12 in IRb and IRa
				push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
				push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
			}elsif ((defined $STRAND2{$_} ne "") and (defined $STRAND3{$_} ne "") and ($GENE1{$_} ne "rps12")) {# ycf3 or clpP
				if (($STRAND1{$_} eq "-") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_})){
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-3"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
				}elsif(($STRAND1{$_} eq "-") and ($START1{$_} > $START2{$_}) and ($START2{$_} > $START3{$_})){
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
				}
			}elsif ((defined $STRAND2{$_} ne "") and (defined $STRAND3{$_} ne "") and ($GENE1{$_} eq "rps12")) {# rps12
				push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
				push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
				push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
			}
		}

		foreach (@output5){
			my @row=split /\t/,$_;
			my $abs1=abs($row[3]-$row[4]);
			my $abs2=abs($row[7]-$row[8]);
			my $abs3=abs($row[11]-$row[12]);
			if (($_=~ /rps12/) and ($abs1 > 50) and ($abs1 < 200)) {
				$_="$row[0]\t$row[1]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\n";
			}
			if (($_=~ /rps12/) and ($abs2 > 50) and ($abs2 < 200)) {
				$_="$row[0]\t$row[1]\t$row[6]\t$row[7]\t$row[8]\t$row[9]\n";
			}
			if (($_=~ /rps12/) and ($abs3 > 50) and ($abs3 < 200)) {
				$_="$row[0]\t$row[1]\t$row[10]\t$row[11]\t$row[12]\t$row[13]\n";
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
		open (my $out_bed,">","$filename_base\_CDS.bed");
		foreach (@output3){
			my @row=split /\t/,$_;
			my $abs=abs($row[3]-$row[4]);
			if ($_!~ /rps12/) {
				print $out_bed $_;
			}
			if (($_=~ /rps12/) and ($abs < 50)) {
				$_=~ s/rps12/rps12+2-2/;
				print $out_bed $_;
			}
			if (($_=~ /rps12/) and ($abs > 50) and ($abs < 200)) {
				$_=~ s/rps12/rps12+1/;
				print $out_bed $_;
			}
			if (($_=~ /rps12/) and ($abs > 200)) {
				$_=~ s/rps12/rps12+2-1/;
				print $out_bed $_;
			}
		}
		close $out_bed;

		#extract_CDS
		my (%SP2,%GENE2,%STRAND4,%START4,%END4,%TYPE4,%STRAND5,%START5,%END5,%TYPE5,%STRAND6,%START6,%END6,%TYPE6,$seq1);
		my $cnt2=0;
		open(my $out_coding,">","$filename_base\_CDS.fasta");

		while (@fasta1){
			my $header=shift @fasta1;
			$seq1=shift @fasta1;
		}
		my %hash_join;
		foreach (@output5){
			$_=~ s/\r|\n//g;
			$cnt2++;
			($SP2{$cnt2},$GENE2{$cnt2},$STRAND4{$cnt2},$START4{$cnt2},$END4{$cnt2},$TYPE4{$cnt2},$STRAND5{$cnt2},$START5{$cnt2},$END5{$cnt2},$TYPE5{$cnt2},$STRAND6{$cnt2},$START6{$cnt2},$END6{$cnt2},$TYPE6{$cnt2})=(split /\s+/,$_)[0,1,2,3,4,5,6,7,8,9,10,11,12,13];
			if (defined $STRAND5{$cnt2} eq "") {
	        	my $str=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
                my $rev_com=reverse $str;
                $rev_com=~ tr/ACGTacgt/TGCAtgca/;
	            if (($STRAND4{$cnt2} eq "-") and ($GENE2{$cnt2} ne "rps12") and (($GENE2{$cnt2}!~ /1_/) and ($GENE2{$cnt2}!~ /2_/))) {
	                print $out_coding ">".$GENE2{$cnt2}."_CDS_".$SP2{$cnt2}."\n".$rev_com."\n";
				}elsif (($STRAND4{$cnt2} eq "-") and ($GENE2{$cnt2} ne "rps12") and ($GENE2{$cnt2}=~ /(^1_(.+))/)) {
					my $header1=$1."_CDS_".$SP2{$cnt2};
					$hash_join{$header1}=$rev_com;
				}elsif (($STRAND4{$cnt2} eq "-") and ($GENE2{$cnt2} ne "rps12") and ($GENE2{$cnt2}=~ /(^2_(.+))/)) {
					my $header2=$1."_CDS_".$SP2{$cnt2};
					$hash_join{$header2}=$rev_com;
	            }elsif (($STRAND4{$cnt2} eq "-") and ($GENE2{$cnt2} eq "rps12")) {
	                print $out_coding ">".$GENE2{$cnt2}."+1"."_CDS_".$SP2{$cnt2}."\n".$rev_com."\n";
	            }elsif(($STRAND4{$cnt2} eq "+") and ($GENE2{$cnt2} ne "rps12") and (($GENE2{$cnt2}!~ /1_/) and ($GENE2{$cnt2}!~ /2_/))){
	                print $out_coding ">".$GENE2{$cnt2}."_CDS_".$SP2{$cnt2}."\n".$str."\n";
				}elsif (($STRAND4{$cnt2} eq "+") and ($GENE2{$cnt2} ne "rps12") and ($GENE2{$cnt2}=~ /(^1_(.+))/)) {
					my $header1=$1."_CDS_".$SP2{$cnt2};
					$hash_join{$header1}=$str;
				}elsif (($STRAND4{$cnt2} eq "+") and ($GENE2{$cnt2} ne "rps12") and ($GENE2{$cnt2}=~ /(^2_(.+))/)) {
					my $header2=$1."_CDS_".$SP2{$cnt2};
					$hash_join{$header2}=$str;
	            }elsif (($STRAND4{$cnt2} eq "+") and ($GENE2{$cnt2} eq "rps12")) {
	                print $out_coding ">".$GENE2{$cnt2}."+1"."_CDS_".$SP2{$cnt2}."\n".$str."\n";
	            }
	        }
		}
		my ($join1,$join2,$header_join);
		foreach my $key (keys %hash_join) {
			$key=~ /^(\d)_((.+)_CDS_(.+))/;
			$header_join=$2;
			if ($1==1) {
				$join1=$hash_join{$key};
			}elsif($1==2) {
				$join2=$hash_join{$key};
			}
		}
		if (defined $header_join) {
			print $out_coding ">".$header_join."\n".$join1.$join2."\n";
		}
		close $out_coding;
		@output5=();

		#generate_intergenic_ranges
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
			next if ((($END7{$_}+1) >= ($START7{$last1}-1)) and (($END7{$_}+1) < ($END7{$last1}-1)));
			next if (($_ > 1) and (($END7{$_}+1) < ($END7{$last0}-1)) and (($END7{$_}+1) < ($START7{$last2}-1)));
			if ((($END7{$_}+1) >= ($START7{$last1}-1)) and (($END7{$_}+1) >= ($END7{$last1}-1))){
	    		push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'.$GENE3{$last2}."\t".$STRAND7{$_}."/".$STRAND7{$last2}."\t".($END7{$_}+1)."\t".($START7{$last2}-1)."\t".$TYPE7{$_}."/".$TYPE7{$last2}."\n";
	        }else{
	    		push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'.$GENE3{$last1}."\t".$STRAND7{$_}."/".$STRAND7{$last1}."\t".($END7{$_}+1)."\t".($START7{$last1}-1)."\t".$TYPE7{$_}."/".$TYPE7{$last1}."\n";
	        }
		}
		foreach (keys %SP3){
			if ($_==$cnt3){
				push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'."end"."\t".$STRAND7{$_}."/"."?"."\t".($END7{$_}+1)."\t".$length."\t".$TYPE7{$_}."/"."?"."\n";
			}
		}
		@output3=();

		#extract_intergenic
		my (%SP4,%GENE4,%STRAND8,%START8,%END8,%TYPE8,$seq2);
		my $cnt4=0;
		open(my $out_noncoding,">","$filename_base\_intergenic.fasta");
		while (@fasta2){
			my $header=shift @fasta2;
			$seq2=shift @fasta2;
		}
		foreach (@output4){
			$_=~ s/\r|\n//g;
			$cnt4++;
			($SP4{$cnt4},$GENE4{$cnt4},$STRAND8{$cnt4},$START8{$cnt4},$END8{$cnt4},$TYPE8{$cnt4})=(split /\s+/,$_)[0,1,2,3,4,5];
			my $str=substr($seq2,($START8{$cnt4}-1),($END8{$cnt4}-$START8{$cnt4}+1));
			print $out_noncoding ">".$GENE4{$cnt4}."_CDS_".$SP4{$cnt4}."\n".$str."\n";
		}
		@output4=();
		close $out_noncoding;
	    unlink "$filename_base\_intergenic.fasta";
    }
	unlink "$filename_base\_temp1";
	unlink "$filename_base\_temp2";
}

my $blast_reference1_random = $reference_path."/".$random."_reference1_blast";
my $blast_reference2_random = $reference_path."/".$random."_reference2_blast";
my $blast_reference3_random = $reference_path."/".$random."_reference3_blast";
my $blast_reference4_random = $reference_path."/".$random."_reference4_blast";

my $screenlog = $output_directory."/"."screen.log";
my $screenperf = $output_directory."/"."screen.perf";


############################################################
## combine_gene_coding_CDS
############################################################
my $dir;
opendir($dir,$reference_directory);
my @input_directory=readdir $dir;
close $dir;
#unlink ("gene.fasta");
#unlink ("coding.fasta");
#unlink ("CDS.fasta");

my (@gene_fasta,@coding_fasta,@CDS_fasta);
foreach my $file (@input_directory){
	if($file=~/_gene.fasta/){
		open (my $input_gene,"<","$reference_directory/$file");
		#open (my $output_gene,">>","gene.fasta");
		while(<$input_gene>){
			#print $output_gene $_;
			push @gene_fasta, $_;
		}
		close $input_gene;
		#close $output_gene;
	}
	if($file=~/_coding.fasta/){
		open (my $input_coding,"<","$reference_directory/$file");
		#open (my $output_coding,">>","coding.fasta");
		while(<$input_coding>){
			#print $output_coding $_;
			push @coding_fasta, $_;
		}
		close $input_coding;
		#close $output_coding;
	}
	if($file=~/_CDS.fasta/){
		open (my $input_CDS,"<","$reference_directory/$file");
		#open (my $output_CDS,">>","CDS.fasta");
		while(<$input_CDS>){
			#print $output_CDS $_;
			push @CDS_fasta, $_;
		}
		close $input_CDS;
		#close $output_CDS;
	}
}
foreach my $reference_filename (@reference_filename) {
	unlink("$reference_filename\_CDS.bed");
	unlink("$reference_filename\_CDS.fasta");
	unlink("$reference_filename\_coding.bed");
	unlink("$reference_filename\_coding.fasta");
	unlink("$reference_filename\_gene.bed");
	unlink("$reference_filename\_gene.fasta");
}


#open (my $in_gene,"<","gene.fasta");
#open (my $in_coding,"<","coding.fasta");
#open (my $in_CDS,"<","CDS.fasta");
#open (my $out_all1,">","all1.fasta");
my @all1;
#while (<$in_gene>){
while (@gene_fasta){
	#print $out_all1 $_;
	$_=shift @gene_fasta;
	push @all1, $_;
}
#while (<$in_coding>){
while (@coding_fasta){
	#print $out_all1 $_;
	$_=shift @coding_fasta;
	push @all1, $_;
}
#while (<$in_CDS>){
while (@CDS_fasta){
	#print $out_all1 $_;
	$_=shift @CDS_fasta;
	push @all1, $_;
}
#close $in_gene;
#close $in_coding;
#close $in_CDS;
#close $out_all1;

#from A-UGC to trnA-UGC
#open (my $in_all1,"<","all1.fasta");
#open (my $out_all2,">","all2.fasta");
my @all2;
my ($header_in_all1,$sequence_in_all1);
#while (defined ($header_in_all1=<$in_all1>) && defined ($sequence_in_all1=<$in_all1>)){
while (@all1){
	$header_in_all1=shift @all1;
	$sequence_in_all1=shift @all1;
	if ($header_in_all1=~ />(.|..)-/){
		$header_in_all1=~ s/>/>trn/g;
		#print $out_all2 $header_in_all1.$sequence_in_all1;
		push @all2, $header_in_all1;
		push @all2, $sequence_in_all1;		
	}else{
		#print $out_all2 $header_in_all1.$sequence_in_all1;
		push @all2, $header_in_all1;
		push @all2, $sequence_in_all1;
	}
}
#close $in_all1;
#close $out_all2;
#unlink("all1.fasta");


#codon_table
my %hash_codon=("---"=>"-","TAA"=>"*","TAG"=>"*","TGA"=>"*","TCA"=>"S","TCC"=>"S","TCG"=>"S","TCT"=>"S","TTC"=>"F","TTT"=>"F","TTA"=>"L","TTG"=>"L","TAC"=>"Y","TAT"=>"Y","TGC"=>"C","TGT"=>"C","TGG"=>"W","CTA"=>"L","CTC"=>"L","CTG"=>"L","CTT"=>"L","CCA"=>"P","CCC"=>"P","CCG"=>"P","CCT"=>"P","CAC"=>"H","CAT"=>"H","CAA"=>"Q","CAG"=>"Q","CGA"=>"R","CGC"=>"R","CGG"=>"R","CGT"=>"R","ATA"=>"I","ATC"=>"I","ATT"=>"I","ATG"=>"M","ACA"=>"T","ACC"=>"T","ACG"=>"T","ACT"=>"T","AAC"=>"N","AAT"=>"N","AAA"=>"K","AAG"=>"K","AGC"=>"S","AGT"=>"S","AGA"=>"R","AGG"=>"R","GTA"=>"V","GTC"=>"V","GTG"=>"V","GTT"=>"V","GCA"=>"A","GCC"=>"A","GCG"=>"A","GCT"=>"A","GAC"=>"D","GAT"=>"D","GAA"=>"E","GAG"=>"E","GGA"=>"G","GGC"=>"G","GGG"=>"G","GGT"=>"G");

#product_of_protein-coding-gene_and_tRNA,e.g., $hash_tRNA{trnQ-UUG}=tRNA-Gln
#open (my $in_product,"<","$Bin/product.txt");
#my %hash_product;
#while (<$in_product>){
#	$_=~ s/\r|\n//g;
#	my @product=split /\t/,$_;
#	$hash_product{$product[0]}=$product[1];
#}
#close $in_product;
#print Dumper \%hash_product;
my %hash_product=("accD"=>"acetyl-CoA carboxylase carboxyltransferase beta subunit","atpA"=>"ATP synthase CF1 alpha subunit","atpB"=>"ATP synthase CF1 beta subunit","atpE"=>"ATP synthase CF1 epsilon subunit","atpF"=>"ATP synthase CF0 subunit I","atpH"=>"ATP synthase CF0 subunit III","atpI"=>"ATP synthase CF0 subunit IV","ccsA"=>"cytochrome c heme attachment protein","cemA"=>"envelope membrane protein","chlB"=>"photochlorophyllide reductase subunit B","chlL"=>"photochlorophyllide reductase subunit L","chlN"=>"photochlorophyllide reductase subunit N","clpP"=>"clp protease proteolytic subunit","infA"=>"translation initiation factor 1","matK"=>"maturase K","ndhA"=>"NADH-plastoquinone oxidoreductase subunit 1","ndhB"=>"NADH-plastoquinone oxidoreductase subunit 2","ndhC"=>"NADH-plastoquinone oxidoreductase subunit 3","ndhD"=>"NADH-plastoquinone oxidoreductase subunit 4","ndhE"=>"NADH-plastoquinone oxidoreductase subunit 4L","ndhF"=>"NADH-plastoquinone oxidoreductase subunit 5","ndhG"=>"NADH-plastoquinone oxidoreductase subunit 6","ndhH"=>"NADH-plastoquinone oxidoreductase subunit 7","ndhI"=>"NADH-plastoquinone oxidoreductase subunit I","ndhJ"=>"NADH-plastoquinone oxidoreductase subunit J","ndhK"=>"NADH-plastoquinone oxidoreductase subunit K","petA"=>"cytochrome f","petB"=>"cytochrome b6","petD"=>"cytochrome b6/f complex subunit IV","petG"=>"cytochrome b6/f complex subunit V","petL"=>"cytochrome b6/f complex subunit VI","petN"=>"cytochrome b6/f complex subunit VIII","psaA"=>"photosystem I P700 apoprotein A1","psaB"=>"photosystem I P700 apoprotein A2","psaC"=>"photosystem I subunit VII","psaI"=>"photosystem I subunit VIII","psaJ"=>"photosystem I subunit IX","psaM"=>"photosystem I protein M","psbA"=>"photosystem II protein D1","psbB"=>"photosystem II CP47 chlorophyll apoprotein","psbC"=>"photosystem II CP43 chlorophyll apoprotein","psbD"=>"photosystem II protein D2","psbE"=>"photosystem II cytochrome b559 alpha subunit","psbF"=>"photosystem II cytochrome b559 beta subunit","psbH"=>"photosystem II phosphoprotein","psbI"=>"photosystem II protein I","psbJ"=>"photosystem II protein J","psbK"=>"photosystem II protein K","psbL"=>"photosystem II protein L","psbM"=>"photosystem II protein M","psbN"=>"photosystem II protein N","psbT"=>"photosystem II protein T","psbZ"=>"photosystem II protein Z","rbcL"=>"ribulose bisophosphate carboxylase","rpl14"=>"ribosomal protein L14","rpl16"=>"50S ribosomal protein L16","rpl2"=>"ribosomal protein L2","rpl20"=>"ribosomal protein L20","rpl22"=>"ribosomal protein L22","rpl23"=>"ribosomal protein L23","rpl32"=>"ribosomal protein L32","rpl33"=>"ribosomal protein L33","rpl36"=>"ribosomal protein L36","rpoA"=>"RNA polymerase alpha subunit","rpoB"=>"RNA polymerase beta subunit","rpoC1"=>"RNA polymerase beta' subunit","rpoC2"=>"RNA polymerase beta'' subunit","rps11"=>"ribosomal protein S11","rps12"=>"ribosomal protein S12","rps12+1"=>"ribosomal protein S12","rps12+2"=>"ribosomal protein S12","rps14"=>"ribosomal protein S14","rps15"=>"ribosomal protein S15","rps16"=>"ribosomal protein S16","rps18"=>"ribosomal protein S18","rps19"=>"ribosomal protein S19","rps2"=>"ribosomal protein S2","rps3"=>"ribosomal protein S3","rps4"=>"ribosomal protein S4","rps7"=>"ribosomal protein S7","rps8"=>"ribosomal protein S8","rrn16"=>"16S ribosomal RNA","rrn23"=>"23S ribosomal RNA","rrn4.5"=>"4.5S ribosomal RNA","rrn5"=>"5S ribosomal RNA","ycf1"=>"hypothetical protein RF1","ycf2"=>"hypothetical protein RF2","ycf3"=>"photosystem I assembly protein Ycf3","ycf4"=>"photosystem I assembly protein Ycf4","ycf15"=>"hypothetical protein RF15","ycf68"=>"hypothetical protein RF68","trnA-UGC"=>"tRNA-Ala","trnC-GCA"=>"tRNA-Cys","trnD-GUC"=>"tRNA-Asp","trnE-UUC"=>"tRNA-Glu","trnF-GAA"=>"tRNA-Phe","trnfM-CAU"=>"tRNA-Met","trnG-GCC"=>"tRNA-Gly","trnG-UCC"=>"tRNA-Gly","trnH-GUG"=>"tRNA-His","trnI-CAU"=>"tRNA-Ile","trnI-GAU"=>"tRNA-Ile","trnK-CUU"=>"tRNA-Asn","trnK-UUU"=>"tRNA-Lys","trnL-CAA"=>"tRNA-Leu","trnL-UAA"=>"tRNA-Leu","trnL-UAG"=>"tRNA-Leu","trnM-CAU"=>"tRNA-Met","trnN-GUU"=>"tRNA-Asn","trnP-GGG"=>"tRNA-Pro","trnP-UGG"=>"tRNA-Pro","trnQ-UUG"=>"tRNA-Gln","trnR-ACG"=>"tRNA-Arg","trnR-CCG"=>"tRNA-Arg","trnR-UCU"=>"tRNA-Arg","trnS-GCU"=>"tRNA-Ser","trnS-GGA"=>"tRNA-Ser","trnS-UGA"=>"tRNA-Ser","trnT-GGU"=>"tRNA-Thr","trnT-UGU"=>"tRNA-Thr","trnV-GAC"=>"tRNA-Val","trnV-UAC"=>"tRNA-Val","trnW-CCA"=>"tRNA-Trp","trnY-GUA"=>"tRNA-Tyr");


#from CDS to CDS_aa(no intron PCG)
#from coding to coding_aa(intron PCG)
#open (my $in_all2,"<","all2.fasta");
my @all22=@all2;
#open (my $out_all3,">","all3.fasta");
my @all3;
my ($header_in_all2,$sequence_in_all2);
#while (defined ($header_in_all2=<$in_all2>) && defined ($sequence_in_all2=<$in_all2>)){
while (@all2){
	$header_in_all2=shift @all2;
	$sequence_in_all2=shift @all2;
	$header_in_all2=~ s/\r|\n//g;
	$sequence_in_all2=~ s/\r|\n//g;

	my $length=length $sequence_in_all2;
	my $aa;

	if ($header_in_all2=~ /_CDS/){
		for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
			my $codon=substr ($sequence_in_all2,$i,3);
			$codon=uc $codon;
			if (exists $hash_codon{$codon}){
				$aa.=$hash_codon{$codon};
			}else{
				$aa.="X";
				my $j=$i+1;
			}
		}
		$header_in_all2=~ s/_CDS/_CDS_aa/g;
		#print $out_all3 "$header_in_all2\n$aa\n";
		push @all3, $header_in_all2;
		push @all3, $aa;
	}elsif ($header_in_all2=~ /^>(?!trn)(.+)-(\d)_coding/){
		if ($2 == 1) {
			for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
				my $codon=substr ($sequence_in_all2,$i,3);
				$codon=uc $codon;
				if (exists $hash_codon{$codon}){
					$aa.=$hash_codon{$codon};
				}else{
					$aa.="X";
					my $j=$i+1;
				}
			}
		}elsif ($2 == 2) {
			if ($length % 3 == 0) {
				for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
					my $codon=substr ($sequence_in_all2,$i,3);
					$codon=uc $codon;
					if (exists $hash_codon{$codon}){
						$aa.=$hash_codon{$codon};
					}else{
						$aa.="X";
						my $j=$i+1;
					}
				}
			}elsif ($length % 3 == 1) {
				for (my $i=1;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
					my $codon=substr ($sequence_in_all2,$i,3);
					$codon=uc $codon;
					if (exists $hash_codon{$codon}){
						$aa.=$hash_codon{$codon};
					}else{
						$aa.="X";
						my $j=$i+1;
					}
				}
			}elsif ($length % 3 == 2) {
				for (my $i=2;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
					my $codon=substr ($sequence_in_all2,$i,3);
					$codon=uc $codon;
					if (exists $hash_codon{$codon}){
						$aa.=$hash_codon{$codon};
					}else{
						$aa.="X";
						my $j=$i+1;
					}
				}
			}
		}elsif ($2 == 3) {
			for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
				my $codon=substr ($sequence_in_all2,$i,3);
				$codon=uc $codon;
				if (exists $hash_codon{$codon}){
					$aa.=$hash_codon{$codon};
				}else{
					$aa.="X";
					my $j=$i+1;
				}
			}
		}
		$header_in_all2=~ s/_coding/_coding_aa/g;
		#print $out_all3 "$header_in_all2\n$aa\n";
		push @all3, $header_in_all2;
		push @all3, $aa;
	}else{
		#print $out_all3 "$header_in_all2\n$sequence_in_all2\n";
		push @all3, $header_in_all2;
		push @all3, $sequence_in_all2;
	}
}
#close $in_all2;
#close $out_all3;

#delete xxx_coding of no intron gene(including PCG and RNA):
#if exists accD_coding and accD_gene,delete accD_coding(no intron PCG);
#if exists trnQ-UUG_coding and trnQ-UUG_gene,delete trnQ-UUG_coding(no intron tRNA);
#if exists rrn16_coding and rrn16_gene,delete rrn16_coding(no intron rRNA);
#open (my $in_all3,"<","all3.fasta");
#open (my $out_reference1,">","reference1.fasta");
#open (my $out_reference2,">","reference2.fasta");
#open (my $out_reference3,">","reference3.fasta");
#open (my $out_reference4,">","reference4.fasta");
my $reference1_random = $reference_path."/".$random."_reference1";
my $reference2_random = $reference_path."/".$random."_reference2";
my $reference3_random = $reference_path."/".$random."_reference3";
my $reference4_random = $reference_path."/".$random."_reference4";
open (my $out_reference1,">",$reference1_random);
open (my $out_reference2,">",$reference2_random);
open (my $out_reference3,">",$reference3_random);
open (my $out_reference4,">",$reference4_random);

my ($header_in_all3,$sequence_in_all3,%hash1_in_all3,%hash2_in_all3,%hash3_in_all3,%hash4_in_all3);
my %gene_number_ref;
#while (defined ($header_in_all3=<$in_all3>) and defined ($sequence_in_all3=<$in_all3>)){
while (@all3){
	$header_in_all3=shift @all3;
	$sequence_in_all3=shift @all3;
	$header_in_all3=~ s/\r|\n//g;
	$sequence_in_all3=~ s/\r|\n//g;

	if ($header_in_all3=~ /^>(trn)(.+)-1_coding/ or $header_in_all3=~ /^>(trn)(.+)-2_coding/ or $header_in_all3=~ /^>(trn)(.+)-3_coding/){
		$hash1_in_all3{$header_in_all3}=$sequence_in_all3;
	}
	if ($header_in_all3=~ /^>((trn)(.+))_gene/ or $header_in_all3=~ /^>((rrn)(.+))_gene/){
		$hash1_in_all3{$header_in_all3}=$sequence_in_all3;
		my $genes=$1;
		$gene_number_ref{$genes}++;
	}
	if ($header_in_all3=~ /^>((?!trn)(.+))_gene/ and $header_in_all3=~ /^>((?!rrn)(.+))_gene/){
		$hash2_in_all3{$header_in_all3}=$sequence_in_all3;
		my $genes=$1;
		$genes=~ s/rps12\+1/rps12/;
		$genes=~ s/rps12\+2/rps12/;
		$gene_number_ref{$genes}++;
	}
	if ($header_in_all3=~ /^>(.+)_CDS_aa/){
		$hash3_in_all3{$header_in_all3}=$sequence_in_all3;
	}
	if ($header_in_all3=~ /^>(.+)_coding_aa/){
		$hash4_in_all3{$header_in_all3}=$sequence_in_all3;
	}
}
foreach my $h1 (sort keys %hash1_in_all3){
	print $out_reference1 "$h1\n$hash1_in_all3{$h1}\n";
}
foreach my $h2 (sort keys %hash2_in_all3){
	print $out_reference2 "$h2\n$hash2_in_all3{$h2}\n";
}
foreach my $h3 (sort keys %hash3_in_all3){
	print $out_reference3 "$h3\n$hash3_in_all3{$h3}\n";
}
foreach my $h4 (sort keys %hash4_in_all3){
	print $out_reference4 "$h4\n$hash4_in_all3{$h4}\n";
}
#close $in_all3;
close $out_reference1;
close $out_reference2;
close $out_reference3;
close $out_reference4;
#unlink("all3.fasta");
#unlink ("gene.fasta");
#unlink ("coding.fasta");
#unlink ("CDS.fasta");
%hash1_in_all3=();
%hash2_in_all3=();
%hash3_in_all3=();
%hash4_in_all3=();

my $now3=&gettime;
print "$now3 || Finish extracting annotations from reference!";
print "\n";


############################################################
## sequence_dealing_in_batch_mode
############################################################
my $pattern4=".fasta";
my $pattern5=".fas";
my $pattern6=".fa";
my $pattern7=".fsa";
my @sequence_filenames;
find(\&target2,$sequence_directory);
sub target2{
	if (/$pattern4/ or /$pattern5/ or /$pattern6/ or /$pattern7/){
		push @sequence_filenames,"$File::Find::name";
	}
	return;
}

my $j=0;


my $target_path;
my $IR_temp_random;
my $blast_IR_temp_random;
while (@sequence_filenames) {
	$j++;
	my $filename_fasta=shift @sequence_filenames;
	my $target_name=substr($filename_fasta,0,rindex($filename_fasta,"\."));
	$target_name=~ s/\s+/_/gg;
	my $filename=substr($target_name,rindex($target_name,"\/")+1);
	$target_path=$target_name;
	$target_path=~ s/$filename//g;
	my $target_name_random = $target_path."/".$random."_".$filename;
	$IR_temp_random = $target_path."/".$random."_IR_temp";
	$blast_IR_temp_random = $target_path."/".$random."_IR_temp_blast";

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

	#fasta_sequence_with_length_cp*2
	open (my $in_fasta_2,"<",$target_name_random);
	#open (my $out_fasta_2,">","fasta_temp");
	open (my $out_fasta_2,">",$IR_temp_random);
	while (<$in_fasta_2>) {
		$_=~ s/\r|\n//g;
		if ($_=~ /^>/) {
			print $out_fasta_2 "$_\n";
		}elsif ($_!~ /^>/) {
			print $out_fasta_2 $_.$_."\n";
		}
	}
	close $in_fasta_2;
	close $out_fasta_2;

	#fasta_sequence
	open (my $in_fasta,"<",$target_name_random);
	my @fasta;
	while (<$in_fasta>){
		$_=~ s/\r|\n//g;
		push @fasta,$_;
	}
	close $in_fasta;
	my ($header,$sequence,$length_cp);
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
	}
	my $rev_coms=reverse $sequence;
	$rev_coms=~ tr/ACGTacgt/TGCAtgca/;


	my $now4=&gettime;
	print "\n";
	print "$now4 || Begin annotating the $j sequence: $filename";
	print "\n";
	print "$now4 || Begin blasting reference to sequence!\n";
	############################################################
	## blast_reference_to_sequence
	############################################################
	#rRNA and no intron tRNA(_gene)
	#intron tRNA(_gene and -1/-2_coding)
	#no intron PCG (_gene and _CDS_aa)
	#intron PCG(_gene and -1/-2/-3_coding_aa)(short exon of rpl16,petB and petD)

	if ($osname eq "MSWin32") {
		system ("makeblastdb.exe -in $target_name_random -hash_index -dbtype nucl -blastdb_version 4 -logfile $screenlog");
		#-max_hsps 1(nucleotide,RNA,including _gene RNA and -1/-2_coding tRNA)
		system ("blastn.exe -task blastn -query $reference1_random -db $target_name_random -outfmt 6 -max_hsps 1 -max_target_seqs 1 -out $blast_reference1_random");# -evalue 0.001 or 0.01
		#-max_hsps 1(nucleotide,PCG,including _gene PCG)
		system ("blastn.exe -task blastn -query $reference2_random -db $target_name_random -outfmt 6 -max_hsps 1 -max_target_seqs 1 -out $blast_reference2_random");# -evalue 0.001 or 0.01 -perc_identity 50?
		#system ("tblastx.exe -query $reference2_random -db $target_name_random -outfmt 6 -max_hsps 1 -max_target_seqs 1 -qcov_hsp_perc 50 -out $blast_reference2_random");# -evalue 0.001 or 0.01
		#-max_hsps 1(amino acid,no intron PCG,including _CDS_aa PCG)
		system ("tblastn.exe -task tblastn -query $reference3_random -db $target_name_random -outfmt 6 -max_hsps 1 -max_target_seqs 1 -out $blast_reference3_random");# -evalue 0.001 or 0.01 -max_intron_length 1000?
		#-max_hsps 1(amino acid,intron PCG,including -1/-2/-3_coding_aa PCG)
		system ("tblastn.exe -task tblastn -query $reference4_random -db $target_name_random -outfmt 6 -max_hsps 1 -max_target_seqs 1 -out $blast_reference4_random");# -evalue 0.001 or 0.01 -qcov_hsp_perc 20?
		#IRb and IRa
		system ("makeblastdb.exe -in $IR_temp_random -hash_index -dbtype nucl -blastdb_version 4 -logfile $screenlog");
		system ("blastn.exe -task blastn -query $IR_temp_random -db $IR_temp_random -outfmt 6 -perc_identity 99 -out $blast_IR_temp_random");
	}elsif ($osname eq "cygwin") {
		system ("makeblastdb -in $target_name_random -hash_index -dbtype nucl -blastdb_version 4 -logfile $screenlog");
		#-max_hsps 1(nucleotide,RNA,including _gene RNA and -1/-2_coding tRNA)
		system ("blastn -task blastn -query $reference1_random -db $target_name_random -outfmt 6 -max_hsps 1 -max_target_seqs 1 -out $blast_reference1_random");# -evalue 0.001 or 0.01
		#-max_hsps 1(nucleotide,PCG,including _gene PCG)
		system ("blastn -task blastn -query $reference2_random -db $target_name_random -outfmt 6 -max_hsps 1 -max_target_seqs 1 -out $blast_reference2_random");# -evalue 0.001 or 0.01 -perc_identity 50?
		#system ("tblastx -query $reference2_random -db $target_name_random -outfmt 6 -max_hsps 1 -max_target_seqs 1 -qcov_hsp_perc 50 -out $blast_reference2_random");# -evalue 0.001 or 0.01
		#-max_hsps 1(amino acid,no intron PCG,including _CDS_aa PCG)
		system ("tblastn -task tblastn -query $reference3_random -db $target_name_random -outfmt 6 -max_hsps 1 -max_target_seqs 1 -out $blast_reference3_random");# -evalue 0.001 or 0.01 -max_intron_length 1000?
		#-max_hsps 1(amino acid,intron PCG,including -1/-2/-3_coding_aa PCG)
		system ("tblastn -task tblastn -query $reference4_random -db $target_name_random -outfmt 6 -max_hsps 1 -max_target_seqs 1 -out $blast_reference4_random");# -evalue 0.001 or 0.01 -qcov_hsp_perc 20?
		#IRb and IRa
		system ("makeblastdb -in $IR_temp_random -hash_index -dbtype nucl -blastdb_version 4 -logfile $screenlog");
		system ("blastn -task blastn -query $IR_temp_random -db $IR_temp_random -outfmt 6 -perc_identity 99 -out $blast_IR_temp_random");
	}elsif ($osname eq "linux") {
		system ("makeblastdb -in $target_name_random -hash_index -dbtype nucl -blastdb_version 4 -logfile $screenlog");
		#-max_hsps 1(nucleotide,RNA,including _gene RNA and -1/-2_coding tRNA)
		system ("blastn -task blastn -query $reference1_random -db $target_name_random -outfmt 6 -max_hsps 1 -max_target_seqs 1 -out $blast_reference1_random");# -evalue 0.001 or 0.01
		#-max_hsps 1(nucleotide,PCG,including _gene PCG)
		system ("blastn -task blastn -query $reference2_random -db $target_name_random -outfmt 6 -max_hsps 1 -max_target_seqs 1 -out $blast_reference2_random");# -evalue 0.001 or 0.01 -perc_identity 50?
		#system ("tblastx -query $reference2_random -db $target_name_random -outfmt 6 -max_hsps 1 -max_target_seqs 1 -qcov_hsp_perc 50 -out $blast_reference2_random");# -evalue 0.001 or 0.01
		#-max_hsps 1(amino acid,no intron PCG,including _CDS_aa PCG)
		system ("tblastn -task tblastn -query $reference3_random -db $target_name_random -outfmt 6 -max_hsps 1 -max_target_seqs 1 -out $blast_reference3_random");# -evalue 0.001 or 0.01 -max_intron_length 1000?
		#-max_hsps 1(amino acid,intron PCG,including -1/-2/-3_coding_aa PCG)
		system ("tblastn -task tblastn -query $reference4_random -db $target_name_random -outfmt 6 -max_hsps 1 -max_target_seqs 1 -out $blast_reference4_random");# -evalue 0.001 or 0.01 -qcov_hsp_perc 20?
		#IRb and IRa
		system ("makeblastdb -in $IR_temp_random -hash_index -dbtype nucl -blastdb_version 4 -logfile $screenlog");
		system ("blastn -task blastn -query $IR_temp_random -db $IR_temp_random -outfmt 6 -perc_identity 99 -out $blast_IR_temp_random");
	}elsif ($osname eq "darwin") {
		system ("makeblastdb -in $target_name_random -hash_index -dbtype nucl -blastdb_version 4 -logfile $screenlog");
		#-max_hsps 1(nucleotide,RNA,including _gene RNA and -1/-2_coding tRNA)
		system ("blastn -task blastn -query $reference1_random -db $target_name_random -outfmt 6 -max_hsps 1 -max_target_seqs 1 -out $blast_reference1_random");# -evalue 0.001 or 0.01
		#-max_hsps 1(nucleotide,PCG,including _gene PCG)
		system ("blastn -task blastn -query $reference2_random -db $target_name_random -outfmt 6 -max_hsps 1 -max_target_seqs 1 -out $blast_reference2_random");# -evalue 0.001 or 0.01 -perc_identity 50?
		#system ("tblastx -query $reference2_random -db $target_name_random -outfmt 6 -max_hsps 1 -max_target_seqs 1 -qcov_hsp_perc 50 -out $blast_reference2_random");# -evalue 0.001 or 0.01
		#-max_hsps 1(amino acid,no intron PCG,including _CDS_aa PCG)
		system ("tblastn -task tblastn -query $reference3_random -db $target_name_random -outfmt 6 -max_hsps 1 -max_target_seqs 1 -out $blast_reference3_random");# -evalue 0.001 or 0.01 -max_intron_length 1000?
		#-max_hsps 1(amino acid,intron PCG,including -1/-2/-3_coding_aa PCG)
		system ("tblastn -task tblastn -query $reference4_random -db $target_name_random -outfmt 6 -max_hsps 1 -max_target_seqs 1 -out $blast_reference4_random");# -evalue 0.001 or 0.01 -qcov_hsp_perc 20?
		#IRb and IRa
		system ("makeblastdb -in $IR_temp_random -hash_index -dbtype nucl -blastdb_version 4 -logfile $screenlog");
		system ("blastn -task blastn -query $IR_temp_random -db $IR_temp_random -outfmt 6 -perc_identity 99 -out $blast_IR_temp_random");
	}


	unlink("$target_name_random");
	unlink("$target_name_random.nhd");
	unlink("$target_name_random.nhi");
	unlink("$target_name_random.nhr");
	unlink("$target_name_random.nin");
	unlink("$target_name_random.nog");
	unlink("$target_name_random.nsd");
	unlink("$target_name_random.nsi");
	unlink("$target_name_random.nsq");
	#unlink ("reference_temp");
	unlink ("$IR_temp_random");
	unlink ("$IR_temp_random.nhd");
	unlink ("$IR_temp_random.nhi");
	unlink ("$IR_temp_random.nhr");
	unlink ("$IR_temp_random.nin");
	unlink ("$IR_temp_random.nog");
	unlink ("$IR_temp_random.nsd");
	unlink ("$IR_temp_random.nsi");
	unlink ("$IR_temp_random.nsq");
	unlink("$screenperf");
	my $now5=&gettime;
	print "$now5 || Finish blasting reference to sequence!";
	print "\n";

	#select IR region
	open (my $input_IR,"<",$blast_IR_temp_random);
	my %IR;
	while (<$input_IR>) {
		$_=~ s/\r|\n//g;
		my ($c1,$c2,$c3,$ir_length,$c5,$c6,$qs,$qe,$ss,$se,$c11,$c12)=split /\t/,$_;
		if (($ir_length != 2 * $length_cp) and ($ir_length != $length_cp) and ($ir_length >= $inverted_repeat) and ($qe <= $length_cp) and ($se <= $length_cp) and ($qs < $ss)) {
			$IR{$ir_length}=$qs."\t".$qe."\t".$ss."\t".$se;
		}
	}
	close $input_IR;
	unlink ("$blast_IR_temp_random");
	my (@IR_length,@boundary);
	foreach my $key (sort {$b <=> $a} keys %IR) {
		push @IR_length,$key;
		push @boundary,$IR{$key};
	}
	my ($AA,$BB,$CC,$DD)=split /\t/,$boundary[0];
	my ($JLB,$JSB,$JLA,$JSA);
	if ($CC <= $length_cp) {
		$JLB=$AA;
		$JSB=$BB;
		$JLA=$CC;
		$JSA=$DD;
	}elsif ($CC > $length_cp) {
		$JLB=$AA;
		$JSB=$BB;
		$JLA=$CC-$length_cp;
		$JSA=$DD;
	}
	#print "$JLB\t$JSB\t$JLA\t$JSA\n";

	#set threshold of similarity value
	open (my $in_bt_ref1,"<",$blast_reference1_random);
	open (my $in_bt_ref2,"<",$blast_reference2_random);
	open (my $in_bt_ref3,"<",$blast_reference3_random);
	open (my $in_bt_ref4,"<",$blast_reference4_random);
	#open (my $out_bt_ref,">>","reference_temp");
	my @reference_temp;
	while(<$in_bt_ref1>){
		$_=~ s/\r|\n//g;
		#print $out_bt_ref "$_\n";
		push @reference_temp, "$_\n";
	}

	my (%gene_remain1,%gene_remain2,%gene_remove1,%gene_remove2,%gene_loss,%gene_loss1,%gene_loss2);
	while(<$in_bt_ref3>){# _CDS_aa
		$_=~ s/\r|\n//g;
		my ($gene,$taxon,$similarity)=split /\t/,$_;
		if (($similarity >= $similarity_value) and ($gene=~ /(.*)_CDS_aa_(.*)/)) {
			$gene_remain1{$1}++;
			#print $out_bt_ref "$_\n";
			push @reference_temp, "$_\n";
		}elsif (($similarity < $similarity_value) and ($gene=~ /(.*)_CDS_aa_(.*)/) and (($1=~ /ycf1/) or ($1=~ /ycf2/))) {
			$gene_remain1{$1}++;
			#print $out_bt_ref "$_\n";
			push @reference_temp, "$_\n";
		}elsif (($similarity < $similarity_value) and ($gene=~ /(.*)_CDS_aa_(.*)/) and (($1!~ /ycf1/) and ($1!~ /ycf2/))) {
			$gene_remove1{$1}++;
		}
	}
	foreach my $key (keys %gene_remove1) {
		if (!exists $gene_remain1{$key}) {
			$gene_loss1{$key}=$gene_remove1{$key};
		}
	}

	while(<$in_bt_ref4>){# -1/2/3_coding_aa
		$_=~ s/\r|\n//g;
		my ($gene,$taxon,$similarity)=split /\t/,$_;
		if (($similarity >= $similarity_value) and ($gene=~ /(.*)-(\d)_coding_aa_(.*)/)) {
			$gene_remain2{$1}++;
		}elsif (($similarity < $similarity_value) and ($gene=~ /(.*)-(\d)_coding_aa_(.*)/) and (($1=~ /ycf1/) or ($1=~ /ycf2/))) {
			$gene_remain2{$1}++;
		}elsif (($similarity < $similarity_value) and ($gene=~ /(.*)-(\d)_coding_aa_(.*)/) and (($1!~ /ycf1/) and ($1!~ /ycf2/))) {
			$gene_remove2{$1}++;
		}
	}
	foreach my $key (keys %gene_remove2) {
		if (!exists $gene_remain2{$key}) {
			$gene_loss2{$key}=$gene_remove2{$key};
		}
	}

	seek($in_bt_ref4,0,0);
	while(<$in_bt_ref4>){# -1/2/3_coding_aa
		$_=~ s/\r|\n//g;
		my ($gene,$taxon,$similarity)=split /\t/,$_;
		my $name;
		if ($gene=~ /(.*)-(\d)_coding_aa_(.*)/) {
			$name=$1;
		}
		if ($similarity >= $similarity_value) {
			#print $out_bt_ref "$_\n";
			push @reference_temp, "$_\n";
		}elsif (($similarity < $similarity_value) and (exists $gene_remain2{$name})) {
			#print $out_bt_ref "$_\n";
			push @reference_temp, "$_\n";
		}
	}

	foreach my $key (keys %gene_loss1) {
		$gene_loss{$key}=$gene_loss1{$key};
	}
	foreach my $key (keys %gene_loss2) {
		$gene_loss{$key}=$gene_loss2{$key};
	}

	while(<$in_bt_ref2>){# _gene
		$_=~ s/\r|\n//g;
		my ($gene,$taxon,$similarity)=split /\t/,$_;
		my $name;
		if ($gene=~ /(.*)_gene/) {
			$name=$1;
		}
		if (!exists $gene_loss{$name}) {
			#print $out_bt_ref "$_\n";
			push @reference_temp, "$_\n";
		}
	}
	close $in_bt_ref1;
	close $in_bt_ref2;
	close $in_bt_ref3;
	close $in_bt_ref4;
	#close $out_bt_ref;


	#short_exon_of_protein-coding-gene_or_tRNA_or_rRNA
	#open (my $in_temp,"<","reference_temp");
	my @reference_temp2=@reference_temp;
	my %hash_temp;
	#while (<$in_temp>) {
	while (@reference_temp) {
		$_=shift @reference_temp;
		$_=~ s/\r|\n//g;
		my ($qseqid,$sseqid,$pident,$length,$mismatch,$gapopen,$qstart,$qend,$sstart,$send,$evalue,$bitscore)=split (/\s+/,$_);
		if ($qseqid=~ /((.+)-1_coding)_(?!aa)(.+)/){# trnK-UUU
			$hash_temp{$1}{$3}=$pident;
		}elsif($qseqid=~ /((.+)-2_coding)_(?!aa)(.+)/){# trnK-UUU
			$hash_temp{$1}{$3}=$pident;
		}elsif($qseqid=~ /((.+)-3_coding)_(?!aa)(.+)/){#
			$hash_temp{$1}{$3}=$pident;
		}elsif ($qseqid=~ /((.+)-1_coding)_aa_(.+)/){# ycf3,atpF
			$hash_temp{$1}{$3}=$pident;
		}elsif($qseqid=~ /((.+)-2_coding)_aa_(.+)/){# ycf3,atpF
			$hash_temp{$1}{$3}=$pident;
		}elsif($qseqid=~ /((.+)-3_coding)_aa_(.+)/){# ycf3
			$hash_temp{$1}{$3}=$pident;
		}elsif($qseqid=~ /((trn)(.+)(?!-\d)_gene)_(.+)/){# trnQ-UUG
			$hash_temp{$1}{$4}=$pident;
		}elsif($qseqid=~ /((rrn)(.+)_gene)_(.+)/) {# rrn16
			$hash_temp{$1}{$4}=$pident;
		}elsif ($qseqid=~ /((.+)_CDS_aa)_(.+)/) {# psbA
			$hash_temp{$1}{$3}=$pident;
		}
	}
	#close $in_temp;
	#print Dumper \%hash_temp;

	my @ref_gene_species;
	foreach my $key1 (keys %hash_temp){
		my @ref_species;
		foreach my $key2 (sort {$hash_temp{$key1}{$b} <=> $hash_temp{$key1}{$a}} keys %{$hash_temp{$key1}}){
			push @ref_species,$key2;
		}
		$key1=~ s/_gene/_coding/g;
		$key1=~ s/_CDS_aa/_coding/g;
		push @ref_gene_species,$key1."_".$ref_species[0];
	}
	%hash_temp=();

	#open (my $input_ref,"<","all2.fasta");
	my @all222;
	@all222=@all22;
	my ($h_ref,$s_ref,%hash_exon,%hash_exon_length,%hash_exon_sequence);
	#while (defined ($h_ref=<$input_ref>) and defined ($s_ref=<$input_ref>)){
	while (@all222){
		$h_ref=shift @all222;
		$s_ref=shift @all222;

		$h_ref=~ s/\r|\n//g;
		$s_ref=~ s/\r|\n//g;

		$h_ref=~ s/^>//g;
		my $len=length $s_ref;
		if ($h_ref=~ /petB-1/) {
			push @ref_gene_species,$h_ref;
		}
		if ($h_ref=~ /petD-1/) {
			push @ref_gene_species,$h_ref;
		}
		if ($h_ref=~ /rpl16-1/) {
			push @ref_gene_species,$h_ref;
		}
		if ($h_ref=~ /rps16-1/) {
			push @ref_gene_species,$h_ref;
		}
		if (grep {$h_ref eq $_} @ref_gene_species) {
			if ($h_ref=~ /((.+)-1_coding)/){
				$hash_exon{$1}{$len}++;
				$hash_exon_length{$1}=$len;
				$hash_exon_sequence{$1}=$s_ref;
			}elsif($h_ref=~ /((.+)-2_coding)/){
				$hash_exon{$1}{$len}++;
				$hash_exon_length{$1}=$len;
				$hash_exon_sequence{$1}=$s_ref;
			}elsif($h_ref=~ /((.+)-3_coding)/){
				$hash_exon{$1}{$len}++;
				$hash_exon_length{$1}=$len;
				$hash_exon_sequence{$1}=$s_ref;
			}elsif($h_ref=~ /((trn)(.+)(?!-\d)_coding)/){
				$hash_exon{$1}{$len}++;
				$hash_exon_length{$1}=$len;
				$hash_exon_sequence{$1}=$s_ref;
			}elsif($h_ref=~ /((rrn)(.+)_coding)/) {
				$hash_exon{$1}{$len}++;
				$hash_exon_length{$1}=$len;
				$hash_exon_sequence{$1}=$s_ref;
			}elsif($h_ref=~ /((.+)_coding)/) {
				$hash_exon{$1}{$len}++;
				$hash_exon_length{$1}=$len;
				$hash_exon_sequence{$1}=$s_ref;
			}
		}
	}
	#close $input_ref;
	#print Dumper \%hash_exon;
	#print Dumper \%hash_exon_length;
	#print Dumper \%hash_exon_sequence;

	#open (my $in_ref,"<","reference_temp");
	#open (my $out_ref,">","reference.tab");
	my @reference_tab;
	my @array_reference;
	#while (<$in_ref>){
	while (@reference_temp2){
		$_=shift @reference_temp2;
		$_=~ s/\r|\n//g;
		my ($qseqid,$sseqid,$pident,$length,$mismatch,$gapopen,$qstart,$qend,$sstart,$send,$evalue,$bitscore)=split (/\s+/,$_);
		push @array_reference,[$qseqid,$sseqid,$pident,$length,$mismatch,$gapopen,$qstart,$qend,$sstart,$send,$evalue,$bitscore];
	}

	foreach my $item (sort {$a->[0] cmp $b->[0] or $a->[1] cmp $b->[1] or $b->[2] <=> $a->[2] or $b->[3] <=> $a->[3]} @array_reference){
		#print $out_ref "$item->[0]\t$item->[1]\t$item->[2]\t$item->[3]\t$item->[4]\t$item->[5]\t$item->[6]\t$item->[7]\t$item->[8]\t$item->[9]\t$item->[10]\t$item->[11]\n";
		push @reference_tab, "$item->[0]\t$item->[1]\t$item->[2]\t$item->[3]\t$item->[4]\t$item->[5]\t$item->[6]\t$item->[7]\t$item->[8]\t$item->[9]\t$item->[10]\t$item->[11]\n";
	}
	#close $in_ref;
	#close $out_ref;
	unlink ("$blast_reference1_random");
	unlink ("$blast_reference2_random");
	unlink ("$blast_reference3_random");
	unlink ("$blast_reference4_random");
	#unlink ("reference_temp");


	############################################################
	## generate_annotation_table
	############################################################
	#open (my $input_reference,"<","reference.tab");
	#open (my $output_annotation_tab,">","$filename.tab");
	my $filename_tab_random = $random."_".$filename.".tab";
	open (my $output_annotation_tab,">",$filename_tab_random);
	my ($row_reference,%hash_tab,%hash1_tab,%hash2_tab,%hash3_tab,%hash4_tab,%hash_gene);
	my $i=0;
	#while (defined ($row_reference=<$input_reference>)){
	while (@reference_tab){
		$row_reference=shift @reference_tab;
		my @array=split /\t/,$row_reference;
		my $contig=$array[1];
		%hash1_tab=();
		if ($array[0]=~ /((.+)_CDS_aa)/){
			$hash_gene{$1}=1;
			$hash1_tab{$2}++;
			$hash_tab{$contig}{$2}{$1}{$array[8]."\t".$array[9]}++ if ((keys %hash1_tab)==1);
		}
		%hash2_tab=();
		if ($array[0]=~ /((.+)_gene)/){
			$hash_gene{$1}=1;
			$hash2_tab{$2}++;
			$hash_tab{$contig}{$2}{$1}{$array[8]."\t".$array[9]}++ if ((keys %hash2_tab)==1);
		}
		%hash3_tab=();
		if ($array[0]=~ /((.+)-(\d)_coding_aa)/){
			$hash_gene{$1}=1;
			$hash3_tab{$2}++;
			$hash_tab{$contig}{$2}{$1}{$array[8]."\t".$array[9]}++ if ((keys %hash3_tab)==1);
		}
		%hash4_tab=();
		if ($array[0]=~ /((.+)-(\d)_coding)(?!_aa)/){
			$hash_gene{$1}=1;
			$hash4_tab{$2}++;
			$hash_tab{$contig}{$2}{$1}{$array[8]."\t".$array[9]}++ if ((keys %hash4_tab)==1);
		}
	}
	#close $input_reference;
	#print Dumper \%hash_tab;

	my @gene_coding_CDS;
	foreach (keys %hash_gene){
		push @gene_coding_CDS,$_;
	}
	my $tick1=0;
	my $tick2=0;
	my $tick3=0;
	if (grep {"clpP-3_coding_aa" eq $_} @gene_coding_CDS) {
		$tick1=1;
	}
	if (grep {"ycf3-3_coding_aa" eq $_} @gene_coding_CDS) {
		$tick2=1;
	}
	if (grep {"rps12+2-2_coding_aa" eq $_} @gene_coding_CDS) {
		$tick3=1;
	}
	foreach my $contig_name (sort keys %hash_tab){
		print $output_annotation_tab "$contig_name\n";
		foreach my $name (sort keys %{$hash_tab{$contig_name}}){
			print $output_annotation_tab "$name\t";
			my $cnt1=0;
			foreach my $gene_coding_CDS (sort {$b cmp $a} keys %{$hash_tab{$contig_name}{$name}}){
				print $output_annotation_tab "$gene_coding_CDS\t" if ($cnt1==0);
				print $output_annotation_tab "\t\t$gene_coding_CDS\t" if ($cnt1>0);
				my @keys;
				my (@array1,@array2,@array3);
				if (($tick1==1) and ($gene_coding_CDS=~ /clpP/)) {
					my ($start,$end);
					foreach my $key (keys %{$hash_tab{$contig_name}{$name}{$gene_coding_CDS}}){
						($start,$end)=split (/\s+/,$key);
						push @array1,[$start,$end];
					}
					if ($start > $end) {
						foreach my $line (sort {$b->[0] <=> $a->[0] or $b->[1] <=> $a->[1]} @array1){
							push @keys,"$line->[0]\t$line->[1]";
						}
					}elsif ($start < $end) {
						foreach my $line (sort {$a->[0] <=> $b->[0] or $a->[1] <=> $b->[1]} @array1){
							push @keys,"$line->[0]\t$line->[1]";
						}
					}
				}
				if (($tick2==1) and ($gene_coding_CDS=~ /ycf3/)) {
					my ($start,$end);
					foreach my $key (keys %{$hash_tab{$contig_name}{$name}{$gene_coding_CDS}}){
						($start,$end)=split (/\s+/,$key);
						push @array2,[$start,$end];
					}
					if ($start > $end) {
						foreach my $line (sort {$b->[0] <=> $a->[0] or $b->[1] <=> $a->[1]} @array2){
							push @keys,"$line->[0]\t$line->[1]";
						}
					}elsif ($start < $end) {
						foreach my $line (sort {$a->[0] <=> $b->[0] or $a->[1] <=> $b->[1]} @array2){
							push @keys,"$line->[0]\t$line->[1]";
						}
					}
				}
				if (($tick3==1) and ($gene_coding_CDS=~ /rps12\+2/)) {
					my ($start,$end);
					foreach my $key (keys %{$hash_tab{$contig_name}{$name}{$gene_coding_CDS}}){
						($start,$end)=split (/\s+/,$key);
						push @array3,[$start,$end];
					}
					if ($start > $end) {
						foreach my $line (sort {$b->[0] <=> $a->[0] or $b->[1] <=> $a->[1]} @array3){
							push @keys,"$line->[0]\t$line->[1]";
						}
					}elsif ($start < $end) {
						foreach my $line (sort {$a->[0] <=> $b->[0] or $a->[1] <=> $b->[1]} @array3){
							push @keys,"$line->[0]\t$line->[1]";
						}
					}
				}
				if (($tick1==0) and ($gene_coding_CDS=~ /clpP/)) {
					foreach my $key (sort {$hash_tab{$contig_name}{$name}{$gene_coding_CDS}{$b} <=> $hash_tab{$contig_name}{$name}{$gene_coding_CDS}{$a}} keys %{$hash_tab{$contig_name}{$name}{$gene_coding_CDS}}){
						push @keys,$key;
					}
				}
				if (($tick2==0) and ($gene_coding_CDS=~ /ycf3/)) {
					foreach my $key (sort {$hash_tab{$contig_name}{$name}{$gene_coding_CDS}{$b} <=> $hash_tab{$contig_name}{$name}{$gene_coding_CDS}{$a}} keys %{$hash_tab{$contig_name}{$name}{$gene_coding_CDS}}){
						push @keys,$key;
					}
				}
				if (($tick3==0) and ($gene_coding_CDS=~ /rps12\+2/)) {
					foreach my $key (sort {$hash_tab{$contig_name}{$name}{$gene_coding_CDS}{$b} <=> $hash_tab{$contig_name}{$name}{$gene_coding_CDS}{$a}} keys %{$hash_tab{$contig_name}{$name}{$gene_coding_CDS}}){
						push @keys,$key;
					}
				}
				if (($gene_coding_CDS!~ /clpP/) and ($gene_coding_CDS!~ /ycf3/) and ($gene_coding_CDS!~ /rps12\+2/)) {
					foreach my $key (sort {$hash_tab{$contig_name}{$name}{$gene_coding_CDS}{$b} <=> $hash_tab{$contig_name}{$name}{$gene_coding_CDS}{$a}} keys %{$hash_tab{$contig_name}{$name}{$gene_coding_CDS}}){
						push @keys,$key;
					}
				}
				my %large_hash;
				if ($gene_coding_CDS=~ /((.+)_CDS_aa)/){
					my $gene1=$2."-1_coding";
					my $gene2=$2."-2_coding";
					my $gene3=$2."-3_coding";
					if ((grep {$gene1 ne $_} @gene_coding_CDS) and (grep {$gene2 ne $_} @gene_coding_CDS)){
						$large_hash{$keys[0]}=$hash_tab{$contig_name}{$name}{$gene_coding_CDS}{$keys[0]};
					}
					if(((grep {$gene1 eq $_} @gene_coding_CDS) or (grep {$gene2 eq $_} @gene_coding_CDS)) and (grep {$gene3 ne $_} @gene_coding_CDS)){
						$large_hash{$keys[0]}=$hash_tab{$contig_name}{$name}{$gene_coding_CDS}{$keys[0]};
						if (defined $keys[1]){
							$large_hash{$keys[1]}=$hash_tab{$contig_name}{$name}{$gene_coding_CDS}{$keys[1]};
						}
					}
					if(grep {$gene3 eq $_} @gene_coding_CDS){
						$large_hash{$keys[0]}=$hash_tab{$contig_name}{$name}{$gene_coding_CDS}{$keys[0]};
						if (defined $keys[1]){
							$large_hash{$keys[1]}=$hash_tab{$contig_name}{$name}{$gene_coding_CDS}{$keys[1]};
						}
						if (defined $keys[2]){
							$large_hash{$keys[2]}=$hash_tab{$contig_name}{$name}{$gene_coding_CDS}{$keys[2]};
						}
					}
				}else{
					$large_hash{$keys[0]}=$hash_tab{$contig_name}{$name}{$gene_coding_CDS}{$keys[0]};
				}
				#print Dumper \%large_hash;
				foreach my $key6 (sort keys %large_hash){
						print $output_annotation_tab "$key6\t$large_hash{$key6}\t";
				}
				$cnt1++;
				print $output_annotation_tab "\n";
				%large_hash=();
			}
		}
	}
	close $output_annotation_tab;
	#unlink("reference.tab");
	%hash_tab=();
	%hash1_tab=();
	%hash2_tab=();
	%hash3_tab=();
	%hash4_tab=();
	%hash_gene=();


	############################################################
	## generate_gb_file_from_annotation_table
	############################################################
	#open (my $in_annotation,"<","$filename.tab");
	open (my $in_annotation,"<",$filename_tab_random);
	my $contig=<$in_annotation>;
	$contig=~ s/\r|\n//g;
	my %hash;
	while (<$in_annotation>){
		my @array=split /\t/,$_;
		if ($_!~ /^\t/){
			my $gene=$array[0];
			if ($array[1]=~ /((.+)_gene)/){
				$hash{$gene}{$1}{$array[2]."\t".$array[3]}=$array[4];
			}
		}elsif($_=~ /^\t/){
			if ($array[2]=~ /((.+)_CDS_aa)/){
				$hash{$2}{$1}{$array[3]."\t".$array[4]}=$array[5];
				#$hash{$2}{$1}{$array[6]."\t".$array[7]}=$array[8] if (defined $array[8]);
			}
			if ($array[2]=~ /((.+)-(\d)_coding_aa)/){
				$hash{$2}{$1}{$array[3]."\t".$array[4]}=$array[5];
			}
			if ($array[2]=~ /((.+)-(\d)_coding(?!_aa))/){
				$hash{$2}{$1}{$array[3]."\t".$array[4]}=$array[5];
			}
		}
	}
	close $in_annotation;
	#unlink("$filename.tab");
	unlink("$filename_tab_random");
	#print Dumper \%hash;

	my $temp=$filename."_temp";
	open (my $out_annotation,">","$output_directory/$temp.gb");
	open (my $logfile,">>","$output_directory/$warning.log");
	my $time=&getdate;
	print $logfile "$filename\n";
	foreach my $key (keys %gene_loss) {
		print $logfile "Warning: $key has not been annotated due to low similarity with reference!\n";
	}
	print $out_annotation "LOCUS       $filename  $length_cp bp    DNA     $type PLN $time"."\n";
	print $out_annotation "FEATURES             Location/Qualifiers"."\n";
	print $out_annotation "     source          "."1..$length_cp"."\n";
	print $out_annotation "                     /organism=\"$filename\""."\n";
	print $out_annotation "                     /organelle=\"plastid:chloroplast\""."\n";
	print $out_annotation "                     /mol_type=\"genomic DNA\""."\n";
	print $out_annotation "                     /note=\"Annotation Method :: PGA-Plastid Genome Annotator\""."\n";

	if (defined $JLB) {
		print $out_annotation "     repeat_region   "."$JLB..$JSB"."\n";
		print $out_annotation "                     /note=\"inverted repeat B\""."\n";
		print $out_annotation "                     /rpt_type=\"inverted\""."\n";
		print $out_annotation "     repeat_region   "."complement($JSA..$JLA)"."\n";
		print $out_annotation "                     /note=\"inverted repeat A\""."\n";
		print $out_annotation "                     /rpt_type=\"inverted\""."\n";
	}

	my %gene_number_seq;
	foreach my $name (sort keys %hash){
		############################################################
		## no-intron RNA
		############################################################
		if ((keys %{$hash{$name}} == 1)){
			foreach my $gene_coding_CDS (keys %{$hash{$name}}){
				my ($start,$end);
				foreach my $position (keys %{$hash{$name}{$gene_coding_CDS}}){
					($start,$end)=(split /\t/,$position)[0,1];
				}

				my $length_ref_rna=$hash_exon_length{"$name\_coding"};
				my $sequence_ref_rna=$hash_exon_sequence{"$name\_coding"};
				$sequence_ref_rna=uc $sequence_ref_rna;
				my $x=31;
				my $y=9;#tRNA
				my $z=9;#rRNA

				if ((($gene_coding_CDS=~ /(((.+)-(\D+))_gene)/) or ($gene_coding_CDS=~ /trnR_gene/) or ($gene_coding_CDS=~ /trnA_gene/)) and ((defined $length_ref_rna) and (defined $sequence_ref_rna))){# tRNA
					if ($start < $end) {
						my $seq_string1;
						if ($start >= 31) {
							$seq_string1=substr($sequence,($start-31),60);
						}elsif ($start < 31) {
							$seq_string1=substr($sequence,($length_cp-(30-$start)-1),(30-$start)).substr($sequence,0,$start).substr($sequence,($start-1),30);
						}
						my $seq_string2;
						if ($end <= $length_cp-30) {
							$seq_string2=substr($sequence,($end-30),60);
						}elsif ($end > $length_cp-30) {
							$seq_string2=substr($sequence,($end-30),($length_cp-$end+30)).substr($sequence,0,(30-($length_cp-$end)));
						}

						my ($start_new,$end_new);
						for (my $i=0;$i<$x;$i+=1) {
							my $match_ref_rna=substr ($sequence_ref_rna,$i,$y);
							my $marks=0;
							for (my $j=$y;$j<60;$j+=1) {
								my $repeat=substr ($seq_string1,-$j,$y);
								if (($repeat eq $match_ref_rna) and ($marks == 0)) {
									$start_new=$start+30-$j-$i if ($start >= 31);
									if ($start < 31) {
										my $start_temp=$start+30-$j-$i;
										if ($start_temp > 0) {
											$start_new=$start_temp;
										}elsif ($start_temp <= 0) {
											$start_new=$length_cp+$start_temp;
										}
									}
									$marks++;
								}
							}
							last if (defined $start_new);
						}
						for (my $i=$y;$i<$x;$i+=1) {
							my $match_ref_rna=substr ($sequence_ref_rna,-$i,$y);
							my $marks=0;
							for (my $j=0;$j<60;$j+=1) {
								my $repeat=substr ($seq_string2,$j,$y);
								if (($repeat eq $match_ref_rna) and ($marks == 0)) {
									$end_new=$end-30+($j+$y)+($i-$y) if ($end <= ($length_cp-30));
									if ($end > ($length_cp-30)) {
										my $end_temp=$end-30+($j+$y)+($i-$y);
										if ($end_temp <= $length_cp) {
											$end_new=$end_temp;
										}elsif ($end_temp > $length_cp) {
											$end_new=$end_temp-$length_cp;
										}
									}
									$marks++;
								}
							}
							last if (defined $end_new);
						}

						if ((defined $start_new) and (defined $end_new)) {
							my $length_trn;
							if ($end_new > $start_new) {
								$length_trn=$end_new-$start_new;
							}elsif ($end_new < $start_new) {
								$length_trn=$end_new+($length_cp-$start_new);
							}

							if ($length_trn >= 70) {
								print $out_annotation "     "."gene"."            ".$start_new."..".$end_new."\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "     "."tRNA"."            ".$start_new."..".$end_new."\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
								$gene_number_seq{$name}++;

								my ($S,$E);
								my $string;
								if ($end_new > $start_new) {
									$string=substr($sequence,($start_new-1),($end_new-$start_new+1));
								}elsif ($end_new < $start_new) {
									$string=substr($sequence,($start_new-1),($length_cp-$start_new+1)).substr($sequence,0,$end_new);
								}
								my $lengths=length $string;
								for (my $i=0;$i<($length_cp-$lengths);$i+=1){
									my $repeat=substr ($sequence,$i,$lengths);
									if (($repeat eq $string) and (($i+1) ne $start_new)) {
										$S=$i+1;
										$E=$i+$lengths-1+1;
										print $out_annotation "     "."gene"."            ".$S."..".$E."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."tRNA"."            ".$S."..".$E."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										$gene_number_seq{$name}++;
									}
								}
								for (my $i=0;$i<($length_cp-$lengths);$i+=1){
									my $repeat=substr ($rev_coms,$i,$lengths);
									if ($repeat eq $string) {
										$S=$length_cp-$i-1+1;
										$E=$length_cp-$i-$lengths+1;
										print $out_annotation "     "."gene"."            "."complement(".$E."..".$S.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."tRNA"."            "."complement(".$E."..".$S.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										$gene_number_seq{$name}++;
									}
								}
							}elsif ($length_trn < 70) {
								print $logfile "Warning: $name (positive no-intron tRNA) has not been annotated due to short length!\n";
							}
						}elsif ((!defined $start_new) or (!defined $end_new)) {
							my $length_trn=$end-$start;
							if ($length_trn >= 70) {
								print $out_annotation "     "."gene"."            ".$start."..".$end."\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "     "."tRNA"."            ".$start."..".$end."\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
								$gene_number_seq{$name}++;
								print $logfile "Warning: $name (positive no-intron tRNA) must be checked due to non-identical boundary with reference!\n";
							}elsif ($length_trn < 70) {
								print $logfile "Warning: $name (positive no-intron tRNA) has not been annotated due to short length!\n";
							}
						}
					}
					if ($start > $end) {
						my $seq_string1;
						if ($start <= ($length_cp-30)) {
							$seq_string1=substr($sequence,($start-30),60);
						}elsif ($start > ($length_cp-30)) {
							$seq_string1=substr($sequence,($start-30),($length_cp-$start+30)).substr($sequence,0,(30-($length_cp-$start)));
						}
						my $seq_string2;
						if ($end >= 31) {
							$seq_string2=substr($sequence,($end-31),60);
						}elsif ($end < 31) {
							$seq_string2=substr($sequence,($length_cp-(30-$end)-1),(30-$end)).substr($sequence,0,$end).substr($sequence,($end-1),30);
						}
						my $rev_seq_string1=reverse $seq_string1;
						$rev_seq_string1=~ tr/ACGTacgt/TGCAtgca/;
						my $rev_seq_string2=reverse $seq_string2;
						$rev_seq_string2=~ tr/ACGTacgt/TGCAtgca/;

						my ($start_new,$end_new);
						for (my $i=0;$i<$x;$i+=1) {
							my $match_ref_rna=substr ($sequence_ref_rna,$i,$y);
							my $marks=0;
							for (my $j=$y;$j<60;$j+=1) {
								my $repeat=substr ($rev_seq_string1,-$j,$y);
								if (($repeat eq $match_ref_rna) and ($marks == 0)) {
									$start_new=$start-30+$j+$i if ($start <= ($length_cp-30));
									if ($start > ($length_cp-30)) {
										my $start_temp=$start-30+$j+$i;
										if ($start_temp <= $length_cp) {
											$start_new=$start_temp;
										}elsif ($start_temp > $length_cp) {
											$start_new=$start_temp-$length_cp;
										}
									}
									$marks++;
								}
							}
							last if (defined $start_new);
						}
						for (my $i=$y;$i<$x;$i+=1) {
							my $match_ref_rna=substr ($sequence_ref_rna,-$i,$y);
							my $marks=0;
							for (my $j=0;$j<60;$j+=1) {
								my $repeat=substr ($rev_seq_string2,$j,$y);
								if (($repeat eq $match_ref_rna) and ($marks == 0)) {
									$end_new=$end+30-($j+$y)-($i-$y) if ($end >= 31);
									if ($end < 31) {
										my $end_temp=$end+30-($j+$y)-($i-$y);
										if ($end_temp > 0) {
											$end_new=$end_temp;
										}elsif ($end_temp <= 0) {
											$end_new=$length_cp+$end_temp;
										}
									}
									$marks++;
								}
							}
							last if (defined $end_new);
						}

						if ((defined $start_new) and (defined $end_new)) {
							my $length_trn;
							if ($end_new < $start_new) {
								$length_trn=$start_new-$end_new;
							}elsif ($end_new > $start_new) {
								$length_trn=$start_new+($length_cp-$end_new);
							}

							if ($length_trn >= 70) {
								print $out_annotation "     "."gene"."            "."complement(".$end_new."..".$start_new.")\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "     "."tRNA"."            "."complement(".$end_new."..".$start_new.")\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
								$gene_number_seq{$name}++;

								my ($S,$E);
								my $string;
								if ($end_new < $start_new) {
									$string=substr($sequence,($end_new-1),($start_new-$end_new+1));
								}elsif ($end_new > $start_new) {
									$string=substr($sequence,($end_new-1),($length_cp-$end_new+1)).substr($sequence,0,$start_new);
								}
								my $lengths=length $string;
								for (my $i=0;$i<($length_cp-$lengths);$i+=1){
									my $repeat=substr ($sequence,$i,$lengths);
									if (($repeat eq $string) and (($i+1) ne $end_new)) {
										$E=$i+1;
										$S=$i+$lengths-1+1;
										print $out_annotation "     "."gene"."            "."complement(".$E."..".$S.")\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."tRNA"."            "."complement(".$E."..".$S.")\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										$gene_number_seq{$name}++;
									}
								}
								for (my $i=0;$i<($length_cp-$lengths);$i+=1){
									my $repeat=substr ($rev_coms,$i,$lengths);
									if ($repeat eq $string) {
										$E=$length_cp-$i-1+1;
										$S=$length_cp-$i-$lengths+1;
										print $out_annotation "     "."gene"."            ".$S."..".$E."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."tRNA"."            ".$S."..".$E."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										$gene_number_seq{$name}++;
									}
								}
							}elsif ($length_trn < 70) {
								print $logfile "Warning: $name (negative no-intron tRNA) has not been annotated due to short length!\n";
							}
						}elsif ((!defined $start_new) or (!defined $end_new)) {
							my $length_trn=$start-$end;
							if ($length_trn >= 70) {
								print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "     "."tRNA"."            "."complement(".$end."..".$start.")\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
								$gene_number_seq{$name}++;
								print $logfile "Warning: $name (negative no-intron tRNA) must be checked due to non-identical boundary with reference!\n";
							}elsif ($length_trn < 70) {
								print $logfile "Warning: $name (negative no-intron tRNA) has not been annotated due to short length!\n";
							}
						}
					}
				}


				if(($gene_coding_CDS=~ /((rrn(\d+\.?\d*))_gene)/) and ((defined $length_ref_rna) and (defined $sequence_ref_rna))){# rRNA
					if ($start < $end) {
						my $seq_string1;
						if ($start >= 801) {
							$seq_string1=substr($sequence,($start-801),900);
						}elsif ($start < 801) {
							$seq_string1=substr($sequence,0,900);
						}
						my $seq_string2=substr($sequence,($end-100),900);

						my ($start_new,$end_new);
						for (my $i=0;$i<$x;$i+=1) {
							my $match_ref_rna=substr ($sequence_ref_rna,$i,$z);
							my $marks=0;
							for (my $j=$z;$j<900;$j+=1) {
								my $repeat=substr ($seq_string1,-$j,$z);
								if (($repeat eq $match_ref_rna) and ($marks == 0)) {
									$start_new=$start+100-$j-$i;
									$marks++;
								}
							}
							last if (defined $start_new);
						}
						for (my $i=$z;$i<$x;$i+=1) {
							my $match_ref_rna=substr ($sequence_ref_rna,-$i,$z);
							my $marks=0;
							for (my $j=0;$j<900;$j+=1) {
								my $repeat=substr ($seq_string2,$j,$z);
								if (($repeat eq $match_ref_rna) and ($marks == 0)) {
									$end_new=$end-100+($j+$z)+($i-$z);
									$marks++;
								}
							}
							last if (defined $end_new);
						}

						if ((defined $start_new) and (defined $end_new)) {
							my $length_rrna=$end_new-$start_new;
							if ($redundancy eq "N") {
								if ((defined $length_ref_rna) and (($length_rrna/$length_ref_rna) < $short)) {
									print $logfile "Warning: $name (positive rRNA) has not been annotated due to query coverage per annotated rRNA less than $short!\n";
								}elsif ((defined $length_ref_rna) and (($length_rrna/$length_ref_rna) > $long)) {
									print $logfile "Warning: $name (positive rRNA) has not been annotated due to query coverage per annotated rRNA more than $long!\n";
								}elsif ((defined $length_ref_rna) and (($length_rrna/$length_ref_rna) > $short) and (($length_rrna/$length_ref_rna) < $long)) {
									print $out_annotation "     "."gene"."            ".$start_new."..".$end_new."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."rRNA"."            ".$start_new."..".$end_new."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									$gene_number_seq{$name}++;
									my ($S,$E);
									my $string=substr($sequence,($start_new-1),($end_new-$start_new+1));
									my $lengths=length $string;
									for (my $i=0;$i<($length_cp-$lengths);$i+=1){
										my $repeat=substr ($sequence,$i,$lengths);
										if (($repeat eq $string) and (($i+1) ne $start_new)) {
											$S=$i+1;
											$E=$i+$lengths-1+1;
											print $out_annotation "     "."gene"."            ".$S."..".$E."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."rRNA"."            ".$S."..".$E."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											$gene_number_seq{$name}++;
										}
									}
									for (my $i=0;$i<($length_cp-$lengths);$i+=1){
										my $repeat=substr ($rev_coms,$i,$lengths);
										if ($repeat eq $string) {
											$S=$length_cp-$i-1+1;
											$E=$length_cp-$i-$lengths+1;
											print $out_annotation "     "."gene"."            "."complement(".$E."..".$S.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."rRNA"."            "."complement(".$E."..".$S.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											$gene_number_seq{$name}++;
										}
									}
								}
							}elsif ($redundancy eq "Y") {
								print $out_annotation "     "."gene"."            ".$start_new."..".$end_new."\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "     "."rRNA"."            ".$start_new."..".$end_new."\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
								$gene_number_seq{$name}++;
								my ($S,$E);
								my $string=substr($sequence,($start_new-1),($end_new-$start_new+1));
								my $lengths=length $string;
								for (my $i=0;$i<($length_cp-$lengths);$i+=1){
									my $repeat=substr ($sequence,$i,$lengths);
									if (($repeat eq $string) and (($i+1) ne $start_new)) {
										$S=$i+1;
										$E=$i+$lengths-1+1;
										print $out_annotation "     "."gene"."            ".$S."..".$E."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."rRNA"."            ".$S."..".$E."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										$gene_number_seq{$name}++;
									}
								}
								for (my $i=0;$i<($length_cp-$lengths);$i+=1){
									my $repeat=substr ($rev_coms,$i,$lengths);
									if ($repeat eq $string) {
										$S=$length_cp-$i-1+1;
										$E=$length_cp-$i-$lengths+1;
										print $out_annotation "     "."gene"."            "."complement(".$E."..".$S.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."rRNA"."            "."complement(".$E."..".$S.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										$gene_number_seq{$name}++;
									}
								}

								if ((defined $length_ref_rna) and (($length_rrna/$length_ref_rna) < $short)) {
									print $logfile "Warning: $name (positive rRNA) need to be checked due to query coverage per annotated rRNA less than $short!\n";
								}elsif ((defined $length_ref_rna) and (($length_rrna/$length_ref_rna) > $long)) {
									print $logfile "Warning: $name (positive rRNA) need to be checked due to query coverage per annotated rRNA more than $long!\n";
								}
							}
						}elsif ((!defined $start_new) or (!defined $end_new)) {
							my $length_rrna=$end-$start;
							if ($redundancy eq "N") {
								if ((defined $length_ref_rna) and (($length_rrna/$length_ref_rna) < $short)) {
									print $logfile "Warning: $name (positive rRNA) has not been annotated due to query coverage per annotated rRNA less than $short!\n";
								}elsif ((defined $length_ref_rna) and (($length_rrna/$length_ref_rna) > $long)) {
									print $logfile "Warning: $name (positive rRNA) has not been annotated due to query coverage per annotated rRNA more than $long!\n";
								}elsif ((defined $length_ref_rna) and (($length_rrna/$length_ref_rna) > $short) and (($length_rrna/$length_ref_rna) < $long)) {
									print $out_annotation "     "."gene"."            ".$start."..".$end."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."rRNA"."            ".$start."..".$end."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									$gene_number_seq{$name}++;
									print $logfile "Warning: $name (positive rRNA) must be checked due to non-identical boundary with reference!\n";
								}
							}elsif ($redundancy eq "Y") {
								print $out_annotation "     "."gene"."            ".$start."..".$end."\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "     "."rRNA"."            ".$start."..".$end."\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
								$gene_number_seq{$name}++;
								print $logfile "Warning: $name (positive rRNA) must be checked due to non-identical boundary with reference!\n";

								if ((defined $length_ref_rna) and (($length_rrna/$length_ref_rna) < $short)) {
									print $logfile "Warning: $name (positive rRNA) need to be checked due to query coverage per annotated rRNA less than $short!\n";
								}elsif ((defined $length_ref_rna) and (($length_rrna/$length_ref_rna) > $long)) {
									print $logfile "Warning: $name (positive rRNA) need to be checked due to query coverage per annotated rRNA more than $long!\n";
								}
							}
						}
					}

					if ($start > $end) {
						my $seq_string1=substr($sequence,($start-100),900);
						my $seq_string2;
						if ($end >= 801) {
							$seq_string2=substr($sequence,($end-801),900);
						}elsif ($end < 801) {
							$seq_string2=substr($sequence,0,900);
						}
						my $rev_seq_string1=reverse $seq_string1;
						$rev_seq_string1=~ tr/ACGTacgt/TGCAtgca/;
						my $rev_seq_string2=reverse $seq_string2;
						$rev_seq_string2=~ tr/ACGTacgt/TGCAtgca/;

						my ($start_new,$end_new);
						for (my $i=0;$i<$x;$i+=1) {
							my $match_ref_rna=substr ($sequence_ref_rna,$i,$z);
							my $marks=0;
							for (my $j=$z;$j<900;$j+=1) {
								my $repeat=substr ($rev_seq_string1,-$j,$z);
								if (($repeat eq $match_ref_rna) and ($marks == 0)) {
									$start_new=$start-100+$j+$i;
									$marks++;
								}
							}
							last if (defined $start_new);
						}
						for (my $i=$z;$i<$x;$i+=1) {
							my $match_ref_rna=substr ($sequence_ref_rna,-$i,$z);
							my $marks=0;
							for (my $j=0;$j<900;$j+=1) {
								my $repeat=substr ($rev_seq_string2,$j,$z);
								if (($repeat eq $match_ref_rna) and ($marks == 0)) {
									$end_new=$end+100-($j+$z)-($i-$z);
									$marks++;
								}
							}
							last if (defined $end_new);
						}

						if ((defined $start_new) and (defined $end_new)) {
							my $length_rrna=$start_new-$end_new;
							if ($redundancy eq "N") {
								if ((defined $length_ref_rna) and (($length_rrna/$length_ref_rna) < $short)) {
									print $logfile "Warning: $name (negative rRNA) has not been annotated due to query coverage per annotated rRNA less than $short!\n";
								}elsif ((defined $length_ref_rna) and (($length_rrna/$length_ref_rna) > $long)) {
									print $logfile "Warning: $name (negative rRNA) has not been annotated due to query coverage per annotated rRNA more than $long!\n";
								}elsif ((defined $length_ref_rna) and (($length_rrna/$length_ref_rna) > $short) and (($length_rrna/$length_ref_rna) < $long)) {
									print $out_annotation "     "."gene"."            "."complement(".$end_new."..".$start_new.")\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."rRNA"."            "."complement(".$end_new."..".$start_new.")\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									$gene_number_seq{$name}++;
									my ($S,$E);
									my $string=substr($sequence,($end_new-1),($start_new-$end_new+1));
									my $lengths=length $string;
									for (my $i=0;$i<($length_cp-$lengths);$i+=1){
										my $repeat=substr ($sequence,$i,$lengths);
										if (($repeat eq $string) and (($i+1) ne $end_new)) {
											$E=$i+1;
											$S=$i+$lengths-1+1;
											print $out_annotation "     "."gene"."            "."complement(".$E."..".$S.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."rRNA"."            "."complement(".$E."..".$S.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											$gene_number_seq{$name}++;
										}
									}
									for (my $i=0;$i<($length_cp-$lengths);$i+=1){
										my $repeat=substr ($rev_coms,$i,$lengths);
										if ($repeat eq $string) {
											$E=$length_cp-$i-1+1;
											$S=$length_cp-$i-$lengths+1;
											print $out_annotation "     "."gene"."            ".$S."..".$E."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."rRNA"."            ".$S."..".$E."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											$gene_number_seq{$name}++;
										}
									}
								}
							}elsif ($redundancy eq "Y") {
								print $out_annotation "     "."gene"."            "."complement(".$end_new."..".$start_new.")\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "     "."rRNA"."            "."complement(".$end_new."..".$start_new.")\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
								$gene_number_seq{$name}++;
								my ($S,$E);
								my $string=substr($sequence,($end_new-1),($start_new-$end_new+1));
								my $lengths=length $string;
								for (my $i=0;$i<($length_cp-$lengths);$i+=1){
									my $repeat=substr ($sequence,$i,$lengths);
									if (($repeat eq $string) and (($i+1) ne $end_new)) {
										$E=$i+1;
										$S=$i+$lengths-1+1;
										print $out_annotation "     "."gene"."            "."complement(".$E."..".$S.")\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."rRNA"."            "."complement(".$E."..".$S.")\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										$gene_number_seq{$name}++;
									}
								}
								for (my $i=0;$i<($length_cp-$lengths);$i+=1){
									my $repeat=substr ($rev_coms,$i,$lengths);
									if ($repeat eq $string) {
										$E=$length_cp-$i-1+1;
										$S=$length_cp-$i-$lengths+1;
										print $out_annotation "     "."gene"."            ".$S."..".$E."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."rRNA"."            ".$S."..".$E."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										$gene_number_seq{$name}++;
									}
								}

								if ((defined $length_ref_rna) and (($length_rrna/$length_ref_rna) < $short)) {
									print $logfile "Warning: $name (negative rRNA) need to be checked due to query coverage per annotated rRNA less than $short!\n";
								}elsif ((defined $length_ref_rna) and (($length_rrna/$length_ref_rna) > $long)) {
									print $logfile "Warning: $name (negative rRNA) need to be checked due to query coverage per annotated rRNA more than $long!\n";
								}
							}
						}elsif ((!defined $start_new) or (!defined $end_new)) {
							my $length_rrna=$start-$end;
							if ($redundancy eq "N") {
								if ((defined $length_ref_rna) and (($length_rrna/$length_ref_rna) < $short)) {
									print $logfile "Warning: $name (negative rRNA) has not been annotated due to query coverage per annotated rRNA less than $short!\n";
								}elsif ((defined $length_ref_rna) and (($length_rrna/$length_ref_rna) > $long)) {
									print $logfile "Warning: $name (negative rRNA) has not been annotated due to query coverage per annotated rRNA more than $long!\n";
								}elsif ((defined $length_ref_rna) and (($length_rrna/$length_ref_rna) > $short) and (($length_rrna/$length_ref_rna) < $long)) {
									print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."rRNA"."            "."complement(".$end."..".$start.")\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									$gene_number_seq{$name}++;
									print $logfile "Warning: $name (negative rRNA) must be checked due to non-identical boundary with reference!\n";
								}
							}elsif ($redundancy eq "Y") {
									print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."rRNA"."            "."complement(".$end."..".$start.")\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									$gene_number_seq{$name}++;
									print $logfile "Warning: $name (negative rRNA) must be checked due to non-identical boundary with reference!\n";

								if ((defined $length_ref_rna) and (($length_rrna/$length_ref_rna) < $short)) {
									print $logfile "Warning: $name (negative rRNA) need to be checked due to query coverage per annotated rRNA less than $short!\n";
								}elsif ((defined $length_ref_rna) and (($length_rrna/$length_ref_rna) > $long)) {
									print $logfile "Warning: $name (negative rRNA) need to be checked due to query coverage per annotated rRNA more than $long!\n";
								}
							}
						}
					}
				}


				if((($gene_coding_CDS=~ /((.+)_gene)/) and $name!~ /rrn/ and $name!~ /trn/ and $name!~ /rps12/) and ((defined $length_ref_rna) and (defined $sequence_ref_rna))){# psaM(CDS_aa without alignment result),pseudogene (e.g., ycf15,ycf68 etal.)
					if ($start < $end){
						#print $out_annotation "     "."gene"."             ".$start."..".$end."\n";
						#print $out_annotation "                     "."/gene=\"$name\""."\n";
						##print $out_annotation "                     "."/translation=\"$name\""."\n";
						#$gene_number_seq{$name}++;
						print $logfile "Warning: $name (positive no-intron PCG) has not been annotated!\n";
					}elsif ($start > $end) {
						#print $out_annotation "     "."gene"."             "."complement(".$end."..".$start.")\n";
						#print $out_annotation "                     "."/gene=\"$name\""."\n";
						##print $out_annotation "                     "."/translation=\"$name\""."\n";
						#$gene_number_seq{$name}++;
						print $logfile "Warning: $name (negative no-intron PCG) has not been annotated!\n";
					}
				}
			}
		}

		############################################################
		## no-intron protein-coding-gene
		############################################################
		if((keys %{$hash{$name}} == 2)){
			my (@position1,@position2);
			foreach my $gene_coding_CDS (sort {$b cmp $a} keys %{$hash{$name}}){
				if ($gene_coding_CDS=~ /((.+)_gene)/){
					foreach (keys %{$hash{$name}{$1}}){
						push @position1,$_."\t".$hash{$name}{$gene_coding_CDS}{$_};
					}
				}elsif($gene_coding_CDS=~ /((.+)_CDS_aa)/){
					foreach (keys %{$hash{$name}{$1}}){
						push @position2,$_."\t".$hash{$name}{$gene_coding_CDS}{$_};
					}
				}
			}

			my $length_ref_pcg=$hash_exon_length{"$name\_coding"};
			my $sequence_ref_pcg=$hash_exon_sequence{"$name\_coding"};
			my $xx=21;
			my $yy=4;
			my $aa_ref_pcg;
			my @aa_ref_pcg;
			if ((defined $length_ref_pcg) and (defined $sequence_ref_pcg)) {
				for (my $i=0;$i<$length_ref_pcg;$i+=3){
					my $codon=substr ($sequence_ref_pcg,$i,3);
					$codon=uc $codon;
					if (exists $hash_codon{$codon}){
						$aa_ref_pcg.=$hash_codon{$codon};
					}else{
						$aa_ref_pcg.="X";
						my $j=$i+1;
						#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
					}
				}
				@aa_ref_pcg=split //,$aa_ref_pcg;
			}

			while (@position1 and @position2){
				my $position1=shift @position1;
				my $position2=shift @position2;
				my ($start1,$end1,$number1)=(split /\t/,$position1)[0,1,2];# _gene
				my ($start2,$end2,$number2)=(split /\t/,$position2)[0,1,2];# _CDS_aa
				if (($start1 == $start2) and ($end1 == $end2)){# identical PCG boundary for _gene and _CDS_aa
					if ($start1 < $end1){# positive
						my $str=substr($sequence,($start1-1),($end1-$start1+1));
						my $length=length $str;
						my $forward_start_codon=substr($str,0,3);
						my $forward_stop_codon=substr($str,($length-3),3);
						my $aa;
						for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
							my $codon=substr ($str,$i,3);
							$codon=uc $codon;
							if (exists $hash_codon{$codon}){
								$aa.=$hash_codon{$codon};
							}else{
								$aa.="X";
								my $j=$i+1;
								#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
							}
						}
						my @aa=split //,$aa;
						if (($length % 3==0) and (!grep {$_=~ /\*/} @aa) and (($forward_start_codon eq "ATG") or ($forward_start_codon eq "GTG")) and ((($forward_stop_codon eq "TAA") or ($forward_stop_codon eq "TAG") or ($forward_stop_codon eq "TGA")) or ($name eq "rps12+1"))){# standard start and stop codon
							my $length_pcg=$end1-$start1;
							if ($redundancy eq "N") {
								if ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) < $short)) {
									print $logfile "Warning: $name (positive no-intron PCG) has not been annotated due to query coverage per annotated PCG less than $short!\n";
								}elsif ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) > $long)) {
									print $logfile "Warning: $name (positive no-intron PCG) has not been annotated due to query coverage per annotated PCG more than $long!\n";
								}elsif ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) > $short) and (($length_pcg/$length_ref_pcg) < $long)) {
									print $out_annotation "     "."gene"."            ".$start1."..".$end1."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."CDS"."             ".$start1."..".$end1."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/codon_start=1"."\n";
									print $out_annotation "                     "."/transl_table=11"."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									#print $out_annotation "                     "."/translation=\"$aa\""."\n";
									$gene_number_seq{$name}++;
									my ($S1,$E1);
									my $string=substr($sequence,($start1-1),($end1-$start1+1));
									my $lengths=length $string;
									for (my $i=0;$i<($length_cp-$lengths);$i+=1){
										my $repeat=substr ($sequence,$i,$lengths);
										if (($repeat eq $string) and (($i+1) ne $start1)) {
											$S1=$i+1;
											$E1=$i+$lengths-1+1;
											print $out_annotation "     "."gene"."            ".$S1."..".$E1."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             ".$S1."..".$E1."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa\""."\n";
											$gene_number_seq{$name}++;
										}
									}
									for (my $i=0;$i<($length_cp-$lengths);$i+=1){
										my $repeat=substr ($rev_coms,$i,$lengths);
										if ($repeat eq $string) {
											$S1=$length_cp-$i-1+1;
											$E1=$length_cp-$i-$lengths+1;
											print $out_annotation "     "."gene"."            "."complement(".$E1."..".$S1.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."complement(".$E1."..".$S1.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa\""."\n";
											$gene_number_seq{$name}++;
										}
									}
								}
							}elsif ($redundancy eq "Y") {
								print $out_annotation "     "."gene"."            ".$start1."..".$end1."\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "     "."CDS"."             ".$start1."..".$end1."\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "                     "."/codon_start=1"."\n";
								print $out_annotation "                     "."/transl_table=11"."\n";
								print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
								#print $out_annotation "                     "."/translation=\"$aa\""."\n";
								$gene_number_seq{$name}++;
								my ($S1,$E1);
								my $string=substr($sequence,($start1-1),($end1-$start1+1));
								my $lengths=length $string;
								for (my $i=0;$i<($length_cp-$lengths);$i+=1){
									my $repeat=substr ($sequence,$i,$lengths);
									if (($repeat eq $string) and (($i+1) ne $start1)) {
										$S1=$i+1;
										$E1=$i+$lengths-1+1;
										print $out_annotation "     "."gene"."            ".$S1."..".$E1."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             ".$S1."..".$E1."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
										$gene_number_seq{$name}++;
									}
								}
								for (my $i=0;$i<($length_cp-$lengths);$i+=1){
									my $repeat=substr ($rev_coms,$i,$lengths);
									if ($repeat eq $string) {
										$S1=$length_cp-$i-1+1;
										$E1=$length_cp-$i-$lengths+1;
										print $out_annotation "     "."gene"."            "."complement(".$E1."..".$S1.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."complement(".$E1."..".$S1.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
										$gene_number_seq{$name}++;
									}
								}

								if ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) < $short)) {
									print $logfile "Warning: $name (positive no-intron PCG) need to be checked due to query coverage per annotated PCG less than $short!\n";
								}elsif ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) > $long)) {
									print $logfile "Warning: $name (positive no-intron PCG) need to be checked due to query coverage per annotated PCG more than $long!\n";
								}
							}
						}elsif((grep {$_=~ /\*/} @aa) or (($forward_start_codon ne "ATG") and ($forward_start_codon ne "GTG")) or ((($forward_stop_codon ne "TAA") and ($forward_stop_codon ne "TAG") and ($forward_stop_codon ne "TGA")) and ($name ne "rps12+1"))){# non-standard start or stop codon
							my $gene_length=($end1-$start1+1);
							my ($start_left,$start_right,$end_left,$end_right);
							if ($start2>=9000) {
								$start_left=$start2-9000;
							}elsif (($start2<9000) and ($start2>30)){
								$start_left=$start2-(int($start2/3)-1)*3;
							}elsif ($start2<=30){
								if ($start2 % 3 == 1) {
									$start_left=1;
								}elsif ($start2 % 3 == 2) {
									$start_left=2;
								}elsif ($start2 % 3 == 0) {
									$start_left=3;
								}
							}
							if ($gene_length >= 60) {
								$start_right=$start2+60;
							}elsif ($gene_length < 60) {
								$start_right=$start2+int($gene_length/3)*3;
							}
							$end_left=$end2-int($gene_length/3)*3;
							if (($length_cp-$end2)>=9000) {
								$end_right=$end2+9000;
							}elsif (($length_cp-$end2)<9000){
								$end_right=$end2+(int(($length_cp-$end2)/3)-1)*3;
							}

							my $str1=substr($sequence,($start_left-1),($start_right-$start_left));
							my $str2=substr($sequence,$end_left,($end_right-$end_left));
							my $length1=length $str1;
							my $length2=length $str2;
							my $aa1;#left range for start codon
							for (my $i=0;$i<($length1-3);$i+=3){# delete stop codon
								my $codon=substr ($str1,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa1.=$hash_codon{$codon};
								}else{
									$aa1.="X";
									my $j=$i+1;
									#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa1=split //,$aa1;
							my (@star1,@start1);
							foreach (0..$#aa1){
								if (($aa1[$_]=~ /\*/) and ($_ < int($length1/3))){
									push @star1,$_;
								}
								if (($aa1[$_]=~ /M/) and ($_ <= int($length1/3))){
									push @start1,$_;
								}
							}
							my $left_star_position=pop @star1;
							my @start2;
							foreach (@start1){
								push @start2,$_ if ((defined $left_star_position) and ($_ > $left_star_position))
							}
							my $left_start_position=$start2[0];
							#my $left_start_position;
							#if ($gene_length > 1000) {
								#$left_start_position=$start2[0];
							#}elsif ($gene_length < 1000) {
								#$left_start_position=pop @start2;
							#}
							my $aa2;#right range for stop codon
							#for (my $i=0;$i<$length2;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
							for (my $i=0;$i<($length2-3);$i+=3){# delete stop codon
								my $codon=substr ($str2,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa2.=$hash_codon{$codon};
								}else{
									$aa2.="X";
									my $j=$i+1;
									#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa2=split //,$aa2;
							my @star2;
							foreach (0..$#aa2){
								if ($aa2[$_]=~ /\*/){
									push @star2,$_;
								}
							}
							my $right_star_position=shift @star2;

							if ((defined $left_start_position) and (defined $right_star_position)){
								my $start=$start_left+$left_start_position * 3;
								my $end=$end_left+($right_star_position * 3+3);
								my $str=substr($sequence,($start-1),($end-$start+1));
								my $length=length $str;
								my $aa;
								for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
									my $codon=substr ($str,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa.=$hash_codon{$codon};
									}else{
										$aa.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}

								my $length_pcg=$end-$start;
								if ($length_pcg >= 60) {
									if ($redundancy eq "N") {
										if ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) < $short)) {
											print $logfile "Warning: $name (positive no-intron PCG) has not been annotated due to query coverage per annotated PCG less than $short!\n";
										}elsif ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) > $long)) {
											print $logfile "Warning: $name (positive no-intron PCG) has not been annotated due to query coverage per annotated PCG more than $long!\n";
										}elsif ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) > $short) and (($length_pcg/$length_ref_pcg) < $long)) {
											print $out_annotation "     "."gene"."            ".$start."..".$end."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             ".$start."..".$end."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa\""."\n";
											$gene_number_seq{$name}++;
											my ($S,$E);
											my $string=substr($sequence,($start-1),($end-$start+1));
											my $lengths=length $string;
											for (my $i=0;$i<($length_cp-$lengths);$i+=1){
												my $repeat=substr ($sequence,$i,$lengths);
												if (($repeat eq $string) and (($i+1) ne $start)) {
													$S=$i+1;
													$E=$i+$lengths-1+1;
													print $out_annotation "     "."gene"."            ".$S."..".$E."\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "     "."CDS"."             ".$S."..".$E."\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "                     "."/codon_start=1"."\n";
													print $out_annotation "                     "."/transl_table=11"."\n";
													print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
													#print $out_annotation "                     "."/translation=\"$aa\""."\n";
													$gene_number_seq{$name}++;
												}
											}
											for (my $i=0;$i<($length_cp-$lengths);$i+=1){
												my $repeat=substr ($rev_coms,$i,$lengths);
												if ($repeat eq $string) {
													$S=$length_cp-$i-1+1;
													$E=$length_cp-$i-$lengths+1;
													print $out_annotation "     "."gene"."            "."complement(".$E."..".$S.")"."\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "     "."CDS"."             "."complement(".$E."..".$S.")"."\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "                     "."/codon_start=1"."\n";
													print $out_annotation "                     "."/transl_table=11"."\n";
													print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
													#print $out_annotation "                     "."/translation=\"$aa\""."\n";
													$gene_number_seq{$name}++;
												}
											}
										}
									}elsif ($redundancy eq "Y") {
										print $out_annotation "     "."gene"."            ".$start."..".$end."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             ".$start."..".$end."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
										$gene_number_seq{$name}++;
										my ($S,$E);
										my $string=substr($sequence,($start-1),($end-$start+1));
										my $lengths=length $string;
										for (my $i=0;$i<($length_cp-$lengths);$i+=1){
											my $repeat=substr ($sequence,$i,$lengths);
											if (($repeat eq $string) and (($i+1) ne $start)) {
												$S=$i+1;
												$E=$i+$lengths-1+1;
												print $out_annotation "     "."gene"."            ".$S."..".$E."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             ".$S."..".$E."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa\""."\n";
												$gene_number_seq{$name}++;
											}
										}
										for (my $i=0;$i<($length_cp-$lengths);$i+=1){
											my $repeat=substr ($rev_coms,$i,$lengths);
											if ($repeat eq $string) {
												$S=$length_cp-$i-1+1;
												$E=$length_cp-$i-$lengths+1;
												print $out_annotation "     "."gene"."            "."complement(".$E."..".$S.")"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."complement(".$E."..".$S.")"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa\""."\n";
												$gene_number_seq{$name}++;
											}
										}

										if ((defined $length_ref_pcg) and ($length_pcg/$length_ref_pcg < $short)) {
											print $logfile "Warning: $name (positive no-intron PCG) need to be checked due to query coverage per annotated PCG less than $short!\n";
										}elsif ((defined $length_ref_pcg) and ($length_pcg/$length_ref_pcg > $long)) {
											print $logfile "Warning: $name (positive no-intron PCG) need to be checked due to query coverage per annotated PCG more than $long!\n";
										}
									}
								}elsif ($length_pcg < 60) {
									print $logfile "Warning: $name (positive no-intron PCG) has not been annotated due to short length!\n";
								}
							}elsif ((!defined $left_start_position) and (defined $right_star_position)){
								my $end=$end_left+($right_star_position * 3+3);

								my $left_star=$start_left+$left_star_position * 3;
								my $seq_string;
								if ($start2 >= 61) {
									#$seq_string=substr($sequence,($start2-61),120);
									$seq_string=substr($sequence,($left_star),(length($start2-$left_star)+60));
								}elsif ($start2 < 61) {
									$seq_string=substr($sequence,0,120);
								}
								my $length_seq_string=length $seq_string;
								my $aa_seq_string;
								for (my $i=0;$i<$length_seq_string;$i+=3){# delete stop codon
									my $codon=substr ($seq_string,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string.=$hash_codon{$codon};
									}else{
										$aa_seq_string.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}

								my $start;
								for (my $i=0;$i<$xx;$i+=1) {
									my $match_aa_ref_pcg=substr ($aa_ref_pcg,$i,$yy);
									my $marks=0;
									for (my $j=$yy;$j<(($length_seq_string)/3);$j+=1) {
										my $repeat=substr ($aa_seq_string,-$j,$yy);
										if (($repeat eq $match_aa_ref_pcg) and ($marks == 0)) {
											$start=$start2+60-$j*3-$i*3 if ($start2 >= 61);
											$start=121-$j*3-$i*3 if ($start2 < 61);
											$marks++;
										}
									}
									if (defined $start){
										last;
									}
								}
								if (!defined $start) {
									$start=$start2;
								}

								if ($start <= ($start_left+$left_star_position*3)) {
									print $logfile "Warning: $name (positive no-intron PCG) has not been annotated due to internal stop codon!\n";
								}elsif ($start > ($start_left+$left_star_position*3)) {
									my $length_pcg=$end-$start;
									if ($length_pcg >= 60) {
										if ($redundancy eq "N") {
											if ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) < $short)) {
												print $logfile "Warning: $name (positive no-intron PCG) has not been annotated due to query coverage per annotated PCG less than $short!\n";
											}elsif ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) > $long)) {
												print $logfile "Warning: $name (positive no-intron PCG) has not been annotated due to query coverage per annotated PCG more than $long!\n";
											}elsif ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) > $short) and (($length_pcg/$length_ref_pcg) < $long)) {
												print $out_annotation "     "."gene"."            ".$start."..".$end."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             ".$start."..".$end."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa2\""."\n";
												$gene_number_seq{$name}++;
												print $logfile "Warning: $name (positive no-intron PCG) has non-canonical start codon!\n";
											}
										}elsif ($redundancy eq "Y") {
											print $out_annotation "     "."gene"."            ".$start."..".$end."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             ".$start."..".$end."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa2\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (positive no-intron PCG) has non-canonical start codon!\n";

											if ((defined $length_ref_pcg) and ($length_pcg/$length_ref_pcg < $short)) {
												print $logfile "Warning: $name (positive no-intron PCG) need to be checked due to query coverage per annotated PCG less than $short!\n";
											}elsif ((defined $length_ref_pcg) and ($length_pcg/$length_ref_pcg > $long)) {
												print $logfile "Warning: $name (positive no-intron PCG) need to be checked due to query coverage per annotated PCG more than $long!\n";
											}
										}
									}elsif ($length_pcg < 60) {
										print $logfile "Warning: $name (positive no-intron PCG) has not been annotated due to short length!\n";
									}
								}
							}
						}
					}elsif($start1 > $end1){# negative
						my $str=substr($sequence,($end1-1),($start1-$end1+1));
						my $rev_com=reverse $str;
						$rev_com=~ tr/ACGTacgt/TGCAtgca/;
						my $length=length $str;
						my $reverse_start_codon=substr($rev_com,0,3);
						my $reverse_stop_codon=substr($rev_com,($length-3),3);
						my $aa;
						for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
							my $codon=substr ($rev_com,$i,3);
							$codon=uc $codon;
							if (exists $hash_codon{$codon}){
								$aa.=$hash_codon{$codon};
							}else{
								$aa.="X";
								my $j=$i+1;
								#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
							}
						}
						my @aa=split //,$aa;
						if (($length % 3==0) and (!grep {$_=~ /\*/} @aa) and (($reverse_start_codon eq "ATG") or ($reverse_start_codon eq "GTG")) and ((($reverse_stop_codon eq "TAA") or ($reverse_stop_codon eq "TAG") or ($reverse_stop_codon eq "TGA")) or ($name eq "rps12+1"))){# standard start and stop codon
							my $length_pcg=$start1-$end1;
							if ($redundancy eq "N") {
								if ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) < $short)) {
									print $logfile "Warning: $name (negative no-intron PCG) has not been annotated due to query coverage per annotated PCG less than $short!\n";
								}elsif ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) > $long)) {
									print $logfile "Warning: $name (negative no-intron PCG) has not been annotated due to query coverage per annotated PCG more than $long!\n";
								}elsif ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) > $short) and (($length_pcg/$length_ref_pcg) < $long)) {
									print $out_annotation "     "."gene"."            "."complement(".$end1."..".$start1.")\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."CDS"."             "."complement(".$end1."..".$start1.")\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/codon_start=1"."\n";
									print $out_annotation "                     "."/transl_table=11"."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									#print $out_annotation "                     "."/translation=\"$aa\""."\n";
									$gene_number_seq{$name}++;
									my ($S1,$E1);
									my $string=substr($sequence,($end1-1),($start1-$end1+1));
									my $lengths=length $string;
									my $rev_coms=reverse $sequence;
									$rev_coms=~ tr/ACGTacgt/TGCAtgca/;
									for (my $i=0;$i<($length_cp-$lengths);$i+=1){
										my $repeat=substr ($sequence,$i,$lengths);
										if (($repeat eq $string) and (($i+1) ne $end1)) {
											$E1=$i+1;
											$S1=$i+$lengths-1+1;
											print $out_annotation "     "."gene"."            "."complement(".$E1."..".$S1.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."complement(".$E1."..".$S1.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa\""."\n";
											$gene_number_seq{$name}++;
										}
									}
									for (my $i=0;$i<($length_cp-$lengths);$i+=1){
										my $repeat=substr ($rev_coms,$i,$lengths);
										if ($repeat eq $string) {
											$E1=$length_cp-$i-1+1;
											$S1=$length_cp-$i-$lengths+1;
											print $out_annotation "     "."gene"."            ".$S1."..".$E1."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             ".$S1."..".$E1."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa\""."\n";
											$gene_number_seq{$name}++;
										}
									}
								}
							}elsif ($redundancy eq "Y") {
								print $out_annotation "     "."gene"."            "."complement(".$end1."..".$start1.")\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "     "."CDS"."             "."complement(".$end1."..".$start1.")\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "                     "."/codon_start=1"."\n";
								print $out_annotation "                     "."/transl_table=11"."\n";
								print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
								#print $out_annotation "                     "."/translation=\"$aa\""."\n";
								$gene_number_seq{$name}++;
								my ($S1,$E1);
								my $string=substr($sequence,($end1-1),($start1-$end1+1));
								my $lengths=length $string;
								my $rev_coms=reverse $sequence;
								$rev_coms=~ tr/ACGTacgt/TGCAtgca/;
								for (my $i=0;$i<($length_cp-$lengths);$i+=1){
									my $repeat=substr ($sequence,$i,$lengths);
									if (($repeat eq $string) and (($i+1) ne $end1)) {
										$E1=$i+1;
										$S1=$i+$lengths-1+1;
										print $out_annotation "     "."gene"."            "."complement(".$E1."..".$S1.")\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."complement(".$E1."..".$S1.")\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
										$gene_number_seq{$name}++;
									}
								}
								for (my $i=0;$i<($length_cp-$lengths);$i+=1){
									my $repeat=substr ($rev_coms,$i,$lengths);
									if ($repeat eq $string) {
										$E1=$length_cp-$i-1+1;
										$S1=$length_cp-$i-$lengths+1;
										print $out_annotation "     "."gene"."            ".$S1."..".$E1."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             ".$S1."..".$E1."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
										$gene_number_seq{$name}++;
									}
								}

								if ((defined $length_ref_pcg) and ($length_pcg/$length_ref_pcg < $short)) {
									print $logfile "Warning: $name (negative no-intron PCG) need to be checked due to query coverage per annotated PCG less than $short!\n";
								}elsif ((defined $length_ref_pcg) and ($length_pcg/$length_ref_pcg > $long)) {
									print $logfile "Warning: $name (negative no-intron PCG) need to be checked due to query coverage per annotated PCG more than $long!\n";
								}
							}
						}elsif((grep {$_=~ /\*/} @aa) or (($reverse_start_codon ne "ATG") and ($reverse_start_codon ne "GTG")) or ((($reverse_stop_codon ne "TAA") and ($reverse_stop_codon ne "TAG") and ($reverse_stop_codon ne "TGA")) and ($name ne "rps12+1"))){# non-standard start or stop codon
							my $gene_length=($start1-$end1+1);
							my ($start_left,$start_right,$end_left,$end_right);
							if ($gene_length >= 60) {
								$start_left=$start2-60;
							}elsif ($gene_length < 60) {
								$start_left=$start2-int($gene_length/3)*3;
							}
							if (($length_cp-$start2)>=9000) {
								$start_right=$start2+9000;
							}elsif (($length_cp-$start2)<9000){
								$start_right=$start2+(int(($length_cp-$start2)/3)-1)*3;
							}
							if ($end2>=9000) {
								$end_left=$end2-9000;
							}elsif (($end2<9000) and ($end2>30)){
								$end_left=$end2-(int($end2/3)-1)*3;
							}elsif ($end2<=30){
								if ($end2 % 3 == 1) {
									$end_left=1;
								}elsif ($end2 % 3 == 2) {
									$end_left=2;
								}elsif ($end2 % 3 == 0) {
									$end_left=3;
								}
							}
							$end_right=$end2+int($gene_length/3)*3;

							my $str1=substr($sequence,($start_left),($start_right-$start_left));
							my $str2=substr($sequence,($end_left-1),($end_right-$end_left));
							my $length1=length $str1;
							my $length2=length $str2;
							my $rev_com1=reverse $str1;
							$rev_com1=~ tr/ACGTacgt/TGCAtgca/;
							my $rev_com2=reverse $str2;
							$rev_com2=~ tr/ACGTacgt/TGCAtgca/;
							my $aa1;#left range for start codon
							for (my $i=0;$i<($length1-3);$i+=3){# delete stop codon
								my $codon=substr ($rev_com1,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa1.=$hash_codon{$codon};
								}else{
									$aa1.="X";
									my $j=$i+1;
									#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa1=split //,$aa1;
							my (@star1,@start1);
							foreach (0..$#aa1){
								if (($aa1[$_]=~ /\*/) and ($_ < int($length1/3))){
									push @star1,$_;
								}
								if (($aa1[$_]=~ /M/) and ($_ <= int($length1/3))){
									push @start1,$_;
								}
							}
							my $right_star_position=pop @star1;
							my @start2;
							foreach (@start1){
								push @start2,$_ if ((defined $right_star_position) and ($_ > $right_star_position))
							}
							my $right_start_position=$start2[0];
							#my $right_start_position;
							#if ($gene_length > 1000) {
								#$right_start_position=$start2[0];
							#}elsif ($gene_length < 1000) {
								#$right_start_position=pop @start2;
							#}
							my $aa2;#right range for stop codon
							#for (my $i=0;$i<$length2;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
							for (my $i=0;$i<($length2-3);$i+=3){# delete stop codon
								my $codon=substr ($rev_com2,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa2.=$hash_codon{$codon};
								}else{
									$aa2.="X";
									my $j=$i+1;
									#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa2=split //,$aa2;
							my @star2;
							foreach (0..$#aa2){
								if ($aa2[$_]=~ /\*/){
									push @star2,$_;
								}
							}
							my $left_star_position=shift @star2;

							if ((defined $right_start_position) and (defined $left_star_position)){
								my $start=$start_right-$right_start_position * 3;
								my $end=$end_right-($left_star_position * 3+3);

								my $str=substr($sequence,($end-1),($start-$end+1));
								my $length=length $str;
								my $rev_com=reverse $str;
								$rev_com=~ tr/ACGTacgt/TGCAtgca/;
								my $aa;
								for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
									my $codon=substr ($rev_com,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa.=$hash_codon{$codon};
									}else{
										$aa.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}

								my $length_pcg=$start-$end;
								if ($length_pcg >= 60) {
									if ($redundancy eq "N") {
										if ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) < $short)) {
											print $logfile "Warning: $name (negative no-intron PCG) has not been annotated due to query coverage per annotated PCG less than $short!\n";
										}elsif ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) > $long)) {
											print $logfile "Warning: $name (negative no-intron PCG) has not been annotated due to query coverage per annotated PCG more than $long!\n";
										}elsif ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) > $short) and (($length_pcg/$length_ref_pcg) < $long)) {
											print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."complement(".$end."..".$start.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa\""."\n";
											$gene_number_seq{$name}++;
											my ($S,$E);
											my $string=substr($sequence,($end-1),($start-$end+1));
											my $lengths=length $string;
											for (my $i=0;$i<($length_cp-$lengths);$i+=1){
												my $repeat=substr ($sequence,$i,$lengths);
												if (($repeat eq $string) and (($i+1) ne $end)) {
													$E=$i+1;
													$S=$i+$lengths-1+1;
													print $out_annotation "     "."gene"."            "."complement(".$E."..".$S.")\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "     "."CDS"."             "."complement(".$E."..".$S.")\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "                     "."/codon_start=1"."\n";
													print $out_annotation "                     "."/transl_table=11"."\n";
													print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
													#print $out_annotation "                     "."/translation=\"$aa\""."\n";
													$gene_number_seq{$name}++;
												}
											}
											for (my $i=0;$i<($length_cp-$lengths);$i+=1){
												my $repeat=substr ($rev_coms,$i,$lengths);
												if ($repeat eq $string) {
													$E=$length_cp-$i-1+1;
													$S=$length_cp-$i-$lengths+1;
													print $out_annotation "     "."gene"."            ".$S."..".$E."\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "     "."CDS"."             ".$S."..".$E."\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "                     "."/codon_start=1"."\n";
													print $out_annotation "                     "."/transl_table=11"."\n";
													print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
													#print $out_annotation "                     "."/translation=\"$aa\""."\n";
													$gene_number_seq{$name}++;
												}
											}
										}
									}elsif ($redundancy eq "Y") {
										print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."complement(".$end."..".$start.")\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
										$gene_number_seq{$name}++;
										my ($S,$E);
										my $string=substr($sequence,($end-1),($start-$end+1));
										my $lengths=length $string;
										for (my $i=0;$i<($length_cp-$lengths);$i+=1){
											my $repeat=substr ($sequence,$i,$lengths);
											if (($repeat eq $string) and (($i+1) ne $end)) {
												$E=$i+1;
												$S=$i+$lengths-1+1;
												print $out_annotation "     "."gene"."            "."complement(".$E."..".$S.")\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."complement(".$E."..".$S.")\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa\""."\n";
												$gene_number_seq{$name}++;
											}
										}
										for (my $i=0;$i<($length_cp-$lengths);$i+=1){
											my $repeat=substr ($rev_coms,$i,$lengths);
											if ($repeat eq $string) {
												$E=$length_cp-$i-1+1;
												$S=$length_cp-$i-$lengths+1;
												print $out_annotation "     "."gene"."            ".$S."..".$E."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             ".$S."..".$E."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa\""."\n";
												$gene_number_seq{$name}++;
											}
										}

										if ((defined $length_ref_pcg) and ($length_pcg/$length_ref_pcg < $short)) {
											print $logfile "Warning: $name (negative no-intron PCG) need to be checked due to query coverage per annotated PCG less than $short!\n";
										}elsif ((defined $length_ref_pcg) and ($length_pcg/$length_ref_pcg > $long)) {
											print $logfile "Warning: $name (negative no-intron PCG) need to be checked due to query coverage per annotated PCG more than $long!\n";
										}
									}
								}elsif ($length_pcg < 60) {
									print $logfile "Warning: $name (negative no-intron PCG) has not been annotated due to short length!\n";
								}
							}elsif ((!defined $right_start_position) and (defined $left_star_position)){
								my $end=$end_right-($left_star_position * 3+3);

								my $right_star=$start_right-$right_star_position * 3;
								#my $seq_string=substr($sequence,($start2-60),120);
								my $seq_string=substr($sequence,($start2-60),(60+length($right_star-$start2)));
								my $length_seq_string=length $seq_string;
								my $rev_seq_string=reverse $seq_string;
								$rev_seq_string=~ tr/ACGTacgt/TGCAtgca/;
								my $aa_seq_string;
								for (my $i=0;$i<$length_seq_string;$i+=3){# delete stop codon
									my $codon=substr ($rev_seq_string,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string.=$hash_codon{$codon};
									}else{
										$aa_seq_string.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}

								my $start;
								for (my $i=0;$i<$xx;$i+=1) {
									my $match_aa_ref_pcg=substr ($aa_ref_pcg,$i,$yy);
									my $marks=0;
									for (my $j=$yy;$j<(($length_seq_string)/3);$j+=1) {
										my $repeat=substr ($aa_seq_string,-$j,$yy);
										if (($repeat eq $match_aa_ref_pcg) and ($marks == 0)) {
											$start=$start2-60+$j*3+$i*3;
											$marks++;
										}
									}
									if (defined $start){
										last;
									}
								}
								if (!defined $start) {
									$start=$start2;
								}

								if ($start >= ($start_right-$right_star_position*3)) {
									print $logfile "Warning: $name (negative no-intron PCG) has not been annotated due to internal stop codon!\n";
								}elsif ($start < ($start_right-$right_star_position*3)) {
									my $length_pcg=$start-$end;
									if ($length_pcg >= 60) {
										if ($redundancy eq "N") {
											if ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) < $short)) {
												print $logfile "Warning: $name (negative no-intron PCG) has not been annotated due to query coverage per annotated PCG less than $short!\n";
											}elsif ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) > $long)) {
												print $logfile "Warning: $name (negative no-intron PCG) has not been annotated due to query coverage per annotated PCG more than $long!\n";
											}elsif ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) > $short) and (($length_pcg/$length_ref_pcg) < $long)) {
												print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."complement(".$end."..".$start.")\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa2\""."\n";
												$gene_number_seq{$name}++;
												print $logfile "Warning: $name (negative no-intron PCG) has non-canonical start codon!\n";
											}
										}elsif ($redundancy eq "Y") {
											print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."complement(".$end."..".$start.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa2\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (negative no-intron PCG) has non-canonical start codon!\n";

											if ((defined $length_ref_pcg) and ($length_pcg/$length_ref_pcg < $short)) {
												print $logfile "Warning: $name (negative no-intron PCG) need to be checked due to query coverage per annotated PCG less than $short!\n";
											}elsif ((defined $length_ref_pcg) and ($length_pcg/$length_ref_pcg > $long)) {
												print $logfile "Warning: $name (negative no-intron PCG) need to be checked due to query coverage per annotated PCG more than $long!\n";
											}
										}
									}elsif ($length_pcg < 60) {
										print $logfile "Warning: $name (negative no-intron PCG) has not been annotated due to short length!\n";
									}
								}
							}
						}
					}
				}elsif(($start1 != $start2) or ($end1 != $end2)){# non-identical PCG boundary for _gene and _CDS_aa
					if ($start2 < $end2){# positive
						my $str1=substr($sequence,($start1-1),($end1-$start1+1));
						my $length1=length $str1;
						my $forward_start_codon1=substr($str1,0,3);
						my $forward_stop_codon1=substr($str1,($length1-3),3);
						my $str2=substr($sequence,($start2-1),($end2-$start2+1));
						my $length2=length $str2;
						my $forward_start_codon2=substr($str2,0,3);
						my $forward_stop_codon2=substr($str2,($length2-3),3);
						my $aa1;
						for (my $i=0;$i<($length1-3);$i+=3){# delete stop codon
							my $codon=substr ($str1,$i,3);
							$codon=uc $codon;
							if (exists $hash_codon{$codon}){
								$aa1.=$hash_codon{$codon};
							}else{
								$aa1.="X";
								my $j=$i+1;
								#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
							}
						}
						my @aa1=split //,$aa1;
						my $aa2;
						#for (my $i=0;$i<$length2;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
						for (my $i=0;$i<($length2-3);$i+=3){# delete stop codon
							my $codon=substr ($str2,$i,3);
							$codon=uc $codon;
							if (exists $hash_codon{$codon}){
								$aa2.=$hash_codon{$codon};
							}else{
								$aa2.="X";
								my $j=$i+1;
								#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
							}
						}
						my @aa2=split //,$aa2;
						if (($length2 % 3==0) and (!grep {$_=~ /\*/} @aa2) and (($forward_start_codon2 eq "ATG") or ($forward_start_codon2 eq "GTG")) and (($forward_stop_codon2 eq "TAA") or ($forward_stop_codon2 eq "TAG") or ($forward_stop_codon2 eq "TGA") or ($name eq "rps12+1"))){# standard start and stop codon for _CDS_aa
							my $length_pcg=$end2-$start2;
							if ($redundancy eq "N") {
								if ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) < $short)) {
									print $logfile "Warning: $name (positive no-intron PCG) has not been annotated due to query coverage per annotated PCG less than $short!\n";
								}elsif ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) > $long)) {
									print $logfile "Warning: $name (positive no-intron PCG) has not been annotated due to query coverage per annotated PCG more than $long!\n";
								}elsif ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) > $short) and (($length_pcg/$length_ref_pcg) < $long)) {
									print $out_annotation "     "."gene"."            ".$start2."..".$end2."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."CDS"."             ".$start2."..".$end2."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/codon_start=1"."\n";
									print $out_annotation "                     "."/transl_table=11"."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									#print $out_annotation "                     "."/translation=\"$aa2\""."\n";
									$gene_number_seq{$name}++;
									my ($S,$E);
									my $string=substr($sequence,($start2-1),($end2-$start2+1));
									my $lengths=length $string;
									my $rev_coms=reverse $sequence;
									$rev_coms=~ tr/ACGTacgt/TGCAtgca/;
									for (my $i=0;$i<($length_cp-$lengths);$i+=1){
										my $repeat=substr ($sequence,$i,$lengths);
										if (($repeat eq $string) and (($i+1) ne $start2)) {
											$S=$i+1;
											$E=$i+$lengths-1+1;
											print $out_annotation "     "."gene"."            ".$S."..".$E."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             ".$S."..".$E."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa2\""."\n";
											$gene_number_seq{$name}++;
										}
									}
									for (my $i=0;$i<($length_cp-$lengths);$i+=1){
										my $repeat=substr ($rev_coms,$i,$lengths);
										if ($repeat eq $string) {
											$S=$length_cp-$i-1+1;
											$E=$length_cp-$i-$lengths+1;
											print $out_annotation "     "."gene"."            "."complement(".$E."..".$S.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."complement(".$E."..".$S.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa2\""."\n";
											$gene_number_seq{$name}++;
										}
									}
								}
							}elsif ($redundancy eq "Y") {
								print $out_annotation "     "."gene"."            ".$start2."..".$end2."\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "     "."CDS"."             ".$start2."..".$end2."\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "                     "."/codon_start=1"."\n";
								print $out_annotation "                     "."/transl_table=11"."\n";
								print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
								#print $out_annotation "                     "."/translation=\"$aa2\""."\n";
								$gene_number_seq{$name}++;
								my ($S,$E);
								my $string=substr($sequence,($start2-1),($end2-$start2+1));
								my $lengths=length $string;
								my $rev_coms=reverse $sequence;
								$rev_coms=~ tr/ACGTacgt/TGCAtgca/;
								for (my $i=0;$i<($length_cp-$lengths);$i+=1){
									my $repeat=substr ($sequence,$i,$lengths);
									if (($repeat eq $string) and (($i+1) ne $start2)) {
										$S=$i+1;
										$E=$i+$lengths-1+1;
										print $out_annotation "     "."gene"."            ".$S."..".$E."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             ".$S."..".$E."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa2\""."\n";
										$gene_number_seq{$name}++;
									}
								}
								for (my $i=0;$i<($length_cp-$lengths);$i+=1){
									my $repeat=substr ($rev_coms,$i,$lengths);
									if ($repeat eq $string) {
										$S=$length_cp-$i-1+1;
										$E=$length_cp-$i-$lengths+1;
										print $out_annotation "     "."gene"."            "."complement(".$E."..".$S.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."complement(".$E."..".$S.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa2\""."\n";
										$gene_number_seq{$name}++;
									}
								}

								if ((defined $length_ref_pcg) and ($length_pcg/$length_ref_pcg < $short)) {
									print $logfile "Warning: $name (positive no-intron PCG) need to be checked due to query coverage per annotated PCG less than $short!\n";
								}elsif ((defined $length_ref_pcg) and ($length_pcg/$length_ref_pcg > $long)) {
									print $logfile "Warning: $name (positive no-intron PCG) need to be checked due to query coverage per annotated PCG more than $long!\n";
								}
							}

#						}elsif(($length1 % 3==0) and (!grep {$_=~ /\*/} @aa1) and (($forward_start_codon1 eq "ATG") or ($forward_start_codon1 eq "GTG")) and (($forward_stop_codon1 eq "TAA") or ($forward_stop_codon1 eq "TAG") or ($forward_stop_codon1 eq "TGA"))){# standard start and stop codon for _gene
#							print $out_annotation "     "."gene"."            ".$start1."..".$end1."\n";
#							print $out_annotation "                     "."/gene=\"$name\""."\n";
#							print $out_annotation "     "."CDS"."             ".$start1."..".$end1."\n";
#							print $out_annotation "                     "."/gene=\"$name\""."\n";
#							print $out_annotation "                     "."/codon_start=1"."\n";
#							print $out_annotation "                     "."/transl_table=11"."\n";
#							print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
#							#print $out_annotation "                     "."/translation=\"$aa1\""."\n";
#							$gene_number_seq{$name}++;
#							my ($S,$E);
#							my $string=substr($sequence,($start1-1),($end1-$start1+1));
#							my $lengths=length $string;
#							my $rev_coms=reverse $sequence;
#							$rev_coms=~ tr/ACGTacgt/TGCAtgca/;
#							for (my $i=0;$i<($length_cp-$lengths);$i+=1){
#								my $repeat=substr ($sequence,$i,$lengths);
#								if (($repeat eq $string) and (($i+1) ne $start1)) {
#									$S=$i+1;
#									$E=$i+$lengths-1+1;
#									print $out_annotation "     "."gene"."            ".$S."..".$E."\n";
#									print $out_annotation "                     "."/gene=\"$name\""."\n";
#									print $out_annotation "     "."CDS"."             ".$S."..".$E."\n";
#									print $out_annotation "                     "."/gene=\"$name\""."\n";
#									print $out_annotation "                     "."/codon_start=1"."\n";
#									print $out_annotation "                     "."/transl_table=11"."\n";
#									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
#									#print $out_annotation "                     "."/translation=\"$aa1\""."\n";
#									$gene_number_seq{$name}++;
#								}
#							}
#							for (my $i=0;$i<($length_cp-$lengths);$i+=1){
#								my $repeat=substr ($rev_coms,$i,$lengths);
#								if ($repeat eq $string) {
#									$S=$length_cp-$i-1+1;
#									$E=$length_cp-$i-$lengths+1;
#									print $out_annotation "     "."gene"."            "."complement(".$E."..".$S.")"."\n";
#									print $out_annotation "                     "."/gene=\"$name\""."\n";
#									print $out_annotation "     "."CDS"."             "."complement(".$E."..".$S.")"."\n";
#									print $out_annotation "                     "."/gene=\"$name\""."\n";
#									print $out_annotation "                     "."/codon_start=1"."\n";
#									print $out_annotation "                     "."/transl_table=11"."\n";
#									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
#									#print $out_annotation "                     "."/translation=\"$aa1\""."\n";
#									$gene_number_seq{$name}++;
#								}
#							}
						}elsif((grep {$_=~ /\*/} @aa2) or (($forward_start_codon2 ne "ATG") and ($forward_start_codon2 ne "GTG")) or (($forward_stop_codon2 ne "TAA") and ($forward_stop_codon2 ne "TAG") and ($forward_stop_codon2 ne "TGA") and ($name ne "rps12+1"))){# non-standard start or stop codon for _CDS_aa
							my $gene_length=($end2-$start2+1);
							my ($start_left,$start_right,$end_left,$end_right);
							if ($start2>=9000) {
								$start_left=$start2-9000;
							}elsif (($start2<9000) and ($start2>30)){
								$start_left=$start2-(int($start2/3)-1)*3;
							}elsif ($start2<=30){
								if ($start2 % 3 == 1) {
									$start_left=1;
								}elsif ($start2 % 3 == 2) {
									$start_left=2;
								}elsif ($start2 % 3 == 0) {
									$start_left=3;
								}
							}
							if ($gene_length >= 60) {
								$start_right=$start2+60;
							}elsif ($gene_length < 60) {
								$start_right=$start2+int($gene_length/3)*3;
							}
							$end_left=$end2-int($gene_length/3)*3;
							if (($length_cp-$end2)>=9000) {
								$end_right=$end2+9000;
							}elsif (($length_cp-$end2)<9000){
								$end_right=$end2+(int(($length_cp-$end2)/3)-1)*3;
							}

							my $str1=substr($sequence,($start_left-1),($start_right-$start_left));
							my $str2=substr($sequence,($end_left),($end_right-$end_left));
							my $length1=length $str1;
							my $length2=length $str2;
							my $aa1;#left range for start codon
							for (my $i=0;$i<($length1-3);$i+=3){# delete stop codon
								my $codon=substr ($str1,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa1.=$hash_codon{$codon};
								}else{
									$aa1.="X";
									my $j=$i+1;
									#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa1=split //,$aa1;
							my (@star1,@start1);
							foreach (0..$#aa1){
								if (($aa1[$_]=~ /\*/) and ($_ < int($length1/3))){
									push @star1,$_;
								}
								if (($aa1[$_]=~ /M/) and ($_ <= int($length1/3))){
									push @start1,$_;
								}
							}
							my $left_star_position=pop @star1;
							my @start2;
							foreach (@start1){
								push @start2,$_ if ((defined $left_star_position) and ($_ > $left_star_position))
							}
							my $left_start_position=$start2[0];
							#my $left_start_position;
							#if ($gene_length > 1000) {
								#$left_start_position=$start2[0];
							#}elsif ($gene_length < 1000) {
								#$left_start_position=pop @start2;
							#}
							my $aa2;#right range for stop codon
							for (my $i=0;$i<($length2-3);$i+=3){# delete stop codon
								my $codon=substr ($str2,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa2.=$hash_codon{$codon};
								}else{
									$aa2.="X";
									my $j=$i+1;
									#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa2=split //,$aa2;
							my @star2;
							foreach (0..$#aa2){
								if ($aa2[$_]=~ /\*/){
									push @star2,$_;
								}
							}
							my $right_star_position=shift @star2;

							if ((defined $left_start_position) and (defined $right_star_position)){
								my $start=$start_left+$left_start_position * 3;
								my $end=$end_left+($right_star_position * 3+3);
								my $str=substr($sequence,($start-1),($end-$start+1));
								my $length=length $str;
								my $aa;
								for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
									my $codon=substr ($str,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa.=$hash_codon{$codon};
									}else{
										$aa.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}

								my $length_pcg=$end-$start;
								if ($length_pcg >= 60) {
									if ($redundancy eq "N") {
										if ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) < $short)) {
											print $logfile "Warning: $name (positive no-intron PCG) has not been annotated due to query coverage per annotated PCG less than $short!\n";
										}elsif ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) > $long)) {
											print $logfile "Warning: $name (positive no-intron PCG) has not been annotated due to query coverage per annotated PCG more than $long!\n";
										}elsif ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) > $short) and (($length_pcg/$length_ref_pcg) < $long)) {
											print $out_annotation "     "."gene"."            ".$start."..".$end."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             ".$start."..".$end."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa\""."\n";
											$gene_number_seq{$name}++;
											my ($S,$E);
											my $string=substr($sequence,($start-1),($end-$start+1));
											my $lengths=length $string;
											for (my $i=0;$i<($length_cp-$lengths);$i+=1){
												my $repeat=substr ($sequence,$i,$lengths);
												if (($repeat eq $string) and (($i+1) ne $start)) {
													$S=$i+1;
													$E=$i+$lengths-1+1;
													print $out_annotation "     "."gene"."            ".$S."..".$E."\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "     "."CDS"."             ".$S."..".$E."\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "                     "."/codon_start=1"."\n";
													print $out_annotation "                     "."/transl_table=11"."\n";
													print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
													#print $out_annotation "                     "."/translation=\"$aa\""."\n";
													$gene_number_seq{$name}++;
												}
											}
											for (my $i=0;$i<($length_cp-$lengths);$i+=1){
												my $repeat=substr ($rev_coms,$i,$lengths);
												if ($repeat eq $string) {
													$S=$length_cp-$i-1+1;
													$E=$length_cp-$i-$lengths+1;
													print $out_annotation "     "."gene"."            "."complement(".$E."..".$S.")"."\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "     "."CDS"."             "."complement(".$E."..".$S.")"."\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "                     "."/codon_start=1"."\n";
													print $out_annotation "                     "."/transl_table=11"."\n";
													print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
													#print $out_annotation "                     "."/translation=\"$aa\""."\n";
													$gene_number_seq{$name}++;
												}
											}
										}
									}elsif ($redundancy eq "Y") {
										print $out_annotation "     "."gene"."            ".$start."..".$end."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             ".$start."..".$end."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
										$gene_number_seq{$name}++;
										my ($S,$E);
										my $string=substr($sequence,($start-1),($end-$start+1));
										my $lengths=length $string;
										for (my $i=0;$i<($length_cp-$lengths);$i+=1){
											my $repeat=substr ($sequence,$i,$lengths);
											if (($repeat eq $string) and (($i+1) ne $start)) {
												$S=$i+1;
												$E=$i+$lengths-1+1;
												print $out_annotation "     "."gene"."            ".$S."..".$E."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             ".$S."..".$E."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa\""."\n";
												$gene_number_seq{$name}++;
											}
										}
										for (my $i=0;$i<($length_cp-$lengths);$i+=1){
											my $repeat=substr ($rev_coms,$i,$lengths);
											if ($repeat eq $string) {
												$S=$length_cp-$i-1+1;
												$E=$length_cp-$i-$lengths+1;
												print $out_annotation "     "."gene"."            "."complement(".$E."..".$S.")"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."complement(".$E."..".$S.")"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa\""."\n";
												$gene_number_seq{$name}++;
											}
										}

										if ((defined $length_ref_pcg) and ($length_pcg/$length_ref_pcg < $short)) {
											print $logfile "Warning: $name (positive no-intron PCG) need to be checked due to query coverage per annotated PCG less than $short!\n";
										}elsif ((defined $length_ref_pcg) and ($length_pcg/$length_ref_pcg > $long)) {
											print $logfile "Warning: $name (positive no-intron PCG) need to be checked due to query coverage per annotated PCG more than $long!\n";
										}
									}
								}elsif ($length_pcg < 60) {
									print $logfile "Warning: $name (positive no-intron PCG) has not been annotated due to short length!\n";
								}
							}elsif ((!defined $left_start_position) and (defined $right_star_position)){
								my $end=$end_left+($right_star_position * 3+3);

								my $left_star=$start_left+$left_star_position * 3;
								my $seq_string;
								if ($start2 >= 61) {
									#$seq_string=substr($sequence,($start2-61),120);
									$seq_string=substr($sequence,($left_star),(length($start2-$left_star)+60));
								}elsif ($start2 < 61) {
									$seq_string=substr($sequence,0,120);
								}
								my $length_seq_string=length $seq_string;
								my $aa_seq_string;
								for (my $i=0;$i<$length_seq_string;$i+=3){# delete stop codon
									my $codon=substr ($seq_string,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string.=$hash_codon{$codon};
									}else{
										$aa_seq_string.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}

								my $start;
								for (my $i=0;$i<$xx;$i+=1) {
									my $match_aa_ref_pcg=substr ($aa_ref_pcg,$i,$yy);
									my $marks=0;
									for (my $j=$yy;$j<(($length_seq_string)/3);$j+=1) {
										my $repeat=substr ($aa_seq_string,-$j,$yy);
										if (($repeat eq $match_aa_ref_pcg) and ($marks == 0)) {
											$start=$start2+60-$j*3-$i*3 if ($start2 >= 61);
											$start=121-$j*3-$i*3 if ($start2 < 61);
											$marks++;
										}
									}
									if (defined $start){
										last;
									}
								}
								if (!defined $start) {
									$start=$start2;
								}

								if ($start <= ($start_left+$left_star_position*3)) {
									print $logfile "Warning: $name (positive no-intron PCG) has not been annotated due to internal stop codon!\n";
								}elsif ($start > ($start_left+$left_star_position*3)) {
									my $length_pcg=$end-$start;
									if ($length_pcg >= 60) {
										if ($redundancy eq "N") {
											if ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) < $short)) {
												print $logfile "Warning: $name (positive no-intron PCG) has not been annotated due to query coverage per annotated PCG less than $short!\n";
											}elsif ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) > $long)) {
												print $logfile "Warning: $name (positive no-intron PCG) has not been annotated due to query coverage per annotated PCG more than $long!\n";
											}elsif ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) > $short) and (($length_pcg/$length_ref_pcg) < $long)) {
												print $out_annotation "     "."gene"."            ".$start."..".$end."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             ".$start."..".$end."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa2\""."\n";
												$gene_number_seq{$name}++;
												print $logfile "Warning: $name (positive no-intron PCG) has non-canonical start codon!\n";
											}
										}elsif ($redundancy eq "Y") {
											print $out_annotation "     "."gene"."            ".$start."..".$end."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             ".$start."..".$end."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa2\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (positive no-intron PCG) has non-canonical start codon!\n";

											if ((defined $length_ref_pcg) and ($length_pcg/$length_ref_pcg < $short)) {
												print $logfile "Warning: $name (positive no-intron PCG) need to be checked due to query coverage per annotated PCG less than $short!\n";
											}elsif ((defined $length_ref_pcg) and ($length_pcg/$length_ref_pcg > $long)) {
												print $logfile "Warning: $name (positive no-intron PCG) need to be checked due to query coverage per annotated PCG more than $long!\n";
											}
										}
									}elsif ($length_pcg < 60) {
										print $logfile "Warning: $name (positive no-intron PCG) has not been annotated due to short length!\n";
									}
								}
							}
						}
					}elsif($start2 > $end2){# negative
						my $str1=substr($sequence,($end1-1),($start1-$end1+1));
						my $rev_com1=reverse $str1;
						$rev_com1=~ tr/ACGTacgt/TGCAtgca/;
						my $length1=length $str1;
						my $reverse_start_codon1=substr($rev_com1,0,3);
						my $reverse_stop_codon1=substr($rev_com1,($length1-3),3);
						my $str2=substr($sequence,($end2-1),($start2-$end2+1));
						my $rev_com2=reverse $str2;
						$rev_com2=~ tr/ACGTacgt/TGCAtgca/;
						my $length2=length $str2;
						my $reverse_start_codon2=substr($rev_com2,0,3);
						my $reverse_stop_codon2=substr($rev_com2,($length2-3),3);
						my $aa1;
						for (my $i=0;$i<($length1-3);$i+=3){# delete stop codon
							my $codon=substr ($rev_com1,$i,3);
							$codon=uc $codon;
							if (exists $hash_codon{$codon}){
								$aa1.=$hash_codon{$codon};
							}else{
								$aa1.="X";
								my $j=$i+1;
								#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
							}
						}
						my @aa1=split //,$aa1;
						my $aa2;
						for (my $i=0;$i<($length2-3);$i+=3){# delete stop codon
							my $codon=substr ($rev_com2,$i,3);
							$codon=uc $codon;
							if (exists $hash_codon{$codon}){
								$aa2.=$hash_codon{$codon};
							}else{
								$aa2.="X";
								my $j=$i+1;
								#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
							}
						}
						my @aa2=split //,$aa2;
						if (($length2 % 3==0) and (!grep {$_=~ /\*/} @aa2) and (($reverse_start_codon2 eq "ATG") or ($reverse_start_codon2 eq "GTG")) and (($reverse_stop_codon2 eq "TAA") or ($reverse_stop_codon2 eq "TAG") or ($reverse_stop_codon2 eq "TGA") or ($name eq "rps12+1"))){# standard start and stop codon for _CDS_aa
							my $length_pcg=$start2-$end2;
							if ($redundancy eq "N") {
								if ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) < $short)) {
									print $logfile "Warning: $name (negative no-intron PCG) has not been annotated due to query coverage per annotated PCG less than $short!\n";
								}elsif ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) > $long)) {
									print $logfile "Warning: $name (negative no-intron PCG) has not been annotated due to query coverage per annotated PCG more than $long!\n";
								}elsif ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) > $short) and (($length_pcg/$length_ref_pcg) < $long)) {
									print $out_annotation "     "."gene"."            "."complement(".$end2."..".$start2.")\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."CDS"."             "."complement(".$end2."..".$start2.")\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/codon_start=1"."\n";
									print $out_annotation "                     "."/transl_table=11"."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									#print $out_annotation "                     "."/translation=\"$aa2\""."\n";
									$gene_number_seq{$name}++;
									my ($S,$E);
									my $string=substr($sequence,($end2-1),($start2-$end2+1));
									my $lengths=length $string;
									my $rev_coms=reverse $sequence;
									$rev_coms=~ tr/ACGTacgt/TGCAtgca/;
									for (my $i=0;$i<($length_cp-$lengths);$i+=1){
										my $repeat=substr ($sequence,$i,$lengths);
										if (($repeat eq $string) and (($i+1) ne $end2)) {
											$E=$i+1;
											$S=$i+$lengths-1+1;
											print $out_annotation "     "."gene"."            "."complement(".$E."..".$S.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."complement(".$E."..".$S.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa2\""."\n";
											$gene_number_seq{$name}++;
										}
									}
									for (my $i=0;$i<($length_cp-$lengths);$i+=1){
										my $repeat=substr ($rev_coms,$i,$lengths);
										if ($repeat eq $string) {
											$E=$length_cp-$i-1+1;
											$S=$length_cp-$i-$lengths+1;
											print $out_annotation "     "."gene"."            ".$S."..".$E."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             ".$S."..".$E."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa2\""."\n";
											$gene_number_seq{$name}++;
										}
									}
								}
							}elsif ($redundancy eq "Y") {
								print $out_annotation "     "."gene"."            "."complement(".$end2."..".$start2.")\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "     "."CDS"."             "."complement(".$end2."..".$start2.")\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "                     "."/codon_start=1"."\n";
								print $out_annotation "                     "."/transl_table=11"."\n";
								print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
								#print $out_annotation "                     "."/translation=\"$aa2\""."\n";
								$gene_number_seq{$name}++;
								my ($S,$E);
								my $string=substr($sequence,($end2-1),($start2-$end2+1));
								my $lengths=length $string;
								my $rev_coms=reverse $sequence;
								$rev_coms=~ tr/ACGTacgt/TGCAtgca/;
								for (my $i=0;$i<($length_cp-$lengths);$i+=1){
									my $repeat=substr ($sequence,$i,$lengths);
									if (($repeat eq $string) and (($i+1) ne $end2)) {
										$E=$i+1;
										$S=$i+$lengths-1+1;
										print $out_annotation "     "."gene"."            "."complement(".$E."..".$S.")\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."complement(".$E."..".$S.")\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa2\""."\n";
										$gene_number_seq{$name}++;
									}
								}
								for (my $i=0;$i<($length_cp-$lengths);$i+=1){
									my $repeat=substr ($rev_coms,$i,$lengths);
									if ($repeat eq $string) {
										$E=$length_cp-$i-1+1;
										$S=$length_cp-$i-$lengths+1;
										print $out_annotation "     "."gene"."            ".$S."..".$E."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             ".$S."..".$E."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa2\""."\n";
										$gene_number_seq{$name}++;
									}
								}

								if ((defined $length_ref_pcg) and ($length_pcg/$length_ref_pcg < $short)) {
									print $logfile "Warning: $name (negative no-intron PCG) need to be checked due to query coverage per annotated PCG less than $short!\n";
								}elsif ((defined $length_ref_pcg) and ($length_pcg/$length_ref_pcg > $long)) {
									print $logfile "Warning: $name (negative no-intron PCG) need to be checked due to query coverage per annotated PCG more than $long!\n";
								}
							}

#						}elsif(($length1 % 3==0) and (!grep {$_=~ /\*/} @aa1) and (($reverse_start_codon1 eq "ATG") or ($reverse_start_codon1 eq "GTG")) and (($reverse_stop_codon1 eq "TAA") or ($reverse_stop_codon1 eq "TAG") or ($reverse_stop_codon1 eq "TGA"))){# standard start and stop codon for _gene
#							print $out_annotation "     "."gene"."            "."complement(".$end1."..".$start1.")\n";
#							print $out_annotation "                     "."/gene=\"$name\""."\n";
#							print $out_annotation "     "."CDS"."             "."complement(".$end1."..".$start1.")\n";
#							print $out_annotation "                     "."/gene=\"$name\""."\n";
#							print $out_annotation "                     "."/codon_start=1"."\n";
#							print $out_annotation "                     "."/transl_table=11"."\n";
#							print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
#							#print $out_annotation "                     "."/translation=\"$aa1\""."\n";
#							$gene_number_seq{$name}++;
#							my ($S,$E);
#							my $string=substr($sequence,($end1-1),($start1-$end1+1));
#							my $lengths=length $string;
#							my $rev_coms=reverse $sequence;
#							$rev_coms=~ tr/ACGTacgt/TGCAtgca/;
#							for (my $i=0;$i<($length_cp-$lengths);$i+=1){
#								my $repeat=substr ($sequence,$i,$lengths);
#								if (($repeat eq $string) and (($i+1) ne $end1)) {
#									$E=$i+1;
#									$S=$i+$lengths-1+1;
#									print $out_annotation "     "."gene"."            "."complement(".$E."..".$S.")\n";
#									print $out_annotation "                     "."/gene=\"$name\""."\n";
#									print $out_annotation "     "."CDS"."             "."complement(".$E."..".$S.")\n";
#									print $out_annotation "                     "."/gene=\"$name\""."\n";
#									print $out_annotation "                     "."/codon_start=1"."\n";
#									print $out_annotation "                     "."/transl_table=11"."\n";
#									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
#									#print $out_annotation "                     "."/translation=\"$aa1\""."\n";
#									$gene_number_seq{$name}++;
#								}
#							}
#							for (my $i=0;$i<($length_cp-$lengths);$i+=1){
#								my $repeat=substr ($rev_coms,$i,$lengths);
#								if ($repeat eq $string) {
#									$E=$length_cp-$i-1+1;
#									$S=$length_cp-$i-$lengths+1;
#									print $out_annotation "     "."gene"."            ".$S."..".$E."\n";
#									print $out_annotation "                     "."/gene=\"$name\""."\n";
#									print $out_annotation "     "."CDS"."             ".$S."..".$E."\n";
#									print $out_annotation "                     "."/gene=\"$name\""."\n";
#									print $out_annotation "                     "."/codon_start=1"."\n";
#									print $out_annotation "                     "."/transl_table=11"."\n";
#									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
#									#print $out_annotation "                     "."/translation=\"$aa1\""."\n";
#									$gene_number_seq{$name}++;
#								}
#							}
						}elsif((grep {$_=~ /\*/} @aa2) or (($reverse_start_codon2 ne "ATG") and ($reverse_start_codon2 ne "GTG")) or (($reverse_stop_codon2 ne "TAA") and ($reverse_stop_codon2 ne "TAG") and ($reverse_stop_codon2 ne "TGA") and ($name ne "rps12+1"))){#non-standard start or stop codon for _CDS_aa
							my $gene_length=($start2-$end2+1);
							my ($start_left,$start_right,$end_left,$end_right);
							if ($gene_length >= 60) {
								$start_left=$start2-60;
							}elsif ($gene_length < 60) {
								$start_left=$start2-int($gene_length/3)*3;
							}
							if (($length_cp-$start2)>=9000) {
								$start_right=$start2+9000;
							}elsif (($length_cp-$start2)<9000){
								$start_right=$start2+(int(($length_cp-$start2)/3)-1)*3;
							}
							if ($end2>=9000) {
								$end_left=$end2-9000;
							}elsif (($end2<9000) and ($end2>30)){
								$end_left=$end2-(int($end2/3)-1)*3;
							}elsif ($end2<=30){
								if ($end2 % 3 == 1) {
									$end_left=1;
								}elsif ($end2 % 3 == 2) {
									$end_left=2;
								}elsif ($end2 % 3 == 0) {
									$end_left=3;
								}
							}
							$end_right=$end2+int($gene_length/3)*3;

							my $str1=substr($sequence,($start_left),($start_right-$start_left));
							my $str2=substr($sequence,($end_left-1),($end_right-$end_left));
							my $length1=length $str1;
							my $length2=length $str2;
							my $rev_com1=reverse $str1;
							$rev_com1=~ tr/ACGTacgt/TGCAtgca/;
							my $rev_com2=reverse $str2;
							$rev_com2=~ tr/ACGTacgt/TGCAtgca/;
							my $aa1;#left range for start codon
							for (my $i=0;$i<($length1-3);$i+=3){# delete stop codon
								my $codon=substr ($rev_com1,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa1.=$hash_codon{$codon};
								}else{
									$aa1.="X";
									my $j=$i+1;
									#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa1=split //,$aa1;
							my (@star1,@start1);
							foreach (0..$#aa1){
								if (($aa1[$_]=~ /\*/) and ($_ < int($length1/3))){
									push @star1,$_;
								}
								if (($aa1[$_]=~ /M/) and ($_ <= int($length1/3))){
									push @start1,$_;
								}
							}
							my $right_star_position=pop @star1;
							my @start2;
							foreach (@start1){
								push @start2,$_ if ((defined $right_star_position) and ($_ > $right_star_position))
							}
							my $right_start_position=$start2[0];
							#my $right_start_position;
							#if ($gene_length > 1000) {
								#$right_start_position=$start2[0];
							#}elsif ($gene_length < 1000) {
								#$right_start_position=pop @start2;
							#}
							my $aa2;#right range for stop codon
							for (my $i=0;$i<$length2;$i+=3){# delete stop codon
								my $codon=substr ($rev_com2,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa2.=$hash_codon{$codon};
								}else{
									$aa2.="X";
									my $j=$i+1;
									#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa2=split //,$aa2;
							my @star2;
							foreach (0..$#aa2){
								if ($aa2[$_]=~ /\*/){
									push @star2,$_;
								}
							}
							my $left_star_position=shift @star2;

							if ((defined $right_start_position) and (defined $left_star_position)){
								my $start=$start_right-$right_start_position * 3;
								my $end=$end_right-($left_star_position * 3+3);
								my $str=substr($sequence,($end-1),($start-$end+1));
								my $length=length $str;
								my $rev_com=reverse $str;
								$rev_com=~ tr/ACGTacgt/TGCAtgca/;
								my $aa;
								for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
									my $codon=substr ($rev_com,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa.=$hash_codon{$codon};
									}else{
										$aa.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}

								my $length_pcg=$start-$end;
								if ($length_pcg >= 60) {
									if ($redundancy eq "N") {
										if ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) < $short)) {
											print $logfile "Warning: $name (negative no-intron PCG) has not been annotated due to query coverage per annotated PCG less than $short!\n";
										}elsif ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) > $long)) {
											print $logfile "Warning: $name (negative no-intron PCG) has not been annotated due to query coverage per annotated PCG more than $long!\n";
										}elsif ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) > $short) and (($length_pcg/$length_ref_pcg) < $long)) {
											print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."complement(".$end."..".$start.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa\""."\n";
											$gene_number_seq{$name}++;
											my ($S,$E);
											my $string=substr($sequence,($end-1),($start-$end+1));
											my $lengths=length $string;
											for (my $i=0;$i<($length_cp-$lengths);$i+=1){
												my $repeat=substr ($sequence,$i,$lengths);
												if (($repeat eq $string) and (($i+1) ne $end)) {
													$E=$i+1;
													$S=$i+$lengths-1+1;
													print $out_annotation "     "."gene"."            "."complement(".$E."..".$S.")\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "     "."CDS"."             "."complement(".$E."..".$S.")\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "                     "."/codon_start=1"."\n";
													print $out_annotation "                     "."/transl_table=11"."\n";
													print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
													#print $out_annotation "                     "."/translation=\"$aa\""."\n";
													$gene_number_seq{$name}++;
												}
											}
											for (my $i=0;$i<($length_cp-$lengths);$i+=1){
												my $repeat=substr ($rev_coms,$i,$lengths);
												if ($repeat eq $string) {
													$E=$length_cp-$i-1+1;
													$S=$length_cp-$i-$lengths+1;
													print $out_annotation "     "."gene"."            ".$S."..".$E."\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "     "."CDS"."             ".$S."..".$E."\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "                     "."/codon_start=1"."\n";
													print $out_annotation "                     "."/transl_table=11"."\n";
													print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
													#print $out_annotation "                     "."/translation=\"$aa\""."\n";
													$gene_number_seq{$name}++;
												}
											}
										}
									}elsif ($redundancy eq "Y") {
										print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."complement(".$end."..".$start.")\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
										$gene_number_seq{$name}++;
										my ($S,$E);
										my $string=substr($sequence,($end-1),($start-$end+1));
										my $lengths=length $string;
										for (my $i=0;$i<($length_cp-$lengths);$i+=1){
											my $repeat=substr ($sequence,$i,$lengths);
											if (($repeat eq $string) and (($i+1) ne $end)) {
												$E=$i+1;
												$S=$i+$lengths-1+1;
												print $out_annotation "     "."gene"."            "."complement(".$E."..".$S.")\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."complement(".$E."..".$S.")\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa\""."\n";
												$gene_number_seq{$name}++;
											}
										}
										for (my $i=0;$i<($length_cp-$lengths);$i+=1){
											my $repeat=substr ($rev_coms,$i,$lengths);
											if ($repeat eq $string) {
												$E=$length_cp-$i-1+1;
												$S=$length_cp-$i-$lengths+1;
												print $out_annotation "     "."gene"."            ".$S."..".$E."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             ".$S."..".$E."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa\""."\n";
												$gene_number_seq{$name}++;
											}
										}

										if ((defined $length_ref_pcg) and ($length_pcg/$length_ref_pcg < $short)) {
											print $logfile "Warning: $name (negative no-intron PCG) need to be checked due to query coverage per annotated PCG less than $short!\n";
										}elsif ((defined $length_ref_pcg) and ($length_pcg/$length_ref_pcg > $long)) {
											print $logfile "Warning: $name (negative no-intron PCG) need to be checked due to query coverage per annotated PCG more than $long!\n";
										}
									}
								}elsif ($length_pcg < 60) {
									print $logfile "Warning: $name (negative no-intron PCG) has not been annotated due to short length!\n";
								}
							}elsif ((!defined $right_start_position) and (defined $left_star_position)){
								my $end=$end_right-($left_star_position * 3+3);

								my $right_star=$start_right-$right_star_position * 3;
								#my $seq_string=substr($sequence,($start2-60),120);
								my $seq_string=substr($sequence,($start2-60),(60+length($right_star-$start2)));
								my $length_seq_string=length $seq_string;
								my $rev_seq_string=reverse $seq_string;
								$rev_seq_string=~ tr/ACGTacgt/TGCAtgca/;
								my $aa_seq_string;
								for (my $i=0;$i<$length_seq_string;$i+=3){# delete stop codon
									my $codon=substr ($rev_seq_string,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string.=$hash_codon{$codon};
									}else{
										$aa_seq_string.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}

								my $start;
								for (my $i=0;$i<$xx;$i+=1) {
									my $match_aa_ref_pcg=substr ($aa_ref_pcg,$i,$yy);
									my $marks=0;
									for (my $j=$yy;$j<(($length_seq_string)/3);$j+=1) {
										my $repeat=substr ($aa_seq_string,-$j,$yy);
										if (($repeat eq $match_aa_ref_pcg) and ($marks == 0)) {
											$start=$start2-60+$j*3+$i*3;
											$marks++;
										}
									}
									if (defined $start){
										last;
									}
								}
								if (!defined $start) {
									$start=$start2;
								}

								if ($start >= ($start_right-$right_star_position*3)) {
									print $logfile "Warning: $name (negative no-intron PCG) has not been annotated due to internal stop codon!\n";
								}elsif ($start < ($start_right-$right_star_position*3)) {
									my $length_pcg=$start-$end;
									if ($length_pcg >= 60) {
										if ($redundancy eq "N") {
											if ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) < $short)) {
												print $logfile "Warning: $name (negative no-intron PCG) has not been annotated due to query coverage per annotated PCG less than $short!\n";
											}elsif ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) > $long)) {
												print $logfile "Warning: $name (negative no-intron PCG) has not been annotated due to query coverage per annotated PCG more than $long!\n";
											}elsif ((defined $length_ref_pcg) and (($length_pcg/$length_ref_pcg) > $short) and (($length_pcg/$length_ref_pcg) < $long)) {
												print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."complement(".$end."..".$start.")\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa\""."\n";
												$gene_number_seq{$name}++;
												print $logfile "Warning: $name (negative no-intron PCG) has non-canonical start codon!\n";
											}
										}elsif ($redundancy eq "Y") {
											print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."complement(".$end."..".$start.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (negative no-intron PCG) has non-canonical start codon!\n";

											if ((defined $length_ref_pcg) and ($length_pcg/$length_ref_pcg < $short)) {
												print $logfile "Warning: $name (negative no-intron PCG) need to be checked due to query coverage per annotated PCG less than $short!\n";
											}elsif ((defined $length_ref_pcg) and ($length_pcg/$length_ref_pcg > $long)) {
												print $logfile "Warning: $name (negative no-intron PCG) need to be checked due to query coverage per annotated PCG more than $long!\n";
											}
										}
									}elsif ($length_pcg < 60) {
										print $logfile "Warning: $name (negative no-intron PCG) has not been annotated due to short length!\n";
									}
								}
							}
						}
					}
				}
			}
		}

		############################################################
		## one-intron or two-intron protein-coding-gene and tRNA gene
		############################################################
		if((keys %{$hash{$name}} >= 2) and (keys %{$hash{$name}} <= 5)){
			my (@position1,@position2,@position3,@position4);
			my $coding;
			foreach my $gene_coding_CDS (sort {$b cmp $a} keys %{$hash{$name}}){
				if ($gene_coding_CDS=~ /^((.+)_gene)$/){
					$coding=$2;
					foreach (keys %{$hash{$name}{$1}}){
						push @position1,$_."\t".$hash{$name}{$gene_coding_CDS}{$_};
					}
				}elsif($gene_coding_CDS=~ /^(((.+)-1)_coding)$/ or $gene_coding_CDS=~ /^(((.+)-1)_coding_aa)$/){
					foreach (keys %{$hash{$name}{$1}}){
						push @position2,$_."\t".$hash{$name}{$gene_coding_CDS}{$_};
					}
				}elsif($gene_coding_CDS=~ /^(((.+)-2)_coding)$/ or $gene_coding_CDS=~ /^(((.+)-2)_coding_aa)$/){
					foreach (keys %{$hash{$name}{$1}}){
						push @position3,$_."\t".$hash{$name}{$gene_coding_CDS}{$_};
					}
				}elsif($gene_coding_CDS=~ /^(((.+)-3)_coding)$/ or $gene_coding_CDS=~ /^(((.+)-3)_coding_aa)$/){
					foreach (keys %{$hash{$name}{$1}}){
						push @position4,$_."\t".$hash{$name}{$gene_coding_CDS}{$_};
					}
				}
			}

			while (@position1 and (@position2 or @position3)){# _gene
				my $coding1=$coding."-1_coding";
				my $coding2=$coding."-2_coding";
				my $coding3=$coding."-3_coding";

				my $position1=shift @position1;
				my ($start1,$end1,$number1)=(split /\t/,$position1)[0,1,2];# _gene
				my $position2=shift @position2;
				my ($start2,$end2,$number2);# -1_coding or -1_coding_aa
				my $position3=shift @position3;
				my ($start3,$end3,$number3);# -2_coding or -2_coding_aa
				my $position4=shift @position4;
				my ($start4,$end4,$number4)=(split /\t/,$position4)[0,1,2] if (defined $position4 ne "");# -3_coding or -3_coding_aa

				my (@coding1,@coding2);
				if (defined $position2 ne ""){
	    			($start2,$end2,$number2)=(split /\t/,$position2)[0,1,2];
				}elsif(defined $position2 eq ""){# very short exon (-1_coding)
					foreach (keys %hash_exon){
						if ($_ eq $coding1){
							foreach (sort {$hash_exon{$_}{$b} <=> $hash_exon{$_}{$a}} keys %{$hash_exon{$_}}){
								push @coding1,$_;
							}
							$start2=$start1;
							if ($start1 < $end1){
								$end2=$start1+$coding1[0]-1;
							}elsif($start1 > $end1){
								$end2=$start1-$coding1[0]+1;
							}
						}
					}
				}
				if (defined $position3 ne ""){
					($start3,$end3,$number3)=(split /\t/,$position3)[0,1,2];
				}elsif(defined $position3 eq ""){# very short exon (-2_coding)
					foreach (keys %hash_exon){
						if ($_ eq $coding2){
							foreach (sort {$hash_exon{$_}{$b} <=> $hash_exon{$_}{$a}} keys %{$hash_exon{$_}}){
								push @coding2,$_;
							}
							$start3=$start1;
							if ($start1 < $end1){
								$end3=$start1+$coding2[0]-1;
							}elsif($start1 > $end1){
								$end3=$start1-$coding2[0]+1;
							}
						}
					}
				}

				############################################################
				## one-intron tRNA
				############################################################
				if($name=~ /(\D+)-(\D){3}/){
					my $length_ref_exon1=$hash_exon_length{"$name-1_coding"};
					my $length_ref_exon2=$hash_exon_length{"$name-2_coding"};
					my $sequence_ref_exon1=$hash_exon_sequence{"$name-1_coding"};
					my $sequence_ref_exon2=$hash_exon_sequence{"$name-2_coding"};
					$sequence_ref_exon1=uc $sequence_ref_exon1;
					$sequence_ref_exon2=uc $sequence_ref_exon2;
					my ($x,$y);
					if (defined $length_ref_exon1 >= 32) {
						$x=31;
					}elsif (defined $length_ref_exon1 < 32) {
						$x=$length_ref_exon1;
					}
					if (defined $length_ref_exon2 >= 32) {
						$y=31;
					}elsif (defined $length_ref_exon2 < 32) {
						$y=$length_ref_exon2;
					}
					my $z=9;

					if ((($start1 == $start2) and ($end1 == $end3)) and ((defined $length_ref_exon1) and (defined $length_ref_exon2) and (defined $sequence_ref_exon1) and (defined $sequence_ref_exon2))){# identical tRNA boundary for _gene and -1_coding, -2_coding
						if ($start1 < $end1){
							my $seq_string1=substr($sequence,($start1-31),60);
							my $seq_string2=substr($sequence,($end2-30),60);
							my $seq_string3=substr($sequence,($start3-31),60);
							my $seq_string4=substr($sequence,($end1-30),60);

							my ($start1_new,$end2_new,$start3_new,$end1_new);
							for (my $i=0;$i<$x;$i+=1) {
								my $match_ref_exon=substr ($sequence_ref_exon1,$i,$z);
								my $marks=0;
								for (my $j=0;$j<60;$j+=1) {
									my $repeat=substr ($seq_string1,$j,$z);
									if (($repeat eq $match_ref_exon) and ($marks == 0)) {
										$start1_new=$start1-30+$j-$i;
										$marks++;
									}
								}
								last if (defined $start1_new);
							}
							for (my $i=$z;$i<$x;$i+=1) {
								my $match_ref_exon=substr ($sequence_ref_exon1,-$i,$z);
								my $marks=0;
								for (my $j=0;$j<60;$j+=1) {
									my $repeat=substr ($seq_string2,$j,$z);
									if (($repeat eq $match_ref_exon) and ($marks == 0)) {
										$end2_new=$end2-30+($j+$z)+($i-$z);
										$marks++;
									}
								}
								last if (defined $end2_new);
							}
							for (my $i=0;$i<$y;$i+=1) {
								my $match_ref_exon=substr ($sequence_ref_exon2,$i,$z);
								my $marks=0;
								for (my $j=0;$j<60;$j+=1) {
									my $repeat=substr ($seq_string3,$j,$z);
									if (($repeat eq $match_ref_exon) and ($marks == 0)) {
										$start3_new=$start3-30+$j-$i;
										$marks++;
									}
								}
								last if (defined $start3_new);
							}
							for (my $i=$z;$i<$y;$i+=1) {
								my $match_ref_exon=substr ($sequence_ref_exon2,-$i,$z);
								my $marks=0;
								for (my $j=0;$j<60;$j+=1) {
									my $repeat=substr ($seq_string4,$j,$z);
									if (($repeat eq $match_ref_exon) and ($marks == 0)) {
										$end1_new=$end1-30+($j+$z)+($i-$z);
										$marks++;
									}
								}
								last if (defined $end1_new);
							}
							my $length_trna=$end1-$start1;

							if ((defined $start1_new) and (defined $end2_new) and (defined $start3_new) and (defined $end1_new)) {
								if (abs($end1_new-$start1_new) < 3000) {
									if ($end2_new < $end1_new) {
										print $out_annotation "     "."gene"."            ".$start1_new."..".$end1_new."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."tRNA"."            "."join(".$start1_new."..".$end2_new.",".$start3_new."..".$end1_new.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										$gene_number_seq{$name}++;
										my ($S1,$E1,$E2,$S3);
										my $string=substr($sequence,($start1_new-1),($end1_new-$start1_new+1));
										my $lengths=length $string;
										my $L1=$end2_new-$start1_new+1;
										my $L2=$end1_new-$start3_new+1;
										my $rev_coms=reverse $sequence;
										$rev_coms=~ tr/ACGTacgt/TGCAtgca/;
										for (my $i=0;$i<($length_cp-$lengths);$i+=1){
											my $repeat=substr ($sequence,$i,$lengths);
											if (($repeat eq $string) and (($i+1) ne $start1_new)) {
												$S1=$i+1;
												$E2=$i+$L1-1+1;
												$S3=$i+$lengths-$L2+1;
												$E1=$i+$lengths-1+1;
												print $out_annotation "     "."gene"."            ".$S1."..".$E1."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."tRNA"."            "."join(".$S1."..".$E2.",".$S3."..".$E1.")"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												$gene_number_seq{$name}++;
											}
										}
										for (my $i=0;$i<($length_cp-$lengths);$i+=1){
											my $repeat=substr ($rev_coms,$i,$lengths);
											if ($repeat eq $string) {
												$S1=$length_cp-$i-1+1;
												$E2=$length_cp-$i-$L1+1;
												$S3=$length_cp-$i-$lengths-1+$L2+1;
												$E1=$length_cp-$i-$lengths+1;
												print $out_annotation "     "."gene"."            "."complement(".$E1."..".$S1.")"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."tRNA"."            "."complement(join(".$E1."..".$S3.",".$E2."..".$S1."))"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												$gene_number_seq{$name}++;
											}
										}
									}elsif ($end2_new > $end1_new) {
										if ($length_trna >= 70) {
											print $out_annotation "     "."gene"."            ".$start1."..".$end1."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."tRNA"."            "."join(".$start1."..".($start1+37).",".($end1-37)."..".$end1.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (positive one-intron tRNA) must be checked due to non-identical boundary with reference!\n";
										}elsif ($length_trna < 70) {
											print $logfile "Warning: $name (positive one-intron tRNA) has not been annotated due to short length!\n";
										}
									}
								}elsif (abs($end1_new-$start1_new) >= 3000) {
									print $logfile "Warning: $name (positive one-intron tRNA) has not been annotated due to two far exons!\n";
								}
							}elsif ((!defined $start1_new) or (!defined $end2_new) or (!defined $start3_new) or (!defined $end1_new)) {
								if ($length_trna >= 70) {
									print $out_annotation "     "."gene"."            ".$start1."..".$end1."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."tRNA"."            "."join(".$start1."..".($start1+37).",".($end1-37)."..".$end1.")"."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									$gene_number_seq{$name}++;
									print $logfile "Warning: $name (positive one-intron tRNA) must be checked due to non-identical boundary with reference!\n";
								}elsif ($length_trna < 70) {
									print $logfile "Warning: $name (positive one-intron tRNA) has not been annotated due to short length!\n";
								}
							}
						}elsif($start1 > $end1){
							my $seq_string1=substr($sequence,($start1-30),60);
							my $seq_string2=substr($sequence,($end2-31),60);
							my $seq_string3=substr($sequence,($start3-30),60);
							my $seq_string4=substr($sequence,($end1-31),60);
							my $rev_seq_string1=reverse $seq_string1;
							$rev_seq_string1=~ tr/ACGTacgt/TGCAtgca/;
							my $rev_seq_string2=reverse $seq_string2;
							$rev_seq_string2=~ tr/ACGTacgt/TGCAtgca/;
							my $rev_seq_string3=reverse $seq_string3;
							$rev_seq_string3=~ tr/ACGTacgt/TGCAtgca/;
							my $rev_seq_string4=reverse $seq_string4;
							$rev_seq_string4=~ tr/ACGTacgt/TGCAtgca/;

							my ($start1_new,$end2_new,$start3_new,$end1_new);
							for (my $i=0;$i<$x;$i+=1) {
								my $match_ref_exon=substr ($sequence_ref_exon1,$i,$z);
								my $marks=0;
								for (my $j=0;$j<60;$j+=1) {
									my $repeat=substr ($rev_seq_string1,$j,$z);
									if (($repeat eq $match_ref_exon) and ($marks == 0)) {
										$start1_new=$start1+30-$j+$i;
										$marks++;
									}
								}
								last if (defined $start1_new);
							}
							for (my $i=$z;$i<$x;$i+=1) {
								my $match_ref_exon=substr ($sequence_ref_exon1,-$i,$z);
								my $marks=0;
								for (my $j=0;$j<60;$j+=1) {
									my $repeat=substr ($rev_seq_string2,$j,$z);
									if (($repeat eq $match_ref_exon) and ($marks == 0)) {
										$end2_new=$end2+30-($j+$z)-($i-$z);
										$marks++;
									}
								}
								last if (defined $end2_new);
							}
							for (my $i=0;$i<$y;$i+=1) {
								my $match_ref_exon=substr ($sequence_ref_exon2,$i,$z);
								my $marks=0;
								for (my $j=0;$j<60;$j+=1) {
									my $repeat=substr ($rev_seq_string3,$j,$z);
									if (($repeat eq $match_ref_exon) and ($marks == 0)) {
										$start3_new=$start3+30-$j+$i;
										$marks++;
									}
								}
								last if (defined $start3_new);
							}
							for (my $i=$z;$i<$y;$i+=1) {
								my $match_ref_exon=substr ($sequence_ref_exon2,-$i,$z);
								my $marks=0;
								for (my $j=0;$j<60;$j+=1) {
									my $repeat=substr ($rev_seq_string4,$j,$z);
									if (($repeat eq $match_ref_exon) and ($marks == 0)) {
										$end1_new=$end1+30-($j+$z)-($i-$z);
										$marks++;
									}
								}
								last if (defined $end1_new);
							}
							my $length_trna=$start1-$end1;

							if ((defined $start1_new) and (defined $end2_new) and (defined $start3_new) and (defined $end1_new)) {
								if (abs($start1_new-$end1_new) < 3000) {
									if ($end2_new > $end1_new) {
										print $out_annotation "     "."gene"."            "."complement(".$end1_new."..".$start1_new.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."tRNA"."            "."complement(join(".$end1_new."..".$start3_new.",".$end2_new."..".$start1_new."))"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										$gene_number_seq{$name}++;
										my ($S1,$E1,$E2,$S3);
										my $string=substr($sequence,($end1_new-1),($start1_new-$end1_new+1));
										my $lengths=length $string;
										my $L1=$start1_new-$end2_new+1;
										my $L2=$start3_new-$end1_new+1;
										my $rev_coms=reverse $sequence;
										$rev_coms=~ tr/ACGTacgt/TGCAtgca/;
										for (my $i=0;$i<($length_cp-$lengths);$i+=1){
											my $repeat=substr ($sequence,$i,$lengths);
											if (($repeat eq $string) and (($i+1) ne $end1_new)) {
												$E1=$i+1;
												$S3=$i+$L2-1+1;
												$E2=$i+$lengths-$L1+1;
												$S1=$i+$lengths-1+1;
												print $out_annotation "     "."gene"."            "."complement(".$E1."..".$S1.")"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."tRNA"."            "."complement(join(".$E1."..".$S3.",".$E2."..".$S1."))"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												$gene_number_seq{$name}++;
											}
										}
										for (my $i=0;$i<($length_cp-$lengths);$i+=1){
											my $repeat=substr ($rev_coms,$i,$lengths);
											if ($repeat eq $string) {
												$E1=$length_cp-$i-1+1;
												$S3=$length_cp-$i-$L2+1;
												$E2=$length_cp-$i-$lengths-1+$L1+1;
												$S1=$length_cp-$i-$lengths+1;
												print $out_annotation "     "."gene"."            ".$S1."..".$E1."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."tRNA"."            "."join(".$S1."..".$E2.",".$S3."..".$E1.")"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												$gene_number_seq{$name}++;
											}
										}
									}elsif ($end2_new < $end1_new) {
										if ($length_trna >= 70) {
											print $out_annotation "     "."gene"."            "."complement(".$end1."..".$start1.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."tRNA"."            "."complement(join(".$end1."..".($end1+37).",".($start1-37)."..".$start1."))"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (negative one-intron tRNA) must be checked due to non-identical boundary with reference!\n";
										}elsif ($length_trna < 70) {
											print $logfile "Warning: $name (negative one-intron tRNA) has not been annotated due to short length!\n";
										}
									}
								}elsif (abs($start1_new-$end1_new) >= 3000) {
									print $logfile "Warning: $name (negative one-intron tRNA) has not been annotated due to two far exons!\n";
								}
							}elsif ((!defined $start1_new) or (!defined $end2_new) or (!defined $start3_new) or (!defined $end1_new)) {
								if ($length_trna >= 70) {
									print $out_annotation "     "."gene"."            "."complement(".$end1."..".$start1.")"."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."tRNA"."            "."complement(join(".$end1."..".($end1+37).",".($start1-37)."..".$start1."))"."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									$gene_number_seq{$name}++;
									print $logfile "Warning: $name (negative one-intron tRNA) must be checked due to non-identical boundary with reference!\n";
								}elsif ($length_trna < 70) {
									print $logfile "Warning: $name (negative one-intron tRNA) has not been annotated due to short length!\n";
								}
							}
						}
					}elsif((($start1 != $start2) or ($end1 != $end3)) and ((defined $length_ref_exon1) and (defined $length_ref_exon2) and (defined $sequence_ref_exon1) and (defined $sequence_ref_exon2))){# non-identical tRNA boundary for _gene and -1_coding, -2_coding
						if ($start1 < $end1){
							my ($seq_string1,$seq_string2);
							if ($length_ref_exon1 > 30) {
								$seq_string1=substr($sequence,($start2-31),60);
								$seq_string2=substr($sequence,($end2-30),60);
							}elsif ($length_ref_exon1 < 30) {
								$seq_string1=substr($sequence,($start1-31),60);
								$seq_string2=substr($sequence,($start1+10),60);
							}
							my $seq_string3=substr($sequence,($start3-31),60);
							my $seq_string4=substr($sequence,($end3-30),60);

							my ($start2_new,$end2_new,$start3_new,$end3_new);
							for (my $i=0;$i<$x;$i+=1) {
								my $match_ref_exon=substr ($sequence_ref_exon1,$i,$z);
								my $marks=0;
								for (my $j=0;$j<60;$j+=1) {
									my $repeat=substr ($seq_string1,$j,$z);
									if (($repeat eq $match_ref_exon) and ($marks == 0)) {
										$start2_new=$start2-30+$j-$i if ($length_ref_exon1 > 30);
										$start2_new=$start1-30+$j-$i if ($length_ref_exon1 < 30);
										$marks++;
									}
								}
								last if (defined $start2_new);
							}
							for (my $i=$z;$i<$x;$i+=1) {
								my $match_ref_exon=substr ($sequence_ref_exon1,-$i,$z);
								my $marks=0;
								for (my $j=0;$j<60;$j+=1) {
									my $repeat=substr ($seq_string2,$j,$z);
									if (($repeat eq $match_ref_exon) and ($marks == 0)) {
										$end2_new=$end2-30+($j+$z)+($i-$z) if ($length_ref_exon1 > 30);
										$end2_new=$start1+10+($j+$z)+($i-$z) if ($length_ref_exon1 < 30);
										$marks++;
									}
								}
								last if (defined $end2_new);
							}
							for (my $i=0;$i<$y;$i+=1) {
								my $match_ref_exon=substr ($sequence_ref_exon2,$i,$z);
								my $marks=0;
								for (my $j=0;$j<60;$j+=1) {
									my $repeat=substr ($seq_string3,$j,$z);
									if (($repeat eq $match_ref_exon) and ($marks == 0)) {
										$start3_new=$start3-30+$j-$i;
										$marks++;
									}
								}
								last if (defined $start3_new);
							}
							for (my $i=$z;$i<$y;$i+=1) {
								my $match_ref_exon=substr ($sequence_ref_exon2,-$i,$z);
								my $marks=0;
								for (my $j=0;$j<60;$j+=1) {
									my $repeat=substr ($seq_string4,$j,$z);
									if (($repeat eq $match_ref_exon) and ($marks == 0)) {
										$end3_new=$end3-30+($j+$z)+($i-$z);
										$marks++;
									}
								}
								last if (defined $end3_new);
							}
							my $length_trna=$end1-$start1;

							if ((defined $start2_new) and (defined $end2_new) and (defined $start3_new) and (defined $end3_new)) {
								if (abs($end3_new-$start2_new) < 3000) {
									if ($end2_new < $end3_new) {
										print $out_annotation "     "."gene"."            ".$start2_new."..".$end3_new."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."tRNA"."            "."join(".$start2_new."..".$end2_new.",".$start3_new."..".$end3_new.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										$gene_number_seq{$name}++;
										my ($S1,$E1,$E2,$S3);
										my $string=substr($sequence,($start2_new-1),($end3_new-$start2_new+1));
										my $lengths=length $string;
										my $L1=$end2_new-$start2_new+1;
										my $L2=$end3_new-$start3_new+1;
										my $rev_coms=reverse $sequence;
										$rev_coms=~ tr/ACGTacgt/TGCAtgca/;
										for (my $i=0;$i<($length_cp-$lengths);$i+=1){
											my $repeat=substr ($sequence,$i,$lengths);
											if (($repeat eq $string) and (($i+1) ne $start2_new)) {
												$S1=$i+1;
												$E2=$i+$L1-1+1;
												$S3=$i+$lengths-$L2+1;
												$E1=$i+$lengths-1+1;
												print $out_annotation "     "."gene"."            ".$S1."..".$E1."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."tRNA"."            "."join(".$S1."..".$E2.",".$S3."..".$E1.")"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												$gene_number_seq{$name}++;
											}
										}
										for (my $i=0;$i<($length_cp-$lengths);$i+=1){
											my $repeat=substr ($rev_coms,$i,$lengths);
											if ($repeat eq $string) {
												$S1=$length_cp-$i-1+1;
												$E2=$length_cp-$i-$L1+1;
												$S3=$length_cp-$i-$lengths-1+$L2+1;
												$E1=$length_cp-$i-$lengths+1;
												print $out_annotation "     "."gene"."            "."complement(".$E1."..".$S1.")"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."tRNA"."            "."complement(join(".$E1."..".$S3.",".$E2."..".$S1."))"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												$gene_number_seq{$name}++;
											}
										}
									}elsif ($end2_new > $end3_new) {
										if ($length_trna >= 70) {
											print $out_annotation "     "."gene"."            ".$start1."..".$end1."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."tRNA"."            "."join(".$start1."..".($start1+37).",".($end1-37)."..".$end1.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (positive one-intron tRNA) must be checked due to non-identical boundary with reference!\n";
										}elsif ($length_trna < 70) {
											print $logfile "Warning: $name (positive one-intron tRNA) has not been annotated due to short length!\n";
										}
									}
								}elsif (abs($end3_new-$start2_new) >= 3000) {
									print $logfile "Warning: $name (positive one-intron tRNA) has not been annotated due to two far exons!\n";
								}
							}elsif ((!defined $start2_new) or (!defined $end2_new) or (!defined $start3_new) or (!defined $end3_new)) {
								if ($length_trna >= 70) {
									print $out_annotation "     "."gene"."            ".$start1."..".$end1."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."tRNA"."            "."join(".$start1."..".($start1+37).",".($end1-37)."..".$end1.")"."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									$gene_number_seq{$name}++;
									print $logfile "Warning: $name (positive one-intron tRNA) must be checked due to non-identical boundary with reference!\n";
								}elsif ($length_trna < 70) {
									print $logfile "Warning: $name (positive one-intron tRNA) has not been annotated due to short length!\n";
								}
							}
						}elsif($start1 > $end1){
							my ($seq_string1,$seq_string2);
							if ($length_ref_exon1 > 30) {
								$seq_string1=substr($sequence,($start2-30),60);
								$seq_string2=substr($sequence,($end2-31),60);
							}elsif ($length_ref_exon1 < 30) {
								$seq_string1=substr($sequence,($start1-30),60);
								$seq_string2=substr($sequence,($start1-71),60);
							}
							my $seq_string3=substr($sequence,($start3-30),60);
							my $seq_string4=substr($sequence,($end3-31),60);
							my $rev_seq_string1=reverse $seq_string1;
							$rev_seq_string1=~ tr/ACGTacgt/TGCAtgca/;
							my $rev_seq_string2=reverse $seq_string2;
							$rev_seq_string2=~ tr/ACGTacgt/TGCAtgca/;
							my $rev_seq_string3=reverse $seq_string3;
							$rev_seq_string3=~ tr/ACGTacgt/TGCAtgca/;
							my $rev_seq_string4=reverse $seq_string4;
							$rev_seq_string4=~ tr/ACGTacgt/TGCAtgca/;

							my ($start2_new,$end2_new,$start3_new,$end3_new);
							for (my $i=0;$i<$x;$i+=1) {
								my $match_ref_exon=substr ($sequence_ref_exon1,$i,$z);
								my $marks=0;
								for (my $j=0;$j<60;$j+=1) {
									my $repeat=substr ($rev_seq_string1,$j,$z);
									if (($repeat eq $match_ref_exon) and ($marks == 0)) {
										$start2_new=$start2+30-$j+$i if ($length_ref_exon1 > 30);
										$start2_new=$start1+30-$j+$i if ($length_ref_exon1 < 30);
										$marks++;
									}
								}
								last if (defined $start2_new);
							}
							for (my $i=$z;$i<$x;$i+=1) {
								my $match_ref_exon=substr ($sequence_ref_exon1,-$i,$z);
								my $marks=0;
								for (my $j=0;$j<60;$j+=1) {
									my $repeat=substr ($rev_seq_string2,$j,$z);
									if (($repeat eq $match_ref_exon) and ($marks == 0)) {
										$end2_new=$end2+30-($j+$z)-($i-$z) if ($length_ref_exon1 > 30);
										$end2_new=$start1-10-($j+$z)-($i-$z) if ($length_ref_exon1 < 30);
										$marks++;
									}
								}
								last if (defined $end2_new);
							}
							for (my $i=0;$i<$y;$i+=1) {
								my $match_ref_exon=substr ($sequence_ref_exon2,$i,$z);
								my $marks=0;
								for (my $j=0;$j<60;$j+=1) {
									my $repeat=substr ($rev_seq_string3,$j,$z);
									if (($repeat eq $match_ref_exon) and ($marks == 0)) {
										$start3_new=$start3+30-$j+$i;
										$marks++;
									}
								}
								last if (defined $start3_new);
							}
							for (my $i=$z;$i<$y;$i+=1) {
								my $match_ref_exon=substr ($sequence_ref_exon2,-$i,$z);
								my $marks=0;
								for (my $j=0;$j<60;$j+=1) {
									my $repeat=substr ($rev_seq_string4,$j,$z);
									if (($repeat eq $match_ref_exon) and ($marks == 0)) {
										$end3_new=$end3+30-($j+$z)-($i-$z);
										$marks++;
									}
								}
								last if (defined $end3_new);
							}
							my $length_trna=$start1-$end1;

							if ((defined $start2_new) and (defined $end2_new) and (defined $start3_new) and (defined $end3_new)) {
								if (abs($start2_new-$end3_new) < 3000) {
									if ($end2_new > $end3_new) {
										print $out_annotation "     "."gene"."            "."complement(".$end3_new."..".$start2_new.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."tRNA"."            "."complement(join(".$end3_new."..".$start3_new.",".$end2_new."..".$start2_new."))"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										$gene_number_seq{$name}++;
										my ($S1,$E1,$E2,$S3);
										my $string=substr($sequence,($end3_new-1),($start2_new-$end3_new+1));
										my $lengths=length $string;
										my $L1=$start2_new-$end2_new+1;
										my $L2=$start3_new-$end3_new+1;
										my $rev_coms=reverse $sequence;
										$rev_coms=~ tr/ACGTacgt/TGCAtgca/;
										for (my $i=0;$i<($length_cp-$lengths);$i+=1){
											my $repeat=substr ($sequence,$i,$lengths);
											if (($repeat eq $string) and (($i+1) ne $end3_new)) {
												$E1=$i+1;
												$S3=$i+$L2-1+1;
												$E2=$i+$lengths-$L1+1;
												$S1=$i+$lengths-1+1;
												print $out_annotation "     "."gene"."            "."complement(".$E1."..".$S1.")"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."tRNA"."            "."complement(join(".$E1."..".$S3.",".$E2."..".$S1."))"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												$gene_number_seq{$name}++;
											}
										}
										for (my $i=0;$i<($length_cp-$lengths);$i+=1){
											my $repeat=substr ($rev_coms,$i,$lengths);
											if ($repeat eq $string) {
												$E1=$length_cp-$i-1+1;
												$S3=$length_cp-$i-$L2+1;
												$E2=$length_cp-$i-$lengths-1+$L1+1;
												$S1=$length_cp-$i-$lengths+1;
												print $out_annotation "     "."gene"."            ".$S1."..".$E1."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."tRNA"."            "."join(".$S1."..".$E2.",".$S3."..".$E1.")"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												$gene_number_seq{$name}++;
											}
										}
									}elsif ($end2_new < $end3_new) {
										if ($length_trna >= 70) {
											print $out_annotation "     "."gene"."            "."complement(".$end1."..".$start1.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."tRNA"."            "."complement(join(".$end1."..".($end1+37).",".($start1-37)."..".$start1."))"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (negative one-intron tRNA) must be checked due to non-identical boundary with reference!\n";
										}elsif ($length_trna < 70) {
											print $logfile "Warning: $name (negative one-intron tRNA) has not been annotated due to short length!\n";
										}
									}
								}elsif (abs($start2_new-$end3_new) >= 3000) {
									print $logfile "Warning: $name (negative one-intron tRNA) has not been annotated due to two far exons!\n";
								}
							}elsif ((!defined $start2_new) or (!defined $end2_new) or (!defined $start3_new) or (!defined $end3_new)) {
								if ($length_trna >= 70) {
									print $out_annotation "     "."gene"."            "."complement(".$end1."..".$start1.")"."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."tRNA"."            "."complement(join(".$end1."..".($end1+37).",".($start1-37)."..".$start1."))"."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									$gene_number_seq{$name}++;
									print $logfile "Warning: $name (negative one-intron tRNA) must be checked due to non-identical boundary with reference!\n";
								}elsif ($length_trna < 70) {
									print $logfile "Warning: $name (negative one-intron tRNA) has not been annotated due to short length!\n";
								}
							}
						}
					}

				############################################################
				## one-intron or two-intron protein-coding-gene
				############################################################
				}elsif ($name!~ /(\D+)-(\D){3}/){
					my $length_ref_exon1=$hash_exon_length{"$name-1_coding"};
					my $length_ref_exon2=$hash_exon_length{"$name-2_coding"};
					my $length_ref_exon3=$hash_exon_length{"$name-3_coding"};
					my $sequence_ref_exon1=$hash_exon_sequence{"$name-1_coding"};
					my $sequence_ref_exon2=$hash_exon_sequence{"$name-2_coding"};
					my $sequence_ref_exon3=$hash_exon_sequence{"$name-3_coding"};

					my ($a,$b,$c);
					if ((defined $length_ref_exon1) and ($length_ref_exon1/3 >= 22)) {
						$a=21;
					}elsif ((defined $length_ref_exon1) and ($length_ref_exon1/3 < 22)) {
						$a=int($length_ref_exon1/3);
					}
					if ((defined $length_ref_exon2) and ($length_ref_exon2/3 >= 22)) {
						$b=21;
					}elsif ((defined $length_ref_exon2) and ($length_ref_exon2/3 < 22)) {
						$b=int($length_ref_exon2/3);
					}
					if ((defined $length_ref_exon3) and ($length_ref_exon3/3 >= 22)) {
						$c=21;
					}elsif ((defined $length_ref_exon3) and ($length_ref_exon3/3 < 22)) {
						$c=int($length_ref_exon3/3);
					}
					my $match=4;

					my $aa_ref_exon1;
					my $aa_ref_exon2;
					my $aa_ref_exon3;
					my @aa_ref_exon1;
					my @aa_ref_exon2;
					my @aa_ref_exon3;
					if (((defined $length_ref_exon1) and (defined $length_ref_exon2) or (defined $length_ref_exon3)) and ((defined $sequence_ref_exon1) and (defined $sequence_ref_exon2) or (defined $sequence_ref_exon3))) {
						if ($length_ref_exon1 % 3==0) {
							for (my $i=0;$i<$length_ref_exon1;$i+=3){
								my $codon=substr ($sequence_ref_exon1,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa_ref_exon1.=$hash_codon{$codon};
								}else{
									$aa_ref_exon1.="X";
									my $j=$i+1;
									#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							@aa_ref_exon1=split //,$aa_ref_exon1;

							if (!defined $start4) {
								for (my $i=0;$i<($length_ref_exon2-3);$i+=3){# delete stop codon
									my $codon=substr ($sequence_ref_exon2,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_ref_exon2.=$hash_codon{$codon};
									}else{
										$aa_ref_exon2.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								@aa_ref_exon2=split //,$aa_ref_exon2;
							}elsif (defined $start4) {
								for (my $i=0;$i<$length_ref_exon2;$i+=3){
									my $codon=substr ($sequence_ref_exon2,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_ref_exon2.=$hash_codon{$codon};
									}else{
										$aa_ref_exon2.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								@aa_ref_exon2=split //,$aa_ref_exon2;

								for (my $i=0;$i<($length_ref_exon3-3);$i+=3){# delete stop codon
									my $codon=substr ($sequence_ref_exon3,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_ref_exon3.=$hash_codon{$codon};
									}else{
										$aa_ref_exon3.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								@aa_ref_exon3=split //,$aa_ref_exon3;
							}
						}elsif ($length_ref_exon1 % 3==1) {
							for (my $i=0;$i<($length_ref_exon1-1);$i+=3){# delete stop codon
								my $codon=substr ($sequence_ref_exon1,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa_ref_exon1.=$hash_codon{$codon};
								}else{
									$aa_ref_exon1.="X";
									my $j=$i+1;
									#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							@aa_ref_exon1=split //,$aa_ref_exon1;

							if (!defined $start4) {
								for (my $i=2;$i<($length_ref_exon2-3);$i+=3){# delete stop codon
									my $codon=substr ($sequence_ref_exon2,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_ref_exon2.=$hash_codon{$codon};
									}else{
										$aa_ref_exon2.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								@aa_ref_exon2=split //,$aa_ref_exon2;
							}elsif (defined $start4) {
								if ($length_ref_exon2 % 3==0) {
									for (my $i=0;$i<$length_ref_exon2;$i+=3){
										my $codon=substr ($sequence_ref_exon2,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_ref_exon2.=$hash_codon{$codon};
										}else{
											$aa_ref_exon2.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
									@aa_ref_exon2=split //,$aa_ref_exon2;

									for (my $i=2;$i<($length_ref_exon3-3);$i+=3){# delete stop codon
										my $codon=substr ($sequence_ref_exon3,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_ref_exon3.=$hash_codon{$codon};
										}else{
											$aa_ref_exon3.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
									@aa_ref_exon3=split //,$aa_ref_exon3;
								}elsif ($length_ref_exon2 % 3==1) {
									for (my $i=1;$i<$length_ref_exon2;$i+=3){
										my $codon=substr ($sequence_ref_exon2,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_ref_exon2.=$hash_codon{$codon};
										}else{
											$aa_ref_exon2.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
									@aa_ref_exon2=split //,$aa_ref_exon2;

									for (my $i=1;$i<($length_ref_exon3-3);$i+=3){# delete stop codon
										my $codon=substr ($sequence_ref_exon3,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_ref_exon3.=$hash_codon{$codon};
										}else{
											$aa_ref_exon3.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
									@aa_ref_exon3=split //,$aa_ref_exon3;
								}elsif ($length_ref_exon2 % 3==2) {
									for (my $i=2;$i<$length_ref_exon2;$i+=3){
										my $codon=substr ($sequence_ref_exon2,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_ref_exon2.=$hash_codon{$codon};
										}else{
											$aa_ref_exon2.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
									@aa_ref_exon2=split //,$aa_ref_exon2;

									for (my $i=0;$i<($length_ref_exon3-3);$i+=3){# delete stop codon
										my $codon=substr ($sequence_ref_exon3,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_ref_exon3.=$hash_codon{$codon};
										}else{
											$aa_ref_exon3.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
									@aa_ref_exon3=split //,$aa_ref_exon3;
								}
							}
						}elsif ($length_ref_exon1 % 3==2) {
							for (my $i=0;$i<($length_ref_exon1-2);$i+=3){# delete stop codon
								my $codon=substr ($sequence_ref_exon1,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa_ref_exon1.=$hash_codon{$codon};
								}else{
									$aa_ref_exon1.="X";
									my $j=$i+1;
									#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							@aa_ref_exon1=split //,$aa_ref_exon1;

							if ($length_ref_exon2 % 3==0) {
								for (my $i=0;$i<$length_ref_exon2;$i+=3){
									my $codon=substr ($sequence_ref_exon2,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_ref_exon2.=$hash_codon{$codon};
									}else{
										$aa_ref_exon2.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								@aa_ref_exon2=split //,$aa_ref_exon2;

								for (my $i=1;$i<($length_ref_exon3-3);$i+=3){# delete stop codon
									my $codon=substr ($sequence_ref_exon3,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_ref_exon3.=$hash_codon{$codon};
									}else{
										$aa_ref_exon3.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								@aa_ref_exon3=split //,$aa_ref_exon3;
							}elsif ($length_ref_exon2 % 3==1) {
								for (my $i=1;$i<$length_ref_exon2;$i+=3){
									my $codon=substr ($sequence_ref_exon2,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_ref_exon2.=$hash_codon{$codon};
									}else{
										$aa_ref_exon2.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								@aa_ref_exon2=split //,$aa_ref_exon2;

								for (my $i=0;$i<($length_ref_exon3-3);$i+=3){# delete stop codon
									my $codon=substr ($sequence_ref_exon3,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_ref_exon3.=$hash_codon{$codon};
									}else{
										$aa_ref_exon3.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								@aa_ref_exon3=split //,$aa_ref_exon3;
							}
						}
					}

					############################################################
					## one-intron protein-coding-gene
					############################################################
					if ((($start1 == $start2) and ($end1 == $end3)) and ((defined $length_ref_exon1) and (defined $length_ref_exon2) and (defined $sequence_ref_exon1) and (defined $sequence_ref_exon2))){# identical PCG boundary for _gene and -1_coding_aa, -2_coding_aa
						if ($start1 < $end1){# positive
							my $str=substr($sequence,($start1-1),($end1-$start1+1));
							my $length=length $str;
							my $forward_start_codon=substr($str,0,3);
							my $forward_stop_codon=substr($str,($length-3),3);
							my ($str1,$str2);
							$str1=substr($sequence,($start1-1),($end2-$start1+1)) if ($start1 < $end2);
							$str1=substr($sequence,($end2-1),($start1-$end2+1)) if ($start1 > $end2);
							$str2=substr($sequence,($start3-1),($end1-$start3+1)) if ($start3 < $end1);
							$str2=substr($sequence,($end1-1),($start3-$end1+1)) if ($start3 > $end1);
							my $str3=$str1.$str2;
							my $length_exon=length $str3;
							my $length_exon2=length $str2;

							my $aa_exon;
							for (my $i=0;$i<($length_exon-3);$i+=3){# delete stop codon
								my $codon=substr ($str3,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa_exon.=$hash_codon{$codon};
								}else{
									$aa_exon.="X";
									my $j=$i+1;
									#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa_exon=split //,$aa_exon;
							if (($name ne "rpl16") and ($name ne "petB") and ($name ne "petD")) {
								if (($length_exon % 3==0) and (!grep {$_=~ /\*/} @aa_exon) and (($forward_start_codon eq "ATG") or ($forward_start_codon eq "GTG") or ($name eq "rps12+2")) and (($forward_stop_codon eq "TAA") or ($forward_stop_codon eq "TAG") or ($forward_stop_codon eq "TGA"))){# standard start and stop codon for _gene
									my $seq_string0=substr($sequence,($start1-31),60);
									my $seq_string1=substr($sequence,($end2-60),((int(($start3-$end2)/3))*3+120));
									my $seq_string2=substr($sequence,($start3-(int(($start3-$end2)/3))*3-61),((int(($start3-$end2)/3))*3+120));
									my $length_seq_string0=length $seq_string0;
									my $length_seq_string1=length $seq_string1;
									my $length_seq_string2=length $seq_string2;
									my $aa_seq_string0;
									my $aa_seq_string1;
									my $aa_seq_string2;
									for (my $i=0;$i<$length_seq_string0;$i+=3){# delete stop codon
										my $codon=substr ($seq_string0,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_seq_string0.=$hash_codon{$codon};
										}else{
											$aa_seq_string0.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
									for (my $i=0;$i<$length_seq_string1;$i+=3){# delete stop codon
										my $codon=substr ($seq_string1,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_seq_string1.=$hash_codon{$codon};
										}else{
											$aa_seq_string1.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
									for (my $i=0;$i<$length_seq_string2;$i+=3){# delete stop codon
										my $codon=substr ($seq_string2,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_seq_string2.=$hash_codon{$codon};
										}else{
											$aa_seq_string2.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}

									my ($str1,$str2,$str,$length,$start1_new,$end2_new,$start3_new);
									my ($string_intron,$intron_length);
									my $mark=0;
									my $ticks=0;
									if ($length_ref_exon1 % 3==0) {
										for (my $i=$match;$i<$a;$i+=1) {
											my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
												my $repeat=substr ($aa_seq_string1,$j,$match);
												if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
													$end2_new=$end2-60+($j+$match)*3+($i-$match)*3;
													$marks++;
												}
											}
											if (defined $end2_new){
												last;
											}
										}
										if (!defined $end2_new) {
											$end2_new=$end2;
											$ticks=1;
										}
										for (my $i=0;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
												my $repeat=substr ($aa_seq_string2,-$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$start3_new=$start3+60-$j*3-$i*3;
													$marks++;
												}
											}
											if (defined $start3_new){
												last;
											}
										}
										if (!defined $start3_new) {
											$start3_new=$start3;
											$ticks=1;
										}
										$intron_length=$start3_new-$end2_new-1;
										if ($intron_length==0) {
											$mark=1;
										}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
											$string_intron=substr ($sequence,$end2_new,$intron_length);
											$intron_length=length $string_intron;
											$mark=2;
										}
									}elsif ($length_ref_exon1 % 3==1) {
										for (my $i=$match;$i<$a;$i+=1) {
											my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
												my $repeat=substr ($aa_seq_string1,$j,$match);
												if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
													$end2_new=$end2-60+($j+$match)*3+($i-$match)*3+1;
													$marks++;
												}
											}
											if (defined $end2_new){
												last;
											}
										}
										if (!defined $end2_new) {
											$end2_new=$end2+1;
											$ticks=1;
										}
										for (my $i=0;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
												my $repeat=substr ($aa_seq_string2,-$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$start3_new=$start3+60-$j*3-$i*3-2;
													$marks++;
												}
											}
											if (defined $start3_new){
												last;
											}
										}
										if (!defined $start3_new) {
											$start3_new=$start3-2;
											$ticks=1;
										}
										$intron_length=$start3_new-$end2_new-1;
										if ($intron_length==0) {
											$mark=1;
										}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
											$string_intron=substr ($sequence,($end2_new-1),$intron_length);
											$intron_length=length $string_intron;
											$mark=2;
										}
									}elsif ($length_ref_exon1 % 3==2) {
										for (my $i=$match;$i<$a;$i+=1) {
											my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
												my $repeat=substr ($aa_seq_string1,$j,$match);
												if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
													$end2_new=$end2-60+($j+$match)*3+($i-$match)*3+2;
													$marks++;
												}
											}
											if (defined $end2_new){
												last;
											}
										}
										if (!defined $end2_new) {
											$end2_new=$end2+2;
											$ticks=1;
										}
										for (my $i=0;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
												my $repeat=substr ($aa_seq_string2,-$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$start3_new=$start3+60-$j*3-$i*3-1;
													$marks++;
												}
											}
											if (defined $start3_new){
												last;
											}
										}
										if (!defined $start3_new) {
											$start3_new=$start3-1;
											$ticks=1;
										}
										$intron_length=$start3_new-$end2_new-1;
										if ($intron_length==0) {
											$mark=1;
										}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
											$string_intron=substr ($sequence,$end2_new+1,$intron_length);
											$intron_length=length $string_intron;
											$mark=2;
										}
									}

									if ($name eq "rps12+2") {
										for (my $i=0;$i<7;$i+=1) {
											my $match_aa_ref_exon1=substr ($aa_ref_exon1,$i,3);
											my $marks=0;
											for (my $j=0;$j<20;$j+=1) {
												my $repeat=substr ($aa_seq_string0,$j,3);
												if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
													$start1_new=$start1-30+$j*3-$i*3;
													$marks++;
												}
											}
										}
									}elsif ($name ne "rps12+2") {
										$start1_new=$start1;
									}
									$str1=substr($sequence,($start1_new-1),($end2_new-($start1_new-1)));
									$str2=substr($sequence,($start3_new-1),($end1-($start3_new-1)));
									$str=$str1.$str2;
									$length=length $str;

									my $aa;
									for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
										my $codon=substr ($str,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa.=$hash_codon{$codon};
										}else{
											$aa.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}

									my $intron_aa;
									if ($mark==2) {
										for (my $i=0;$i<$intron_length;$i+=3){
											my $codon=substr ($string_intron,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$intron_aa.=$hash_codon{$codon};
											}else{
												$intron_aa.="X";
												my $j=$i+1;
												#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}
									}

									if ((abs($end3-$start2) < 4000) and ($start2 < $end2) and ($start3 < $end3) and ($start2 < $end3)) {
										if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/))) {# spliced exon,not 3x or have stop codon
											print $out_annotation "     "."gene"."            ".$start1_new."..".$end1."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."join(".$start1_new."..".$end2_new.",".$start3_new."..".$end1.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
											$gene_number_seq{$name}++;
										}elsif (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/))) {# joined exon,no intron or no stop codon
											print $out_annotation "     "."gene"."            ".$start1_new."..".$end1."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             ".$start1_new."..".$end1."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (positive one-intron PCG) lost intron!\n";
										}

										my ($S1,$E1,$E2,$S3);
										my $string=substr($sequence,($start1_new-1),($end1-$start1_new+1));
										my $lengths=length $string;
										my $L1=$end2_new-$start1_new+1;
										my $L2=$end1-$start3_new+1;
										for (my $i=0;$i<($length_cp-$lengths);$i+=1){
											my $repeat=substr ($sequence,$i,$lengths);
											if (($repeat eq $string) and (($i+1) ne $start1_new)) {
												$S1=$i+1;
												$E2=$i+$L1-1+1;
												$S3=$i+$lengths-$L2+1;
												$E1=$i+$lengths-1+1;
												print $out_annotation "     "."gene"."            ".$S1."..".$E1."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."join(".$S1."..".$E2.",".$S3."..".$E1.")"."\n" if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)));
												print $out_annotation "     "."CDS"."             ".$S1."..".$E1."\n" if (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/)));
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
												$gene_number_seq{$name}++;
											}
										}
										for (my $i=0;$i<($length_cp-$lengths);$i+=1){
											my $repeat=substr ($rev_coms,$i,$lengths);
											if ($repeat eq $string) {
												$S1=$length_cp-$i-1+1;
												$E2=$length_cp-$i-$L1+1;
												$S3=$length_cp-$i-$lengths-1+$L2+1;
												$E1=$length_cp-$i-$lengths+1;
												print $out_annotation "     "."gene"."            "."complement(".$E1."..".$S1.")"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."complement(join(".$E1."..".$S3.",".$E2."..".$S1."))"."\n" if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)));
												print $out_annotation "     "."CDS"."             "."complement(".$E1."..".$S1.")"."\n" if (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/)));
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
												$gene_number_seq{$name}++;
											}
										}
										if (($ticks == 1) and (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)))) {
											print $logfile "Warning: $name (positive one-intron PCG) intron-exon boundary need to be checked!\n";
										}
									}elsif ((abs($end3-$start2) < 4000) and (($start2 > $end2) or ($start3 > $end3) or ($start2 > $end3))) {
										print $logfile "Warning: $name (positive one-intron PCG) has not been annotated due to two different direction exons!\n";
									}elsif (abs($end3-$start2) >= 4000) {
										print $logfile "Warning: $name (positive one-intron PCG) has not been annotated due to two far exons!\n";
									}
								}elsif((grep {$_=~ /\*/} @aa_exon) or (($forward_start_codon ne "ATG") and ($forward_start_codon ne "GTG")) or (($forward_stop_codon ne "TAA") and ($forward_stop_codon ne "TAG") and ($forward_stop_codon ne "TGA"))){# non-standard start or stop codon for _gene
									my $start_left;
									if ($start1>=9000) {
										$start_left=$start1-9000;
									}elsif ($start1<9000){
										$start_left=$start1-(int($start1/3)-1)*3;
									}
									my $start_right=$start1+60;
									my $length_exon1=($end2-$start1+1);
									my $end_right;
									if (($length_cp-$start3)>=9000) {
										$end_right=$start3+9000;
									}elsif (($length_cp-$start3)<9000){
										$end_right=$start3+((int(($length_cp-$start3)/3)-1)*3);
									}
									my $str1=substr($sequence,($start_left-1),($start_right-$start_left));
									my $str2=substr($sequence,($start3-1),($end_right-$start3));
									my $length1=length $str1;
									my $length2=length $str2;
									my $aa1;#left range for start codon
									for (my $i=0;$i<($length1-3);$i+=3){# delete stop codon
										my $codon=substr ($str1,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa1.=$hash_codon{$codon};
										}else{
											$aa1.="X";
											my $m=$i+1;
											#print "Bad codon $codon in position $m of gene $name in species $contig!\n";
										}
									}
									my @aa1=split //,$aa1;
									my (@star1,@start1);
									foreach (0..$#aa1){
										if (($aa1[$_]=~ /\*/) and ($_ < 3000)){
											push @star1,$_;
										}
										if (($aa1[$_]=~ /M/) and ($_ <= 3060)){
											push @start1,$_;
										}
									}
									my $left_star_position=pop @star1;
									until ($left_star_position<(3000+$length_exon1/3)){
										$left_star_position=pop @star1;
									}
									my @start2;
									foreach (@start1){
										push @start2,$_ if ((defined $left_star_position) and ($_ > $left_star_position))
									}
									my $left_start_position=$start2[0];
									#my $left_start_position=pop @start2;
									my $aa2;#right range for stop codon
									for (my $k=0;$k<($length2-3);$k+=3){# delete stop codon
										my $codon=substr ($str2,$k,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa2.=$hash_codon{$codon};
										}else{
											$aa2.="X";
											my $n=$k+1;
											#print "Bad codon $codon in position $n of gene $name in species $contig!\n";
										}
									}
									my @aa2=split //,$aa2;
									my @star2;
									foreach (0..$#aa2){
										if ($aa2[$_]=~ /\*/){
											push @star2,$_;
										}
									}
									my $right_star_position=shift @star2;

									my $seq_string0=substr($sequence,($start1-31),60);
									my $seq_string1=substr($sequence,($end2-60),((int(($start3-$end2)/3))*3+120));
									my $seq_string2=substr($sequence,($start3-(int(($start3-$end2)/3))*3-61),((int(($start3-$end2)/3))*3+120));
									my $length_seq_string0=length $seq_string0;
									my $length_seq_string1=length $seq_string1;
									my $length_seq_string2=length $seq_string2;
									my $aa_seq_string0;
									my $aa_seq_string1;
									my $aa_seq_string2;
									for (my $i=0;$i<$length_seq_string0;$i+=3){# delete stop codon
										my $codon=substr ($seq_string0,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_seq_string0.=$hash_codon{$codon};
										}else{
											$aa_seq_string0.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
									for (my $i=0;$i<$length_seq_string1;$i+=3){# delete stop codon
										my $codon=substr ($seq_string1,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_seq_string1.=$hash_codon{$codon};
										}else{
											$aa_seq_string1.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
									for (my $i=0;$i<$length_seq_string2;$i+=3){# delete stop codon
										my $codon=substr ($seq_string2,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_seq_string2.=$hash_codon{$codon};
										}else{
											$aa_seq_string2.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}

									if ((defined $left_start_position) and (defined $right_star_position)){
										my $start=$start_left+$left_start_position * 3;
										my $end=($start3-1)+($right_star_position * 3+3);

										my ($str1,$str2,$str,$length,$start_new,$end2_new,$start3_new);
										my ($string_intron,$intron_length);
										my $mark=0;
										my $ticks=0;
										if ($length_ref_exon1 % 3==0) {
											for (my $i=$match;$i<$a;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string1,$j,$match);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$end2_new=$end2-60+($j+$match)*3+($i-$match)*3;
														$marks++;
													}
												}
												if (defined $end2_new){
													last;
												}
											}
											if (!defined $end2_new) {
												$end2_new=$end2;
												$ticks=1;
											}
											for (my $i=0;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string2,-$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$start3_new=$start3+60-$j*3-$i*3;
														$marks++;
													}
												}
												if (defined $start3_new){
													last;
												}
											}
											if (!defined $start3_new) {
												$start3_new=$start3;
												$ticks=1;
											}
											$intron_length=$start3_new-$end2_new-1;
											if ($intron_length==0) {
												$mark=1;
											}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
												$string_intron=substr ($sequence,$end2_new,$intron_length);
												$intron_length=length $string_intron;
												$mark=2;
											}
										}elsif ($length_ref_exon1 % 3==1) {
											for (my $i=$match;$i<$a;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string1,$j,$match);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$end2_new=$end2-60+($j+$match)*3+($i-$match)*3+1;
														$marks++;
													}
												}
												if (defined $end2_new){
													last;
												}
											}
											if (!defined $end2_new) {
												$end2_new=$end2+1;
												$ticks=1;
											}
											for (my $i=0;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string2,-$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$start3_new=$start3+60-$j*3-$i*3-2;
														$marks++;
													}
												}
												if (defined $start3_new){
													last;
												}
											}
											if (!defined $start3_new) {
												$start3_new=$start3-2;
												$ticks=1;
											}
											$intron_length=$start3_new-$end2_new-1;
											if ($intron_length==0) {
												$mark=1;
											}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
												$string_intron=substr ($sequence,($end2_new-1),$intron_length);
												$intron_length=length $string_intron;
												$mark=2;
											}
										}elsif ($length_ref_exon1 % 3==2) {
											for (my $i=$match;$i<$a;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string1,$j,$match);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$end2_new=$end2-60+($j+$match)*3+($i-$match)*3+2;
														$marks++;
													}
												}
												if (defined $end2_new){
													last;
												}
											}
											if (!defined $end2_new) {
												$end2_new=$end2+2;
												$ticks=1;
											}
											for (my $i=0;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string2,-$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$start3_new=$start3+60-$j*3-$i*3-1;
														$marks++;
													}
												}
												if (defined $start3_new){
													last;
												}
											}
											if (!defined $start3_new) {
												$start3_new=$start3-1;
												$ticks=1;
											}
											$intron_length=$start3_new-$end2_new-1;
											if ($intron_length==0) {
												$mark=1;
											}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
												$string_intron=substr ($sequence,$end2_new+1,$intron_length);
												$intron_length=length $string_intron;
												$mark=2;
											}
										}


										if ($name eq "rps12+2") {
											for (my $i=0;$i<7;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,$i,3);
												my $marks=0;
												for (my $j=0;$j<20;$j+=1) {
													my $repeat=substr ($aa_seq_string0,$j,3);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$start_new=$start1-30+$j*3-$i*3;
														$marks++;
													}
												}
											}
										}elsif ($name ne "rps12+2") {
											$start_new=$start;
										}
										$str1=substr($sequence,($start_new-1),($end2_new-($start_new-1)));
										$str2=substr($sequence,($start3_new-1),($end-($start3_new-1)));
										$str=$str1.$str2;
										$length=length $str;

										my $aa;
										for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
											my $codon=substr ($str,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$aa.=$hash_codon{$codon};
											}else{
												$aa.="X";
												my $j=$i+1;
												#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}

										my $intron_aa;
										if ($mark==2) {
											for (my $i=0;$i<$intron_length;$i+=3){
												my $codon=substr ($string_intron,$i,3);
												$codon=uc $codon;
												if (exists $hash_codon{$codon}){
													$intron_aa.=$hash_codon{$codon};
												}else{
													$intron_aa.="X";
													my $j=$i+1;
													#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
												}
											}
										}

										if ((abs($end3-$start2) < 4000) and ($start2 < $end2) and ($start3 < $end3) and ($start2 < $end3)) {
											if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/))) {# spliced exon,not 3x or have stop codon
												print $out_annotation "     "."gene"."            ".$start_new."..".$end."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."join(".$start_new."..".$end2_new.",".$start3_new."..".$end.")"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa\""."\n";
												$gene_number_seq{$name}++;
											}elsif (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/))) {# joined exon,no intron or no stop codon
												print $out_annotation "     "."gene"."            ".$start_new."..".$end."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             ".$start_new."..".$end."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa\""."\n";
												$gene_number_seq{$name}++;
												print $logfile "Warning: $name (positive one-intron PCG) lost intron!\n";
											}

											my ($S1,$E1,$E2,$S3);
											my $string=substr($sequence,($start_new-1),($end-$start_new+1));
											my $lengths=length $string;
											my $L1=$end2_new-$start_new+1;
											my $L2=$end-$start3_new+1;
											for (my $i=0;$i<($length_cp-$lengths);$i+=1){
												my $repeat=substr ($sequence,$i,$lengths);
												if (($repeat eq $string) and (($i+1) ne $start_new)) {
													$S1=$i+1;
													$E2=$i+$L1-1+1;
													$S3=$i+$lengths-$L2+1;
													$E1=$i+$lengths-1+1;
													print $out_annotation "     "."gene"."            ".$S1."..".$E1."\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "     "."CDS"."             "."join(".$S1."..".$E2.",".$S3."..".$E1.")"."\n" if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)));
													print $out_annotation "     "."CDS"."             ".$S1."..".$E1."\n" if (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/)));
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "                     "."/codon_start=1"."\n";
													print $out_annotation "                     "."/transl_table=11"."\n";
													print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
													#print $out_annotation "                     "."/translation=\"$aa\""."\n";
													$gene_number_seq{$name}++;
												}
											}
											for (my $i=0;$i<($length_cp-$lengths);$i+=1){
												my $repeat=substr ($rev_coms,$i,$lengths);
												if ($repeat eq $string) {
													$S1=$length_cp-$i-1+1;
													$E2=$length_cp-$i-$L1+1;
													$S3=$length_cp-$i-$lengths-1+$L2+1;
													$E1=$length_cp-$i-$lengths+1;
													print $out_annotation "     "."gene"."            "."complement(".$E1."..".$S1.")"."\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "     "."CDS"."             "."complement(join(".$E1."..".$S3.",".$E2."..".$S1."))"."\n" if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)));
													print $out_annotation "     "."CDS"."             "."complement(".$E1."..".$S1.")"."\n" if (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/)));
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "                     "."/codon_start=1"."\n";
													print $out_annotation "                     "."/transl_table=11"."\n";
													print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
													#print $out_annotation "                     "."/translation=\"$aa\""."\n";
													$gene_number_seq{$name}++;
												}
											}
											if (($ticks == 1) and (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)))) {
												print $logfile "Warning: $name (positive one-intron PCG) intron-exon boundary need to be checked!\n";
											}
										}elsif ((abs($end3-$start2) < 4000) and (($start2 > $end2) or ($start3 > $end3) or ($start2 > $end3))) {
											print $logfile "Warning: $name (positive one-intron PCG) has not been annotated due to two different direction exons!\n";
										}elsif (abs($end3-$start2) >= 4000) {
											print $logfile "Warning: $name (positive one-intron PCG) has not been annotated due to two far exons!\n";
										}
									}elsif ((!defined $left_start_position) and (defined $right_star_position)){
										my $end=($start3-1)+($right_star_position * 3+3);

										my ($str1,$str2,$str,$length,$start1_new,$end2_new,$start3_new);
										my ($string_intron,$intron_length);
										my $mark=0;
										my $ticks=0;
										if ($length_ref_exon1 % 3==0) {
											for (my $i=$match;$i<$a;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string1,$j,$match);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$end2_new=$end2-60+($j+$match)*3+($i-$match)*3;
														$marks++;
													}
												}
												if (defined $end2_new){
													last;
												}
											}
											if (!defined $end2_new) {
												$end2_new=$end2;
												$ticks=1;
											}
											for (my $i=0;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string2,-$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$start3_new=$start3+60-$j*3-$i*3;
														$marks++;
													}
												}
												if (defined $start3_new){
													last;
												}
											}
											if (!defined $start3_new) {
												$start3_new=$start3;
												$ticks=1;
											}
											$intron_length=$start3_new-$end2_new-1;
											if ($intron_length==0) {
												$mark=1;
											}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
												$string_intron=substr ($sequence,$end2_new,$intron_length);
												$intron_length=length $string_intron;
												$mark=2;
											}
										}elsif ($length_ref_exon1 % 3==1) {
											for (my $i=$match;$i<$a;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string1,$j,$match);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$end2_new=$end2-60+($j+$match)*3+($i-$match)*3+1;
														$marks++;
													}
												}
												if (defined $end2_new){
													last;
												}
											}
											if (!defined $end2_new) {
												$end2_new=$end2+1;
												$ticks=1;
											}
											for (my $i=0;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string2,-$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$start3_new=$start3+60-$j*3-$i*3-2;
														$marks++;
													}
												}
												if (defined $start3_new){
													last;
												}
											}
											if (!defined $start3_new) {
												$start3_new=$start3-2;
												$ticks=1;
											}
											$intron_length=$start3_new-$end2_new-1;
											if ($intron_length==0) {
												$mark=1;
											}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
												$string_intron=substr ($sequence,($end2_new-1),$intron_length);
												$intron_length=length $string_intron;
												$mark=2;
											}
										}elsif ($length_ref_exon1 % 3==2) {
											for (my $i=$match;$i<$a;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string1,$j,$match);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$end2_new=$end2-60+($j+$match)*3+($i-$match)*3+2;
														$marks++;
													}
												}
												if (defined $end2_new){
													last;
												}
											}
											if (!defined $end2_new) {
												$end2_new=$end2+2;
												$ticks=1;
											}
											for (my $i=0;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string2,-$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$start3_new=$start3+60-$j*3-$i*3-1;
														$marks++;
													}
												}
												if (defined $start3_new){
													last;
												}
											}
											if (!defined $start3_new) {
												$start3_new=$start3-1;
												$ticks=1;
											}
											$intron_length=$start3_new-$end2_new-1;
											if ($intron_length==0) {
												$mark=1;
											}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
												$string_intron=substr ($sequence,$end2_new+1,$intron_length);
												$intron_length=length $string_intron;
												$mark=2;
											}
										}

										if ($name eq "rps12+2") {
											for (my $i=0;$i<7;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,$i,3);
												my $marks=0;
												for (my $j=0;$j<20;$j+=1) {
													my $repeat=substr ($aa_seq_string0,$j,3);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$start1_new=$start1-30+$j*3-$i*3;
														$marks++;
													}
												}
											}
										}elsif ($name ne "rps12+2") {
											$start1_new=$start1;
										}
										$str1=substr($sequence,($start1_new-1),($end2_new-($start1_new-1)));
										$str2=substr($sequence,($start3_new-1),($end-($start3_new-1)));
										$str=$str1.$str2;
										$length=length $str;

										my $aa;
										for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
											my $codon=substr ($str,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$aa.=$hash_codon{$codon};
											}else{
												$aa.="X";
												my $j=$i+1;
												#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}

										my $intron_aa;
										if ($mark==2) {
											for (my $i=0;$i<$intron_length;$i+=3){
												my $codon=substr ($string_intron,$i,3);
												$codon=uc $codon;
												if (exists $hash_codon{$codon}){
													$intron_aa.=$hash_codon{$codon};
												}else{
													$intron_aa.="X";
													my $j=$i+1;
													#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
												}
											}
										}

										if ((abs($end3-$start2) < 4000) and ($start2 < $end2) and ($start3 < $end3) and ($start2 < $end3)) {
											if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/))) {# spliced exon,not 3x or have stop codon
												print $out_annotation "     "."gene"."            ".$start1_new."..".$end."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."join(".$start1_new."..".$end2_new.",".$start3_new."..".$end.")"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa\""."\n";
												$gene_number_seq{$name}++;
											}elsif (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/))) {# joined exon,no intron or no stop codon
												print $out_annotation "     "."gene"."            ".$start1_new."..".$end."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             ".$start1_new."..".$end."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa\""."\n";
												$gene_number_seq{$name}++;
												print $logfile "Warning: $name (positive one-intron PCG) lost intron!\n";
											}

											my ($S1,$E1,$E2,$S3);
											my $string=substr($sequence,($start1_new-1),($end-$start1_new+1));
											my $lengths=length $string;
											my $L1=$end2_new-$start1_new+1;
											my $L2=$end-$start3_new+1;
											for (my $i=0;$i<($length_cp-$lengths);$i+=1){
												my $repeat=substr ($sequence,$i,$lengths);
												if (($repeat eq $string) and (($i+1) ne $start1_new)) {
													$S1=$i+1;
													$E2=$i+$L1-1+1;
													$S3=$i+$lengths-$L2+1;
													$E1=$i+$lengths-1+1;
													print $out_annotation "     "."gene"."            ".$S1."..".$E1."\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "     "."CDS"."             "."join(".$S1."..".$E2.",".$S3."..".$E1.")"."\n" if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)));
													print $out_annotation "     "."CDS"."             ".$S1."..".$E1."\n" if (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/)));
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "                     "."/codon_start=1"."\n";
													print $out_annotation "                     "."/transl_table=11"."\n";
													print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
													#print $out_annotation "                     "."/translation=\"$aa\""."\n";
													$gene_number_seq{$name}++;
												}
											}
											for (my $i=0;$i<($length_cp-$lengths);$i+=1){
												my $repeat=substr ($rev_coms,$i,$lengths);
												if ($repeat eq $string) {
													$S1=$length_cp-$i-1+1;
													$E2=$length_cp-$i-$L1+1;
													$S3=$length_cp-$i-$lengths-1+$L2+1;
													$E1=$length_cp-$i-$lengths+1;
													print $out_annotation "     "."gene"."            "."complement(".$E1."..".$S1.")"."\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "     "."CDS"."             "."complement(join(".$E1."..".$S3.",".$E2."..".$S1."))"."\n" if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)));
													print $out_annotation "     "."CDS"."             "."complement(".$E1."..".$S1.")"."\n" if (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/)));
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "                     "."/codon_start=1"."\n";
													print $out_annotation "                     "."/transl_table=11"."\n";
													print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
													#print $out_annotation "                     "."/translation=\"$aa\""."\n";
													$gene_number_seq{$name}++;
												}
											}
											if (($ticks == 1) and (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)))) {
												print $logfile "Warning: $name (positive one-intron PCG) has non-canonical start codon and intron-exon boundary need to be checked!\n" if ($name ne "rps12+2");
											}elsif ($ticks == 0) {
												print $logfile "Warning: $name (positive one-intron PCG) has non-canonical start codon!\n" if ($name ne "rps12+2");
											}
										}elsif ((abs($end3-$start2) < 4000) and (($start2 > $end2) or ($start3 > $end3) or ($start2 > $end3))) {
											print $logfile "Warning: $name (positive one-intron PCG) has not been annotated due to two different direction exons!\n";
										}elsif (abs($end3-$start2) >= 4000) {
											print $logfile "Warning: $name (positive one-intron PCG) has not been annotated due to two far exons!\n";
										}
									}else{
										#print $out_annotation "     "."gene"."            ".$start1."..".$end1."\n";
										#print $out_annotation "                     "."/gene=\"$name\""."\n";
										#print $out_annotation "     "."CDS"."             "."join(".$start1."..".$end2.",".$start3."..".$end1.")"."\n";
										#print $out_annotation "                     "."/gene=\"$name\""."\n";
										#print $out_annotation "                     "."/codon_start=1"."\n";
										#print $out_annotation "                     "."/transl_table=11"."\n";
										#print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										##print $out_annotation "                     "."/translation=\"$aa\""."\n";
										#$gene_number_seq{$name}++;
										#print $logfile "Warning: $name (positive one-intron PCG) maybe pseudogene!\n";
										print $logfile "Warning: $name (positive one-intron PCG) has not been annotated!\n";
									}
								}
							}elsif (($name eq "rpl16") or ($name eq "petB") or ($name eq "petD")) {
								if ($length_exon2 < 150) {
									print $logfile "Warning: $name (positive one-intron PCG) has not been annotated due to short length!\n";
								}elsif ($length_exon2 >= 150) {
									if (($length_ref_exon1 % 3==0) and ($length_exon % 3==0) and (!grep {$_=~ /\*/} @aa_exon) and (($forward_start_codon eq "ATG") or ($forward_start_codon eq "GTG")) and (($forward_stop_codon eq "TAA") or ($forward_stop_codon eq "TAG") or ($forward_stop_codon eq "TGA"))){# standard start and stop codon for _gene
										print $out_annotation "     "."gene"."            ".$start1."..".$end1."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."join(".$start1."..".$end2.",".$start3."..".$end1.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
										$gene_number_seq{$name}++;
									}elsif(($length_ref_exon1 % 3!=0) or (grep {$_=~ /\*/} @aa_exon) or (($forward_start_codon ne "ATG") and ($forward_start_codon ne "GTG")) or (($forward_stop_codon ne "TAA") and ($forward_stop_codon ne "TAG") and ($forward_stop_codon ne "TGA"))){# non-standard start or stop codon for _gene
										my $start_left;
										if ($start1>=9000) {
											$start_left=$start1-9000;
										}elsif ($start1<9000){
											$start_left=$start1-(int($start1/3)-1)*3;
										}
										my $start_right=$start1+60;
										my $length_exon1=($end2-$start1+1);
										my $end_right;
										if (($length_cp-$start3)>=9000) {
											$end_right=$start3+9000;
										}elsif (($length_cp-$start3)<9000){
											$end_right=$start3+((int(($length_cp-$start3)/3)-1)*3);
										}
										my $str1=substr($sequence,($start_left-1),($start_right-$start_left));
										my $str2=substr($sequence,($start3-1),($end_right-$start3));
										my $length1=length $str1;
										my $length2=length $str2;
										my $coding_1_length=$hash_exon_length{"$name-1_coding"};
										my $coding_1_sequence=$hash_exon_sequence{"$name-1_coding"};
										my $left_start_position;
										for (my $i=0;$i<($length1-3);$i+=1){
											my $exon=substr ($str1,$i,$coding_1_length);
											$exon=lc $exon;
											if ($exon eq $coding_1_sequence) {
												$left_start_position=$i;
											}
										}

										my $aa2;#right range for stop codon
										for (my $k=0;$k<($length2-3);$k+=3){# delete stop codon
											my $codon=substr ($str2,$k,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$aa2.=$hash_codon{$codon};
											}else{
												$aa2.="X";
												my $n=$k+1;
												#print "Bad codon $codon in position $n of gene $name in species $contig!\n";
											}
										}
										my @aa2=split //,$aa2;
										my @star2;
										foreach (0..$#aa2){
											if ($aa2[$_]=~ /\*/){
												push @star2,$_;
											}
										}
										my $right_star_position=shift @star2;
										if ((defined $left_start_position) and (defined $right_star_position)){
											my $start=$start_left+$left_start_position;
											my $end=($start3-1)+($right_star_position * 3+3);

											my ($str1,$str2,$str,$length,$end2_new,$start3_new);
											$end2_new=$start+$coding_1_length-1;
											if ($length_ref_exon1 % 3==0) {
												$start3_new=$start3;
											}elsif ($length_ref_exon1 % 3==1) {
												$start3_new=$start3-2;
											}elsif ($length_ref_exon1 % 3==2) {
												$start3_new=$start3-1;
											}
											$str1=substr($sequence,($start-1),($end2_new-($start-1)));
											$str2=substr($sequence,($start3_new-1),($end-($start3_new-1)));
											$str=$str1.$str2;
											$length=length $str;

											my $aa;
											for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
												my $codon=substr ($str,$i,3);
												$codon=uc $codon;
												if (exists $hash_codon{$codon}){
													$aa.=$hash_codon{$codon};
												}else{
													$aa.="X";
													my $j=$i+1;
													#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
												}
											}
											print $out_annotation "     "."gene"."            ".$start."..".$end."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."join(".$start."..".$end2_new.",".$start3_new."..".$end.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa\""."\n";
											$gene_number_seq{$name}++;
										}elsif ((!defined $left_start_position) and (defined $right_star_position)){
											my $end=($start3-1)+($right_star_position * 3+3);

											my ($str1,$str2,$str,$length,$end2_new,$start3_new);
											$end2_new=$start1+$coding_1_length-1;
											if ($length_ref_exon1 % 3==0) {
												$start3_new=$start3;
											}elsif ($length_ref_exon1 % 3==1) {
												$start3_new=$start3-2;
											}elsif ($length_ref_exon1 % 3==2) {
												$start3_new=$start3-1;
											}
											$str1=substr($sequence,($start1-1),($end2_new-($start1-1)));
											$str2=substr($sequence,($start3_new-1),($end-($start3_new-1)));
											$str=$str1.$str2;
											$length=length $str;

											my $aa;
											#for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
											for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
												my $codon=substr ($str,$i,3);
												$codon=uc $codon;
												if (exists $hash_codon{$codon}){
													$aa.=$hash_codon{$codon};
												}else{
													$aa.="X";
													my $j=$i+1;
													#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
												}
											}
											#print $out_annotation "     "."gene"."            ".$start1."..".$end."\n";
											print $out_annotation "     "."gene"."            ".$start3."..".$end."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											#print $out_annotation "     "."CDS"."             "."join(".$start1."..".$end2_new.",".$start3_new."..".$end.")"."\n";
											print $out_annotation "     "."CDS"."             ".$start3."..".$end."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (positive one-intron PCG) must be checked due to exon 1 not found!\n";
										}
									}
								}
							}
						}elsif($start1 > $end1){# negative
							my $str=substr($sequence,($end1-1),($start1-$end1+1));
							my $rev_com=reverse $str;
							$rev_com=~ tr/ACGTacgt/TGCAtgca/;
							my $length=length $str;
							my $reverse_start_codon=substr($rev_com,0,3);
							my $reverse_stop_codon=substr($rev_com,($length-3),3);
							my ($str1,$str2);
							$str1=substr($sequence,($end2-1),($start1-$end2+1)) if ($end2 < $start1);
							$str1=substr($sequence,($start1-1),($end2-$start1+1)) if ($end2 > $start1);
							$str2=substr($sequence,($end1-1),($start3-$end1+1)) if ($end1 < $start3);
							$str2=substr($sequence,($start3-1),($end1-$start3+1)) if ($end1 > $start3);
							my $rev_com1=reverse $str1;
							$rev_com1=~ tr/ACGTacgt/TGCAtgca/;
							my $rev_com2=reverse $str2;
							$rev_com2=~ tr/ACGTacgt/TGCAtgca/;
							my $str3=$rev_com1.$rev_com2;
							my $length_exon=length $str3;
							my $length_exon2=length $str2;
							my $aa_exon;
							for (my $i=0;$i<($length_exon-3);$i+=3){# delete stop codon
								my $codon=substr ($str3,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa_exon.=$hash_codon{$codon};
								}else{
									$aa_exon.="X";
									my $j=$i+1;
									#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa_exon=split //,$aa_exon;
							if (($name ne "rpl16") and ($name ne "petB") and ($name ne "petD")) {
								if (($length_exon % 3==0) and (!grep {$_=~ /\*/} @aa_exon) and (($reverse_start_codon eq "ATG") or ($reverse_start_codon eq "GTG") or ($name eq "rps12+2")) and (($reverse_stop_codon eq "TAA") or ($reverse_stop_codon eq "TAG") or ($reverse_stop_codon eq "TGA"))){# standard start and stop codon for _gene
									my $seq_string0=substr($sequence,($start1-30),60);
									my $seq_string1=substr($sequence,($end2-(int(($end2-$start3)/3))*3-61),((int(($end2-$start3)/3))*3+120));
									my $seq_string2=substr($sequence,($start3-60),((int(($end2-$start3)/3))*3+120));
									my $length_seq_string0=length $seq_string0;
									my $length_seq_string1=length $seq_string1;
									my $length_seq_string2=length $seq_string2;
									my $rev_seq_string0=reverse $seq_string0;
									$rev_seq_string0=~ tr/ACGTacgt/TGCAtgca/;
									my $rev_seq_string1=reverse $seq_string1;
									$rev_seq_string1=~ tr/ACGTacgt/TGCAtgca/;
									my $rev_seq_string2=reverse $seq_string2;
									$rev_seq_string2=~ tr/ACGTacgt/TGCAtgca/;
									my $aa_seq_string0;
									my $aa_seq_string1;
									my $aa_seq_string2;
									for (my $i=0;$i<$length_seq_string0;$i+=3){# delete stop codon
										my $codon=substr ($rev_seq_string0,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_seq_string0.=$hash_codon{$codon};
										}else{
											$aa_seq_string0.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
									for (my $i=0;$i<$length_seq_string1;$i+=3){# delete stop codon
										my $codon=substr ($rev_seq_string1,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_seq_string1.=$hash_codon{$codon};
										}else{
											$aa_seq_string1.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
									for (my $i=0;$i<$length_seq_string2;$i+=3){# delete stop codon
										my $codon=substr ($rev_seq_string2,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_seq_string2.=$hash_codon{$codon};
										}else{
											$aa_seq_string2.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}

									my ($str1,$str2,$str,$length,$rev_com,$start1_new,$end2_new,$start3_new);
									my ($string_intron,$intron_length,$string_intron_rev_com);
									my $mark=0;
									my $ticks=0;
									if ($length_ref_exon1 % 3==0) {
										for (my $i=$match;$i<$a;$i+=1) {
											my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
												my $repeat=substr ($aa_seq_string1,$j,$match);
												if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
													$end2_new=$end2+60-($j+$match)*3-($i-$match)*3;
													$marks++;
												}
											}
											if (defined $end2_new){
												last;
											}
										}
										if (!defined $end2_new) {
											$end2_new=$end2;
											$ticks=1;
										}
										for (my $i=0;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
												my $repeat=substr ($aa_seq_string2,-$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$start3_new=$start3-60+$j*3+$i*3;
													$marks++;
												}
											}
											if (defined $start3_new){
												last;
											}
										}
										if (!defined $start3_new) {
											$start3_new=$start3;
											$ticks=1;
										}
										$intron_length=$end2_new-$start3_new-1;
										if ($intron_length==0) {
											$mark=1;
										}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
											$string_intron=substr ($sequence,$start3_new,$intron_length);
											$intron_length=length $string_intron;
											$string_intron_rev_com=reverse $string_intron;
											$string_intron_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark=2;
										}
									}elsif ($length_ref_exon1 % 3==1) {
										for (my $i=3;$i<$a;$i+=1) {
											my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
												my $repeat=substr ($aa_seq_string1,$j,$match);
												if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
													$end2_new=$end2+60-($j+$match)*3-($i-$match)*3-1;
													$marks++;
												}
											}
											last if (defined $end2_new);
										}
										if (!defined $end2_new) {
											$end2_new=$end2-1;
											$ticks=1;
										}
										for (my $i=0;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
												my $repeat=substr ($aa_seq_string2,-$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$start3_new=$start3-60+$j*3+$i*3+2;
													$marks++;
												}
											}
											if (defined $start3_new){
												last;
											}
										}
										if (!defined $start3_new) {
											$start3_new=$start3+2;
											$ticks=1;
										}
										$intron_length=$end2_new-$start3_new-1;
										if ($intron_length==0) {
											$mark=1;
										}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
											$string_intron=substr ($sequence,($start3_new+1),$intron_length);
											$intron_length=length $string_intron;
											$string_intron_rev_com=reverse $string_intron;
											$string_intron_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark=2;
										}
									}elsif ($length_ref_exon1 % 3==2) {
										for (my $i=3;$i<$a;$i+=1) {
											my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
												my $repeat=substr ($aa_seq_string1,$j,$match);
												if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
													$end2_new=$end2+60-($j+$match)*3-($i-$match)*3-2;
													$marks++;
												}
											}
											last if (defined $end2_new);
										}
										if (!defined $end2_new) {
											$end2_new=$end2-2;
											$ticks=1;
										}
										for (my $i=0;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
												my $repeat=substr ($aa_seq_string2,-$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$start3_new=$start3-60+$j*3+$i*3+1;
													$marks++;
												}
											}
											if (defined $start3_new){
												last;
											}
										}
										if (!defined $start3_new) {
											$start3_new=$start3+1;
											$ticks=1;
										}
										$intron_length=$end2_new-$start3_new-1;
										if ($intron_length==0) {
											$mark=1;
										}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
											$string_intron=substr ($sequence,($start3_new-1),$intron_length);
											$intron_length=length $string_intron;
											$string_intron_rev_com=reverse $string_intron;
											$string_intron_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark=2;
										}
									}

									if ($name eq "rps12+2") {
										for (my $i=0;$i<7;$i+=1) {
											my $match_aa_ref_exon1=substr ($aa_ref_exon1,$i,3);
											my $marks=0;
											for (my $j=0;$j<20;$j+=1) {
												my $repeat=substr ($aa_seq_string0,$j,3);
												if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
													$start1_new=$start1+30-$j*3+$i*3;
													$marks++;
												}
											}
										}
									}elsif ($name ne "rps12+2") {
										$start1_new=$start1;
									}
									$str1=substr($sequence,($end2_new-1),($start1_new-($end2_new-1)));
									$str2=substr($sequence,($end1-1),($start3_new-($end1-1)));
									$str=$str1.$str2;
									$length=length $str;
									$rev_com=reverse $str;
									$rev_com=~ tr/ACGTacgt/TGCAtgca/;

									my $aa;
									for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
										my $codon=substr ($rev_com,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa.=$hash_codon{$codon};
										}else{
											$aa.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}

									my $intron_aa;
									if ($mark==2) {
										for (my $i=0;$i<$intron_length;$i+=3){
											my $codon=substr ($string_intron_rev_com,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$intron_aa.=$hash_codon{$codon};
											}else{
												$intron_aa.="X";
												my $j=$i+1;
												#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}
									}

									if ((abs($start2-$end3) < 4000) and ($start2 > $end2) and ($start3 > $end3) and ($start2 > $end3)) {
										if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/))) {# spliced exon,not 3x or have stop codon
											print $out_annotation "     "."gene"."            "."complement(".$end1."..".$start1_new.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."complement(join(".$end1."..".$start3_new.",".$end2_new."..".$start1_new."))"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
											$gene_number_seq{$name}++;
										}elsif (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/))) {# joined exon,no intron or no stop codon
											print $out_annotation "     "."gene"."            "."complement(".$end1."..".$start1_new.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."complement(".$end1."..".$start1_new.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (negative one-intron PCG) lost intron!\n";
										}

										my ($S1,$E1,$E2,$S3);
										my $string=substr($sequence,($end1-1),($start1_new-$end1+1));
										my $lengths=length $string;
										my $L1=$start1_new-$end2_new+1;
										my $L2=$start3_new-$end1+1;
										my $rev_coms=reverse $sequence;
										$rev_coms=~ tr/ACGTacgt/TGCAtgca/;
										for (my $i=0;$i<($length_cp-$lengths);$i+=1){
											my $repeat=substr ($sequence,$i,$lengths);
											if (($repeat eq $string) and (($i+1) ne $end1)) {
												$E1=$i+1;
												$S3=$i+$L2-1+1;
												$E2=$i+$lengths-$L1+1;
												$S1=$i+$lengths-1+1;
												print $out_annotation "     "."gene"."            "."complement(".$E1."..".$S1.")"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."complement(join(".$E1."..".$S3.",".$E2."..".$S1."))"."\n" if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)));
												print $out_annotation "     "."CDS"."             "."complement(".$E1."..".$S1.")"."\n" if (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/)));
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
												$gene_number_seq{$name}++;
											}
										}
										for (my $i=0;$i<($length_cp-$lengths);$i+=1){
											my $repeat=substr ($rev_coms,$i,$lengths);
											if ($repeat eq $string) {
												$E1=$length_cp-$i-1+1;
												$S3=$length_cp-$i-$L2+1;
												$E2=$length_cp-$i-$lengths-1+$L1+1;
												$S1=$length_cp-$i-$lengths+1;
												print $out_annotation "     "."gene"."            ".$S1."..".$E1."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."join(".$S1."..".$E2.",".$S3."..".$E1.")"."\n" if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)));
												print $out_annotation "     "."CDS"."             ".$S1."..".$E1."\n" if (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/)));
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
												$gene_number_seq{$name}++;
											}
										}
										if (($ticks == 1) and (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)))) {
											print $logfile "Warning: $name (negative one-intron PCG) intron-exon boundary need to be checked!\n";
										}
									}elsif ((abs($start2-$end3) < 4000) and (($start2 < $end2) or ($start3 < $end3) or ($start2 < $end3))) {
										print $logfile "Warning: $name (negative one-intron PCG) has not been annotated due to two different direction exons!\n";
									}elsif (abs($start2-$end3) >= 4000) {
										print $logfile "Warning: $name (negative one-intron PCG) has not been annotated due to two far exons!\n";
									}
								}elsif((grep {$_=~ /\*/} @aa_exon) or (($reverse_start_codon ne "ATG") and ($reverse_start_codon ne "GTG")) or (($reverse_stop_codon ne "TAA") and ($reverse_stop_codon ne "TAG") and ($reverse_stop_codon ne "TGA"))){# non-standard start or stop codon for _gene
									my $start_left=$start1-60;
									my $start_right;
									if (($length_cp-$start1)>=9000) {
										$start_right=$start1+9000;
									}elsif (($length_cp-$start1)<9000){
										$start_right=$start1+(int(($length_cp-$start1)/3)-1)*3;
									}
									my $length_exon1=($start1-$end2+1);
									my $end_left;
									if ($start3>=9000) {
										$end_left=$start3-9000;
									}elsif ($start3<9000){
										$end_left=$start3-((int($start3/3)-1)*3);
									}
									my $str1=substr($sequence,($start_left),($start_right-$start_left));
									my $str2=substr($sequence,($end_left),($start3-$end_left));
									my $length1=length $str1;
									my $length2=length $str2;
									my $rev_com1=reverse $str1;
									$rev_com1=~ tr/ACGTacgt/TGCAtgca/;
									my $rev_com2=reverse $str2;
									$rev_com2=~ tr/ACGTacgt/TGCAtgca/;
									my $aa1;#left range for start codon
									for (my $i=0;$i<($length1-3);$i+=3){# delete stop codon
										my $codon=substr ($rev_com1,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa1.=$hash_codon{$codon};
										}else{
											$aa1.="X";
											my $m=$i+1;
											#print "Bad codon $codon in position $m of gene $name in species $contig!\n";
										}
									}
									my @aa1=split //,$aa1;
									my (@star1,@start1);
									foreach (0..$#aa1){
										if (($aa1[$_]=~ /\*/) and ($_ < 3000)){
											push @star1,$_;
										}
										if (($aa1[$_]=~ /M/) and ($_ <= 3060)){
											push @start1,$_;
										}
									}
									my $right_star_position=pop @star1;
									until ($right_star_position<(3000+$length_exon1/3)){
										$right_star_position=pop @star1;
									}
									my @start2;
									foreach (@start1){
										push @start2,$_ if ((defined $right_star_position) and ($_ > $right_star_position))
									}
									my $right_start_position=$start2[0];
									#my $right_start_position=pop @start2;
									my $aa2;#right range for stop codon
									for (my $k=0;$k<($length2-3);$k+=3){# delete stop codon
										my $codon=substr ($rev_com2,$k,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa2.=$hash_codon{$codon};
										}else{
											$aa2.="X";
											my $n=$k+1;
											#print "Bad codon $codon in position $n of gene $name in species $contig!\n";
										}
									}
									my @aa2=split //,$aa2;
									my @star2;
									foreach (0..$#aa2){
										if ($aa2[$_]=~ /\*/){
											push @star2,$_;
										}
									}
									my $left_star_position=shift @star2;

									my $seq_string0=substr($sequence,($start1-30),60);
									my $seq_string1=substr($sequence,($end2-(int(($end2-$start3)/3))*3-61),((int(($end2-$start3)/3))*3+120));
									my $seq_string2=substr($sequence,($start3-60),((int(($end2-$start3)/3))*3+120));
									my $length_seq_string0=length $seq_string0;
									my $length_seq_string1=length $seq_string1;
									my $length_seq_string2=length $seq_string2;
									my $rev_seq_string0=reverse $seq_string0;
									$rev_seq_string0=~ tr/ACGTacgt/TGCAtgca/;
									my $rev_seq_string1=reverse $seq_string1;
									$rev_seq_string1=~ tr/ACGTacgt/TGCAtgca/;
									my $rev_seq_string2=reverse $seq_string2;
									$rev_seq_string2=~ tr/ACGTacgt/TGCAtgca/;
									my $aa_seq_string0;
									my $aa_seq_string1;
									my $aa_seq_string2;
									for (my $i=0;$i<$length_seq_string0;$i+=3){# delete stop codon
										my $codon=substr ($rev_seq_string0,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_seq_string0.=$hash_codon{$codon};
										}else{
											$aa_seq_string0.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
									for (my $i=0;$i<$length_seq_string1;$i+=3){# delete stop codon
										my $codon=substr ($rev_seq_string1,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_seq_string1.=$hash_codon{$codon};
										}else{
											$aa_seq_string1.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
									for (my $i=0;$i<$length_seq_string2;$i+=3){# delete stop codon
										my $codon=substr ($rev_seq_string2,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_seq_string2.=$hash_codon{$codon};
										}else{
											$aa_seq_string2.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}

									if ((defined $right_start_position) and (defined $left_star_position)){
										my $start=$start_right-$right_start_position * 3;
										my $end=($start3+1)-($left_star_position * 3+3);

										my ($str1,$str2,$str,$length,$rev_com,$start_new,$end2_new,$start3_new);
										my ($string_intron,$intron_length,$string_intron_rev_com);
										my $mark=0;
										my $ticks=0;
										if ($length_ref_exon1 % 3==0) {
											for (my $i=$match;$i<$a;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string1,$j,$match);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$end2_new=$end2+60-($j+$match)*3-($i-$match)*3;
														$marks++;
													}
												}
												if (defined $end2_new){
													last;
												}
											}
											if (!defined $end2_new) {
												$end2_new=$end2;
												$ticks=1;
											}
											for (my $i=0;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string2,-$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$start3_new=$start3-60+$j*3+$i*3;
														$marks++;
													}
												}
												if (defined $start3_new){
													last;
												}
											}
											if (!defined $start3_new) {
												$start3_new=$start3;
												$ticks=1;
											}
											$intron_length=$end2_new-$start3_new-1;
											if ($intron_length==0) {
												$mark=1;
											}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
												$string_intron=substr ($sequence,$start3_new,$intron_length);
												$intron_length=length $string_intron;
												$string_intron_rev_com=reverse $string_intron;
												$string_intron_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark=2;
											}
										}elsif ($length_ref_exon1 % 3==1) {
											for (my $i=$match;$i<$a;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string1,$j,$match);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$end2_new=$end2+60-($j+$match)*3-($i-$match)*3-1;
														$marks++;
													}
												}
												if (defined $end2_new){
													last;
												}
											}
											if (!defined $end2_new) {
												$end2_new=$end2-1;
												$ticks=1;
											}
											for (my $i=0;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string2,-$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$start3_new=$start3-60+$j*3+$i*3+2;
														$marks++;
													}
												}
												if (defined $start3_new){
													last;
												}
											}
											if (!defined $start3_new) {
												$start3_new=$start3+2;
												$ticks=1;
											}
											$intron_length=$end2_new-$start3_new-1;
											if ($intron_length==0) {
												$mark=1;
											}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
												$string_intron=substr ($sequence,($start3_new+1),$intron_length);
												$intron_length=length $string_intron;
												$string_intron_rev_com=reverse $string_intron;
												$string_intron_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark=2;
											}
										}elsif ($length_ref_exon1 % 3==2) {
											for (my $i=$match;$i<$a;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string1,$j,$match);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$end2_new=$end2+60-($j+$match)*3-($i-$match)*3-2;
														$marks++;
													}
												}
												if (defined $end2_new){
													last;
												}
											}
											if (!defined $end2_new) {
												$end2_new=$end2-2;
												$ticks=1;
											}
											for (my $i=0;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string2,-$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$start3_new=$start3-60+$j*3+$i*3+1;
														$marks++;
													}
												}
												if (defined $start3_new){
													last;
												}
											}
											if (!defined $start3_new) {
												$start3_new=$start3+1;
												$ticks=1;
											}
											$intron_length=$end2_new-$start3_new-1;
											if ($intron_length==0) {
												$mark=1;
											}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
												$string_intron=substr ($sequence,($start3_new-1),$intron_length);
												$intron_length=length $string_intron;
												$string_intron_rev_com=reverse $string_intron;
												$string_intron_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark=2;
											}
										}

										if ($name eq "rps12+2") {
											for (my $i=0;$i<7;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,$i,3);
												my $marks=0;
												for (my $j=0;$j<20;$j+=1) {
													my $repeat=substr ($aa_seq_string0,$j,3);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$start_new=$start1+30-$j*3+$i*3;
														$marks++;
													}
												}
											}
										}elsif ($name ne "rps12+2") {
											$start_new=$start;
										}
										$str1=substr($sequence,($end2_new-1),($start_new-($end2_new-1)));
										$str2=substr($sequence,($end-1),($start3_new-($end-1)));
										$str=$str1.$str2;
										$length=length $str;
										$rev_com=reverse $str;
										$rev_com=~ tr/ACGTacgt/TGCAtgca/;

										my $aa;
										for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
											my $codon=substr ($rev_com,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$aa.=$hash_codon{$codon};
											}else{
												$aa.="X";
												my $j=$i+1;
												#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}

										my $intron_aa;
										if ($mark==2) {
											for (my $i=0;$i<$intron_length;$i+=3){
												my $codon=substr ($string_intron_rev_com,$i,3);
												$codon=uc $codon;
												if (exists $hash_codon{$codon}){
													$intron_aa.=$hash_codon{$codon};
												}else{
													$intron_aa.="X";
													my $j=$i+1;
													#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
												}
											}
										}

										if ((abs($start2-$end3) < 4000) and ($start2 > $end2) and ($start3 > $end3) and ($start2 > $end3)) {
											if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/))) {# spliced exon,not 3x or have stop codon
												print $out_annotation "     "."gene"."            "."complement(".$end."..".$start_new.")\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."complement(join(".$end."..".$start3_new.",".$end2_new."..".$start_new."))"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa\""."\n";
												$gene_number_seq{$name}++;
											}elsif (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/))) {# joined exon,no intron or no stop codon
												print $out_annotation "     "."gene"."            "."complement(".$end."..".$start_new.")\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."complement(".$end."..".$start_new.")\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa\""."\n";
												$gene_number_seq{$name}++;
												print $logfile "Warning: $name (negative one-intron PCG) lost intron!\n";
											}

											my ($S1,$E1,$E2,$S3);
											my $string=substr($sequence,($end-1),($start_new-$end+1));
											my $lengths=length $string;
											my $L1=$start_new-$end2_new+1;
											my $L2=$start3_new-$end+1;
											my $rev_coms=reverse $sequence;
											$rev_coms=~ tr/ACGTacgt/TGCAtgca/;
											for (my $i=0;$i<($length_cp-$lengths);$i+=1){
												my $repeat=substr ($sequence,$i,$lengths);
												if (($repeat eq $string) and (($i+1) ne $end)) {
													$E1=$i+1;
													$S3=$i+$L2-1+1;
													$E2=$i+$lengths-$L1+1;
													$S1=$i+$lengths-1+1;
													print $out_annotation "     "."gene"."            "."complement(".$E1."..".$S1.")\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "     "."CDS"."             "."complement(join(".$E1."..".$S3.",".$E2."..".$S1."))"."\n" if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)));
													print $out_annotation "     "."CDS"."             "."complement(".$E1."..".$S1.")"."\n" if (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/)));
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "                     "."/codon_start=1"."\n";
													print $out_annotation "                     "."/transl_table=11"."\n";
													print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
													#print $out_annotation "                     "."/translation=\"$aa\""."\n";
													$gene_number_seq{$name}++;
												}
											}
											for (my $i=0;$i<($length_cp-$lengths);$i+=1){
												my $repeat=substr ($rev_coms,$i,$lengths);
												if ($repeat eq $string) {
													$E1=$length_cp-$i-1+1;
													$S3=$length_cp-$i-$L2+1;
													$E2=$length_cp-$i-$lengths-1+$L1+1;
													$S1=$length_cp-$i-$lengths+1;
													print $out_annotation "     "."gene"."            ".$S1."..".$E1."\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "     "."CDS"."             "."join(".$S1."..".$E2.",".$S3."..".$E1.")"."\n" if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)));
													print $out_annotation "     "."CDS"."             ".$S1."..".$E1."\n" if (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/)));
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "                     "."/codon_start=1"."\n";
													print $out_annotation "                     "."/transl_table=11"."\n";
													print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
													#print $out_annotation "                     "."/translation=\"$aa\""."\n";
													$gene_number_seq{$name}++;
												}
											}
											if (($ticks == 1) and (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)))) {
												print $logfile "Warning: $name (negative one-intron PCG) intron-exon boundary need to be checked!\n";
											}
										}elsif ((abs($start2-$end3) < 4000) and (($start2 < $end2) or ($start3 < $end3) or ($start2 < $end3))) {
											print $logfile "Warning: $name (negative one-intron PCG) has not been annotated due to two different direction exons!\n";
										}elsif (abs($start2-$end3) >= 4000) {
											print $logfile "Warning: $name (negative one-intron PCG) has not been annotated due to two far exons!\n";
										}
									}elsif ((!defined $right_start_position) and (defined $left_star_position)){
										my $end=($start3+1)-($left_star_position * 3+3);

										my ($str1,$str2,$str,$length,$rev_com,$start1_new,$end2_new,$start3_new);
										my ($string_intron,$intron_length,$string_intron_rev_com);
										my $mark=0;
										my $ticks=1;
										if ($length_ref_exon1 % 3==0) {
											for (my $i=$match;$i<$a;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string1,$j,$match);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$end2_new=$end2+60-($j+$match)*3-($i-$match)*3;
														$marks++;
													}
												}
												if (defined $end2_new){
													last;
												}
											}
											if (!defined $end2_new) {
												$end2_new=$end2;
												$ticks=1;
											}
											for (my $i=0;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string2,-$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$start3_new=$start3-60+$j*3+$i*3;
														$marks++;
													}
												}
												if (defined $start3_new){
													last;
												}
											}
											if (!defined $start3_new) {
												$start3_new=$start3;
												$ticks=1;
											}
											$intron_length=$end2_new-$start3_new-1;
											if ($intron_length==0) {
												$mark=1;
											}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
												$string_intron=substr ($sequence,$start3_new,$intron_length);
												$intron_length=length $string_intron;
												$string_intron_rev_com=reverse $string_intron;
												$string_intron_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark=2;
											}
										}elsif ($length_ref_exon1 % 3==1) {
											for (my $i=$match;$i<$a;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string1,$j,$match);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$end2_new=$end2+60-($j+$match)*3-($i-$match)*3-1;
														$marks++;
													}
												}
												if (defined $end2_new){
													last;
												}
											}
											if (!defined $end2_new) {
												$end2_new=$end2-1;
												$ticks=1;
											}
											for (my $i=0;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string2,-$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$start3_new=$start3-60+$j*3+$i*3+2;
														$marks++;
													}
												}
												if (defined $start3_new){
													last;
												}
											}
											if (!defined $start3_new) {
												$start3_new=$start3+2;
												$ticks=1;
											}
											$intron_length=$end2_new-$start3_new-1;
											if ($intron_length==0) {
												$mark=1;
											}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
												$string_intron=substr ($sequence,($start3_new+1),$intron_length);
												$intron_length=length $string_intron;
												$string_intron_rev_com=reverse $string_intron;
												$string_intron_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark=2;
											}
										}elsif ($length_ref_exon1 % 3==2) {
											for (my $i=$match;$i<$a;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string1,$j,$match);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$end2_new=$end2+60-($j+$match)*3-($i-$match)*3-2;
														$marks++;
													}
												}
												if (defined $end2_new){
													last;
												}
											}
											if (!defined $end2_new) {
												$end2_new=$end2-2;
												$ticks=1;
											}
											for (my $i=0;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string2,-$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$start3_new=$start3-60+$j*3+$i*3+1;
														$marks++;
													}
												}
												if (defined $start3_new){
													last;
												}
											}
											if (!defined $start3_new) {
												$start3_new=$start3+1;
												$ticks=1;
											}
											$intron_length=$end2_new-$start3_new-1;
											if ($intron_length==0) {
												$mark=1;
											}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
												$string_intron=substr ($sequence,($start3_new-1),$intron_length);
												$intron_length=length $string_intron;
												$string_intron_rev_com=reverse $string_intron;
												$string_intron_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark=2;
											}
										}

										if ($name eq "rps12+2") {
											for (my $i=0;$i<7;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,$i,3);
												my $marks=0;
												for (my $j=0;$j<20;$j+=1) {
													my $repeat=substr ($aa_seq_string0,$j,3);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$start1_new=$start1+30-$j*3+$i*3;
														$marks++;
													}
												}
											}
										}elsif ($name ne "rps12+2") {
											$start1_new=$start1;
										}
										$str1=substr($sequence,($end2_new-1),($start1_new-($end2_new-1)));
										$str2=substr($sequence,($end-1),($start3_new-($end-1)));
										$str=$str1.$str2;
										$length=length $str;
										$rev_com=reverse $str;
										$rev_com=~ tr/ACGTacgt/TGCAtgca/;

										my $aa;
										for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
											my $codon=substr ($rev_com,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$aa.=$hash_codon{$codon};
											}else{
												$aa.="X";
												my $j=$i+1;
												#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}

										my $intron_aa;
										if ($mark==2) {
											for (my $i=0;$i<$intron_length;$i+=3){
												my $codon=substr ($string_intron_rev_com,$i,3);
												$codon=uc $codon;
												if (exists $hash_codon{$codon}){
													$intron_aa.=$hash_codon{$codon};
												}else{
													$intron_aa.="X";
													my $j=$i+1;
													#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
												}
											}
										}

										if ((abs($start2-$end3) < 4000) and ($start2 > $end2) and ($start3 > $end3) and ($start2 > $end3)) {
											if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/))) {# spliced exon,not 3x or have stop codon
												print $out_annotation "     "."gene"."            "."complement(".$end."..".$start1_new.")\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."complement(join(".$end."..".$start3_new.",".$end2_new."..".$start1_new."))"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa\""."\n";
												$gene_number_seq{$name}++;
											}elsif (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/))) {# joined exon,no intron or no stop codon
												print $out_annotation "     "."gene"."            "."complement(".$end."..".$start1_new.")\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."complement(".$end."..".$start1_new.")\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa\""."\n";
												$gene_number_seq{$name}++;
												print $logfile "Warning: $name (negative one-intron PCG) lost intron!\n";
											}

											my ($S1,$E1,$E2,$S3);
											my $string=substr($sequence,($end-1),($start1_new-$end+1));
											my $lengths=length $string;
											my $L1=$start1_new-$end2_new+1;
											my $L2=$start3_new-$end+1;
											my $rev_coms=reverse $sequence;
											$rev_coms=~ tr/ACGTacgt/TGCAtgca/;
											for (my $i=0;$i<($length_cp-$lengths);$i+=1){
												my $repeat=substr ($sequence,$i,$lengths);
												if (($repeat eq $string) and (($i+1) ne $end)) {
													$E1=$i+1;
													$S3=$i+$L2-1+1;
													$E2=$i+$lengths-$L1+1;
													$S1=$i+$lengths-1+1;
													print $out_annotation "     "."gene"."            "."complement(".$E1."..".$S1.")\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "     "."CDS"."             "."complement(join(".$E1."..".$S3.",".$E2."..".$S1."))"."\n" if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)));
													print $out_annotation "     "."CDS"."             "."complement(".$E1."..".$S1.")"."\n" if (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/)));
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "                     "."/codon_start=1"."\n";
													print $out_annotation "                     "."/transl_table=11"."\n";
													print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
													#print $out_annotation "                     "."/translation=\"$aa\""."\n";
													$gene_number_seq{$name}++;
												}
											}
											for (my $i=0;$i<($length_cp-$lengths);$i+=1){
												my $repeat=substr ($rev_coms,$i,$lengths);
												if ($repeat eq $string) {
													$E1=$length_cp-$i-1+1;
													$S3=$length_cp-$i-$L2+1;
													$E2=$length_cp-$i-$lengths-1+$L1+1;
													$S1=$length_cp-$i-$lengths+1;
													print $out_annotation "     "."gene"."            ".$S1."..".$E1."\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "     "."CDS"."             "."join(".$S1."..".$E2.",".$S3."..".$E1.")"."\n" if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)));
													print $out_annotation "     "."CDS"."             ".$S1."..".$E1."\n" if (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/)));
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "                     "."/codon_start=1"."\n";
													print $out_annotation "                     "."/transl_table=11"."\n";
													print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
													#print $out_annotation "                     "."/translation=\"$aa\""."\n";
													$gene_number_seq{$name}++;
												}
											}
											if (($ticks == 1) and (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)))) {
												print $logfile "Warning: $name (negative one-intron PCG) has non-canonical start codon and intron-exon boundary need to be checked!\n" if ($name ne "rps12+2");
											}elsif ($ticks == 0) {
												print $logfile "Warning: $name (negative one-intron PCG) has non-canonical start codon!\n" if ($name ne "rps12+2");
											}
										}elsif ((abs($start2-$end3) < 4000) and (($start2 < $end2) or ($start3 < $end3) or ($start2 < $end3))) {
											print $logfile "Warning: $name (negative one-intron PCG) has not been annotated due to two different direction exons!\n";
										}elsif (abs($start2-$end3) >= 4000) {
											print $logfile "Warning: $name (negative one-intron PCG) has not been annotated due to two far exons!\n";
										}
									}else{
										#print $out_annotation "     "."gene"."            "."complement(".$end1."..".$start1.")\n";
										#print $out_annotation "                     "."/gene=\"$name\""."\n";
										#print $out_annotation "     "."CDS"."             "."complement(join(".$end1."..".$start3.",".$end2."..".$start1."))"."\n";
										#print $out_annotation "                     "."/gene=\"$name\""."\n";
										#print $out_annotation "                     "."/codon_start=1"."\n";
										#print $out_annotation "                     "."/transl_table=11"."\n";
										#print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										##print $out_annotation "                     "."/translation=\"$aa\""."\n";
										#$gene_number_seq{$name}++;
										#print $logfile "Warning: $name (negative one-intron PCG) maybe pseudogene!\n";
										print $logfile "Warning: $name (negative one-intron PCG) has not been annotated!\n";
									}
								}
							}elsif (($name eq "rpl16") or ($name eq "petB") or ($name eq "petD")) {
								if ($length_exon2 < 150) {
									print $logfile "Warning: $name (negative one-intron PCG) has not been annotated due to short length!\n";
								}elsif ($length_exon2 >= 150) {
									if (($length_ref_exon1 % 3==0) and ($length_exon % 3==0) and (!grep {$_=~ /\*/} @aa_exon) and (($reverse_start_codon eq "ATG") or ($reverse_start_codon eq "GTG")) and (($reverse_stop_codon eq "TAA") or ($reverse_stop_codon eq "TAG") or ($reverse_stop_codon eq "TGA"))){# standard start and stop codon for _gene
										print $out_annotation "     "."gene"."            "."complement(".$end1."..".$start1.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."complement(join(".$end1."..".$start3.",".$end2."..".$start1."))"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
										$gene_number_seq{$name}++;
									}elsif(($length_ref_exon1 % 3!=0) or (grep {$_=~ /\*/} @aa_exon) or (($reverse_start_codon ne "ATG") and ($reverse_start_codon ne "GTG")) or (($reverse_stop_codon ne "TAA") and ($reverse_stop_codon ne "TAG") and ($reverse_stop_codon ne "TGA"))){# non-standard start or stop codon for _gene
										my $start_left=$start1-60;
										my $start_right;
										if (($length_cp-$start1)>=9000) {
											$start_right=$start1+9000;
										}elsif (($length_cp-$start1)<9000){
											$start_right=$start1+(int(($length_cp-$start1)/3)-1)*3;
										}
										my $length_exon1=($start1-$end2+1);
										my $end_left;
										if ($start3>=9000) {
											$end_left=$start3-9000;
										}elsif ($start3<9000){
											$end_left=$start3-((int($start3/3)-1)*3);
										}
										my $str1=substr($sequence,($start_left),($start_right-$start_left));
										my $str2=substr($sequence,($end_left),($start3-$end_left));
										my $length1=length $str1;
										my $length2=length $str2;
										my $rev_com1=reverse $str1;
										$rev_com1=~ tr/ACGTacgt/TGCAtgca/;
										my $rev_com2=reverse $str2;
										$rev_com2=~ tr/ACGTacgt/TGCAtgca/;
										my $coding_1_length=$hash_exon_length{"$name-1_coding"};
										my $coding_1_sequence=$hash_exon_sequence{"$name-1_coding"};
										my $right_start_position;
										for (my $i=0;$i<($length1-3);$i+=1){
											my $exon=substr ($rev_com1,$i,$coding_1_length);
											$exon=lc $exon;
											if ($exon eq $coding_1_sequence) {
												$right_start_position=$i;
											}
										}

										my $aa2;#right range for stop codon
										for (my $k=0;$k<($length2-3);$k+=3){# delete stop codon
											my $codon=substr ($rev_com2,$k,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$aa2.=$hash_codon{$codon};
											}else{
												$aa2.="X";
												my $n=$k+1;
												#print "Bad codon $codon in position $n of gene $name in species $contig!\n";
											}
										}
										my @aa2=split //,$aa2;
										my @star2;
										foreach (0..$#aa2){
											if ($aa2[$_]=~ /\*/){
												push @star2,$_;
											}
										}
										my $left_star_position=shift @star2;
										if ((defined $right_start_position) and (defined $left_star_position)){
											my $start=$start_right-$right_start_position;
											my $end=($start3+1)-($left_star_position * 3+3);

											my ($str1,$str2,$str,$length,$rev_com,$end2_new,$start3_new);
											$end2_new=$start-$coding_1_length+1;
											if ($length_ref_exon1 % 3==0) {
												$start3_new=$start3;
											}elsif ($length_ref_exon1 % 3==1) {
												$start3_new=$start3+2;
											}elsif ($length_ref_exon1 % 3==2) {
												$start3_new=$start3+1;
											}
											$str1=substr($sequence,($end2_new-1),($start-($end2_new-1)));
											$str2=substr($sequence,($end-1),($start3_new-($end-1)));
											$str=$str1.$str2;
											$length=length $str;
											$rev_com=reverse $str;
											$rev_com=~ tr/ACGTacgt/TGCAtgca/;

											my $aa;
											for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
												my $codon=substr ($rev_com,$i,3);
												$codon=uc $codon;
												if (exists $hash_codon{$codon}){
													$aa.=$hash_codon{$codon};
												}else{
													$aa.="X";
													my $j=$i+1;
													#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
												}
											}
											print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."complement(join(".$end."..".$start3_new.",".$end2_new."..".$start."))"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa\""."\n";
											$gene_number_seq{$name}++;
										}elsif ((!defined $right_start_position) and (defined $left_star_position)){
											my $end=($start3+1)-($left_star_position * 3+3);

											my ($str1,$str2,$str,$length,$rev_com,$end2_new,$start3_new);
											$end2_new=$start1-$coding_1_length+1;
											if ($length_ref_exon1 % 3==0) {
												$start3_new=$start3;
											}elsif ($length_ref_exon1 % 3==1) {
												$start3_new=$start3+2;
											}elsif ($length_ref_exon1 % 3==2) {
												$start3_new=$start3+1;
											}
											$str1=substr($sequence,($end2_new-1),($start1-($end2_new-1)));
											$str2=substr($sequence,($end-1),($start3_new-($end-1)));
											$str=$str1.$str2;
											$length=length $str;
											$rev_com=reverse $str;
											$rev_com=~ tr/ACGTacgt/TGCAtgca/;

											my $aa;
											for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
												my $codon=substr ($rev_com,$i,3);
												$codon=uc $codon;
												if (exists $hash_codon{$codon}){
													$aa.=$hash_codon{$codon};
												}else{
													$aa.="X";
													my $j=$i+1;
													#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
												}
											}
											#print $out_annotation "     "."gene"."            "."complement(".$end."..".$start1.")\n";
											print $out_annotation "     "."gene"."            "."complement(".$end."..".$start3.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											#print $out_annotation "     "."CDS"."             "."complement(join(".$end."..".$start3_new.",".$end2_new."..".$start1."))"."\n";
											print $out_annotation "     "."CDS"."             "."complement(".$end."..".$start3.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (negative one-intron PCG) must be checked due to exon 1 not found!\n";
										}
									}
								}
							}
						}
					}elsif (((($start1 != $start2) and (!defined $start4)) or ((!defined $start4) and ($end1 != $end3))) and ((defined $length_ref_exon1) and (defined $length_ref_exon2) and (defined $sequence_ref_exon1) and (defined $sequence_ref_exon2))){# non-identical PCG boundary for _gene and -1_coding_aa, -2_coding_aa
						if ($start2 < $end3){# positive
							my ($str1,$str2);
							$str1=substr($sequence,($start2-1),($end2-$start2+1)) if ($end2 > $start2);
							$str1=substr($sequence,($end2-1),($start2-$end2+1)) if ($end2 < $start2);
							$str2=substr($sequence,($start3-1),($end3-$start3+1)) if ($end3 > $start3);
							$str2=substr($sequence,($end3-1),($start3-$end3+1)) if ($end3 < $start3);
							my $str3=$str1.$str2;
							my $length_exon=length $str3;
							my $length_exon2=length $str2;
							my $forward_start_codon=substr($str3,0,3);
							my $forward_stop_codon=substr($str3,($length_exon-3),3);
							my $aa_exon;
							for (my $i=0;$i<($length_exon-3);$i+=3){# delete stop codon
								my $codon=substr ($str3,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa_exon.=$hash_codon{$codon};
								}else{
									$aa_exon.="X";
									my $j=$i+1;
									#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa_exon=split //,$aa_exon;
							if (($name ne "rpl16") and ($name ne "petB") and ($name ne "petD")) {
								if (($length_exon % 3==0) and (!grep {$_=~ /\*/} @aa_exon) and (($forward_start_codon eq "ATG") or ($forward_start_codon eq "GTG") or ($name eq "rps12+2")) and (($forward_stop_codon eq "TAA") or ($forward_stop_codon eq "TAG") or ($forward_stop_codon eq "TGA"))){# standard start and stop codon for _gene
									my $seq_string0=substr($sequence,($start2-31),60);
									my $seq_string1=substr($sequence,($end2-60),((int(($start3-$end2)/3))*3+120));
									my $seq_string2=substr($sequence,($start3-(int(($start3-$end2)/3))*3-61),((int(($start3-$end2)/3))*3+120));
									my $length_seq_string0=length $seq_string0;
									my $length_seq_string1=length $seq_string1;
									my $length_seq_string2=length $seq_string2;
									my $aa_seq_string0;
									my $aa_seq_string1;
									my $aa_seq_string2;
									for (my $i=0;$i<$length_seq_string0;$i+=3){# delete stop codon
										my $codon=substr ($seq_string0,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_seq_string0.=$hash_codon{$codon};
										}else{
											$aa_seq_string0.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
									for (my $i=0;$i<$length_seq_string1;$i+=3){# delete stop codon
										my $codon=substr ($seq_string1,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_seq_string1.=$hash_codon{$codon};
										}else{
											$aa_seq_string1.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
									for (my $i=0;$i<$length_seq_string2;$i+=3){# delete stop codon
										my $codon=substr ($seq_string2,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_seq_string2.=$hash_codon{$codon};
										}else{
											$aa_seq_string2.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}

									my ($str1,$str2,$str,$length,$start2_new,$end2_new,$start3_new);
									my ($string_intron,$intron_length);
									my $mark=0;
									my $ticks=0;
									if ($length_ref_exon1 % 3==0) {
										for (my $i=$match;$i<$a;$i+=1) {
											my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
												my $repeat=substr ($aa_seq_string1,$j,$match);
												if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
													$end2_new=$end2-60+($j+$match)*3+($i-$match)*3;
													$marks++;
												}
											}
											if (defined $end2_new){
												last;
											}
										}
										if (!defined $end2_new) {
											$end2_new=$end2;
											$ticks=1;
										}
										for (my $i=0;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
												my $repeat=substr ($aa_seq_string2,-$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$start3_new=$start3+60-$j*3-$i*3;
													$marks++;
												}
											}
											if (defined $start3_new){
												last;
											}
										}
										if (!defined $start3_new) {
											$start3_new=$start3;
											$ticks=1;
										}
										$intron_length=$start3_new-$end2_new-1;
										if ($intron_length==0) {
											$mark=1;
										}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
											$string_intron=substr ($sequence,$end2_new,$intron_length);
											$intron_length=length $string_intron;
											$mark=2;
										}
									}elsif ($length_ref_exon1 % 3==1) {
										for (my $i=$match;$i<$a;$i+=1) {
											my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
												my $repeat=substr ($aa_seq_string1,$j,$match);
												if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
													$end2_new=$end2-60+($j+$match)*3+($i-$match)*3+1;
													$marks++;
												}
											}
											if (defined $end2_new){
												last;
											}
										}
										if (!defined $end2_new) {
											$end2_new=$end2+1;
											$ticks=1;
										}
										for (my $i=0;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
												my $repeat=substr ($aa_seq_string2,-$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$start3_new=$start3+60-$j*3-$i*3-2;
													$marks++;
												}
											}
											if (defined $start3_new){
												last;
											}
										}
										if (!defined $start3_new) {
											$start3_new=$start3-2;
											$ticks=1;
										}
										$intron_length=$start3_new-$end2_new-1;
										if ($intron_length==0) {
											$mark=1;
										}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
											$string_intron=substr ($sequence,($end2_new-1),$intron_length);
											$intron_length=length $string_intron;
											$mark=2;
										}
									}elsif ($length_ref_exon1 % 3==2) {
										for (my $i=$match;$i<$a;$i+=1) {
											my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
												my $repeat=substr ($aa_seq_string1,$j,$match);
												if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
													$end2_new=$end2-60+($j+$match)*3+($i-$match)*3+2;
													$marks++;
												}
											}
											if (defined $end2_new){
												last;
											}
										}
										if (!defined $end2_new) {
											$end2_new=$end2+2;
											$ticks=1;
										}
										for (my $i=0;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
												my $repeat=substr ($aa_seq_string2,-$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$start3_new=$start3+60-$j*3-$i*3-1;
													$marks++;
												}
											}
											if (defined $start3_new){
												last;
											}
										}
										if (!defined $start3_new) {
											$start3_new=$start3-1;
											$ticks=1;
										}
										$intron_length=$start3_new-$end2_new-1;
										if ($intron_length==0) {
											$mark=1;
										}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
											$string_intron=substr ($sequence,$end2_new+1,$intron_length);
											$intron_length=length $string_intron;
											$mark=2;
										}
									}

									if ($name eq "rps12+2") {
										for (my $i=0;$i<7;$i+=1) {
											my $match_aa_ref_exon1=substr ($aa_ref_exon1,$i,3);
											my $marks=0;
											for (my $j=0;$j<20;$j+=1) {
												my $repeat=substr ($aa_seq_string0,$j,3);
												if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
													$start2_new=$start2-30+$j*3-$i*3;
													$marks++;
												}
											}
										}
									}elsif ($name ne "rps12+2") {
										$start2_new=$start2;
									}
									$str1=substr($sequence,($start2_new-1),($end2_new-($start2_new-1)));
									$str2=substr($sequence,($start3_new-1),($end3-($start3_new-1)));
									$str=$str1.$str2;
									$length=length $str;

									my $aa;
									for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
										my $codon=substr ($str,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa.=$hash_codon{$codon};
										}else{
											$aa.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}

									my $intron_aa;
									if ($mark==2) {
										for (my $i=0;$i<$intron_length;$i+=3){
											my $codon=substr ($string_intron,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$intron_aa.=$hash_codon{$codon};
											}else{
												$intron_aa.="X";
												my $j=$i+1;
												#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}
									}

									if ((abs($end3-$start2) < 4000) and ($start2 < $end2) and ($start3 < $end3) and ($start2 < $end3)) {
										if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/))) {# spliced exon,not 3x or have stop codon
											print $out_annotation "     "."gene"."            ".$start2_new."..".$end3."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."join(".$start2_new."..".$end2_new.",".$start3_new."..".$end3.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
											$gene_number_seq{$name}++;
										}elsif (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/))) {# joined exon,no intron or no stop codon
											print $out_annotation "     "."gene"."            ".$start2_new."..".$end3."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             ".$start2_new."..".$end3."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (positive one-intron PCG) lost intron!\n";
										}

										my ($S1,$E1,$E2,$S3);
										my $string=substr($sequence,($start2_new-1),($end3-$start2_new+1));
										my $lengths=length $string;
										my $L1=$end2_new-$start2_new+1;
										my $L2=$end3-$start3_new+1;
										for (my $i=0;$i<($length_cp-$lengths);$i+=1){
											my $repeat=substr ($sequence,$i,$lengths);
											if (($repeat eq $string) and (($i+1) ne $start2_new)) {
												$S1=$i+1;
												$E2=$i+$L1-1+1;
												$S3=$i+$lengths-$L2+1;
												$E1=$i+$lengths-1+1;
												print $out_annotation "     "."gene"."            ".$S1."..".$E1."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."join(".$S1."..".$E2.",".$S3."..".$E1.")"."\n" if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)));
												print $out_annotation "     "."CDS"."             ".$S1."..".$E1."\n" if (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/)));
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
												$gene_number_seq{$name}++;
											}
										}
										for (my $i=0;$i<($length_cp-$lengths);$i+=1){
											my $repeat=substr ($rev_coms,$i,$lengths);
											if ($repeat eq $string) {
												$S1=$length_cp-$i-1+1;
												$E2=$length_cp-$i-$L1+1;
												$S3=$length_cp-$i-$lengths-1+$L2+1;
												$E1=$length_cp-$i-$lengths+1;
												print $out_annotation "     "."gene"."            "."complement(".$E1."..".$S1.")"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."complement(join(".$E1."..".$S3.",".$E2."..".$S1."))"."\n" if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)));
												print $out_annotation "     "."CDS"."             "."complement(".$E1."..".$S1.")"."\n" if (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/)));
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
												$gene_number_seq{$name}++;
											}
										}
										if (($ticks == 1) and (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)))) {
											print $logfile "Warning: $name (positive one-intron PCG) intron-exon boundary need to be checked!\n";
										}
									}elsif ((abs($end3-$start2) < 4000) and (($start2 > $end2) or ($start3 > $end3) or ($start2 > $end3))) {
										print $logfile "Warning: $name (positive one-intron PCG) has not been annotated due to two different direction exons!\n";
									}elsif (abs($end3-$start2) >= 4000) {
										print $logfile "Warning: $name (positive one-intron PCG) has not been annotated due to two far exons!\n";
									}
								}elsif((grep {$_=~ /\*/} @aa_exon) or (($forward_start_codon ne "ATG") and ($forward_start_codon ne "GTG")) or (($forward_stop_codon ne "TAA") and ($forward_stop_codon ne "TAG") and ($forward_stop_codon ne "TGA"))){# non-standard start or stop codon for _gene
									my $start_left;
									if ($start2>=9000) {
										$start_left=$start2-9000;
									}elsif ($start2<9000){
										$start_left=$start2-(int($start2/3)-1)*3;
									}
									my $start_right=$start2+60;
									my $length_exon1=($end2-$start2+1);
									my $end_right;
									if (($length_cp-$start3)>=9000) {
										$end_right=$start3+9000;
									}elsif (($length_cp-$start3)<9000){
										$end_right=$start3+((int(($length_cp-$start3)/3)-1)*3);
									}
									my $str1=substr($sequence,($start_left-1),($start_right-$start_left));
									my $str2=substr($sequence,($start3-1),($end_right-$start3));
									my $length1=length $str1;
									my $length2=length $str2;
									my $aa1;#left range for start codon
									for (my $i=0;$i<($length1-3);$i+=3){# delete stop codon
										my $codon=substr ($str1,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa1.=$hash_codon{$codon};
										}else{
											$aa1.="X";
											my $m=$i+1;
											#print "Bad codon $codon in position $m of gene $name in species $contig!\n";
										}
									}
									my @aa1=split //,$aa1;
									my (@star1,@start1);
									foreach (0..$#aa1){
										if (($aa1[$_]=~ /\*/) and ($_ < 3000)){
											push @star1,$_;
										}
										if (($aa1[$_]=~ /M/) and ($_ <= 3060)){
											push @start1,$_;
											}
									}
									my $left_star_position=pop @star1;
									until ($left_star_position<(3000+$length_exon1/3)){
										$left_star_position=pop @star1;
									}
									my @start2;
									foreach (@start1){
										push @start2,$_ if ((defined $left_star_position) and ($_ > $left_star_position))
									}
									my $left_start_position=$start2[0];
									#my $left_start_position=pop @start2;
									my $aa2;#right range for stop codon
									for (my $k=0;$k<($length2-3);$k+=3){# delete stop codon
										my $codon=substr ($str2,$k,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa2.=$hash_codon{$codon};
										}else{
											$aa2.="X";
											my $n=$k+1;
											#print "Bad codon $codon in position $n of gene $name in species $contig!\n";
										}
									}
									my @aa2=split //,$aa2;
									my @star2;
									foreach (0..$#aa2){
										if ($aa2[$_]=~ /\*/){
											push @star2,$_;
										}
									}
									my $right_star_position=shift @star2;

									my $seq_string0=substr($sequence,($start2-31),60);
									my $seq_string1=substr($sequence,($end2-60),((int(($start3-$end2)/3))*3+120));
									my $seq_string2=substr($sequence,($start3-(int(($start3-$end2)/3))*3-61),((int(($start3-$end2)/3))*3+120));
									my $length_seq_string0=length $seq_string0;
									my $length_seq_string1=length $seq_string1;
									my $length_seq_string2=length $seq_string2;
									my $aa_seq_string0;
									my $aa_seq_string1;
									my $aa_seq_string2;
									for (my $i=0;$i<$length_seq_string0;$i+=3){# delete stop codon
										my $codon=substr ($seq_string0,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_seq_string0.=$hash_codon{$codon};
										}else{
											$aa_seq_string0.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
									for (my $i=0;$i<$length_seq_string1;$i+=3){# delete stop codon
										my $codon=substr ($seq_string1,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_seq_string1.=$hash_codon{$codon};
										}else{
											$aa_seq_string1.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
									for (my $i=0;$i<$length_seq_string2;$i+=3){# delete stop codon
										my $codon=substr ($seq_string2,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_seq_string2.=$hash_codon{$codon};
										}else{
											$aa_seq_string2.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}

									if ((defined $left_start_position) and (defined $right_star_position)){
										my $start=$start_left+$left_start_position * 3;
										my $end=($start3-1)+($right_star_position * 3+3);

										my ($str1,$str2,$str,$length,$start_new,$end2_new,$start3_new);
										my ($string_intron,$intron_length);
										my $mark=0;
										my $ticks=0;
										if ($length_ref_exon1 % 3==0) {
											for (my $i=$match;$i<$a;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string1,$j,$match);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$end2_new=$end2-60+($j+$match)*3+($i-$match)*3;
														$marks++;
													}
												}
												if (defined $end2_new){
													last;
												}
											}
											if (!defined $end2_new) {
												$end2_new=$end2;
												$ticks=1;
											}
											for (my $i=0;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string2,-$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$start3_new=$start3+60-$j*3-$i*3;
														$marks++;
													}
												}
												if (defined $start3_new){
													last;
												}
											}
											if (!defined $start3_new) {
												$start3_new=$start3;
												$ticks=1;
											}
											$intron_length=$start3_new-$end2_new-1;
											if ($intron_length==0) {
												$mark=1;
											}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
												$string_intron=substr ($sequence,$end2_new,$intron_length);
												$intron_length=length $string_intron;
												$mark=2;
											}
										}elsif ($length_ref_exon1 % 3==1) {
											for (my $i=$match;$i<$a;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string1,$j,$match);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$end2_new=$end2-60+($j+$match)*3+($i-$match)*3+1;
														$marks++;
													}
												}
												if (defined $end2_new){
													last;
												}
											}
											if (!defined $end2_new) {
												$end2_new=$end2+1;
												$ticks=1;
											}
											for (my $i=0;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string2,-$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$start3_new=$start3+60-$j*3-$i*3-2;
														$marks++;
													}
												}
												if (defined $start3_new){
													last;
												}
											}
											if (!defined $start3_new) {
												$start3_new=$start3-2;
												$ticks=1;
											}
											$intron_length=$start3_new-$end2_new-1;
											if ($intron_length==0) {
												$mark=1;
											}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
												$string_intron=substr ($sequence,($end2_new-1),$intron_length);
												$intron_length=length $string_intron;
												$mark=2;
											}
										}elsif ($length_ref_exon1 % 3==2) {
											for (my $i=$match;$i<$a;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string1,$j,$match);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$end2_new=$end2-60+($j+$match)*3+($i-$match)*3+2;
														$marks++;
													}
												}
												if (defined $end2_new){
													last;
												}
											}
											if (!defined $end2_new) {
												$end2_new=$end2+2;
												$ticks=1;
											}
											for (my $i=0;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string2,-$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$start3_new=$start3+60-$j*3-$i*3-1;
														$marks++;
													}
												}
												if (defined $start3_new){
													last;
												}
											}
											if (!defined $start3_new) {
												$start3_new=$start3-1;
												$ticks=1;
											}
											$intron_length=$start3_new-$end2_new-1;
											if ($intron_length==0) {
												$mark=1;
											}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
												$string_intron=substr ($sequence,$end2_new+1,$intron_length);
												$intron_length=length $string_intron;
												$mark=2;
											}
										}

										if ($name eq "rps12+2") {
											for (my $i=0;$i<7;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,$i,3);
												my $marks=0;
												for (my $j=0;$j<20;$j+=1) {
													my $repeat=substr ($aa_seq_string0,$j,3);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$start_new=$start2-30+$j*3-$i*3;
														$marks++;
													}
												}
											}
										}elsif ($name ne "rps12+2") {
											$start_new=$start;
										}
										$str1=substr($sequence,($start_new-1),($end2_new-($start_new-1)));
										$str2=substr($sequence,($start3_new-1),($end-($start3_new-1)));
										$str=$str1.$str2;
										$length=length $str;

										my $aa;
										for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
											my $codon=substr ($str,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$aa.=$hash_codon{$codon};
											}else{
												$aa.="X";
												my $j=$i+1;
												#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}

										my $intron_aa;
										if ($mark==2) {
											for (my $i=0;$i<$intron_length;$i+=3){
												my $codon=substr ($string_intron,$i,3);
												$codon=uc $codon;
												if (exists $hash_codon{$codon}){
													$intron_aa.=$hash_codon{$codon};
												}else{
													$intron_aa.="X";
													my $j=$i+1;
													#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
												}
											}
										}

										if ((abs($end3-$start2) < 4000) and ($start2 < $end2) and ($start3 < $end3) and ($start2 < $end3)) {
											if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/))) {# spliced exon,not 3x or have stop codon
												print $out_annotation "     "."gene"."            ".$start_new."..".$end."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."join(".$start_new."..".$end2_new.",".$start3_new."..".$end.")"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa\""."\n";
												$gene_number_seq{$name}++;
											}elsif (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/))) {# joined exon,no intron or no stop codon
												print $out_annotation "     "."gene"."            ".$start_new."..".$end."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             ".$start_new."..".$end."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa\""."\n";
												$gene_number_seq{$name}++;
												print $logfile "Warning: $name (positive one-intron PCG) lost intron!\n";
											}

											my ($S1,$E1,$E2,$S3);
											my $string=substr($sequence,($start_new-1),($end-$start_new+1));
											my $lengths=length $string;
											my $L1=$end2_new-$start_new+1;
											my $L2=$end-$start3_new+1;
											for (my $i=0;$i<($length_cp-$lengths);$i+=1){
												my $repeat=substr ($sequence,$i,$lengths);
												if (($repeat eq $string) and (($i+1) ne $start_new)) {
													$S1=$i+1;
													$E2=$i+$L1-1+1;
													$S3=$i+$lengths-$L2+1;
													$E1=$i+$lengths-1+1;
													print $out_annotation "     "."gene"."            ".$S1."..".$E1."\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "     "."CDS"."             "."join(".$S1."..".$E2.",".$S3."..".$E1.")"."\n" if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)));
													print $out_annotation "     "."CDS"."             ".$S1."..".$E1."\n" if (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/)));
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "                     "."/codon_start=1"."\n";
													print $out_annotation "                     "."/transl_table=11"."\n";
													print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
													#print $out_annotation "                     "."/translation=\"$aa\""."\n";
													$gene_number_seq{$name}++;
												}
											}
											for (my $i=0;$i<($length_cp-$lengths);$i+=1){
												my $repeat=substr ($rev_coms,$i,$lengths);
												if ($repeat eq $string) {
													$S1=$length_cp-$i-1+1;
													$E2=$length_cp-$i-$L1+1;
													$S3=$length_cp-$i-$lengths-1+$L2+1;
													$E1=$length_cp-$i-$lengths+1;
													print $out_annotation "     "."gene"."            "."complement(".$E1."..".$S1.")"."\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "     "."CDS"."             "."complement(join(".$E1."..".$S3.",".$E2."..".$S1."))"."\n" if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)));
													print $out_annotation "     "."CDS"."             "."complement(".$E1."..".$S1.")"."\n" if (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/)));
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "                     "."/codon_start=1"."\n";
													print $out_annotation "                     "."/transl_table=11"."\n";
													print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
													#print $out_annotation "                     "."/translation=\"$aa\""."\n";
													$gene_number_seq{$name}++;
												}
											}
											if (($ticks == 1) and (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)))) {
												print $logfile "Warning: $name (positive one-intron PCG) intron-exon boundary need to be checked!\n";
											}
										}elsif ((abs($end3-$start2) < 4000) and (($start2 > $end2) or ($start3 > $end3) or ($start2 > $end3))) {
											print $logfile "Warning: $name (positive one-intron PCG) has not been annotated due to two different direction exons!\n";
										}elsif (abs($end3-$start2) >= 4000) {
											print $logfile "Warning: $name (positive one-intron PCG) has not been annotated due to two far exons!\n";
										}
									}elsif ((!defined $left_start_position) and (defined $right_star_position)){
										my $end=($start3-1)+($right_star_position * 3+3);

										my ($str1,$str2,$str,$length,$start2_new,$end2_new,$start3_new);
										my ($string_intron,$intron_length);
										my $mark=0;
										my $ticks=0;
										if ($length_ref_exon1 % 3==0) {
											for (my $i=$match;$i<$a;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string1,$j,$match);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$end2_new=$end2-60+($j+$match)*3+($i-$match)*3;
														$marks++;
													}
												}
												if (defined $end2_new){
													last;
												}
											}
											if (!defined $end2_new) {
												$end2_new=$end2;
												$ticks=1;
											}
											for (my $i=0;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string2,-$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$start3_new=$start3+60-$j*3-$i*3;
														$marks++;
													}
												}
												if (defined $start3_new){
													last;
												}
											}
											if (!defined $start3_new) {
												$start3_new=$start3;
												$ticks=1;
											}
											$intron_length=$start3_new-$end2_new-1;
											if ($intron_length==0) {
												$mark=1;
											}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
												$string_intron=substr ($sequence,$end2_new,$intron_length);
												$intron_length=length $string_intron;
												$mark=2;
											}
										}elsif ($length_ref_exon1 % 3==1) {
											for (my $i=$match;$i<$a;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string1,$j,$match);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$end2_new=$end2-60+($j+$match)*3+($i-$match)*3+1;
														$marks++;
													}
												}
												if (defined $end2_new){
													last;
												}
											}
											if (!defined $end2_new) {
												$end2_new=$end2+1;
												$ticks=1;
											}
											for (my $i=0;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string2,-$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$start3_new=$start3+60-$j*3-$i*3-2;
														$marks++;
													}
												}
												if (defined $start3_new){
													last;
												}
											}
											if (!defined $start3_new) {
												$start3_new=$start3-2;
												$ticks=1;
											}
											$intron_length=$start3_new-$end2_new-1;
											if ($intron_length==0) {
												$mark=1;
											}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
												$string_intron=substr ($sequence,($end2_new-1),$intron_length);
												$intron_length=length $string_intron;
												$mark=2;
											}
										}elsif ($length_ref_exon1 % 3==2) {
											for (my $i=$match;$i<$a;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string1,$j,$match);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$end2_new=$end2-60+($j+$match)*3+($i-$match)*3+2;
														$marks++;
													}
												}
												if (defined $end2_new){
													last;
												}
											}
											if (!defined $end2_new) {
												$end2_new=$end2+2;
												$ticks=1;
											}
											for (my $i=0;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string2,-$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$start3_new=$start3+60-$j*3-$i*3-1;
														$marks++;
													}
												}
												if (defined $start3_new){
													last;
												}
											}
											if (!defined $start3_new) {
												$start3_new=$start3-1;
												$ticks=1;
											}
											$intron_length=$start3_new-$end2_new-1;
											if ($intron_length==0) {
												$mark=1;
											}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
												$string_intron=substr ($sequence,$end2_new+1,$intron_length);
												$intron_length=length $string_intron;
												$mark=2;
											}
										}

										if ($name eq "rps12+2") {
											for (my $i=0;$i<7;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,$i,3);
												my $marks=0;
												for (my $j=0;$j<20;$j+=1) {
													my $repeat=substr ($aa_seq_string0,$j,3);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$start2_new=$start2-30+$j*3-$i*3;
														$marks++;
													}
												}
											}
										}elsif ($name ne "rps12+2") {
											$start2_new=$start2;
										}
										$str1=substr($sequence,($start2_new-1),($end2_new-($start2_new-1)));
										$str2=substr($sequence,($start3_new-1),($end-($start3_new-1)));
										$str=$str1.$str2;
										$length=length $str;

										my $aa;
										for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
											my $codon=substr ($str,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$aa.=$hash_codon{$codon};
											}else{
												$aa.="X";
												my $j=$i+1;
												#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}

										my $intron_aa;
										if ($mark==2) {
											for (my $i=0;$i<$intron_length;$i+=3){
												my $codon=substr ($string_intron,$i,3);
												$codon=uc $codon;
												if (exists $hash_codon{$codon}){
													$intron_aa.=$hash_codon{$codon};
												}else{
													$intron_aa.="X";
													my $j=$i+1;
													#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
												}
											}
										}

										if ((abs($end3-$start2) < 4000) and ($start2 < $end2) and ($start3 < $end3)) {
											if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/))) {# spliced exon,not 3x or have stop codon
												print $out_annotation "     "."gene"."            ".$start2_new."..".$end."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."join(".$start2_new."..".$end2_new.",".$start3_new."..".$end.")"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa\""."\n";
												$gene_number_seq{$name}++;
											}elsif (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/))) {# joined exon,no intron or no stop codon
												print $out_annotation "     "."gene"."            ".$start2_new."..".$end."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             ".$start2_new."..".$end."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa\""."\n";
												$gene_number_seq{$name}++;
												print $logfile "Warning: $name (positive one-intron PCG) lost intron!\n";
											}

											my ($S1,$E1,$E2,$S3);
											my $string=substr($sequence,($start2_new-1),($end-$start2_new+1));
											my $lengths=length $string;
											my $L1=$end2_new-$start2_new+1;
											my $L2=$end-$start3_new+1;
											for (my $i=0;$i<($length_cp-$lengths);$i+=1){
												my $repeat=substr ($sequence,$i,$lengths);
												if (($repeat eq $string) and (($i+1) ne $start2_new)) {
													$S1=$i+1;
													$E2=$i+$L1-1+1;
													$S3=$i+$lengths-$L2+1;
													$E1=$i+$lengths-1+1;
													print $out_annotation "     "."gene"."            ".$S1."..".$E1."\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "     "."CDS"."             "."join(".$S1."..".$E2.",".$S3."..".$E1.")"."\n" if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)));
													print $out_annotation "     "."CDS"."             ".$S1."..".$E1."\n" if (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/)));
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "                     "."/codon_start=1"."\n";
													print $out_annotation "                     "."/transl_table=11"."\n";
													print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
													#print $out_annotation "                     "."/translation=\"$aa\""."\n";
													$gene_number_seq{$name}++;
												}
											}
											for (my $i=0;$i<($length_cp-$lengths);$i+=1){
												my $repeat=substr ($rev_coms,$i,$lengths);
												if ($repeat eq $string) {
													$S1=$length_cp-$i-1+1;
													$E2=$length_cp-$i-$L1+1;
													$S3=$length_cp-$i-$lengths-1+$L2+1;
													$E1=$length_cp-$i-$lengths+1;
													print $out_annotation "     "."gene"."            "."complement(".$E1."..".$S1.")"."\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "     "."CDS"."             "."complement(join(".$E1."..".$S3.",".$E2."..".$S1."))"."\n" if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)));
													print $out_annotation "     "."CDS"."             "."complement(".$E1."..".$S1.")"."\n" if (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/)));
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "                     "."/codon_start=1"."\n";
													print $out_annotation "                     "."/transl_table=11"."\n";
													print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
													#print $out_annotation "                     "."/translation=\"$aa\""."\n";
													$gene_number_seq{$name}++;
												}
											}
											if (($ticks == 1) and (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)))) {
												print $logfile "Warning: $name (positive one-intron PCG) has non-canonical start codon and intron-exon boundary need to be checked!\n" if ($name ne "rps12+2");
											}elsif ($ticks == 0) {
												print $logfile "Warning: $name (positive one-intron PCG) has non-canonical start codon!\n" if ($name ne "rps12+2");
											}
										}elsif ((abs($end3-$start2) < 4000) and (($start2 > $end2) or ($start3 > $end3))) {
											print $logfile "Warning: $name (positive one-intron PCG) has not been annotated due to two different direction exons!\n";
										}elsif (abs($end3-$start2) >= 4000) {
											print $logfile "Warning: $name (positive one-intron PCG) has not been annotated due to two far exons!\n";
										}
									}else{
										#print $out_annotation "     "."gene"."            ".$start1."..".$end1."\n";
										#print $out_annotation "                     "."/gene=\"$name\""."\n";
										#print $out_annotation "     "."CDS"."             "."join(".$start1."..".$end2.",".$start3."..".$end1.")"."\n";
										#print $out_annotation "                     "."/gene=\"$name\""."\n";
										#print $out_annotation "                     "."/codon_start=1"."\n";
										#print $out_annotation "                     "."/transl_table=11"."\n";
										#print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										##print $out_annotation "                     "."/translation=\"$aa\""."\n";
										#$gene_number_seq{$name}++;
										#print $logfile "Warning: $name (positive one-intron PCG) maybe pseudogene!\n";
										print $logfile "Warning: $name (positive one-intron PCG) has not been annotated!\n";
									}
								}
							}elsif (($name eq "rpl16") or ($name eq "petB") or ($name eq "petD")) {
								if ($length_exon2 < 150) {
									print $logfile "Warning: $name (positive one-intron PCG) has not been annotated due to short length!\n";
								}elsif ($length_exon2 >= 150) {
									if (($length_ref_exon1 % 3==0) and ($length_exon % 3==0) and (!grep {$_=~ /\*/} @aa_exon) and (($forward_start_codon eq "ATG") or ($forward_start_codon eq "GTG")) and (($forward_stop_codon eq "TAA") or ($forward_stop_codon eq "TAG") or ($forward_stop_codon eq "TGA"))){# standard start and stop codon for _gene
										print $out_annotation "     "."gene"."            ".$start2."..".$end3."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."join(".$start2."..".$end2.",".$start3."..".$end3.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
										$gene_number_seq{$name}++;
									}elsif(($length_ref_exon1 % 3!=0) or (grep {$_=~ /\*/} @aa_exon) or (($forward_start_codon ne "ATG") and ($forward_start_codon ne "GTG")) or (($forward_stop_codon ne "TAA") and ($forward_stop_codon ne "TAG") and ($forward_stop_codon ne "TGA"))){# non-standard start or stop codon for _gene
										my $start_left;
										if ($start2>=9000) {
											$start_left=$start2-9000;
										}elsif ($start2<9000){
											$start_left=$start2-(int($start2/3)-1)*3;
										}
										my $start_right=$start2+60;
										my $length_exon1=($end2-$start2+1);
										my $end_right;
										if (($length_cp-$start3)>=9000) {
											$end_right=$start3+9000;
										}elsif (($length_cp-$start3)<9000){
											$end_right=$start3+((int(($length_cp-$start3)/3)-1)*3);
										}
										my $str1=substr($sequence,($start_left-1),($start_right-$start_left+1));
										my $str2=substr($sequence,($start3-1),($end_right-$start3));
										my $length1=length $str1;
										my $length2=length $str2;
										my $coding_1_length=$hash_exon_length{"$name-1_coding"};
										my $coding_1_sequence=$hash_exon_sequence{"$name-1_coding"};
										my $left_start_position;
										for (my $i=0;$i<($length1-3);$i+=1){
											my $exon=substr ($str1,$i,$coding_1_length);
											$exon=lc $exon;
											if ($exon eq $coding_1_sequence) {
												$left_start_position=$i;
											}
										}

										my $aa2;#right range for stop codon
										for (my $k=0;$k<($length2-3);$k+=3){# delete stop codon
											my $codon=substr ($str2,$k,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$aa2.=$hash_codon{$codon};
											}else{
												$aa2.="X";
												my $n=$k+1;
												#print "Bad codon $codon in position $n of gene $name in species $contig!\n";
											}
										}
										my @aa2=split //,$aa2;
										my @star2;
										foreach (0..$#aa2){
											if ($aa2[$_]=~ /\*/){
												push @star2,$_;
											}
										}
										my $right_star_position=shift @star2;

										if ((defined $left_start_position) and (defined $right_star_position)){
											my $start=$start_left+$left_start_position;
											my $end=($start3-1)+($right_star_position * 3+3);

											my ($str1,$str2,$str,$length,$end2_new,$start3_new);
											$end2_new=$start+$coding_1_length-1;
											if ($length_ref_exon1 % 3==0) {
												$start3_new=$start3;
											}elsif ($length_ref_exon1 % 3==1) {
												$start3_new=$start3-2;
											}elsif ($length_ref_exon1 % 3==2) {
												$start3_new=$start3-1;
											}
											$str1=substr($sequence,($start-1),($end2_new-($start-1)));
											$str2=substr($sequence,($start3_new-1),($end-($start3_new-1)));
											$str=$str1.$str2;
											$length=length $str;

											my $aa;
											for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
												my $codon=substr ($str,$i,3);
												$codon=uc $codon;
												if (exists $hash_codon{$codon}){
													$aa.=$hash_codon{$codon};
												}else{
													$aa.="X";
													my $j=$i+1;
													#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
												}
											}
											print $out_annotation "     "."gene"."            ".$start."..".$end."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."join(".$start."..".$end2_new.",".$start3_new."..".$end.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa\""."\n";
											$gene_number_seq{$name}++;
										}elsif ((!defined $left_start_position) and (defined $right_star_position)){
											my $end=($start3-1)+($right_star_position * 3+3);

											my ($str1,$str2,$str,$length,$end2_new,$start3_new);
											$end2_new=$start1+$coding_1_length-1;
											if ($length_ref_exon1 % 3==0) {
												$start3_new=$start3;
											}elsif ($length_ref_exon1 % 3==1) {
												$start3_new=$start3-2;
											}elsif ($length_ref_exon1 % 3==2) {
												$start3_new=$start3-1;
											}
											$str1=substr($sequence,($start1-1),($end2_new-($start1-1)));
											$str2=substr($sequence,($start3_new-1),($end-($start3_new-1)));
											$str=$str1.$str2;
											$length=length $str;

											my $aa;
											for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
												my $codon=substr ($str,$i,3);
												$codon=uc $codon;
												if (exists $hash_codon{$codon}){
													$aa.=$hash_codon{$codon};
												}else{
													$aa.="X";
													my $j=$i+1;
													#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
												}
											}
											#print $out_annotation "     "."gene"."            ".$start2."..".$end."\n";
											print $out_annotation "     "."gene"."            ".$start3."..".$end."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											#print $out_annotation "     "."CDS"."             "."join(".$start2."..".$end2_new.",".$start3_new."..".$end.")"."\n";
											print $out_annotation "     "."CDS"."             ".$start3."..".$end."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (positive one-intron PCG) must be checked due to exon 1 not found!\n";
										}
									}
								}
							}
						}elsif($start2 > $end3){# negative
							my ($str1,$str2);
							$str1=substr($sequence,($end2-1),($start2-$end2+1)) if ($end2 < $start2);
							$str1=substr($sequence,($start2-1),($end2-$start2+1)) if ($end2 > $start2);
							$str2=substr($sequence,($end3-1),($start3-$end3+1)) if ($end3 < $start3);
							$str2=substr($sequence,($start3-1),($end3-$start3+1)) if ($end3 > $start3);
							my $rev_com1=reverse $str1;
							$rev_com1=~ tr/ACGTacgt/TGCAtgca/;
							my $rev_com2=reverse $str2;
							$rev_com2=~ tr/ACGTacgt/TGCAtgca/;
							my $str3=$rev_com1.$rev_com2;
							my $length_exon=length $str3;
							my $length_exon2=length $str2;
							my $reverse_start_codon=substr($str3,0,3);
							my $reverse_stop_codon=substr($str3,($length_exon-3),3);
							my $aa_exon;
							for (my $i=0;$i<($length_exon-3);$i+=3){# delete stop codon
								my $codon=substr ($str3,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa_exon.=$hash_codon{$codon};
								}else{
									$aa_exon.="X";
									my $j=$i+1;
									#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa_exon=split //,$aa_exon;
							if (($name ne "rpl16") and ($name ne "petB") and ($name ne "petD")) {
								if (($length_exon % 3==0) and (!grep {$_=~ /\*/} @aa_exon) and (($reverse_start_codon eq "ATG") or ($reverse_start_codon eq "GTG") or ($name eq "rps12+2")) and (($reverse_stop_codon eq "TAA") or ($reverse_stop_codon eq "TAG") or ($reverse_stop_codon eq "TGA"))){# standard start and stop codon for _gene
									my $seq_string0=substr($sequence,($start2-30),60);
									my $seq_string1=substr($sequence,($end2-(int(($end2-$start3)/3))*3-61),((int(($end2-$start3)/3))*3+120));
									my $seq_string2=substr($sequence,($start3-60),((int(($end2-$start3)/3))*3+120));
									my $length_seq_string0=length $seq_string0;
									my $length_seq_string1=length $seq_string1;
									my $length_seq_string2=length $seq_string2;
									my $rev_seq_string0=reverse $seq_string0;
									$rev_seq_string0=~ tr/ACGTacgt/TGCAtgca/;
									my $rev_seq_string1=reverse $seq_string1;
									$rev_seq_string1=~ tr/ACGTacgt/TGCAtgca/;
									my $rev_seq_string2=reverse $seq_string2;
									$rev_seq_string2=~ tr/ACGTacgt/TGCAtgca/;
									my $aa_seq_string0;
									my $aa_seq_string1;
									my $aa_seq_string2;
									for (my $i=0;$i<$length_seq_string0;$i+=3){# delete stop codon
										my $codon=substr ($rev_seq_string0,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_seq_string0.=$hash_codon{$codon};
										}else{
											$aa_seq_string0.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
									for (my $i=0;$i<$length_seq_string1;$i+=3){# delete stop codon
										my $codon=substr ($rev_seq_string1,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_seq_string1.=$hash_codon{$codon};
										}else{
											$aa_seq_string1.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
									for (my $i=0;$i<$length_seq_string2;$i+=3){# delete stop codon
										my $codon=substr ($rev_seq_string2,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_seq_string2.=$hash_codon{$codon};
										}else{
											$aa_seq_string2.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}

									my ($str1,$str2,$str,$length,$rev_com,$start2_new,$end2_new,$start3_new);
									my ($string_intron,$intron_length,$string_intron_rev_com);
									my $mark=0;
									my $ticks=0;
									if ($length_ref_exon1 % 3==0) {
										for (my $i=$match;$i<$a;$i+=1) {
											my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
												my $repeat=substr ($aa_seq_string1,$j,$match);
												if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
													$end2_new=$end2+60-($j+$match)*3-($i-$match)*3;
													$marks++;
												}
											}
											if (defined $end2_new){
												last;
											}
										}
										if (!defined $end2_new) {
											$end2_new=$end2;
											$ticks=1;
										}
										for (my $i=0;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
												my $repeat=substr ($aa_seq_string2,-$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$start3_new=$start3-60+$j*3+$i*3;
													$marks++;
												}
											}
											if (defined $start3_new){
												last;
											}
										}
										if (!defined $start3_new) {
											$start3_new=$start3;
											$ticks=1;
										}
										$intron_length=$end2_new-$start3_new-1;
										if ($intron_length==0) {
											$mark=1;
										}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
											$string_intron=substr ($sequence,$start3_new,$intron_length);
											$intron_length=length $string_intron;
											$string_intron_rev_com=reverse $string_intron;
											$string_intron_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark=2;
										}
									}elsif ($length_ref_exon1 % 3==1) {
										for (my $i=$match;$i<$a;$i+=1) {
											my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
												my $repeat=substr ($aa_seq_string1,$j,$match);
												if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
													$end2_new=$end2+60-($j+$match)*3-($i-$match)*3-1;
													$marks++;
												}
											}
											if (defined $end2_new){
												last;
											}
										}
										if (!defined $end2_new) {
											$end2_new=$end2-1;
											$ticks=1;
										}
										for (my $i=0;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
												my $repeat=substr ($aa_seq_string2,-$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$start3_new=$start3-60+$j*3+$i*3+2;
													$marks++;
												}
											}
											if (defined $start3_new){
												last;
											}
										}
										if (!defined $start3_new) {
											$start3_new=$start3+2;
											$ticks=1;
										}
										$intron_length=$end2_new-$start3_new-1;
										if ($intron_length==0) {
											$mark=1;
										}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
											$string_intron=substr ($sequence,($start3_new+1),$intron_length);
											$intron_length=length $string_intron;
											$string_intron_rev_com=reverse $string_intron;
											$string_intron_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark=2;
										}
									}elsif ($length_ref_exon1 % 3==2) {
										for (my $i=$match;$i<$a;$i+=1) {
											my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
												my $repeat=substr ($aa_seq_string1,$j,$match);
												if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
													$end2_new=$end2+60-($j+$match)*3-($i-$match)*3-2;
													$marks++;
												}
											}
											if (defined $end2_new){
												last;
											}
										}
										if (!defined $end2_new) {
											$end2_new=$end2-2;
											$ticks=1;
										}
										for (my $i=0;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
												my $repeat=substr ($aa_seq_string2,-$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$start3_new=$start3-60+$j*3+$i*3+1;
													$marks++;
												}
											}
											if (defined $start3_new){
												last;
											}
										}
										if (!defined $start3_new) {
											$start3_new=$start3+1;
											$ticks=1;
										}
										$intron_length=$end2_new-$start3_new-1;
										if ($intron_length==0) {
											$mark=1;
										}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
											$string_intron=substr ($sequence,($start3_new-1),$intron_length);
											$intron_length=length $string_intron;
											$string_intron_rev_com=reverse $string_intron;
											$string_intron_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark=2;
										}
									}

									if ($name eq "rps12+2") {
										for (my $i=0;$i<7;$i+=1) {
											my $match_aa_ref_exon1=substr ($aa_ref_exon1,$i,3);
											my $marks=0;
											for (my $j=0;$j<20;$j+=1) {
												my $repeat=substr ($aa_seq_string0,$j,3);
												if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
													$start2_new=$start2+30-$j*3+$i*3;
													$marks++;
												}
											}
										}
									}elsif ($name ne "rps12+2") {
										$start2_new=$start2;
									}
									$str1=substr($sequence,($end2_new-1),($start2_new-($end2_new-1)));
									$str2=substr($sequence,($end3-1),($start3_new-($end3-1)));
									$str=$str1.$str2;
									$length=length $str;
									$rev_com=reverse $str;
									$rev_com=~ tr/ACGTacgt/TGCAtgca/;

									my $aa;
									for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
										my $codon=substr ($rev_com,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa.=$hash_codon{$codon};
										}else{
											$aa.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}

									my $intron_aa;
									if ($mark==2) {
										for (my $i=0;$i<$intron_length;$i+=3){
											my $codon=substr ($string_intron_rev_com,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$intron_aa.=$hash_codon{$codon};
											}else{
												$intron_aa.="X";
												my $j=$i+1;
												#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}
									}

									if ((abs($start2-$end3) < 4000) and ($start2 > $end2) and ($start3 > $end3) and ($start2 > $end3)) {
										if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/))) {# spliced exon,not 3x or have stop codon
											print $out_annotation "     "."gene"."            "."complement(".$end3."..".$start2_new.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."complement(join(".$end3."..".$start3_new.",".$end2_new."..".$start2_new."))"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
											$gene_number_seq{$name}++;
										}elsif (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/))) {# joined exon,no intron or no stop codon
											print $out_annotation "     "."gene"."            "."complement(".$end3."..".$start2_new.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."complement(".$end3."..".$start2_new.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (negative one-intron PCG) lost intron!\n";
										}

										my ($S1,$E1,$E2,$S3);
										my $string=substr($sequence,($end3-1),($start2_new-$end3+1));
										my $lengths=length $string;
										my $L1=$start2_new-$end2_new+1;
										my $L2=$start3_new-$end3+1;
										for (my $i=0;$i<($length_cp-$lengths);$i+=1){
											my $repeat=substr ($sequence,$i,$lengths);
											if (($repeat eq $string) and (($i+1) ne $end3)) {
												$E1=$i+1;
												$S3=$i+$L2-1+1;
												$E2=$i+$lengths-$L1+1;
												$S1=$i+$lengths-1+1;
												print $out_annotation "     "."gene"."            "."complement(".$E1."..".$S1.")"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."complement(join(".$E1."..".$S3.",".$E2."..".$S1."))"."\n" if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)));
												print $out_annotation "     "."CDS"."             "."complement(".$E1."..".$S1.")"."\n" if (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/)));
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
												$gene_number_seq{$name}++;
											}
										}
										for (my $i=0;$i<($length_cp-$lengths);$i+=1){
											my $repeat=substr ($rev_coms,$i,$lengths);
											if ($repeat eq $string) {
												$E1=$length_cp-$i-1+1;
												$S3=$length_cp-$i-$L2+1;
												$E2=$length_cp-$i-$lengths-1+$L1+1;
												$S1=$length_cp-$i-$lengths+1;
												print $out_annotation "     "."gene"."            ".$S1."..".$E1."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."join(".$S1."..".$E2.",".$S3."..".$E1.")"."\n" if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)));
												print $out_annotation "     "."CDS"."             ".$S1."..".$E1."\n" if (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/)));
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
												$gene_number_seq{$name}++;
											}
										}
										if (($ticks == 1) and (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)))) {
											print $logfile "Warning: $name (negative one-intron PCG) need to be checked!\n";
										}
									}elsif ((abs($start2-$end3) < 4000) and (($start2 < $end2) or ($start3 < $end3) or ($start2 < $end3))) {
										print $logfile "Warning: $name (negative one-intron PCG) has not been annotated due to two different direction exons!\n";
									}elsif (abs($start2-$end3) >= 4000) {
										print $logfile "Warning: $name (negative one-intron PCG) has not been annotated due to two far exons!\n";
									}
								}elsif((grep {$_=~ /\*/} @aa_exon) or (($reverse_start_codon ne "ATG") and ($reverse_start_codon ne "GTG")) or (($reverse_stop_codon ne "TAA") and ($reverse_stop_codon ne "TAG") and ($reverse_stop_codon ne "TGA"))){# non-standard start or stop codon for _gene
									my $start_left=$start2-60;
									my $start_right;
									if (($length_cp-$start2)>=9000) {
										$start_right=$start2+9000;
									}elsif (($length_cp-$start2)<9000){
										$start_right=$start2+(int(($length_cp-$start2)/3)-1)*3;
									}
									my $length_exon1=($start2-$end2+1);
									my $end_left;
									if ($start3>=9000) {
										$end_left=$start3-9000;
									}elsif ($start3<9000){
										$end_left=$start3-(int($start3/3)-1)*3;
									}
									my $str1=substr($sequence,($start_left),($start_right-$start_left));
									my $str2=substr($sequence,($end_left),($start3-$end_left));
									my $length1=length $str1;
									my $length2=length $str2;
									my $rev_com1=reverse $str1;
									$rev_com1=~ tr/ACGTacgt/TGCAtgca/;
									my $rev_com2=reverse $str2;
									$rev_com2=~ tr/ACGTacgt/TGCAtgca/;
									my $aa1;#left range for start codon
									for (my $i=0;$i<($length1-3);$i+=3){# delete stop codon
										my $codon=substr ($rev_com1,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa1.=$hash_codon{$codon};
										}else{
											$aa1.="X";
											my $m=$i+1;
											#print "Bad codon $codon in position $m of gene $name in species $contig!\n";
										}
									}
									my @aa1=split //,$aa1;
									my (@star1,@start1);
									foreach (0..$#aa1){
										if (($aa1[$_]=~ /\*/) and ($_ < 3000)){
											push @star1,$_;
										}
										if (($aa1[$_]=~ /M/) and ($_ <= 3060)){
											push @start1,$_;
										}
									}
									my $right_star_position=pop @star1;
									until ($right_star_position<(3000+$length_exon1/3)){
										$right_star_position=pop @star1;
									}
									my @start2;
									foreach (@start1){
										push @start2,$_ if ((defined $right_star_position) and ($_ > $right_star_position))
									}
									my $right_start_position=$start2[0];
									#my $right_start_position=pop @start2;
									my $aa2;#right range for stop codon
									for (my $k=0;$k<($length2-3);$k+=3){# delete stop codon
										my $codon=substr ($rev_com2,$k,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa2.=$hash_codon{$codon};
										}else{
											$aa2.="X";
											my $n=$k+1;
											#print "Bad codon $codon in position $n of gene $name in species $contig!\n";
										}
									}
									my @aa2=split //,$aa2;
									my @star2;
									foreach (0..$#aa2){
										if ($aa2[$_]=~ /\*/){
											push @star2,$_;
										}
									}
									my $left_star_position=shift @star2;

									my $seq_string0=substr($sequence,($start2-30),60);
									my $seq_string1=substr($sequence,($end2-(int(($end2-$start3)/3))*3-61),((int(($end2-$start3)/3))*3+120));
									my $seq_string2=substr($sequence,($start3-60),((int(($end2-$start3)/3))*3+120));
									my $length_seq_string0=length $seq_string0;
									my $length_seq_string1=length $seq_string1;
									my $length_seq_string2=length $seq_string2;
									my $rev_seq_string0=reverse $seq_string0;
									$rev_seq_string0=~ tr/ACGTacgt/TGCAtgca/;
									my $rev_seq_string1=reverse $seq_string1;
									$rev_seq_string1=~ tr/ACGTacgt/TGCAtgca/;
									my $rev_seq_string2=reverse $seq_string2;
									$rev_seq_string2=~ tr/ACGTacgt/TGCAtgca/;
									my $aa_seq_string0;
									my $aa_seq_string1;
									my $aa_seq_string2;
									for (my $i=0;$i<$length_seq_string0;$i+=3){# delete stop codon
										my $codon=substr ($rev_seq_string0,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_seq_string0.=$hash_codon{$codon};
										}else{
											$aa_seq_string0.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
									for (my $i=0;$i<$length_seq_string1;$i+=3){# delete stop codon
										my $codon=substr ($rev_seq_string1,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_seq_string1.=$hash_codon{$codon};
										}else{
											$aa_seq_string1.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
									for (my $i=0;$i<$length_seq_string2;$i+=3){# delete stop codon
										my $codon=substr ($rev_seq_string2,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa_seq_string2.=$hash_codon{$codon};
										}else{
											$aa_seq_string2.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}

									if ((defined $right_start_position) and (defined $left_star_position)){
										my $start=$start_right-$right_start_position * 3;
										my $end=($start3+1)-($left_star_position * 3+3);

										my ($str1,$str2,$str,$length,$rev_com,$start_new,$end2_new,$start3_new);
										my ($string_intron,$intron_length,$string_intron_rev_com);
										my $mark=0;
										my $ticks=0;
										if ($length_ref_exon1 % 3==0) {
											for (my $i=$match;$i<$a;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string1,$j,$match);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$end2_new=$end2+60-($j+$match)*3-($i-$match)*3;
														$marks++;
													}
												}
												if (defined $end2_new){
													last;
												}
											}
											if (!defined $end2_new) {
												$end2_new=$end2;
												$ticks=1;
											}
											for (my $i=0;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string2,-$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$start3_new=$start3-60+$j*3+$i*3;
														$marks++;
													}
												}
												if (defined $start3_new){
													last;
												}
											}
											if (!defined $start3_new) {
												$start3_new=$start3;
												$ticks=1;
											}
											$intron_length=$end2_new-$start3_new-1;
											if ($intron_length==0) {
												$mark=1;
											}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
												$string_intron=substr ($sequence,$start3_new,$intron_length);
												$intron_length=length $string_intron;
												$string_intron_rev_com=reverse $string_intron;
												$string_intron_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark=2;
											}
										}elsif ($length_ref_exon1 % 3==1) {
											for (my $i=$match;$i<$a;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string1,$j,$match);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$end2_new=$end2+60-($j+$match)*3-($i-$match)*3-1;
														$marks++;
													}
												}
												if (defined $end2_new){
													last;
												}
											}
											if (!defined $end2_new) {
												$end2_new=$end2-1;
												$ticks=1;
											}
											for (my $i=0;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string2,-$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$start3_new=$start3-60+$j*3+$i*3+2;
														$marks++;
													}
												}
												if (defined $start3_new){
													last;
												}
											}
											if (!defined $start3_new) {
												$start3_new=$start3+2;
												$ticks=1;
											}
											$intron_length=$end2_new-$start3_new-1;
											if ($intron_length==0) {
												$mark=1;
											}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
												$string_intron=substr ($sequence,($start3_new+1),$intron_length);
												$intron_length=length $string_intron;
												$string_intron_rev_com=reverse $string_intron;
												$string_intron_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark=2;
											}
										}elsif ($length_ref_exon1 % 3==2) {
											for (my $i=$match;$i<$a;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string1,$j,$match);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$end2_new=$end2+60-($j+$match)*3-($i-$match)*3-2;
														$marks++;
													}
												}
												if (defined $end2_new){
													last;
												}
											}
											if (!defined $end2_new) {
												$end2_new=$end2-2;
												$ticks=1;
											}
											for (my $i=0;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string2,-$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$start3_new=$start3-60+$j*3+$i*3+1;
														$marks++;
													}
												}
												if (defined $start3_new){
													last;
												}
											}
											if (!defined $start3_new) {
												$start3_new=$start3+1;
												$ticks=1;
											}
											$intron_length=$end2_new-$start3_new-1;
											if ($intron_length==0) {
												$mark=1;
											}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
												$string_intron=substr ($sequence,($start3_new-1),$intron_length);
												$intron_length=length $string_intron;
												$string_intron_rev_com=reverse $string_intron;
												$string_intron_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark=2;
											}
										}

										if ($name eq "rps12+2") {
											for (my $i=0;$i<7;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,$i,3);
												my $marks=0;
												for (my $j=0;$j<20;$j+=1) {
													my $repeat=substr ($aa_seq_string0,$j,3);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$start_new=$start2+30-$j*3+$i*3;
														$marks++;
													}
												}
											}
										}elsif ($name ne "rps12+2") {
											$start_new=$start;
										}
										$str1=substr($sequence,($end2_new-1),($start_new-($end2_new-1)));
										$str2=substr($sequence,($end-1),($start3_new-($end-1)));
										$str=$str1.$str2;
										$length=length $str;
										$rev_com=reverse $str;
										$rev_com=~ tr/ACGTacgt/TGCAtgca/;

										my $aa;
										for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
											my $codon=substr ($rev_com,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$aa.=$hash_codon{$codon};
											}else{
												$aa.="X";
												my $j=$i+1;
												#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}

										my $intron_aa;
										if ($mark==2) {
											for (my $i=0;$i<$intron_length;$i+=3){
												my $codon=substr ($string_intron_rev_com,$i,3);
												$codon=uc $codon;
												if (exists $hash_codon{$codon}){
													$intron_aa.=$hash_codon{$codon};
												}else{
													$intron_aa.="X";
													my $j=$i+1;
													#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
												}
											}
										}

										if ((abs($start2-$end3) < 4000) and ($start2 > $end2) and ($start3 > $end3) and ($start2 > $end3)) {
											if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/))) {# spliced exon,not 3x or have stop codon
												print $out_annotation "     "."gene"."            "."complement(".$end."..".$start_new.")\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."complement(join(".$end."..".$start3_new.",".$end2_new."..".$start_new."))"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa\""."\n";
												$gene_number_seq{$name}++;
											}elsif (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/))) {# joined exon,no intron or no stop codon
												print $out_annotation "     "."gene"."            "."complement(".$end."..".$start_new.")\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."complement(".$end."..".$start_new.")"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa\""."\n";
												$gene_number_seq{$name}++;
												print $logfile "Warning: $name (negative one-intron PCG) lost intron!\n";
											}

											my ($S1,$E1,$E2,$S3);
											my $string=substr($sequence,($end-1),($start_new-$end+1));
											my $lengths=length $string;
											my $L1=$start_new-$end2_new+1;
											my $L2=$start3_new-$end+1;
											for (my $i=0;$i<($length_cp-$lengths);$i+=1){
												my $repeat=substr ($sequence,$i,$lengths);
												if (($repeat eq $string) and (($i+1) ne $end)) {
													$E1=$i+1;
													$S3=$i+$L2-1+1;
													$E2=$i+$lengths-$L1+1;
													$S1=$i+$lengths-1+1;
													print $out_annotation "     "."gene"."            "."complement(".$E1."..".$S1.")\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "     "."CDS"."             "."complement(join(".$E1."..".$S3.",".$E2."..".$S1."))"."\n" if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)));
													print $out_annotation "     "."CDS"."             "."complement(".$E1."..".$S1.")"."\n" if (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/)));
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "                     "."/codon_start=1"."\n";
													print $out_annotation "                     "."/transl_table=11"."\n";
													print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
													#print $out_annotation "                     "."/translation=\"$aa\""."\n";
													$gene_number_seq{$name}++;
												}
											}
											for (my $i=0;$i<($length_cp-$lengths);$i+=1){
												my $repeat=substr ($rev_coms,$i,$lengths);
												if ($repeat eq $string) {
													$E1=$length_cp-$i-1+1;
													$S3=$length_cp-$i-$L2+1;
													$E2=$length_cp-$i-$lengths-1+$L1+1;
													$S1=$length_cp-$i-$lengths+1;
													print $out_annotation "     "."gene"."            ".$S1."..".$E1."\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "     "."CDS"."             "."join(".$S1."..".$E2.",".$S3."..".$E1.")"."\n" if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)));
													print $out_annotation "     "."CDS"."             ".$S1."..".$E1."\n" if (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/)));
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "                     "."/codon_start=1"."\n";
													print $out_annotation "                     "."/transl_table=11"."\n";
													print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
													#print $out_annotation "                     "."/translation=\"$aa\""."\n";
													$gene_number_seq{$name}++;
												}
											}
											if (($ticks == 1) and (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)))) {
												print $logfile "Warning: $name (negative one-intron PCG) intron-exon boundary need to be checked!\n";
											}
										}elsif ((abs($start2-$end3) < 4000) and (($start2 < $end2) or ($start3 < $end3) or ($start2 < $end3))) {
											print $logfile "Warning: $name (negative one-intron PCG) has not been annotated due to two different direction exons!\n";
										}elsif (abs($start2-$end3) >= 4000) {
											print $logfile "Warning: $name (negative one-intron PCG) has not been annotated due to two far exons!\n";
										}
									}elsif ((!defined $right_start_position) and (defined $left_star_position)){
										my $end=($start3+1)-($left_star_position * 3+3);

										my ($str1,$str2,$str,$length,$rev_com,$start2_new,$end2_new,$start3_new);
										my ($string_intron,$intron_length,$string_intron_rev_com);
										my $mark=0;
										my $ticks=1;
										if ($length_ref_exon1 % 3==0) {
											for (my $i=$match;$i<$a;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string1,$j,$match);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$end2_new=$end2+60-($j+$match)*3-($i-$match)*3;
														$marks++;
													}
												}
												if (defined $end2_new){
													last;
												}
											}
											if (!defined $end2_new) {
												$end2_new=$end2;
												$ticks=1;
											}
											for (my $i=0;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string2,-$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$start3_new=$start3-60+$j*3+$i*3;
														$marks++;
													}
												}
												if (defined $start3_new){
													last;
												}
											}
											if (!defined $start3_new) {
												$start3_new=$start3;
												$ticks=1;
											}
											$intron_length=$end2_new-$start3_new-1;
											if ($intron_length==0) {
												$mark=1;
											}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
												$string_intron=substr ($sequence,$start3_new,$intron_length);
												$intron_length=length $string_intron;
												$string_intron_rev_com=reverse $string_intron;
												$string_intron_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark=2;
											}
										}elsif ($length_ref_exon1 % 3==1) {
											for (my $i=$match;$i<$a;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string1,$j,$match);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$end2_new=$end2+60-($j+$match)*3-($i-$match)*3-1;
														$marks++;
													}
												}
												if (defined $end2_new){
													last;
												}
											}
											if (!defined $end2_new) {
												$end2_new=$end2-1;
												$ticks=1;
											}
											for (my $i=0;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string2,-$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$start3_new=$start3-60+$j*3+$i*3+2;
														$marks++;
													}
												}
												if (defined $start3_new){
													last;
												}
											}
											if (!defined $start3_new) {
												$start3_new=$start3+2;
												$ticks=1;
											}
											$intron_length=$end2_new-$start3_new-1;
											if ($intron_length==0) {
												$mark=1;
											}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
												$string_intron=substr ($sequence,($start3_new+1),$intron_length);
												$intron_length=length $string_intron;
												$string_intron_rev_com=reverse $string_intron;
												$string_intron_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark=2;
											}
										}elsif ($length_ref_exon1 % 3==2) {
											for (my $i=$match;$i<$a;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<(($length_seq_string1)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string1,$j,$match);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$end2_new=$end2+60-($j+$match)*3-($i-$match)*3-2;
														$marks++;
													}
												}
												if (defined $end2_new){
													last;
												}
											}
											if (!defined $end2_new) {
												$end2_new=$end2-2;
												$ticks=1;
											}
											for (my $i=0;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<(($length_seq_string2)/3);$j+=1) {
													my $repeat=substr ($aa_seq_string2,-$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$start3_new=$start3-60+$j*3+$i*3+1;
														$marks++;
													}
												}
												if (defined $start3_new){
													last;
												}
											}
											if (!defined $start3_new) {
												$start3_new=$start3+1;
												$ticks=1;
											}
											$intron_length=$end2_new-$start3_new-1;
											if ($intron_length==0) {
												$mark=1;
											}elsif (($intron_length!=0) and ($intron_length % 3==0)) {
												$string_intron=substr ($sequence,($start3_new-1),$intron_length);
												$intron_length=length $string_intron;
												$string_intron_rev_com=reverse $string_intron;
												$string_intron_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark=2;
											}
										}

										if ($name eq "rps12+2") {
											for (my $i=0;$i<7;$i+=1) {
												my $match_aa_ref_exon1=substr ($aa_ref_exon1,$i,3);
												my $marks=0;
												for (my $j=0;$j<20;$j+=1) {
													my $repeat=substr ($aa_seq_string0,$j,3);
													if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
														$start2_new=$start2+30-$j*3+$i*3;
														$marks++;
													}
												}
											}
										}elsif ($name ne "rps12+2") {
											$start2_new=$start2;
										}
										$str1=substr($sequence,($end2_new-1),($start2_new-($end2_new-1)));
										$str2=substr($sequence,($end-1),($start3_new-($end-1)));
										$str=$str1.$str2;
										$length=length $str;
										$rev_com=reverse $str;
										$rev_com=~ tr/ACGTacgt/TGCAtgca/;

										my $aa;
										for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
											my $codon=substr ($rev_com,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$aa.=$hash_codon{$codon};
											}else{
												$aa.="X";
												my $j=$i+1;
												#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}

										my $intron_aa;
										if ($mark==2) {
											for (my $i=0;$i<$intron_length;$i+=3){
												my $codon=substr ($string_intron_rev_com,$i,3);
												$codon=uc $codon;
												if (exists $hash_codon{$codon}){
													$intron_aa.=$hash_codon{$codon};
												}else{
													$intron_aa.="X";
													my $j=$i+1;
													#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
												}
											}
										}

										if ((abs($start2-$end3) < 4000) and ($start2 > $end2) and ($start3 > $end3) and ($start2 > $end3)) {
											if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/))) {# spliced exon,not 3x or have stop codon
												print $out_annotation "     "."gene"."            "."complement(".$end."..".$start2_new.")\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."complement(join(".$end."..".$start3_new.",".$end2_new."..".$start2_new."))"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa\""."\n";
												$gene_number_seq{$name}++;
											}elsif (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/))) {# joined exon,no intron or no stop codon
												print $out_annotation "     "."gene"."            "."complement(".$end."..".$start2_new.")\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."complement(".$end."..".$start2_new.")\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa\""."\n";
												$gene_number_seq{$name}++;
												print $logfile "Warning: $name (negative one-intron PCG) lost intron!\n";
											}

											my ($S1,$E1,$E2,$S3);
											my $string=substr($sequence,($end-1),($start2_new-$end+1));
											my $lengths=length $string;
											my $L1=$start2_new-$end2_new+1;
											my $L2=$start3_new-$end+1;
											my $rev_coms=reverse $sequence;
											$rev_coms=~ tr/ACGTacgt/TGCAtgca/;
											for (my $i=0;$i<($length_cp-$lengths);$i+=1){
												my $repeat=substr ($sequence,$i,$lengths);
												if (($repeat eq $string) and (($i+1) ne $end)) {
													$E1=$i+1;
													$S3=$i+$L2-1+1;
													$E2=$i+$lengths-$L1+1;
													$S1=$i+$lengths-1+1;
													print $out_annotation "     "."gene"."            "."complement(".$E1."..".$S1.")\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "     "."CDS"."             "."complement(join(".$E1."..".$S3.",".$E2."..".$S1."))"."\n" if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)));
													print $out_annotation "     "."CDS"."             "."complement(".$E1."..".$S1.")"."\n" if (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/)));
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "                     "."/codon_start=1"."\n";
													print $out_annotation "                     "."/transl_table=11"."\n";
													print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
													#print $out_annotation "                     "."/translation=\"$aa\""."\n";
													$gene_number_seq{$name}++;
												}
											}
											for (my $i=0;$i<($length_cp-$lengths);$i+=1){
												my $repeat=substr ($rev_coms,$i,$lengths);
												if ($repeat eq $string) {
													$E1=$length_cp-$i-1+1;
													$S3=$length_cp-$i-$L2+1;
													$E2=$length_cp-$i-$lengths-1+$L1+1;
													$S1=$length_cp-$i-$lengths+1;
													print $out_annotation "     "."gene"."            ".$S1."..".$E1."\n";
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "     "."CDS"."             "."join(".$S1."..".$E2.",".$S3."..".$E1.")"."\n" if (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)));
													print $out_annotation "     "."CDS"."             ".$S1."..".$E1."\n" if (($mark==1) or ((defined $intron_aa) and ($intron_aa!~ /\*/)));
													print $out_annotation "                     "."/gene=\"$name\""."\n";
													print $out_annotation "                     "."/codon_start=1"."\n";
													print $out_annotation "                     "."/transl_table=11"."\n";
													print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
													#print $out_annotation "                     "."/translation=\"$aa\""."\n";
													$gene_number_seq{$name}++;
												}
											}
											if (($ticks == 1) and (($mark==0) or ((defined $intron_aa) and ($intron_aa=~ /\*/)))) {
												print $logfile "Warning: $name (negative one-intron PCG) has non-canonical start codon and intron-exon boundary need to be checked!\n" if ($name ne "rps12+2");
											}elsif ($ticks == 0) {
												print $logfile "Warning: $name (negative one-intron PCG) has non-canonical start codon!\n" if ($name ne "rps12+2");
											}
										}elsif ((abs($start2-$end3) < 4000) and (($start2 < $end2) or ($start3 < $end3) or ($start2 < $end3))) {
											print $logfile "Warning: $name (negative one-intron PCG) has not been annotated due to two different direction exons!\n";
										}elsif (abs($start2-$end3) >= 4000) {
											print $logfile "Warning: $name (negative one-intron PCG) has not been annotated due to two far exons!\n";
										}
									}else{
										#print $out_annotation "     "."gene"."            "."complement(".$end1."..".$start1.")\n";
										#print $out_annotation "                     "."/gene=\"$name\""."\n";
										#print $out_annotation "     "."CDS"."             "."complement(join(".$end1."..".$start3.",".$end2."..".$start1."))"."\n";
										#print $out_annotation "                     "."/gene=\"$name\""."\n";
										#print $out_annotation "                     "."/codon_start=1"."\n";
										#print $out_annotation "                     "."/transl_table=11"."\n";
										#print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										##print $out_annotation "                     "."/translation=\"$aa\""."\n";
										#$gene_number_seq{$name}++;
										#print $logfile "Warning: $name (negative one-intron PCG) maybe pseudogene!\n";
										print $logfile "Warning: $name (negative one-intron PCG) has not been annotated!\n";
									}
								}
							}elsif (($name eq "rpl16") or ($name eq "petB") or ($name eq "petD")) {
								if ($length_exon2 < 150) {
									print $logfile "Warning: $name (negative one-intron PCG) has not been annotated due to short length!\n";
								}elsif ($length_exon2 >= 150) {
									if (($length_ref_exon1 % 3==0) and ($length_exon % 3==0) and (!grep {$_=~ /\*/} @aa_exon) and (($reverse_start_codon eq "ATG") or ($reverse_start_codon eq "GTG")) and (($reverse_stop_codon eq "TAA") or ($reverse_stop_codon eq "TAG") or ($reverse_stop_codon eq "TGA"))){# standard start and stop codon for _gene
										print $out_annotation "     "."gene"."            "."complement(".$end3."..".$start2.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."complement(join(".$end3."..".$start3.",".$end2."..".$start2."))"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
										$gene_number_seq{$name}++;
									}elsif(($length_ref_exon1 % 3!=0) or (grep {$_=~ /\*/} @aa_exon) or (($reverse_start_codon ne "ATG") and ($reverse_start_codon ne "GTG")) or (($reverse_stop_codon ne "TAA") and ($reverse_stop_codon ne "TAG") and ($reverse_stop_codon ne "TGA"))){# non-standard start or stop codon for _gene
										my $start_left=$start2-60;
										my $start_right;
										if (($length_cp-$start2)>=9000) {
											$start_right=$start2+9000;
										}elsif (($length_cp-$start2)<9000){
											$start_right=$start2+(int(($length_cp-$start2)/3)-1)*3;
										}
										my $length_exon1=($start2-$end2+1);
										my $end_left;
										if ($start3>=9000) {
											$end_left=$start3-9000;
										}elsif ($start3<9000){
											$end_left=$start3-(int($start3/3)-1)*3;
										}
										my $str1=substr($sequence,($start_left),($start_right-$start_left));
										my $str2=substr($sequence,($end_left),($start3-$end_left));
										my $length1=length $str1;
										my $length2=length $str2;
										my $rev_com1=reverse $str1;
										$rev_com1=~ tr/ACGTacgt/TGCAtgca/;
										my $rev_com2=reverse $str2;
										$rev_com2=~ tr/ACGTacgt/TGCAtgca/;
										my $coding_1_length=$hash_exon_length{"$name-1_coding"};
										my $coding_1_sequence=$hash_exon_sequence{"$name-1_coding"};
										my $right_start_position;
										for (my $i=0;$i<($length1-3);$i+=1){
											my $exon=substr ($rev_com1,$i,$coding_1_length);
											$exon=lc $exon;
											if ($exon eq $coding_1_sequence) {
												$right_start_position=$i;
											}
										}

										my $aa2;#right range for stop codon
										for (my $k=0;$k<($length2-3);$k+=3){# delete stop codon
											my $codon=substr ($rev_com2,$k,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$aa2.=$hash_codon{$codon};
											}else{
												$aa2.="X";
												my $n=$k+1;
												#print "Bad codon $codon in position $n of gene $name in species $contig!\n";
											}
										}
										my @aa2=split //,$aa2;
										my @star2;
										foreach (0..$#aa2){
											if ($aa2[$_]=~ /\*/){
												push @star2,$_;
											}
										}
										my $left_star_position=shift @star2;
										if ((defined $right_start_position) and (defined $left_star_position)){
											my $start=$start_right-$right_start_position;
											my $end=($start3+1)-($left_star_position * 3+3);

											my ($str1,$str2,$str,$length,$rev_com,$end2_new,$start3_new);
											$end2_new=$start-$coding_1_length+1;
											if ($length_ref_exon1 % 3==0) {
												$start3_new=$start3;
											}elsif ($length_ref_exon1 % 3==1) {
												$start3_new=$start3+2;
											}elsif ($length_ref_exon1 % 3==2) {
												$start3_new=$start3+1;
											}
											$str1=substr($sequence,($end2_new-1),($start-($end2_new-1)));
											$str2=substr($sequence,($end-1),($start3_new-($end-1)));
											$str=$str1.$str2;
											$length=length $str;
											$rev_com=reverse $str;
											$rev_com=~ tr/ACGTacgt/TGCAtgca/;

											my $aa;
											for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
												my $codon=substr ($rev_com,$i,3);
												$codon=uc $codon;
												if (exists $hash_codon{$codon}){
													$aa.=$hash_codon{$codon};
												}else{
													$aa.="X";
													my $j=$i+1;
													#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
												}
											}
											print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."complement(join(".$end."..".$start3_new.",".$end2_new."..".$start."))"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa\""."\n";
											$gene_number_seq{$name}++;
										}elsif ((!defined $right_start_position) and (defined $left_star_position)){
											my $end=($start3+1)-($left_star_position * 3+3);

											my ($str1,$str2,$str,$length,$rev_com,$end2_new,$start3_new);
											$end2_new=$start1-$coding_1_length+1;
											if ($length_ref_exon1 % 3==0) {
												$start3_new=$start3;
											}elsif ($length_ref_exon1 % 3==1) {
												$start3_new=$start3+2;
											}elsif ($length_ref_exon1 % 3==2) {
												$start3_new=$start3+1;
											}
											$str1=substr($sequence,($end2_new-1),($start1-($end2_new-1)));
											$str2=substr($sequence,($end-1),($start3_new-($end-1)));
											$str=$str1.$str2;
											$length=length $str;
											$rev_com=reverse $str;
											$rev_com=~ tr/ACGTacgt/TGCAtgca/;

											my $aa;
											for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
												my $codon=substr ($rev_com,$i,3);
												$codon=uc $codon;
												if (exists $hash_codon{$codon}){
													$aa.=$hash_codon{$codon};
												}else{
													$aa.="X";
													my $j=$i+1;
													#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
												}
											}
											#print $out_annotation "     "."gene"."            "."complement(".$end."..".$start2.")\n";
											print $out_annotation "     "."gene"."            "."complement(".$end."..".$start3.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											#print $out_annotation "     "."CDS"."             "."complement(join(".$end."..".$start3_new.",".$end2_new."..".$start2."))"."\n";
											print $out_annotation "     "."CDS"."             "."complement(".$end."..".$start3.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (negative one-intron PCG) must be checked due to exon 1 not found!\n";
										}
									}
								}
							}
						}
					}

					############################################################
					## two-intron protein-coding-gene
					############################################################
					if((($start1 == $start2) and (defined $end4) and ($end1 == $end4)) and ((defined $length_ref_exon1) and (defined $length_ref_exon2) and (defined $length_ref_exon3) and (defined $sequence_ref_exon1) and (defined $sequence_ref_exon2) and (defined $sequence_ref_exon3))){# identical PCG boundary for _gene and -1_coding, -3_coding
						if ($start3 < $end3){# positive
							my $str=substr($sequence,($start1-1),($end1-$start1+1));
							my $length_3=length $str;
							my $forward_start_codon=substr($str,0,3);
							my $forward_stop_codon=substr($str,($length_3-3),3);
							my $str1=substr($sequence,($start1-1),($end2-$start1+1));
							my $str2=substr($sequence,($start3-1),($end3-$start3+1));
							my $str3=substr($sequence,($start4-1),($end1-$start4+1));
							my $str4=$str1.$str2.$str3;
							my $length_exon=length $str4;
							my $aa_exon;
							for (my $i=0;$i<($length_exon-3);$i+=3){# delete stop codon
								my $codon=substr ($str4,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa_exon.=$hash_codon{$codon};
								}else{
									$aa_exon.="X";
									my $j=$i+1;
									#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa_exon=split //,$aa_exon;
							if (($length_exon % 3==0) and (!grep {$_=~ /\*/} @aa_exon) and (($forward_start_codon eq "ATG") or ($forward_start_codon eq "GTG")) and (($forward_stop_codon eq "TAA") or ($forward_stop_codon eq "TAG") or ($forward_stop_codon eq "TGA"))){# standard start and stop codon for _gene
								my $seq_string1=substr($sequence,($end2-60),120);
								my $seq_string2=substr($sequence,($start3-61),120);
								my $seq_string3=substr($sequence,($end3-60),120);
								my $seq_string4=substr($sequence,($start4-61),120);
								my $aa_seq_string1;
								my $aa_seq_string2;
								my $aa_seq_string3;
								my $aa_seq_string4;
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($seq_string1,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string1.=$hash_codon{$codon};
									}else{
										$aa_seq_string1.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($seq_string2,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string2.=$hash_codon{$codon};
									}else{
										$aa_seq_string2.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($seq_string3,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string3.=$hash_codon{$codon};
									}else{
										$aa_seq_string3.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($seq_string4,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string4.=$hash_codon{$codon};
									}else{
										$aa_seq_string4.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}

								my ($str1,$str2,$str3,$str,$length,$end2_new,$start3_new,$end3_new,$start4_new);
								my ($string_intron1,$string_intron2,$intron1_length,$intron2_length);
								my $mark1=0;
								my $mark2=0;
								if ($length_ref_exon1 % 3==0) {
									for (my $i=$match;$i<$a;$i+=1) {
										my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
										my $marks=0;
										for (my $j=0;$j<40;$j+=1) {
											my $repeat=substr ($aa_seq_string1,$j,$match);
											if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
												$end2_new=$end2-60+($j+$match)*3+($i-$match)*3;
												$marks++;
											}
										}
										if (defined $end2_new){
											last;
										}
									}
									if (!defined $end2_new) {
										$end2_new=$end2;
									}
									for (my $i=0;$i<$b;$i+=1) {
										my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
										my $marks=0;
										for (my $j=$match;$j<40;$j+=1) {
											my $repeat=substr ($aa_seq_string2,-$j,$match);
											if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
												$start3_new=$start3+60-$j*3-$i*3;
												$marks++;
											}
										}
										if (defined $start3_new){
											last;
										}
									}
									if (!defined $start3_new) {
										$start3_new=$start3;
									}
									$intron1_length=$start3_new-$end2_new-1;
									if ($intron1_length==0) {
										$mark1=1;
									}
									if (($intron1_length!=0) and ($intron1_length % 3==0)) {
										$string_intron1=substr ($sequence,$end2_new,$intron1_length);
										$intron1_length=length $string_intron1;
										$mark1=2;
									}

									if ($length_ref_exon2 % 3==0) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3-60+($j+$match)*3+($i-$match)*3;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4+60-$j*3-$i*3;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4;
										}
										$intron2_length=$start4_new-$end3_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,$end3_new,$intron2_length);
											$intron2_length=length $string_intron2;
											$mark2=2;
										}
									}elsif ($length_ref_exon2 % 3==1) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3-60+($j+$match)*3+($i-$match)*3+1;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3+1;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4+60-$j*3-$i*3-2;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4-2;
										}
										$intron2_length=$start4_new-$end3_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,($end3_new-1),$intron2_length);
											$intron2_length=length $string_intron2;
											$mark2=2;
										}
									}elsif ($length_ref_exon2 % 3==2) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3-60+($j+$match)*3+($i-$match)*3+2;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3+2;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4+60-$j*3-$i*3-1;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4-1;
										}
										$intron2_length=$start4_new-$end3_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,($end3_new+1),$intron2_length);
											$intron2_length=length $string_intron2;
											$mark2=2;
										}
									}
								}elsif ($length_ref_exon1 % 3==1) {
									for (my $i=$match;$i<$a;$i+=1) {
										my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
										my $marks=0;
										for (my $j=0;$j<40;$j+=1) {
											my $repeat=substr ($aa_seq_string1,$j,$match);
											if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
												$end2_new=$end2-60+($j+$match)*3+($i-$match)*3+1;
												$marks++;
											}
										}
										if (defined $end2_new){
											last;
										}
									}
									if (!defined $end2_new) {
										$end2_new=$end2+1;
									}
									for (my $i=0;$i<$b;$i+=1) {
										my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
										my $marks=0;
										for (my $j=$match;$j<40;$j+=1) {
											my $repeat=substr ($aa_seq_string2,-$j,$match);
											if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
												$start3_new=$start3+60-$j*3-$i*3-2;
												$marks++;
											}
										}
										if (defined $start3_new){
											last;
										}
									}
									if (!defined $start3_new) {
										$start3_new=$start3-2;
									}
									$intron1_length=$start3_new-$end2_new-1;
									if ($intron1_length==0) {
										$mark1=1;
									}
									if (($intron1_length!=0) and ($intron1_length % 3==0)) {
										$string_intron1=substr ($sequence,($end2_new-1),$intron1_length);
										$intron1_length=length $string_intron1;
										$mark1=2;
									}

									if ($length_ref_exon2 % 3==0) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3-60+($j+$match)*3+($i-$match)*3+1;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3+1;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4+60-$j*3-$i*3-2;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4-2;
										}
										$intron2_length=$start4_new-$end3_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,($end3_new-1),$intron2_length);
											$intron2_length=length $string_intron2;
											$mark2=2;
										}
									}elsif ($length_ref_exon2 % 3==1) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3-60+($j+$match)*3+($i-$match)*3+2;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3+2;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4+60-$j*3-$i*3-1;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4-1;
										}
										$intron2_length=$start4_new-$end3_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,($end3_new+1),$intron2_length);
											$intron2_length=length $string_intron2;
											$mark2=2;
										}
									}elsif ($length_ref_exon2 % 3==2) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3-60+($j+$match)*3+($i-$match)*3;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4+60-$j*3-$i*3;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4;
										}
										$intron2_length=$start4_new-$end3_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,$end3_new,$intron2_length);
											$intron2_length=length $string_intron2;
											$mark2=2;
										}
									}
								}elsif ($length_ref_exon1 % 3==2) {
									for (my $i=$match;$i<$a;$i+=1) {
										my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
										my $marks=0;
										for (my $j=0;$j<40;$j+=1) {
											my $repeat=substr ($aa_seq_string1,$j,$match);
											if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
												$end2_new=$end2-60+($j+$match)*3+($i-$match)*3+2;
												$marks++;
											}
										}
										if (defined $end2_new){
											last;
										}
									}
									if (!defined $end2_new) {
										$end2_new=$end2+2;
									}
									for (my $i=0;$i<$b;$i+=1) {
										my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
										my $marks=0;
										for (my $j=$match;$j<40;$j+=1) {
											my $repeat=substr ($aa_seq_string2,-$j,$match);
											if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
												$start3_new=$start3+60-$j*3-$i*3-1;
												$marks++;
											}
										}
										if (defined $start3_new){
											last;
										}
									}
									if (!defined $start3_new) {
										$start3_new=$start3-1;
									}
									$intron1_length=$start3_new-$end2_new-1;
									if ($intron1_length==0) {
										$mark1=1;
									}elsif (($intron1_length!=0) and ($intron1_length % 3==0)) {
										$string_intron1=substr ($sequence,($end2_new+1),$intron1_length);
										$intron1_length=length $string_intron1;
										$mark1=2;
									}

									if ($length_ref_exon2 % 3==0) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3-60+($j+$match)*3+($i-$match)*3+2;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3+2;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4+60-$j*3-$i*3-1;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4-1;
										}
										$intron2_length=$start4_new-$end3_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,($end3_new+1),$intron2_length);
											$intron2_length=length $string_intron2;
											$mark2=2;
										}
									}elsif ($length_ref_exon2 % 3==1) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3-60+($j+$match)*3+($i-$match)*3;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4+60-$j*3-$i*3;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4;
										}
										$intron2_length=$start4_new-$end3_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,$end3_new,$intron2_length);
											$intron2_length=length $string_intron2;
											$mark2=2;
										}
									}elsif ($length_ref_exon2 % 3==2) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3-60+($j+$match)*3+($i-$match)*3+1;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3+1;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4+60-$j*3-$i*3-2;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4-2;
										}
										$intron2_length=$start4_new-$end3_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,($end3_new-1),$intron2_length);
											$intron2_length=length $string_intron2;
											$mark2=2;
										}
									}
								}

								$str1=substr($sequence,($start1-1),($end2_new-($start1-1)));
								$str2=substr($sequence,($start3_new-1),($end3_new-($start3_new-1)));
								$str3=substr($sequence,($start4_new-1),($end1-($start4_new-1)));
								$str=$str1.$str2.$str3;
								$length=length $str;

								my $aa;
								for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
									my $codon=substr ($str,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa.=$hash_codon{$codon};
									}else{
										$aa.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}

								my $intron1_aa;
								if ($mark1==2) {
									for (my $i=0;$i<$intron1_length;$i+=3){
										my $codon=substr ($string_intron1,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$intron1_aa.=$hash_codon{$codon};
										}else{
											$intron1_aa.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
								}
								my $intron2_aa;
								if ($mark2==2) {
									for (my $i=0;$i<$intron2_length;$i+=3){
										my $codon=substr ($string_intron2,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$intron2_aa.=$hash_codon{$codon};
										}else{
											$intron2_aa.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
								}

								#if ((defined $end2_new) and (defined $start3_new) and (defined $end3_new) and (defined $start4_new)) {
								if (abs($end4-$start2) < 3000) {
									if (($start3 > $end2) and ($start4 > $end3)) {
										if ((($mark1==0) or ((defined $intron1_aa) and ($intron1_aa=~ /\*/))) and (($mark2==0) or ((defined $intron2_aa) and ($intron2_aa=~ /\*/))) and (($intron1_length > 30) and ($intron2_length > 30))) {# (two introns) spliced exon,not 3x or have stop codon
											print $out_annotation "     "."gene"."            ".$start1."..".$end1."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."join(".$start1."..".$end2_new.",".$start3_new."..".$end3_new.",".$start4_new."..".$end1.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
											$gene_number_seq{$name}++;
										}elsif ((($mark1==1) or ((defined $intron1_aa) and ($intron1_aa!~ /\*/))) and (($mark2==1) or ((defined $intron2_aa) and ($intron2_aa!~ /\*/))) and (($intron1_length < 30) and ($intron2_length < 30))) {# (no intron) joined exon,no intron or no stop codon
											print $out_annotation "     "."gene"."            ".$start1."..".$end1."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             ".$start1."..".$end1."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (positive two-intron PCG) lost two introns!\n";
										}elsif ((($mark1==0) or ((defined $intron1_aa) and ($intron1_aa=~ /\*/))) and (($mark2==1) or ((defined $intron2_aa) and ($intron2_aa!~ /\*/))) and (($intron1_length > 30) and ($intron2_length < 30))) {# (loss intron2)
											print $out_annotation "     "."gene"."            ".$start1."..".$end1."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."join(".$start1."..".$end2_new.",".$start3_new."..".$end1.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (positive two-intron PCG) lost intron 2!\n";
										}elsif ((($mark1==1) or ((defined $intron1_aa) and ($intron1_aa!~ /\*/))) and (($mark2==0) or ((defined $intron2_aa) and ($intron2_aa=~ /\*/))) and (($intron1_length < 30) and ($intron2_length > 30))) {# (loss intron1)
											print $out_annotation "     "."gene"."            ".$start1."..".$end1."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."join(".$start1."..".$end3_new.",".$start4_new."..".$end1.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (positive two-intron PCG) lost intron 1!\n";
										}
									}elsif (($start3 < $end2) or ($start4 < $end3)) {
										print $logfile "Warning: $name (positive two-intron PCG) has not been annotated!\n";
									}
								#}elsif ((!defined $end2_new) or (!defined $start3_new) or (!defined $end3_new) or (!defined $start4_new)) {
								}elsif (abs($end4-$start2) > 3000) {
									print $out_annotation "     "."gene"."            ".$start3_new."..".$end3_new."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."CDS"."             ".$start3_new."..".$end3_new."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/codon_start=1"."\n";
									print $out_annotation "                     "."/transl_table=11"."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
									$gene_number_seq{$name}++;
									print $logfile "Warning: $name (positive two-intron PCG) must be checked!\n";
								}
							}elsif((grep {$_=~ /\*/} @aa_exon) or (($forward_start_codon ne "ATG") and ($forward_start_codon ne "GTG")) or (($forward_stop_codon ne "TAA") and ($forward_stop_codon ne "TAG") and ($forward_stop_codon ne "TGA"))){# non-standard start or stop codon for _gene
								my $start_left;
								if ($start1>=9000) {
									$start_left=$start1-9000;
								}elsif ($start1<9000){
									$start_left=$start1-int($start1/3)*3;
								}
								my $start_right=$start1+60;
								my $length_exon1=($end2-$start1+1);
								my $length_exon2=($end3-$start3+1);
								my $length_exon12=$length_exon1+$length_exon2;
								my $end_right;
								if (($length_cp-$start4)>=9002) {
									if ($length_exon12 % 3==0){
										$end_right=$start4+9000;
									}
									if ($length_exon12 % 3==1){
										$end_right=$start4+9001;
									}
									if ($length_exon12 % 3==2){
										$end_right=$start4+9002;
									}
								}elsif (($length_cp-$start4)<9002){
									if ($length_exon12 % 3==0){
										$end_right=$start4+(int(($length_cp-$start4)/3)*3);
									}
									if ($length_exon12 % 3==1){
										$end_right=$start4+(int(($length_cp-$start4)/3)*3+1);
									}
									if ($length_exon12 % 3==2){
										$end_right=$start4+(int(($length_cp-$start4)/3)*3+2);
									}
								}
								my $str1=substr($sequence,($start_left-1),($start_right-$start_left+1));
								my $str2=substr($sequence,($start4-1),($end_right-$start4+1));
								my $length1=length $str1;
								my $length2=length $str2;
								my $aa1;#left range for start codon
								for (my $i=0;$i<($length1-3);$i+=3){# delete stop codon
									my $codon=substr ($str1,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa1.=$hash_codon{$codon};
									}else{
										$aa1.="X";
										my $m=$i+1;
										#print "Bad codon $codon in position $m of gene $name in species $contig!\n";
									}
								}
								my @aa1=split //,$aa1;
								my (@star1,@start1);
								foreach (0..$#aa1){
									if (($aa1[$_]=~ /\*/) and ($_ < 3000)){
										push @star1,$_;
									}
									if (($aa1[$_]=~ /M/) and ($_ <= 3060)){
										push @start1,$_;
										}
								}
								my $left_star_position=pop @star1;
								until ($left_star_position<(3000+$length_exon1/3)){
									$left_star_position=pop @star1;
								}
								my @start2;
								foreach (@start1){
									push @start2,$_ if ((defined $left_star_position) and ($_ > $left_star_position))
								}
								#my $left_start_position=$start2[0];
								my $left_start_position=pop @start2;
								my $aa2;#right range for stop codon
								my $j;
								if ($length_exon1 % 3==0){
									$j=0;
								}
								if ($length_exon1 % 3==1){
									$j=2;
								}
								if ($length_exon1 % 3==2){
									$j=1;
								}
								for (my $k=$j;$k<($length2-3);$k+=3){# delete stop codon
									my $codon=substr ($str2,$k,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa2.=$hash_codon{$codon};
									}else{
										$aa2.="X";
										my $n=$k+1;
										#print "Bad codon $codon in position $n of gene $name in species $contig!\n";
									}
								}
								my @aa2=split //,$aa2;
								my @star2;
								foreach (0..$#aa2){
									if ($aa2[$_]=~ /\*/){
										push @star2,$_;
									}
								}
								my $right_star_position=shift @star2;

								my $seq_string1=substr($sequence,($end2-60),120);
								my $seq_string2=substr($sequence,($start3-61),120);
								my $seq_string3=substr($sequence,($end3-60),120);
								my $seq_string4=substr($sequence,($start4-61),120);
								my $aa_seq_string1;
								my $aa_seq_string2;
								my $aa_seq_string3;
								my $aa_seq_string4;
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($seq_string1,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string1.=$hash_codon{$codon};
									}else{
										$aa_seq_string1.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($seq_string2,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string2.=$hash_codon{$codon};
									}else{
										$aa_seq_string2.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($seq_string3,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string3.=$hash_codon{$codon};
									}else{
										$aa_seq_string3.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($seq_string4,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string4.=$hash_codon{$codon};
									}else{
										$aa_seq_string4.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}

								if ((defined $left_start_position) and (defined $right_star_position)){
									my $start=$start_left+$left_start_position * 3;
									my $end;
									if ($length_exon12 % 3==0){
										$end=($start4-1)+($right_star_position * 3+3);
									}
									if ($length_exon12 % 3==1){
										$end=($start4+1)+($right_star_position * 3+3);
									}
									if ($length_exon12 % 3==2){
										$end=($start4)+($right_star_position * 3+3);
									}

									my ($str1,$str2,$str3,$str,$length,$end2_new,$start3_new,$end3_new,$start4_new);
									my ($string_intron1,$string_intron2,$intron1_length,$intron2_length);
									my $mark1=0;
									my $mark2=0;
									if ($length_ref_exon1 % 3==0) {
										for (my $i=$match;$i<$a;$i+=1) {
											my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string1,$j,$match);
												if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
													$end2_new=$end2-60+($j+$match)*3+($i-$match)*3;
													$marks++;
												}
											}
											if (defined $end2_new){
												last;
											}
										}
										if (!defined $end2_new) {
											$end2_new=$end2;
										}
										for (my $i=0;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string2,-$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$start3_new=$start3+60-$j*3-$i*3;
													$marks++;
												}
											}
											if (defined $start3_new){
												last;
											}
										}
										if (!defined $start3_new) {
											$start3_new=$start3;
										}
										$intron1_length=$start3_new-$end2_new-1;
										if ($intron1_length==0) {
											$mark1=1;
										}
										if (($intron1_length!=0) and ($intron1_length % 3==0)) {
											$string_intron1=substr ($sequence,$end2_new,$intron1_length);
											$intron1_length=length $string_intron1;
											$mark1=2;
										}

										if ($length_ref_exon2 % 3==0) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3-60+($j+$match)*3+($i-$match)*3;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4+60-$j*3-$i*3;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4;
											}
											$intron2_length=$start4_new-$end3_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,$end3_new,$intron2_length);
												$intron2_length=length $string_intron2;
												$mark2=2;
											}
										}elsif ($length_ref_exon2 % 3==1) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3-60+($j+$match)*3+($i-$match)*3+1;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3+1;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4+60-$j*3-$i*3-2;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4-2;
											}
											$intron2_length=$start4_new-$end3_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,($end3_new-1),$intron2_length);
												$intron2_length=length $string_intron2;
												$mark2=2;
											}
										}elsif ($length_ref_exon2 % 3==2) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3-60+($j+$match)*3+($i-$match)*3+2;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3+2;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4+60-$j*3-$i*3-1;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4-1;
											}
											$intron2_length=$start4_new-$end3_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,($end3_new+1),$intron2_length);
												$intron2_length=length $string_intron2;
												$mark2=2;
											}
										}
									}elsif ($length_ref_exon1 % 3==1) {
										for (my $i=$match;$i<$a;$i+=1) {
											my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string1,$j,$match);
												if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
													$end2_new=$end2-60+($j+$match)*3+($i-$match)*3+1;
													$marks++;
												}
											}
											if (defined $end2_new){
												last;
											}
										}
										if (!defined $end2_new) {
											$end2_new=$end2+1;
										}
										for (my $i=0;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string2,-$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$start3_new=$start3+60-$j*3-$i*3-2;
													$marks++;
												}
											}
											if (defined $start3_new){
												last;
											}
										}
										if (!defined $start3_new) {
											$start3_new=$start3-2;
										}
										$intron1_length=$start3_new-$end2_new-1;
										if ($intron1_length==0) {
											$mark1=1;
										}
										if (($intron1_length!=0) and ($intron1_length % 3==0)) {
											$string_intron1=substr ($sequence,($end2_new-1),$intron1_length);
											$intron1_length=length $string_intron1;
											$mark1=2;
										}

										if ($length_ref_exon2 % 3==0) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3-60+($j+$match)*3+($i-$match)*3+1;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3+1;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4+60-$j*3-$i*3-2;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4-2;
											}
											$intron2_length=$start4_new-$end3_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,($end3_new-1),$intron2_length);
												$intron2_length=length $string_intron2;
												$mark2=2;
											}
										}elsif ($length_ref_exon2 % 3==1) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3-60+($j+$match)*3+($i-$match)*3+2;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3+2;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4+60-$j*3-$i*3-1;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4-1;
											}
											$intron2_length=$start4_new-$end3_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,($end3_new+1),$intron2_length);
												$intron2_length=length $string_intron2;
												$mark2=2;
											}
										}elsif ($length_ref_exon2 % 3==2) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3-60+($j+$match)*3+($i-$match)*3;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4+60-$j*3-$i*3;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4;
											}
											$intron2_length=$start4_new-$end3_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,$end3_new,$intron2_length);
												$intron2_length=length $string_intron2;
												$mark2=2;
											}
										}
									}elsif ($length_ref_exon1 % 3==2) {
										for (my $i=$match;$i<$a;$i+=1) {
											my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string1,$j,$match);
												if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
													$end2_new=$end2-60+($j+$match)*3+($i-$match)*3+2;
													$marks++;
												}
											}
											if (defined $end2_new){
												last;
											}
										}
										if (!defined $end2_new) {
											$end2_new=$end2+2;
										}
										for (my $i=0;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string2,-$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$start3_new=$start3+60-$j*3-$i*3-1;
													$marks++;
												}
											}
											if (defined $start3_new){
												last;
											}
										}
										if (!defined $start3_new) {
											$start3_new=$start3-1;
										}
										$intron1_length=$start3_new-$end2_new-1;
										if ($intron1_length==0) {
											$mark1=1;
										}elsif (($intron1_length!=0) and ($intron1_length % 3==0)) {
											$string_intron1=substr ($sequence,($end2_new+1),$intron1_length);
											$intron1_length=length $string_intron1;
											$mark1=2;
										}

										if ($length_ref_exon2 % 3==0) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3-60+($j+$match)*3+($i-$match)*3+2;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3+2;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4+60-$j*3-$i*3-1;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4-1;
											}
											$intron2_length=$start4_new-$end3_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,($end3_new+1),$intron2_length);
												$intron2_length=length $string_intron2;
												$mark2=2;
											}
										}elsif ($length_ref_exon2 % 3==1) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3-60+($j+$match)*3+($i-$match)*3;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4+60-$j*3-$i*3;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4;
											}
											$intron2_length=$start4_new-$end3_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,$end3_new,$intron2_length);
												$intron2_length=length $string_intron2;
												$mark2=2;
											}
										}elsif ($length_ref_exon2 % 3==2) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3-60+($j+$match)*3+($i-$match)*3+1;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3+1;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4+60-$j*3-$i*3-2;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4-2;
											}
											$intron2_length=$start4_new-$end3_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,($end3_new-1),$intron2_length);
												$intron2_length=length $string_intron2;
												$mark2=2;
											}
										}
									}

									$str1=substr($sequence,($start-1),($end2_new-($start-1)));
									$str2=substr($sequence,($start3_new-1),($end3_new-($start3_new-1)));
									$str3=substr($sequence,($start4_new-1),($end-($start4_new-1)));
									$str=$str1.$str2.$str3;
									$length=length $str;
	
									my $aa;
									for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
										my $codon=substr ($str,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa.=$hash_codon{$codon};
										}else{
											$aa.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
	
									my $intron1_aa;
									if ($mark1==2) {
										for (my $i=0;$i<$intron1_length;$i+=3){
											my $codon=substr ($string_intron1,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$intron1_aa.=$hash_codon{$codon};
											}else{
												$intron1_aa.="X";
												my $j=$i+1;
												#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}
									}
									my $intron2_aa;
									if ($mark2==2) {
										for (my $i=0;$i<$intron2_length;$i+=3){
											my $codon=substr ($string_intron2,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$intron2_aa.=$hash_codon{$codon};
											}else{
												$intron2_aa.="X";
												my $j=$i+1;
												#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}
									}

									#if ((defined $end2_new) and (defined $start3_new) and (defined $end3_new) and (defined $start4_new)) {
									if (abs($end4-$start2) < 3000) {
										if (($start3 > $end2) and ($start4 > $end3)) {
											if ((($mark1==0) or ((defined $intron1_aa) and ($intron1_aa=~ /\*/))) and (($mark2==0) or ((defined $intron2_aa) and ($intron2_aa=~ /\*/))) and (($intron1_length > 30) and ($intron2_length > 30))) {# (two introns) spliced exon,not 3x or have stop codon
												print $out_annotation "     "."gene"."            ".$start."..".$end."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."join(".$start."..".$end2_new.",".$start3_new."..".$end3_new.",".$start4_new."..".$end.")"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
												$gene_number_seq{$name}++;
											}elsif ((($mark1==1) or ((defined $intron1_aa) and ($intron1_aa!~ /\*/))) and (($mark2==1) or ((defined $intron2_aa) and ($intron2_aa!~ /\*/))) and (($intron1_length < 30) and ($intron2_length < 30))) {# (no intron) joined exon,no intron or no stop codon
												print $out_annotation "     "."gene"."            ".$start."..".$end."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             ".$start."..".$end."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
												$gene_number_seq{$name}++;
												print $logfile "Warning: $name (positive two-intron PCG) lost two introns!\n";
											}elsif ((($mark1==0) or ((defined $intron1_aa) and ($intron1_aa=~ /\*/))) and (($mark2==1) or ((defined $intron2_aa) and ($intron2_aa!~ /\*/))) and (($intron1_length > 30) and ($intron2_length < 30))) {# (loss intron2)
												print $out_annotation "     "."gene"."            ".$start."..".$end."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."join(".$start."..".$end2_new.",".$start3_new."..".$end.")"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
												$gene_number_seq{$name}++;
												print $logfile "Warning: $name (positive two-intron PCG) lost intron 2!\n";
											}elsif ((($mark1==1) or ((defined $intron1_aa) and ($intron1_aa!~ /\*/))) and (($mark2==0) or ((defined $intron2_aa) and ($intron2_aa=~ /\*/))) and (($intron1_length < 30) and ($intron2_length > 30))) {# (loss intron1)
												print $out_annotation "     "."gene"."            ".$start."..".$end."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."join(".$start."..".$end3_new.",".$start4_new."..".$end.")"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
												$gene_number_seq{$name}++;
												print $logfile "Warning: $name (positive two-intron PCG) lost intron 1!\n";
											}
										}elsif (($start3 < $end2) or ($start4 < $end3)) {
											print $logfile "Warning: $name (positive two-intron PCG) has not been annotated!\n";
										}
									#}elsif ((!defined $end2_new) or (!defined $start3_new) or (!defined $end3_new) or (!defined $start4_new)) {
									}elsif (abs($end4-$start2) > 3000) {
										print $out_annotation "     "."gene"."            ".$start3_new."..".$end3_new."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             ".$start3_new."..".$end3_new."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
										$gene_number_seq{$name}++;
										print $logfile "Warning: $name (positive two-intron PCG) must be checked!\n";
									}
								}elsif ((!defined $left_start_position) and (defined $right_star_position)){
									my $end;
									if ($length_exon12 % 3==0){
										$end=($start4-1)+($left_star_position * 3+3);
									}
									if ($length_exon12 % 3==1){
										$end=($start4+1)+($left_star_position * 3+3);
									}
									if ($length_exon12 % 3==2){
										$end=($start4)+($left_star_position * 3+3);
									}
									if (abs($end4-$start2) < 3000) {
										if (($start3 > $end2) and ($start4 > $end3)) {
											print $out_annotation "     "."gene"."            ".$start1."..".$end."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."join(".$start1."..".$end2.",".$start3."..".$end3.",".$start4."..".$end.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (positive two-intron PCG) has non-canonical start codon!\n";
										}elsif (($start3 < $end2) or ($start4 < $end3)) {
											print $logfile "Warning: $name (positive two-intron PCG) has not been annotated!\n";
										}
									}elsif (abs($end4-$start2) > 3000) {
										print $out_annotation "     "."gene"."            ".$start3."..".$end3."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             ".$start3."..".$end3."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
										$gene_number_seq{$name}++;
										print $logfile "Warning: $name (positive two-intron PCG) must be checked!\n";
									}
								}
							}else{
								print "The exon length of $name (positive two-intron PCG) is not an intergral multiple of 3!\n";
							}
						}elsif($start3 > $end3){# negative
							my $str=substr($sequence,($end1-1),($start1-$end1+1));
							my $rev_com=reverse $str;
							$rev_com=~ tr/ACGTacgt/TGCAtgca/;
							my $length_3=length $str;
							my $reverse_start_codon=substr($rev_com,0,3);
							my $reverse_stop_codon=substr($rev_com,($length_3-3),3);
							my $str1=substr($sequence,($end2-1),($start1-$end2+1));
							my $str2=substr($sequence,($end3-1),($start3-$end3+1));
							my $str3=substr($sequence,($end1-1),($start4-$end1+1));
							my $rev_com1=reverse $str1;
							$rev_com1=~ tr/ACGTacgt/TGCAtgca/;
							my $rev_com2=reverse $str2;
							$rev_com2=~ tr/ACGTacgt/TGCAtgca/;
							my $rev_com3=reverse $str3;
							$rev_com3=~ tr/ACGTacgt/TGCAtgca/;
							my $str4=$rev_com1.$rev_com2.$rev_com3;
							my $length_exon=length $str4;
							my $aa_exon;
							for (my $i=0;$i<($length_exon-3);$i+=3){# delete stop codon
								my $codon=substr ($str4,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa_exon.=$hash_codon{$codon};
								}else{
									$aa_exon.="X";
									my $j=$i+1;
									#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa_exon=split //,$aa_exon;
							if (($length_exon % 3==0) and (!grep {$_=~ /\*/} @aa_exon) and (($reverse_start_codon eq "ATG") or ($reverse_start_codon eq "GTG")) and (($reverse_stop_codon eq "TAA") or ($reverse_stop_codon eq "TAG") or ($reverse_stop_codon eq "TGA"))){# standard start and stop codon
								my $seq_string1=substr($sequence,($end2-61),120);
								my $seq_string2=substr($sequence,($start3-60),120);
								my $seq_string3=substr($sequence,($end3-61),120);
								my $seq_string4=substr($sequence,($start4-60),120);
								my $rev_seq_string1=reverse $seq_string1;
								$rev_seq_string1=~ tr/ACGTacgt/TGCAtgca/;
								my $rev_seq_string2=reverse $seq_string2;
								$rev_seq_string2=~ tr/ACGTacgt/TGCAtgca/;
								my $rev_seq_string3=reverse $seq_string3;
								$rev_seq_string3=~ tr/ACGTacgt/TGCAtgca/;
								my $rev_seq_string4=reverse $seq_string4;
								$rev_seq_string4=~ tr/ACGTacgt/TGCAtgca/;
								my $aa_seq_string1;
								my $aa_seq_string2;
								my $aa_seq_string3;
								my $aa_seq_string4;
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($rev_seq_string1,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string1.=$hash_codon{$codon};
									}else{
										$aa_seq_string1.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($rev_seq_string2,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string2.=$hash_codon{$codon};
									}else{
										$aa_seq_string2.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($rev_seq_string3,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string3.=$hash_codon{$codon};
									}else{
										$aa_seq_string3.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($rev_seq_string4,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string4.=$hash_codon{$codon};
									}else{
										$aa_seq_string4.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}

								my ($str1,$str2,$str3,$str,$rev_com,$length,$end2_new,$start3_new,$end3_new,$start4_new);
								my ($string_intron1,$string_intron2,$intron1_length,$intron2_length,$string_intron1_rev_com,$string_intron2_rev_com);
								my $mark1=0;
								my $mark2=0;
								if ($length_ref_exon1 % 3==0) {
									for (my $i=$match;$i<$a;$i+=1) {
										my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
										my $marks=0;
										for (my $j=0;$j<40;$j+=1) {
											my $repeat=substr ($aa_seq_string1,$j,$match);
											if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
												$end2_new=$end2+60-($j+$match)*3-($i-$match)*3;
												$marks++;
											}
										}
										if (defined $end2_new){
											last;
										}
									}
									if (!defined $end2_new) {
										$end2_new=$end2;
									}
									for (my $i=0;$i<$b;$i+=1) {
										my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
										my $marks=0;
										for (my $j=$match;$j<40;$j+=1) {
											my $repeat=substr ($aa_seq_string2,-$j,$match);
											if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
												$start3_new=$start3-60+$j*3+$i*3;
												$marks++;
											}
										}
										if (defined $start3_new){
											last;
										}
									}
									if (!defined $start3_new) {
										$start3_new=$start3;
									}
									$intron1_length=$end2_new-$start3_new-1;
									if ($intron1_length==0) {
										$mark1=1;
									}
									if (($intron1_length!=0) and ($intron1_length % 3==0)) {
										$string_intron1=substr ($sequence,$start3_new,$intron1_length);
										$intron1_length=length $string_intron1;
										$string_intron1_rev_com=reverse $string_intron1;
										$string_intron1_rev_com=~ tr/ACGTacgt/TGCAtgca/;
										$mark1=2;
									}

									if ($length_ref_exon2 % 3==0) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3+60-($j+$match)*3-($i-$match)*3;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4-60+$j*3+$i*3;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4;
										}
										$intron2_length=$end3_new-$start4_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,$start4_new,$intron2_length);
											$intron2_length=length $string_intron2;
											$string_intron2_rev_com=reverse $string_intron2;
											$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark2=2;
										}
									}elsif ($length_ref_exon2 % 3==1) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3+60-($j+$match)*3-($i-$match)*3-1;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3-1;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4-60+$j*3+$i*3+2;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4+2;
										}
										$intron2_length=$end3_new-$start4_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,($start4_new+1),$intron2_length);
											$intron2_length=length $string_intron2;
											$string_intron2_rev_com=reverse $string_intron2;
											$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark2=2;
										}
									}elsif ($length_ref_exon2 % 3==2) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3+60-($j+$match)*3-($i-$match)*3-2;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3-2;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4-60+$j*3+$i*3+1;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4+1;
										}
										$intron2_length=$end3_new-$start4_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,($start4_new-1),$intron2_length);
											$intron2_length=length $string_intron2;
											$string_intron2_rev_com=reverse $string_intron2;
											$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark2=2;
										}
									}
								}elsif ($length_ref_exon1 % 3==1) {
									for (my $i=$match;$i<$a;$i+=1) {
										my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
										my $marks=0;
										for (my $j=0;$j<40;$j+=1) {
											my $repeat=substr ($aa_seq_string1,$j,$match);
											if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
												$end2_new=$end2+60-($j+$match)*3-($i-$match)*3-1;
												$marks++;
											}
										}
										if (defined $end2_new){
											last;
										}
									}
									if (!defined $end2_new) {
										$end2_new=$end2-1;
									}
									for (my $i=0;$i<$b;$i+=1) {
										my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
										my $marks=0;
										for (my $j=$match;$j<40;$j+=1) {
											my $repeat=substr ($aa_seq_string2,-$j,$match);
											if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
												$start3_new=$start3-60+$j*3+$i*3+2;
												$marks++;
											}
										}
										if (defined $start3_new){
											last;
										}
									}
									if (!defined $start3_new) {
										$start3_new=$start3+2;
									}
									$intron1_length=$end2_new-$start3_new-1;
									if ($intron1_length==0) {
										$mark1=1;
									}
									if (($intron1_length!=0) and ($intron1_length % 3==0)) {
										$string_intron1=substr ($sequence,($start3_new+1),$intron1_length);
										$intron1_length=length $string_intron1;
										$string_intron1_rev_com=reverse $string_intron1;
										$string_intron1_rev_com=~ tr/ACGTacgt/TGCAtgca/;
										$mark1=2;
									}

									if ($length_ref_exon2 % 3==0) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3+60-($j+$match)*3-($i-$match)*3-1;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3-1;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4-60+$j*3+$i*3+2;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4+2;
										}
										$intron2_length=$end3_new-$start4_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,($start4_new+1),$intron2_length);
											$intron2_length=length $string_intron2;
											$string_intron2_rev_com=reverse $string_intron2;
											$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark2=2;
										}
									}elsif ($length_ref_exon2 % 3==1) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3+60-($j+$match)*3-($i-$match)*3-2;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3-2;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4-60+$j*3+$i*3+1;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4+1;
										}
										$intron2_length=$end3_new-$start4_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,($start4_new-1),$intron2_length);
											$intron2_length=length $string_intron2;
											$string_intron2_rev_com=reverse $string_intron2;
											$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark2=2;
										}
									}elsif ($length_ref_exon2 % 3==2) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3+60-($j+$match)*3-($i-$match)*3;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4-60+$j*3+$i*3;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4;
										}
										$intron2_length=$end3_new-$start4_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,$start4_new,$intron2_length);
											$intron2_length=length $string_intron2;
											$string_intron2_rev_com=reverse $string_intron2;
											$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark2=2;
										}
									}
								}elsif ($length_ref_exon1 % 3==2) {
									for (my $i=$match;$i<$a;$i+=1) {
										my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
										my $marks=0;
										for (my $j=0;$j<40;$j+=1) {
											my $repeat=substr ($aa_seq_string1,$j,$match);
											if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
												$end2_new=$end2+60-($j+$match)*3-($i-$match)*3-2;
												$marks++;
											}
										}
										if (defined $end2_new){
											last;
										}
									}
									if (!defined $end2_new) {
										$end2_new=$end2-2;
									}
									for (my $i=0;$i<$b;$i+=1) {
										my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
										my $marks=0;
										for (my $j=$match;$j<40;$j+=1) {
											my $repeat=substr ($aa_seq_string2,-$j,$match);
											if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
												$start3_new=$start3-60+$j*3+$i*3+1;
												$marks++;
											}
										}
										if (defined $start3_new){
											last;
										}
									}
									if (!defined $start3_new) {
										$start3_new=$start3+1;
									}
									$intron1_length=$end2_new-$start3_new-1;
									if ($intron1_length==0) {
										$mark1=1;
									}elsif (($intron1_length!=0) and ($intron1_length % 3==0)) {
										$string_intron1=substr ($sequence,($start3_new-1),$intron1_length);
										$intron1_length=length $string_intron1;
										$string_intron1_rev_com=reverse $string_intron1;
										$string_intron1_rev_com=~ tr/ACGTacgt/TGCAtgca/;
										$mark1=2;
									}

									if ($length_ref_exon2 % 3==0) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3+60-($j+$match)*3-($i-$match)*3-2;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3-2;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4-60+$j*3+$i*3+1;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4+1;
										}
										$intron2_length=$end3_new-$start4_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,($start4_new-1),$intron2_length);
											$intron2_length=length $string_intron2;
											$string_intron2_rev_com=reverse $string_intron2;
											$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark2=2;
										}
									}elsif ($length_ref_exon2 % 3==1) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3+60-($j+$match)*3-($i-$match)*3;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4-60+$j*3+$i*3;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4;
										}
										$intron2_length=$end3_new-$start4_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,$start4_new,$intron2_length);
											$intron2_length=length $string_intron2;
											$string_intron2_rev_com=reverse $string_intron2;
											$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark2=2;
										}
									}elsif ($length_ref_exon2 % 3==2) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3+60-($j+$match)*3-($i-$match)*3-1;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3-1;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4-60+$j*3+$i*3+2;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4+2;
										}
										$intron2_length=$end3_new-$start4_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,($start4_new+1),$intron2_length);
											$intron2_length=length $string_intron2;
											$string_intron2_rev_com=reverse $string_intron2;
											$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark2=2;
										}
									}
								}

								$str1=substr($sequence,($end2_new-1),($start1-($end2_new-1)));
								$str2=substr($sequence,($end3_new-1),($start3_new-($end3_new-1)));
								$str3=substr($sequence,($end1-1),($start4_new-($end1-1)));
								$str=$str1.$str2.$str3;
								$length=length $str;
								$rev_com=reverse $str;
								$rev_com=~ tr/ACGTacgt/TGCAtgca/;

								my $aa;
								for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
									my $codon=substr ($rev_com,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa.=$hash_codon{$codon};
									}else{
										$aa.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}

								my $intron1_aa;
								if ($mark1==2) {
									for (my $i=0;$i<$intron1_length;$i+=3){
										my $codon=substr ($string_intron1_rev_com,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$intron1_aa.=$hash_codon{$codon};
										}else{
											$intron1_aa.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
								}
								my $intron2_aa;
								if ($mark2==2) {
									for (my $i=0;$i<$intron2_length;$i+=3){
										my $codon=substr ($string_intron2_rev_com,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$intron2_aa.=$hash_codon{$codon};
										}else{
											$intron2_aa.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
								}

								#if ((defined $end2_new) and (defined $start3_new) and (defined $end3_new) and (defined $start4_new)) {
								if (abs($start2-$end4) < 3000) {
									if (($start3 < $end2) and ($start4 < $end3)) {
										if ((($mark1==0) or ((defined $intron1_aa) and ($intron1_aa=~ /\*/))) and (($mark2==0) or ((defined $intron2_aa) and ($intron2_aa=~ /\*/))) and (($intron1_length > 30) and ($intron2_length > 30))) {# (two introns) spliced exon,not 3x or have stop codon
											print $out_annotation "     "."gene"."            "."complement(".$end1."..".$start1.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."complement(join(".$end1."..".$start4_new.",".$end3_new."..".$start3_new.",".$end2_new."..".$start1."))"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
											$gene_number_seq{$name}++;
										}elsif ((($mark1==1) or ((defined $intron1_aa) and ($intron1_aa!~ /\*/))) and (($mark2==1) or ((defined $intron2_aa) and ($intron2_aa!~ /\*/))) and (($intron1_length < 30) and ($intron2_length < 30))) {# (no intron) joined exon,no intron or no stop codon
											print $out_annotation "     "."gene"."            "."complement(".$end1."..".$start1.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."complement(".$end1."..".$start1.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (negative two-intron PCG) lost two introns!\n";
										}elsif ((($mark1==0) or ((defined $intron1_aa) and ($intron1_aa=~ /\*/))) and (($mark2==1) or ((defined $intron2_aa) and ($intron2_aa!~ /\*/))) and (($intron1_length > 30) and ($intron2_length < 30))) {# (loss intron2)
											print $out_annotation "     "."gene"."            "."complement(".$end1."..".$start1.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."complement(join(".$end1."..".$start3_new.",".$end2_new."..".$start1."))"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (negative two-intron PCG) lost intron 2!\n";
										}elsif ((($mark1==1) or ((defined $intron1_aa) and ($intron1_aa!~ /\*/))) and (($mark2==0) or ((defined $intron2_aa) and ($intron2_aa=~ /\*/))) and (($intron1_length < 30) and ($intron2_length > 30))) {# (loss intron1)
											print $out_annotation "     "."gene"."            "."complement(".$end1."..".$start1.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."complement(join(".$end1."..".$start4_new.",".$end3_new."..".$start1."))"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (negative two-intron PCG) lost intron 1!\n";
										}
									}elsif (($start3 > $end2) or ($start4 > $end3)) {
										print $logfile "Warning: $name (negative two-intron PCG) has not been annotated!\n";
									}
								#}elsif ((!defined $end2_new) or (!defined $start3_new) or (!defined $end3_new) or (!defined $start4_new)) {
								}elsif (abs($start2-$end4) > 3000) {
									print $out_annotation "     "."gene"."            "."complement(".$end3_new."..".$start3_new.")\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."CDS"."             "."complement(".$end3_new."..".$start3_new.")"."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/codon_start=1"."\n";
									print $out_annotation "                     "."/transl_table=11"."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
									$gene_number_seq{$name}++;
									print $logfile "Warning: $name (negative two-intron PCG) must be checked!\n";
								}
							}elsif((grep {$_=~ /\*/} @aa_exon) or (($reverse_start_codon ne "ATG") and ($reverse_start_codon ne "GTG")) or (($reverse_stop_codon ne "TAA") and ($reverse_stop_codon ne "TAG") and ($reverse_stop_codon ne "TGA"))){# non-standard start or stop codon for _gene
								my $start_left=$start1-60;
								my $start_right;
								if (($length_cp-$start1)>=9000) {
									$start_right=$start1+9000;
								}elsif (($length_cp-$start1)<9000){
									$start_right=$start1+int(($length_cp-$start1)/3)*3;
								}
								my $length_exon1=($start1-$end2+1);
								my $length_exon2=($start3-$end3+1);
								my $length_exon12=$length_exon1+$length_exon2;
								my $end_left;
								if ($start4>=9002) {
									if ($length_exon12 % 3==0){
										$end_left=$start4-9000;
									}
									if ($length_exon12 % 3==1){
										$end_left=$start4-9002;
									}
									if ($length_exon12 % 3==2){
										$end_left=$start4-9001;
									}
								}elsif ($start4<9002){
									if ($length_exon12 % 3==0){
										$end_left=$start4-(int($start4/3)*3);
									}
									if ($length_exon12 % 3==1){
										$end_left=$start4-(int($start4/3)*3+2);
									}
									if ($length_exon12 % 3==2){
										$end_left=$start4-(int($start4/3)*3+1);
									}
								}
								my $str1=substr($sequence,($start_left-1),($start_right-$start_left+1));
								my $str2=substr($sequence,($end_left-1),($start4-$end_left+1));
								my $length1=length $str1;
								my $length2=length $str2;
								my $rev_com1=reverse $str1;
								$rev_com1=~ tr/ACGTacgt/TGCAtgca/;
								my $rev_com2=reverse $str2;
								$rev_com2=~ tr/ACGTacgt/TGCAtgca/;
								my $aa1;#left range for start codon
								for (my $i=0;$i<($length1-3);$i+=3){# delete stop codon
									my $codon=substr ($rev_com1,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa1.=$hash_codon{$codon};
									}else{
										$aa1.="X";
										my $m=$i+1;
										#print "Bad codon $codon in position $m of gene $name in species $contig!\n";
									}
								}
								my @aa1=split //,$aa1;
								my (@star1,@start1);
								foreach (0..$#aa1){
									if (($aa1[$_]=~ /\*/) and ($_ < 3000)){
										push @star1,$_;
									}
									if (($aa1[$_]=~ /M/) and ($_ <= 3060)){
										push @start1,$_;
									}
								}
								my $right_star_position=pop @star1;
								until ($right_star_position<(3000+$length_exon1/3)){
									$right_star_position=pop @star1;
								}
								my @start2;
								foreach (@start1){
									push @start2,$_ if ((defined $right_star_position) and ($_ > $right_star_position))
								}
								#my $right_start_position=$start2[0];
								my $right_start_position=pop @start2;
								my $aa2;#right range for stop codon
								my $j;
								if ($length_exon1 % 3==0){
									$j=0;
								}
								if ($length_exon1 % 3==1){
									$j=2;
								}
								if ($length_exon1 % 3==2){
									$j=1;
								}
								for (my $k=$j;$k<($length2-3);$k+=3){# delete stop codon
									my $codon=substr ($rev_com2,$k,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa2.=$hash_codon{$codon};
									}else{
										$aa2.="X";
										my $n=$k+1;
										#print "Bad codon $codon in position $n of gene $name in species $contig!\n";
									}
								}
								my @aa2=split //,$aa2;
								my @star2;
								foreach (0..$#aa2){
									if ($aa2[$_]=~ /\*/){
										push @star2,$_;
									}
								}
								my $left_star_position=shift @star2;

								my $seq_string1=substr($sequence,($end2-61),120);
								my $seq_string2=substr($sequence,($start3-60),120);
								my $seq_string3=substr($sequence,($end3-61),120);
								my $seq_string4=substr($sequence,($start4-60),120);
								my $rev_seq_string1=reverse $seq_string1;
								$rev_seq_string1=~ tr/ACGTacgt/TGCAtgca/;
								my $rev_seq_string2=reverse $seq_string2;
								$rev_seq_string2=~ tr/ACGTacgt/TGCAtgca/;
								my $rev_seq_string3=reverse $seq_string3;
								$rev_seq_string3=~ tr/ACGTacgt/TGCAtgca/;
								my $rev_seq_string4=reverse $seq_string4;
								$rev_seq_string4=~ tr/ACGTacgt/TGCAtgca/;
								my $aa_seq_string1;
								my $aa_seq_string2;
								my $aa_seq_string3;
								my $aa_seq_string4;
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($rev_seq_string1,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string1.=$hash_codon{$codon};
									}else{
										$aa_seq_string1.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($rev_seq_string2,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string2.=$hash_codon{$codon};
									}else{
										$aa_seq_string2.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($rev_seq_string3,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string3.=$hash_codon{$codon};
									}else{
										$aa_seq_string3.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($rev_seq_string4,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string4.=$hash_codon{$codon};
									}else{
										$aa_seq_string4.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}

								if ((defined $right_start_position) and (defined $left_star_position)){
									my $start=$start_right-$right_start_position * 3;
									my $end;
									if ($length_exon12 % 3==0){
										$end=($start4+1)-($left_star_position * 3+3);
									}
									if ($length_exon12 % 3==1){
										$end=($start4-1)-($left_star_position * 3+3);
									}
									if ($length_exon12 % 3==2){
										$end=($start4)-($left_star_position * 3+3);
									}

									my ($str1,$str2,$str3,$str,$rev_com,$length,$end2_new,$start3_new,$end3_new,$start4_new);
									my ($string_intron1,$string_intron2,$intron1_length,$intron2_length,$string_intron1_rev_com,$string_intron2_rev_com);
									my $mark1=0;
									my $mark2=0;
									if ($length_ref_exon1 % 3==0) {
										for (my $i=$match;$i<$a;$i+=1) {
											my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string1,$j,$match);
												if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
													$end2_new=$end2+60-($j+$match)*3-($i-$match)*3;
													$marks++;
												}
											}
											if (defined $end2_new){
												last;
											}
										}
										if (!defined $end2_new) {
											$end2_new=$end2;
										}
										for (my $i=0;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string2,-$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$start3_new=$start3-60+$j*3+$i*3;
													$marks++;
												}
											}
											if (defined $start3_new){
												last;
											}
										}
										if (!defined $start3_new) {
											$start3_new=$start3;
										}
										$intron1_length=$end2_new-$start3_new-1;
										if ($intron1_length==0) {
											$mark1=1;
										}
										if (($intron1_length!=0) and ($intron1_length % 3==0)) {
											$string_intron1=substr ($sequence,$start3_new,$intron1_length);
											$intron1_length=length $string_intron1;
											$string_intron1_rev_com=reverse $string_intron1;
											$string_intron1_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark1=2;
										}

										if ($length_ref_exon2 % 3==0) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3+60-($j+$match)*3-($i-$match)*3;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4-60+$j*3+$i*3;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4;
											}
											$intron2_length=$end3_new-$start4_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,$start4_new,$intron2_length);
												$intron2_length=length $string_intron2;
												$string_intron2_rev_com=reverse $string_intron2;
												$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark2=2;
											}
										}elsif ($length_ref_exon2 % 3==1) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3+60-($j+$match)*3-($i-$match)*3-1;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3-1;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4-60+$j*3+$i*3+2;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4+2;
											}
											$intron2_length=$end3_new-$start4_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,($start4_new+1),$intron2_length);
												$intron2_length=length $string_intron2;
												$string_intron2_rev_com=reverse $string_intron2;
												$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark2=2;
											}
										}elsif ($length_ref_exon2 % 3==2) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3+60-($j+$match)*3-($i-$match)*3-2;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3-2;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4-60+$j*3+$i*3+1;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4+1;
											}
											$intron2_length=$end3_new-$start4_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,($start4_new-1),$intron2_length);
												$intron2_length=length $string_intron2;
												$string_intron2_rev_com=reverse $string_intron2;
												$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark2=2;
											}
										}
									}elsif ($length_ref_exon1 % 3==1) {
										for (my $i=$match;$i<$a;$i+=1) {
											my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string1,$j,$match);
												if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
													$end2_new=$end2+60-($j+$match)*3-($i-$match)*3-1;
													$marks++;
												}
											}
											if (defined $end2_new){
												last;
											}
										}
										if (!defined $end2_new) {
											$end2_new=$end2-1;
										}
										for (my $i=0;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string2,-$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$start3_new=$start3-60+$j*3+$i*3+2;
													$marks++;
												}
											}
											if (defined $start3_new){
												last;
											}
										}
										if (!defined $start3_new) {
											$start3_new=$start3+2;
										}
										$intron1_length=$end2_new-$start3_new-1;
										if ($intron1_length==0) {
											$mark1=1;
										}
										if (($intron1_length!=0) and ($intron1_length % 3==0)) {
											$string_intron1=substr ($sequence,($start3_new+1),$intron1_length);
											$intron1_length=length $string_intron1;
											$string_intron1_rev_com=reverse $string_intron1;
											$string_intron1_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark1=2;
										}

										if ($length_ref_exon2 % 3==0) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3+60-($j+$match)*3-($i-$match)*3-1;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3-1;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4-60+$j*3+$i*3+2;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4+2;
											}
											$intron2_length=$end3_new-$start4_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,($start4_new+1),$intron2_length);
												$intron2_length=length $string_intron2;
												$string_intron2_rev_com=reverse $string_intron2;
												$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark2=2;
											}
										}elsif ($length_ref_exon2 % 3==1) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3+60-($j+$match)*3-($i-$match)*3-2;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3-2;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4-60+$j*3+$i*3+1;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4+1;
											}
											$intron2_length=$end3_new-$start4_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,($start4_new-1),$intron2_length);
												$intron2_length=length $string_intron2;
												$string_intron2_rev_com=reverse $string_intron2;
												$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark2=2;
											}
										}elsif ($length_ref_exon2 % 3==2) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3+60-($j+$match)*3-($i-$match)*3;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4-60+$j*3+$i*3;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4;
											}
											$intron2_length=$end3_new-$start4_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,$start4_new,$intron2_length);
												$intron2_length=length $string_intron2;
												$string_intron2_rev_com=reverse $string_intron2;
												$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark2=2;
											}
										}
									}elsif ($length_ref_exon1 % 3==2) {
										for (my $i=$match;$i<$a;$i+=1) {
											my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string1,$j,$match);
												if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
													$end2_new=$end2+60-($j+$match)*3-($i-$match)*3-2;
													$marks++;
												}
											}
											if (defined $end2_new){
												last;
											}
										}
										if (!defined $end2_new) {
											$end2_new=$end2-2;
										}
										for (my $i=0;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string2,-$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$start3_new=$start3-60+$j*3+$i*3+1;
													$marks++;
												}
											}
											if (defined $start3_new){
												last;
											}
										}
										if (!defined $start3_new) {
											$start3_new=$start3+1;
										}
										$intron1_length=$end2_new-$start3_new-1;
										if ($intron1_length==0) {
											$mark1=1;
										}elsif (($intron1_length!=0) and ($intron1_length % 3==0)) {
											$string_intron1=substr ($sequence,($start3_new-1),$intron1_length);
											$intron1_length=length $string_intron1;
											$string_intron1_rev_com=reverse $string_intron1;
											$string_intron1_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark1=2;
										}

										if ($length_ref_exon2 % 3==0) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3+60-($j+$match)*3-($i-$match)*3-2;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3-2;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4-60+$j*3+$i*3+1;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4+1;
											}
											$intron2_length=$end3_new-$start4_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,($start4_new-1),$intron2_length);
												$intron2_length=length $string_intron2;
												$string_intron2_rev_com=reverse $string_intron2;
												$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark2=2;
											}
										}elsif ($length_ref_exon2 % 3==1) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3+60-($j+$match)*3-($i-$match)*3;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4-60+$j*3+$i*3;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4;
											}
											$intron2_length=$end3_new-$start4_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,$start4_new,$intron2_length);
												$intron2_length=length $string_intron2;
												$string_intron2_rev_com=reverse $string_intron2;
												$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark2=2;
											}
										}elsif ($length_ref_exon2 % 3==2) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3+60-($j+$match)*3-($i-$match)*3-1;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3-1;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4-60+$j*3+$i*3+2;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4+2;
											}
											$intron2_length=$end3_new-$start4_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,($start4_new+1),$intron2_length);
												$intron2_length=length $string_intron2;
												$string_intron2_rev_com=reverse $string_intron2;
												$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark2=2;
											}
										}
									}

									$str1=substr($sequence,($end2_new-1),($start-($end2_new-1)));
									$str2=substr($sequence,($end3_new-1),($start3_new-($end3_new-1)));
									$str3=substr($sequence,($end-1),($start4_new-($end-1)));
									$str=$str1.$str2.$str3;
									$length=length $str;
									$rev_com=reverse $str;
									$rev_com=~ tr/ACGTacgt/TGCAtgca/;

									my $aa;
									for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
										my $codon=substr ($rev_com,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa.=$hash_codon{$codon};
										}else{
											$aa.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}

									my $intron1_aa;
									if ($mark1==2) {
										for (my $i=0;$i<$intron1_length;$i+=3){
											my $codon=substr ($string_intron1_rev_com,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$intron1_aa.=$hash_codon{$codon};
											}else{
												$intron1_aa.="X";
												my $j=$i+1;
												#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}
									}
									my $intron2_aa;
									if ($mark2==2) {
										for (my $i=0;$i<$intron2_length;$i+=3){
											my $codon=substr ($string_intron2_rev_com,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$intron2_aa.=$hash_codon{$codon};
											}else{
												$intron2_aa.="X";
												my $j=$i+1;
												#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}
									}

									#if ((defined $end2_new) and (defined $start3_new) and (defined $end3_new) and (defined $start4_new)) {
									if (abs($start2-$end4) < 3000) {
										if (($start3 < $end2) and ($start4 < $end3)) {
											if ((($mark1==0) or ((defined $intron1_aa) and ($intron1_aa=~ /\*/))) and (($mark2==0) or ((defined $intron2_aa) and ($intron2_aa=~ /\*/))) and (($intron1_length > 30) and ($intron2_length > 30))) {# (two introns) spliced exon,not 3x or have stop codon
												print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."complement(join(".$end."..".$start4_new.",".$end3_new."..".$start3_new.",".$end2_new."..".$start."))"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
												$gene_number_seq{$name}++;
											}elsif ((($mark1==1) or ((defined $intron1_aa) and ($intron1_aa!~ /\*/))) and (($mark2==1) or ((defined $intron2_aa) and ($intron2_aa!~ /\*/))) and (($intron1_length < 30) and ($intron2_length < 30))) {# (no intron) joined exon,no intron or no stop codon
												print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."complement(".$end."..".$start.")\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
												$gene_number_seq{$name}++;
												print $logfile "Warning: $name (negative two-intron PCG) lost two introns!\n";
											}elsif ((($mark1==0) or ((defined $intron1_aa) and ($intron1_aa=~ /\*/))) and (($mark2==1) or ((defined $intron2_aa) and ($intron2_aa!~ /\*/))) and (($intron1_length > 30) and ($intron2_length < 30))) {# (loss intron2)
												print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."complement(join(".$end."..".$start3_new.",".$end2_new."..".$start."))"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
												$gene_number_seq{$name}++;
												print $logfile "Warning: $name (negative two-intron PCG) lost intron 2!\n";
											}elsif ((($mark1==1) or ((defined $intron1_aa) and ($intron1_aa!~ /\*/))) and (($mark2==0) or ((defined $intron2_aa) and ($intron2_aa=~ /\*/))) and (($intron1_length < 30) and ($intron2_length > 30))) {# (loss intron1)
												print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."complement(join(".$end."..".$start4_new.",".$end3_new."..".$start."))"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
												$gene_number_seq{$name}++;
												print $logfile "Warning: $name (negative two-intron PCG) lost intron 1!\n";
											}
										}elsif (($start3 > $end2) or ($start4 > $end3)) {
											print $logfile "Warning: $name (negative two-intron PCG) has not been annotated!\n";
										}
									#}elsif ((!defined $end2_new) or (!defined $start3_new) or (!defined $end3_new) or (!defined $start4_new)) {
									}elsif (abs($start2-$end4) > 3000) {
										print $out_annotation "     "."gene"."            "."complement(".$end3_new."..".$start3_new.")\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."complement(".$end3_new."..".$start3_new.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
										$gene_number_seq{$name}++;
										print $logfile "Warning: $name (negative two-intron PCG) must be checked!\n";
									}
								}elsif ((!defined $right_start_position) and (defined $left_star_position)){
									my $end;
									if ($length_exon12 % 3==0){
										$end=($start4+1)-($left_star_position * 3+3);
									}
									if ($length_exon12 % 3==1){
										$end=($start4-1)-($left_star_position * 3+3);
									}
									if ($length_exon12 % 3==2){
										$end=($start4)-($left_star_position * 3+3);
									}
									if (abs($start2-$end4) < 3000) {
										if (($start3 < $end2) and ($start4 < $end3)) {
											print $out_annotation "     "."gene"."            "."complement(".$end."..".$start1.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."complement(join(".$end."..".$start4.",".$end3."..".$start3.",".$end2."..".$start1."))"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (negative two-intron PCG) has non-canonical start codon!\n";
										}elsif (($start3 > $end2) or ($start4 > $end3)) {
											print $logfile "Warning: $name (negative two-intron PCG) has not been annotated!\n";
										}
									}elsif (abs($start2-$end4) > 3000) {
										print $out_annotation "     "."gene"."            "."complement(".$end3."..".$start3.")\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."complement(".$end3."..".$start3.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
										$gene_number_seq{$name}++;
										print $logfile "Warning: $name (negative two-intron PCG) must be checked!\n";
									}
								}
							}else{
								print "The exon length of $name (negative two-intron PCG) is not an intergral multiple of 3!\n";
							}
						}
					}elsif(((($start1 != $start2) and (defined $end4)) or ((defined $end4) and ($end1 != $end4))) and ((defined $length_ref_exon1) and (defined $length_ref_exon2) and (defined $length_ref_exon3) and (defined $sequence_ref_exon1) and (defined $sequence_ref_exon2) and (defined $sequence_ref_exon3))){# non-identical PCG boundary for _gene and -1_coding, -3_coding
						if ($start3 < $end3){# positive
							my $str1=substr($sequence,($start2-1),($end2-$start2+1));
							my $str2=substr($sequence,($start3-1),($end3-$start3+1));
							my $str3=substr($sequence,($start4-1),($end4-$start4+1));
							my $str4=$str1.$str2.$str3;
							my $length_exon=length $str4;
							my $forward_start_codon=substr($str4,0,3);
							my $forward_stop_codon=substr($str4,($length_exon-3),3);
							my $aa_exon;
							for (my $i=0;$i<($length_exon-3);$i+=3){# delete stop codon
								my $codon=substr ($str4,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa_exon.=$hash_codon{$codon};
								}else{
									$aa_exon.="X";
									my $j=$i+1;
									#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa_exon=split //,$aa_exon;
							if (($length_exon % 3==0) and (!grep {$_=~ /\*/} @aa_exon) and (($forward_start_codon eq "ATG") or ($forward_start_codon eq "GTG")) and (($forward_stop_codon eq "TAA") or ($forward_stop_codon eq "TAG") or ($forward_stop_codon eq "TGA"))){# standard start and stop codon for _gene
								my $seq_string1=substr($sequence,($end2-60),120);
								my $seq_string2=substr($sequence,($start3-61),120);
								my $seq_string3=substr($sequence,($end3-60),120);
								my $seq_string4=substr($sequence,($start4-61),120);
								my $aa_seq_string1;
								my $aa_seq_string2;
								my $aa_seq_string3;
								my $aa_seq_string4;
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($seq_string1,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string1.=$hash_codon{$codon};
									}else{
										$aa_seq_string1.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($seq_string2,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string2.=$hash_codon{$codon};
									}else{
										$aa_seq_string2.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($seq_string3,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string3.=$hash_codon{$codon};
									}else{
										$aa_seq_string3.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($seq_string4,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string4.=$hash_codon{$codon};
									}else{
										$aa_seq_string4.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}

								my ($str1,$str2,$str3,$str,$length,$end2_new,$start3_new,$end3_new,$start4_new);
								my ($string_intron1,$string_intron2,$intron1_length,$intron2_length);
								my $mark1=0;
								my $mark2=0;
								if ($length_ref_exon1 % 3==0) {
									for (my $i=$match;$i<$a;$i+=1) {
										my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
										my $marks=0;
										for (my $j=0;$j<40;$j+=1) {
											my $repeat=substr ($aa_seq_string1,$j,$match);
											if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
												$end2_new=$end2-60+($j+$match)*3+($i-$match)*3;
												$marks++;
											}
										}
										if (defined $end2_new){
											last;
										}
									}
									if (!defined $end2_new) {
										$end2_new=$end2;
									}
									for (my $i=0;$i<$b;$i+=1) {
										my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
										my $marks=0;
										for (my $j=$match;$j<40;$j+=1) {
											my $repeat=substr ($aa_seq_string2,-$j,$match);
											if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
												$start3_new=$start3+60-$j*3-$i*3;
												$marks++;
											}
										}
										if (defined $start3_new){
											last;
										}
									}
									if (!defined $start3_new) {
										$start3_new=$start3;
									}
									$intron1_length=$start3_new-$end2_new-1;
									if ($intron1_length==0) {
										$mark1=1;
									}
									if (($intron1_length!=0) and ($intron1_length % 3==0)) {
										$string_intron1=substr ($sequence,$end2_new,$intron1_length);
										$intron1_length=length $string_intron1;
										$mark1=2;
									}

									if ($length_ref_exon2 % 3==0) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3-60+($j+$match)*3+($i-$match)*3;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4+60-$j*3-$i*3;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4;
										}
										$intron2_length=$start4_new-$end3_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,$end3_new,$intron2_length);
											$intron2_length=length $string_intron2;
											$mark2=2;
										}
									}elsif ($length_ref_exon2 % 3==1) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3-60+($j+$match)*3+($i-$match)*3+1;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3+1;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4+60-$j*3-$i*3-2;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4-2;
										}
										$intron2_length=$start4_new-$end3_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,($end3_new-1),$intron2_length);
											$intron2_length=length $string_intron2;
											$mark2=2;
										}
									}elsif ($length_ref_exon2 % 3==2) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3-60+($j+$match)*3+($i-$match)*3+2;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3+2;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4+60-$j*3-$i*3-1;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4-1;
										}
										$intron2_length=$start4_new-$end3_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,($end3_new+1),$intron2_length);
											$intron2_length=length $string_intron2;
											$mark2=2;
										}
									}
								}elsif ($length_ref_exon1 % 3==1) {
									for (my $i=$match;$i<$a;$i+=1) {
										my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
										my $marks=0;
										for (my $j=0;$j<40;$j+=1) {
											my $repeat=substr ($aa_seq_string1,$j,$match);
											if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
												$end2_new=$end2-60+($j+$match)*3+($i-$match)*3+1;
												$marks++;
											}
										}
										if (defined $end2_new){
											last;
										}
									}
									if (!defined $end2_new) {
										$end2_new=$end2+1;
									}
									for (my $i=0;$i<$b;$i+=1) {
										my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
										my $marks=0;
										for (my $j=$match;$j<40;$j+=1) {
											my $repeat=substr ($aa_seq_string2,-$j,$match);
											if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
												$start3_new=$start3+60-$j*3-$i*3-2;
												$marks++;
											}
										}
										if (defined $start3_new){
											last;
										}
									}
									if (!defined $start3_new) {
										$start3_new=$start3-2;
									}
									$intron1_length=$start3_new-$end2_new-1;
									if ($intron1_length==0) {
										$mark1=1;
									}
									if (($intron1_length!=0) and ($intron1_length % 3==0)) {
										$string_intron1=substr ($sequence,($end2_new-1),$intron1_length);
										$intron1_length=length $string_intron1;
										$mark1=2;
									}

									if ($length_ref_exon2 % 3==0) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3-60+($j+$match)*3+($i-$match)*3+1;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3+1;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4+60-$j*3-$i*3-2;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4-2;
										}
										$intron2_length=$start4_new-$end3_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,($end3_new-1),$intron2_length);
											$intron2_length=length $string_intron2;
											$mark2=2;
										}
									}elsif ($length_ref_exon2 % 3==1) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3-60+($j+$match)*3+($i-$match)*3+2;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3+2;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4+60-$j*3-$i*3-1;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4-1;
										}
										$intron2_length=$start4_new-$end3_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,($end3_new+1),$intron2_length);
											$intron2_length=length $string_intron2;
											$mark2=2;
										}
									}elsif ($length_ref_exon2 % 3==2) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3-60+($j+$match)*3+($i-$match)*3;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4+60-$j*3-$i*3;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4;
										}
										$intron2_length=$start4_new-$end3_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,$end3_new,$intron2_length);
											$intron2_length=length $string_intron2;
											$mark2=2;
										}
									}
								}elsif ($length_ref_exon1 % 3==2) {
									for (my $i=$match;$i<$a;$i+=1) {
										my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
										my $marks=0;
										for (my $j=0;$j<40;$j+=1) {
											my $repeat=substr ($aa_seq_string1,$j,$match);
											if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
												$end2_new=$end2-60+($j+$match)*3+($i-$match)*3+2;
												$marks++;
											}
										}
										if (defined $end2_new){
											last;
										}
									}
									if (!defined $end2_new) {
										$end2_new=$end2+2;
									}
									for (my $i=0;$i<$b;$i+=1) {
										my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
										my $marks=0;
										for (my $j=$match;$j<40;$j+=1) {
											my $repeat=substr ($aa_seq_string2,-$j,$match);
											if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
												$start3_new=$start3+60-$j*3-$i*3-1;
												$marks++;
											}
										}
										if (defined $start3_new){
											last;
										}
									}
									if (!defined $start3_new) {
										$start3_new=$start3-1;
									}
									$intron1_length=$start3_new-$end2_new-1;
									if ($intron1_length==0) {
										$mark1=1;
									}elsif (($intron1_length!=0) and ($intron1_length % 3==0)) {
										$string_intron1=substr ($sequence,($end2_new+1),$intron1_length);
										$intron1_length=length $string_intron1;
										$mark1=2;
									}

									if ($length_ref_exon2 % 3==0) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3-60+($j+$match)*3+($i-$match)*3+2;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3+2;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4+60-$j*3-$i*3-1;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4-1;
										}
										$intron2_length=$start4_new-$end3_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,($end3_new+1),$intron2_length);
											$intron2_length=length $string_intron2;
											$mark2=2;
										}
									}elsif ($length_ref_exon2 % 3==1) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3-60+($j+$match)*3+($i-$match)*3;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4+60-$j*3-$i*3;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4;
										}
										$intron2_length=$start4_new-$end3_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,$end3_new,$intron2_length);
											$intron2_length=length $string_intron2;
											$mark2=2;
										}
									}elsif ($length_ref_exon2 % 3==2) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3-60+($j+$match)*3+($i-$match)*3+1;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3+1;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4+60-$j*3-$i*3-2;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4-2;
										}
										$intron2_length=$start4_new-$end3_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,($end3_new-1),$intron2_length);
											$intron2_length=length $string_intron2;
											$mark2=2;
										}
									}
								}

								$str1=substr($sequence,($start2-1),($end2_new-($start2-1)));
								$str2=substr($sequence,($start3_new-1),($end3_new-($start3_new-1)));
								$str3=substr($sequence,($start4_new-1),($end4-($start4_new-1)));
								$str=$str1.$str2.$str3;
								$length=length $str;

								my $aa;
								for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
									my $codon=substr ($str,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa.=$hash_codon{$codon};
									}else{
										$aa.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}

								my $intron1_aa;
								if ($mark1==2) {
									for (my $i=0;$i<$intron1_length;$i+=3){
										my $codon=substr ($string_intron1,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$intron1_aa.=$hash_codon{$codon};
										}else{
											$intron1_aa.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
								}
								my $intron2_aa;
								if ($mark2==2) {
									for (my $i=0;$i<$intron2_length;$i+=3){
										my $codon=substr ($string_intron2,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$intron2_aa.=$hash_codon{$codon};
										}else{
											$intron2_aa.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
								}

								#if ((defined $end2_new) and (defined $start3_new) and (defined $end3_new) and (defined $start4_new)) {
								if (abs($end4-$start2) < 3000) {
									if (($start3 > $end2) and ($start4 > $end3)) {
										if ((($mark1==0) or ((defined $intron1_aa) and ($intron1_aa=~ /\*/))) and (($mark2==0) or ((defined $intron2_aa) and ($intron2_aa=~ /\*/))) and (($intron1_length > 30) and ($intron2_length > 30))) {# (two introns) spliced exon,not 3x or have stop codon
											print $out_annotation "     "."gene"."            ".$start2."..".$end4."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."join(".$start2."..".$end2_new.",".$start3_new."..".$end3_new.",".$start4_new."..".$end4.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
											$gene_number_seq{$name}++;
										}elsif ((($mark1==1) or ((defined $intron1_aa) and ($intron1_aa!~ /\*/))) and (($mark2==1) or ((defined $intron2_aa) and ($intron2_aa!~ /\*/))) and (($intron1_length < 30) and ($intron2_length < 30))) {# (no intron) joined exon,no intron or no stop codon
											print $out_annotation "     "."gene"."            ".$start2."..".$end4."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             ".$start2."..".$end4."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (positive two-intron PCG) lost two introns!\n";
										}elsif ((($mark1==0) or ((defined $intron1_aa) and ($intron1_aa=~ /\*/))) and (($mark2==1) or ((defined $intron2_aa) and ($intron2_aa!~ /\*/))) and (($intron1_length > 30) and ($intron2_length < 30))) {# (loss intron2)
											print $out_annotation "     "."gene"."            ".$start2."..".$end4."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."join(".$start2."..".$end2_new.",".$start3_new."..".$end4.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (positive two-intron PCG) lost intron 2!\n";
										}elsif ((($mark1==1) or ((defined $intron1_aa) and ($intron1_aa!~ /\*/))) and (($mark2==0) or ((defined $intron2_aa) and ($intron2_aa=~ /\*/))) and (($intron1_length < 30) and ($intron2_length > 30))) {# (loss intron1)
											print $out_annotation "     "."gene"."            ".$start2."..".$end4."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."join(".$start2."..".$end3_new.",".$start4_new."..".$end4.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (positive two-intron PCG) lost intron 1!\n";
										}
									}elsif (($start3 < $end2) or ($start4 < $end3)) {
										print $logfile "Warning: $name (positive two-intron PCG) has not been annotated!\n";
									}
								#}elsif ((!defined $end2_new) or (!defined $start3_new) or (!defined $end3_new) or (!defined $start4_new)) {
								}elsif (abs($end4-$start2) > 3000) {
									print $out_annotation "     "."gene"."            ".$start3_new."..".$end3_new."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."CDS"."             ".$start3_new."..".$end3_new."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/codon_start=1"."\n";
									print $out_annotation "                     "."/transl_table=11"."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
									$gene_number_seq{$name}++;
									print $logfile "Warning: $name (positive two-intron PCG) must be checked!\n";
								}
							}elsif((grep {$_=~ /\*/} @aa_exon) or (($forward_start_codon ne "ATG") and ($forward_start_codon ne "GTG")) or (($forward_stop_codon ne "TAA") and ($forward_stop_codon ne "TAG") and ($forward_stop_codon ne "TGA"))){# non-standard start or stop codon for _gene
								my $start_left;
								if ($start2>=9000) {
									$start_left=$start2-9000;
								}elsif ($start2<9000){
									$start_left=$start2-int($start2/3)*3;
								}
								my $start_right=$start2+60;
								my $length_exon1=($end2-$start2+1);
								my $length_exon2=($end3-$start3+1);
								my $length_exon12=$length_exon1+$length_exon2;
								my $end_right;
								if (($length_cp-$start4)>=9002) {
									if ($length_exon12 % 3==0){
										$end_right=$start4+9000;
									}
									if ($length_exon12 % 3==1){
										$end_right=$start4+9001;
									}
									if ($length_exon12 % 3==2){
										$end_right=$start4+9002;
									}
								}elsif (($length_cp-$start4)<9002){
									if ($length_exon12 % 3==0){
										$end_right=$start4+(int(($length_cp-$start4)/3)*3);
									}
									if ($length_exon12 % 3==1){
										$end_right=$start4+(int(($length_cp-$start4)/3)*3+1);
									}
									if ($length_exon12 % 3==2){
										$end_right=$start4+(int(($length_cp-$start4)/3)*3+2);
									}
								}
								my $str1=substr($sequence,($start_left-1),($start_right-$start_left+1));
								my $str2=substr($sequence,($start4-1),($end_right-$start4+1));
								my $length1=length $str1;
								my $length2=length $str2;
								my $aa1;#left range for start codon
								for (my $i=0;$i<($length1-3);$i+=3){# delete stop codon
									my $codon=substr ($str1,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa1.=$hash_codon{$codon};
									}else{
										$aa1.="X";
										my $m=$i+1;
										#print "Bad codon $codon in position $m of gene $name in species $contig!\n";
									}
								}
								my @aa1=split //,$aa1;
								my (@star1,@start1);
								foreach (0..$#aa1){
									if (($aa1[$_]=~ /\*/) and ($_ < 3000)){
										push @star1,$_;
									}
									if (($aa1[$_]=~ /M/) and ($_ <= 3060)){
										push @start1,$_;
										}
								}
								my $left_star_position=pop @star1;
								until ($left_star_position<(3000+$length_exon1/3)){
									$left_star_position=pop @star1;
								}
								my @start2;
								foreach (@start1){
									push @start2,$_ if ((defined $left_star_position) and ($_ > $left_star_position))
								}
								#my $left_start_position=$start2[0];
								my $left_start_position=pop @start2;
								my $aa2;#right range for stop codon
								my $j;
								if ($length_exon1 % 3==0){
									$j=0;
								}
								if ($length_exon1 % 3==1){
									$j=2;
								}
								if ($length_exon1 % 3==2){
									$j=1;
								}
								for (my $k=$j;$k<($length2-3);$k+=3){# delete stop codon
									my $codon=substr ($str2,$k,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa2.=$hash_codon{$codon};
									}else{
										$aa2.="X";
										my $n=$k+1;
										#print "Bad codon $codon in position $n of gene $name in species $contig!\n";
									}
								}
								my @aa2=split //,$aa2;
								my @star2;
								foreach (0..$#aa2){
									if ($aa2[$_]=~ /\*/){
										push @star2,$_;
									}
								}
								my $right_star_position=shift @star2;

								my $seq_string1=substr($sequence,($end2-60),120);
								my $seq_string2=substr($sequence,($start3-61),120);
								my $seq_string3=substr($sequence,($end3-60),120);
								my $seq_string4=substr($sequence,($start4-61),120);
								my $aa_seq_string1;
								my $aa_seq_string2;
								my $aa_seq_string3;
								my $aa_seq_string4;
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($seq_string1,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string1.=$hash_codon{$codon};
									}else{
										$aa_seq_string1.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($seq_string2,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string2.=$hash_codon{$codon};
									}else{
										$aa_seq_string2.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($seq_string3,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string3.=$hash_codon{$codon};
									}else{
										$aa_seq_string3.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($seq_string4,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string4.=$hash_codon{$codon};
									}else{
										$aa_seq_string4.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}

								if ((defined $left_start_position) and (defined $right_star_position)){
									my $start=$start_left+$left_start_position * 3;
									my $end;
									if ($length_exon12 % 3==0){
										$end=($start4-1)+($right_star_position * 3+3);
									}
									if ($length_exon12 % 3==1){
										$end=($start4+1)+($right_star_position * 3+3);
									}
									if ($length_exon12 % 3==2){
										$end=($start4)+($right_star_position * 3+3);
									}

									my ($str1,$str2,$str3,$str,$length,$end2_new,$start3_new,$end3_new,$start4_new);
									my ($string_intron1,$string_intron2,$intron1_length,$intron2_length);
									my $mark1=0;
									my $mark2=0;
									if ($length_ref_exon1 % 3==0) {
										for (my $i=$match;$i<$a;$i+=1) {
											my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string1,$j,$match);
												if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
													$end2_new=$end2-60+($j+$match)*3+($i-$match)*3;
													$marks++;
												}
											}
											if (defined $end2_new){
												last;
											}
										}
										if (!defined $end2_new) {
											$end2_new=$end2;
										}
										for (my $i=0;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string2,-$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$start3_new=$start3+60-$j*3-$i*3;
													$marks++;
												}
											}
											if (defined $start3_new){
												last;
											}
										}
										if (!defined $start3_new) {
											$start3_new=$start3;
										}
										$intron1_length=$start3_new-$end2_new-1;
										if ($intron1_length==0) {
											$mark1=1;
										}
										if (($intron1_length!=0) and ($intron1_length % 3==0)) {
											$string_intron1=substr ($sequence,$end2_new,$intron1_length);
											$intron1_length=length $string_intron1;
											$mark1=2;
										}

										if ($length_ref_exon2 % 3==0) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3-60+($j+$match)*3+($i-$match)*3;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4+60-$j*3-$i*3;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4;
											}
											$intron2_length=$start4_new-$end3_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,$end3_new,$intron2_length);
												$intron2_length=length $string_intron2;
												$mark2=2;
											}
										}elsif ($length_ref_exon2 % 3==1) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3-60+($j+$match)*3+($i-$match)*3+1;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3+1;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4+60-$j*3-$i*3-2;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4-2;
											}
											$intron2_length=$start4_new-$end3_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,($end3_new-1),$intron2_length);
												$intron2_length=length $string_intron2;
												$mark2=2;
											}
										}elsif ($length_ref_exon2 % 3==2) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3-60+($j+$match)*3+($i-$match)*3+2;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3+2;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4+60-$j*3-$i*3-1;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4-1;
											}
											$intron2_length=$start4_new-$end3_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,($end3_new+1),$intron2_length);
												$intron2_length=length $string_intron2;
												$mark2=2;
											}
										}
									}elsif ($length_ref_exon1 % 3==1) {
										for (my $i=$match;$i<$a;$i+=1) {
											my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string1,$j,$match);
												if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
													$end2_new=$end2-60+($j+$match)*3+($i-$match)*3+1;
													$marks++;
												}
											}
											if (defined $end2_new){
												last;
											}
										}
										if (!defined $end2_new) {
											$end2_new=$end2+1;
										}
										for (my $i=0;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string2,-$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$start3_new=$start3+60-$j*3-$i*3-2;
													$marks++;
												}
											}
											if (defined $start3_new){
												last;
											}
										}
										if (!defined $start3_new) {
											$start3_new=$start3-2;
										}
										$intron1_length=$start3_new-$end2_new-1;
										if ($intron1_length==0) {
											$mark1=1;
										}
										if (($intron1_length!=0) and ($intron1_length % 3==0)) {
											$string_intron1=substr ($sequence,($end2_new-1),$intron1_length);
											$intron1_length=length $string_intron1;
											$mark1=2;
										}

										if ($length_ref_exon2 % 3==0) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3-60+($j+$match)*3+($i-$match)*3+1;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3+1;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4+60-$j*3-$i*3-2;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4-2;
											}
											$intron2_length=$start4_new-$end3_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,($end3_new-1),$intron2_length);
												$intron2_length=length $string_intron2;
												$mark2=2;
											}
										}elsif ($length_ref_exon2 % 3==1) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3-60+($j+$match)*3+($i-$match)*3+2;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3+2;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4+60-$j*3-$i*3-1;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4-1;
											}
											$intron2_length=$start4_new-$end3_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,($end3_new+1),$intron2_length);
												$intron2_length=length $string_intron2;
												$mark2=2;
											}
										}elsif ($length_ref_exon2 % 3==2) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3-60+($j+$match)*3+($i-$match)*3;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4+60-$j*3-$i*3;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4;
											}
											$intron2_length=$start4_new-$end3_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,$end3_new,$intron2_length);
												$intron2_length=length $string_intron2;
												$mark2=2;
											}
										}
									}elsif ($length_ref_exon1 % 3==2) {
										for (my $i=$match;$i<$a;$i+=1) {
											my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string1,$j,$match);
												if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
													$end2_new=$end2-60+($j+$match)*3+($i-$match)*3+2;
													$marks++;
												}
											}
											if (defined $end2_new){
												last;
											}
										}
										if (!defined $end2_new) {
											$end2_new=$end2+2;
										}
										for (my $i=0;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string2,-$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$start3_new=$start3+60-$j*3-$i*3-1;
													$marks++;
												}
											}
											if (defined $start3_new){
												last;
											}
										}
										if (!defined $start3_new) {
											$start3_new=$start3-1;
										}
										$intron1_length=$start3_new-$end2_new-1;
										if ($intron1_length==0) {
											$mark1=1;
										}elsif (($intron1_length!=0) and ($intron1_length % 3==0)) {
											$string_intron1=substr ($sequence,($end2_new+1),$intron1_length);
											$intron1_length=length $string_intron1;
											$mark1=2;
										}

										if ($length_ref_exon2 % 3==0) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3-60+($j+$match)*3+($i-$match)*3+2;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3+2;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4+60-$j*3-$i*3-1;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4-1;
											}
											$intron2_length=$start4_new-$end3_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,($end3_new+1),$intron2_length);
												$intron2_length=length $string_intron2;
												$mark2=2;
											}
										}elsif ($length_ref_exon2 % 3==1) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3-60+($j+$match)*3+($i-$match)*3;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4+60-$j*3-$i*3;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4;
											}
											$intron2_length=$start4_new-$end3_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,$end3_new,$intron2_length);
												$intron2_length=length $string_intron2;
												$mark2=2;
											}
										}elsif ($length_ref_exon2 % 3==2) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3-60+($j+$match)*3+($i-$match)*3+1;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3+1;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4+60-$j*3-$i*3-2;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4-2;
											}
											$intron2_length=$start4_new-$end3_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,($end3_new-1),$intron2_length);
												$intron2_length=length $string_intron2;
												$mark2=2;
											}
										}
									}

									$str1=substr($sequence,($start-1),($end2_new-($start-1)));
									$str2=substr($sequence,($start3_new-1),($end3_new-($start3_new-1)));
									$str3=substr($sequence,($start4_new-1),($end-($start4_new-1)));
									$str=$str1.$str2.$str3;
									$length=length $str;
	
									my $aa;
									for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
										my $codon=substr ($str,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa.=$hash_codon{$codon};
										}else{
											$aa.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
	
									my $intron1_aa;
									if ($mark1==2) {
										for (my $i=0;$i<$intron1_length;$i+=3){
											my $codon=substr ($string_intron1,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$intron1_aa.=$hash_codon{$codon};
											}else{
												$intron1_aa.="X";
												my $j=$i+1;
												#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}
									}
									my $intron2_aa;
									if ($mark2==2) {
										for (my $i=0;$i<$intron2_length;$i+=3){
											my $codon=substr ($string_intron2,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$intron2_aa.=$hash_codon{$codon};
											}else{
												$intron2_aa.="X";
												my $j=$i+1;
												#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}
									}
	
									#if ((defined $end2_new) and (defined $start3_new) and (defined $end3_new) and (defined $start4_new)) {
									if (abs($end4-$start2) < 3000) {
										if (($start3 > $end2) and ($start4 > $end3)) {
											if ((($mark1==0) or ((defined $intron1_aa) and ($intron1_aa=~ /\*/))) and (($mark2==0) or ((defined $intron2_aa) and ($intron2_aa=~ /\*/))) and (($intron1_length > 30) and ($intron2_length > 30))) {# (two introns) spliced exon,not 3x or have stop codon
												print $out_annotation "     "."gene"."            ".$start."..".$end."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."join(".$start."..".$end2_new.",".$start3_new."..".$end3_new.",".$start4_new."..".$end.")"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
												$gene_number_seq{$name}++;
											}elsif ((($mark1==1) or ((defined $intron1_aa) and ($intron1_aa!~ /\*/))) and (($mark2==1) or ((defined $intron2_aa) and ($intron2_aa!~ /\*/))) and (($intron1_length < 30) and ($intron2_length < 30))) {# (no intron) joined exon,no intron or no stop codon
												print $out_annotation "     "."gene"."            ".$start."..".$end."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             ".$start."..".$end."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
												$gene_number_seq{$name}++;
												print $logfile "Warning: $name (positive two-intron PCG) lost two introns!\n";
											}elsif ((($mark1==0) or ((defined $intron1_aa) and ($intron1_aa=~ /\*/))) and (($mark2==1) or ((defined $intron2_aa) and ($intron2_aa!~ /\*/))) and (($intron1_length > 30) and ($intron2_length < 30))) {# (loss intron2)
												print $out_annotation "     "."gene"."            ".$start."..".$end."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."join(".$start."..".$end2_new.",".$start3_new."..".$end.")"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
												$gene_number_seq{$name}++;
												print $logfile "Warning: $name (positive two-intron PCG) lost intron 2!\n";
											}elsif ((($mark1==1) or ((defined $intron1_aa) and ($intron1_aa!~ /\*/))) and (($mark2==0) or ((defined $intron2_aa) and ($intron2_aa=~ /\*/))) and (($intron1_length < 30) and ($intron2_length > 30))) {# (loss intron1)
												print $out_annotation "     "."gene"."            ".$start."..".$end."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."join(".$start."..".$end3_new.",".$start4_new."..".$end.")"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
												$gene_number_seq{$name}++;
												print $logfile "Warning: $name (positive two-intron PCG) lost intron 1!\n";
											}
										}elsif (($start3 < $end2) or ($start4 < $end3)) {
											print $logfile "Warning: $name (positive two-intron PCG) has not been annotated!\n";
										}
									#}elsif ((!defined $end2_new) or (!defined $start3_new) or (!defined $end3_new) or (!defined $start4_new)) {
									}elsif (abs($end4-$start2) > 3000) {
										print $out_annotation "     "."gene"."            ".$start3_new."..".$end3_new."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             ".$start3_new."..".$end3_new."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
										$gene_number_seq{$name}++;
										print $logfile "Warning: $name (positive two-intron PCG) must be checked!\n";
									}
								}elsif ((!defined $left_start_position) and (defined $right_star_position)){
									my $end;
									if ($length_exon12 % 3==0){
										$end=($start4-1)+($left_star_position * 3+3);
									}
									if ($length_exon12 % 3==1){
										$end=($start4+1)+($left_star_position * 3+3);
									}
									if ($length_exon12 % 3==2){
										$end=($start4)+($left_star_position * 3+3);
									}
									if (abs($end4-$start2) < 3000) {
										if (($start3 > $end2) and ($start4 > $end3)) {
											print $out_annotation "     "."gene"."            ".$start2."..".$end."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."join(".$start2."..".$end2.",".$start3."..".$end3.",".$start4."..".$end.")"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (positive two-intron PCG) has non-canonical start codon!\n";
										}elsif (($start3 < $end2) or ($start4 < $end3)) {
											print $logfile "Warning: $name (positive two-intron PCG) has not been annotated!\n";
										}
									}elsif (abs($end4-$start2) > 3000) {
										print $out_annotation "     "."gene"."            ".$start3."..".$end3."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             ".$start3."..".$end3."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
										$gene_number_seq{$name}++;
										print $logfile "Warning: $name (positive two-intron PCG) must be checked!\n";
									}
								}
							}else{
								print "The exon length of $name (positive two-intron PCG) is not an intergral multiple of 3!\n";
							}
						}elsif($start3 > $end3){# negative
							my $str=substr($sequence,($end4-1),($start2-$end4+1));
							my $rev_com=reverse $str;
							$rev_com=~ tr/ACGTacgt/TGCAtgca/;
							my $length_3=length $str;
							my $reverse_start_codon=substr($rev_com,0,3);
							my $reverse_stop_codon=substr($rev_com,($length_3-3),3);
							my $str1=substr($sequence,($end2-1),($start2-$end2+1));
							my $str2=substr($sequence,($end3-1),($start3-$end3+1));
							my $str3=substr($sequence,($end4-1),($start4-$end4+1));
							my $rev_com1=reverse $str1;
							$rev_com1=~ tr/ACGTacgt/TGCAtgca/;
							my $rev_com2=reverse $str2;
							$rev_com2=~ tr/ACGTacgt/TGCAtgca/;
							my $rev_com3=reverse $str3;
							$rev_com3=~ tr/ACGTacgt/TGCAtgca/;
							my $str4=$rev_com1.$rev_com2.$rev_com3;
							my $length_exon=length $str4;
							my $aa_exon;
							for (my $i=0;$i<($length_exon-3);$i+=3){# delete stop codon
								my $codon=substr ($str4,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa_exon.=$hash_codon{$codon};
								}else{
									$aa_exon.="X";
									my $j=$i+1;
									#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa_exon=split //,$aa_exon;
							if (($length_exon % 3==0) and (!grep {$_=~ /\*/} @aa_exon) and (($reverse_start_codon eq "ATG") or ($reverse_start_codon eq "GTG")) and (($reverse_stop_codon eq "TAA") or ($reverse_stop_codon eq "TAG") or ($reverse_stop_codon eq "TGA"))){# standard start and stop codon
								my $seq_string1=substr($sequence,($end2-61),120);
								my $seq_string2=substr($sequence,($start3-60),120);
								my $seq_string3=substr($sequence,($end3-61),120);
								my $seq_string4=substr($sequence,($start4-60),120);
								my $rev_seq_string1=reverse $seq_string1;
								$rev_seq_string1=~ tr/ACGTacgt/TGCAtgca/;
								my $rev_seq_string2=reverse $seq_string2;
								$rev_seq_string2=~ tr/ACGTacgt/TGCAtgca/;
								my $rev_seq_string3=reverse $seq_string3;
								$rev_seq_string3=~ tr/ACGTacgt/TGCAtgca/;
								my $rev_seq_string4=reverse $seq_string4;
								$rev_seq_string4=~ tr/ACGTacgt/TGCAtgca/;
								my $aa_seq_string1;
								my $aa_seq_string2;
								my $aa_seq_string3;
								my $aa_seq_string4;
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($rev_seq_string1,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string1.=$hash_codon{$codon};
									}else{
										$aa_seq_string1.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($rev_seq_string2,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string2.=$hash_codon{$codon};
									}else{
										$aa_seq_string2.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($rev_seq_string3,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string3.=$hash_codon{$codon};
									}else{
										$aa_seq_string3.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($rev_seq_string4,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string4.=$hash_codon{$codon};
									}else{
										$aa_seq_string4.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}

								my ($str1,$str2,$str3,$str,$rev_com,$length,$end2_new,$start3_new,$end3_new,$start4_new);
								my ($string_intron1,$string_intron2,$intron1_length,$intron2_length,$string_intron1_rev_com,$string_intron2_rev_com);
								my $mark1=0;
								my $mark2=0;
								if ($length_ref_exon1 % 3==0) {
									for (my $i=$match;$i<$a;$i+=1) {
										my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
										my $marks=0;
										for (my $j=0;$j<40;$j+=1) {
											my $repeat=substr ($aa_seq_string1,$j,$match);
											if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
												$end2_new=$end2+60-($j+$match)*3-($i-$match)*3;
												$marks++;
											}
										}
										if (defined $end2_new){
											last;
										}
									}
									if (!defined $end2_new) {
										$end2_new=$end2;
									}
									for (my $i=0;$i<$b;$i+=1) {
										my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
										my $marks=0;
										for (my $j=$match;$j<40;$j+=1) {
											my $repeat=substr ($aa_seq_string2,-$j,$match);
											if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
												$start3_new=$start3-60+$j*3+$i*3;
												$marks++;
											}
										}
										if (defined $start3_new){
											last;
										}
									}
									if (!defined $start3_new) {
										$start3_new=$start3;
									}
									$intron1_length=$end2_new-$start3_new-1;
									if ($intron1_length==0) {
										$mark1=1;
									}
									if (($intron1_length!=0) and ($intron1_length % 3==0)) {
										$string_intron1=substr ($sequence,$start3_new,$intron1_length);
										$intron1_length=length $string_intron1;
										$string_intron1_rev_com=reverse $string_intron1;
										$string_intron1_rev_com=~ tr/ACGTacgt/TGCAtgca/;
										$mark1=2;
									}

									if ($length_ref_exon2 % 3==0) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3+60-($j+$match)*3-($i-$match)*3;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4-60+$j*3+$i*3;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4;
										}
										$intron2_length=$end3_new-$start4_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,$start4_new,$intron2_length);
											$intron2_length=length $string_intron2;
											$string_intron2_rev_com=reverse $string_intron2;
											$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark2=2;
										}
									}elsif ($length_ref_exon2 % 3==1) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3+60-($j+$match)*3-($i-$match)*3-1;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3-1;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4-60+$j*3+$i*3+2;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4+2;
										}
										$intron2_length=$end3_new-$start4_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,($start4_new+1),$intron2_length);
											$intron2_length=length $string_intron2;
											$string_intron2_rev_com=reverse $string_intron2;
											$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark2=2;
										}
									}elsif ($length_ref_exon2 % 3==2) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3+60-($j+$match)*3-($i-$match)*3-2;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3-2;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4-60+$j*3+$i*3+1;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4+1;
										}
										$intron2_length=$end3_new-$start4_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,($start4_new-1),$intron2_length);
											$intron2_length=length $string_intron2;
											$string_intron2_rev_com=reverse $string_intron2;
											$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark2=2;
										}
									}
								}elsif ($length_ref_exon1 % 3==1) {
									for (my $i=$match;$i<$a;$i+=1) {
										my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
										my $marks=0;
										for (my $j=0;$j<40;$j+=1) {
											my $repeat=substr ($aa_seq_string1,$j,$match);
											if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
												$end2_new=$end2+60-($j+$match)*3-($i-$match)*3-1;
												$marks++;
											}
										}
										if (defined $end2_new){
											last;
										}
									}
									if (!defined $end2_new) {
										$end2_new=$end2-1;
									}
									for (my $i=0;$i<$b;$i+=1) {
										my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
										my $marks=0;
										for (my $j=$match;$j<40;$j+=1) {
											my $repeat=substr ($aa_seq_string2,-$j,$match);
											if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
												$start3_new=$start3-60+$j*3+$i*3+2;
												$marks++;
											}
										}
										if (defined $start3_new){
											last;
										}
									}
									if (!defined $start3_new) {
										$start3_new=$start3+2;
									}
									$intron1_length=$end2_new-$start3_new-1;
									if ($intron1_length==0) {
										$mark1=1;
									}
									if (($intron1_length!=0) and ($intron1_length % 3==0)) {
										$string_intron1=substr ($sequence,($start3_new+1),$intron1_length);
										$intron1_length=length $string_intron1;
										$string_intron1_rev_com=reverse $string_intron1;
										$string_intron1_rev_com=~ tr/ACGTacgt/TGCAtgca/;
										$mark1=2;
									}

									if ($length_ref_exon2 % 3==0) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3+60-($j+$match)*3-($i-$match)*3-1;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3-1;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4-60+$j*3+$i*3+2;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4+2;
										}
										$intron2_length=$end3_new-$start4_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,($start4_new+1),$intron2_length);
											$intron2_length=length $string_intron2;
											$string_intron2_rev_com=reverse $string_intron2;
											$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark2=2;
										}
									}elsif ($length_ref_exon2 % 3==1) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3+60-($j+$match)*3-($i-$match)*3-2;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3-2;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4-60+$j*3+$i*3+1;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4+1;
										}
										$intron2_length=$end3_new-$start4_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,($start4_new-1),$intron2_length);
											$intron2_length=length $string_intron2;
											$string_intron2_rev_com=reverse $string_intron2;
											$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark2=2;
										}
									}elsif ($length_ref_exon2 % 3==2) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3+60-($j+$match)*3-($i-$match)*3;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4-60+$j*3+$i*3;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4;
										}
										$intron2_length=$end3_new-$start4_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,$start4_new,$intron2_length);
											$intron2_length=length $string_intron2;
											$string_intron2_rev_com=reverse $string_intron2;
											$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark2=2;
										}
									}
								}elsif ($length_ref_exon1 % 3==2) {
									for (my $i=$match;$i<$a;$i+=1) {
										my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
										my $marks=0;
										for (my $j=0;$j<40;$j+=1) {
											my $repeat=substr ($aa_seq_string1,$j,$match);
											if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
												$end2_new=$end2+60-($j+$match)*3-($i-$match)*3-2;
												$marks++;
											}
										}
										if (defined $end2_new){
											last;
										}
									}
									if (!defined $end2_new) {
										$end2_new=$end2-2;
									}
									for (my $i=0;$i<$b;$i+=1) {
										my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
										my $marks=0;
										for (my $j=$match;$j<40;$j+=1) {
											my $repeat=substr ($aa_seq_string2,-$j,$match);
											if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
												$start3_new=$start3-60+$j*3+$i*3+1;
												$marks++;
											}
										}
										if (defined $start3_new){
											last;
										}
									}
									if (!defined $start3_new) {
										$start3_new=$start3+1;
									}
									$intron1_length=$end2_new-$start3_new-1;
									if ($intron1_length==0) {
										$mark1=1;
									}elsif (($intron1_length!=0) and ($intron1_length % 3==0)) {
										$string_intron1=substr ($sequence,($start3_new-1),$intron1_length);
										$intron1_length=length $string_intron1;
										$string_intron1_rev_com=reverse $string_intron1;
										$string_intron1_rev_com=~ tr/ACGTacgt/TGCAtgca/;
										$mark1=2;
									}

									if ($length_ref_exon2 % 3==0) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3+60-($j+$match)*3-($i-$match)*3-2;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3-2;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4-60+$j*3+$i*3+1;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4+1;
										}
										$intron2_length=$end3_new-$start4_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,($start4_new-1),$intron2_length);
											$intron2_length=length $string_intron2;
											$string_intron2_rev_com=reverse $string_intron2;
											$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark2=2;
										}
									}elsif ($length_ref_exon2 % 3==1) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3+60-($j+$match)*3-($i-$match)*3;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4-60+$j*3+$i*3;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4;
										}
										$intron2_length=$end3_new-$start4_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,$start4_new,$intron2_length);
											$intron2_length=length $string_intron2;
											$string_intron2_rev_com=reverse $string_intron2;
											$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark2=2;
										}
									}elsif ($length_ref_exon2 % 3==2) {
										for (my $i=$match;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string3,$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$end3_new=$end3+60-($j+$match)*3-($i-$match)*3-1;
													$marks++;
												}
											}
											if (defined $end3_new){
												last;
											}
										}
										if (!defined $end3_new) {
											$end3_new=$end3-1;
										}
										for (my $i=0;$i<$c;$i+=1) {
											my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string4,-$j,$match);
												if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
													$start4_new=$start4-60+$j*3+$i*3+2;
													$marks++;
												}
											}
											if (defined $start4_new){
												last;
											}
										}
										if (!defined $start4_new) {
											$start4_new=$start4+2;
										}
										$intron2_length=$end3_new-$start4_new-1;
										if ($intron2_length==0) {
											$mark2=1;
										}
										if (($intron2_length!=0) and ($intron2_length % 3==0)) {
											$string_intron2=substr ($sequence,($start4_new+1),$intron2_length);
											$intron2_length=length $string_intron2;
											$string_intron2_rev_com=reverse $string_intron2;
											$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark2=2;
										}
									}
								}

								$str1=substr($sequence,($end2_new-1),($start2-($end2_new-1)));
								$str2=substr($sequence,($end3_new-1),($start3_new-($end3_new-1)));
								$str3=substr($sequence,($end4-1),($start4_new-($end4-1)));
								$str=$str1.$str2.$str3;
								$length=length $str;
								$rev_com=reverse $str;
								$rev_com=~ tr/ACGTacgt/TGCAtgca/;

								my $aa;
								for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
									my $codon=substr ($rev_com,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa.=$hash_codon{$codon};
									}else{
										$aa.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}

								my $intron1_aa;
								if ($mark1==2) {
									for (my $i=0;$i<$intron1_length;$i+=3){
										my $codon=substr ($string_intron1_rev_com,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$intron1_aa.=$hash_codon{$codon};
										}else{
											$intron1_aa.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
								}
								my $intron2_aa;
								if ($mark2==2) {
									for (my $i=0;$i<$intron2_length;$i+=3){
										my $codon=substr ($string_intron2_rev_com,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$intron2_aa.=$hash_codon{$codon};
										}else{
											$intron2_aa.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
								}

								#if ((defined $end2_new) and (defined $start3_new) and (defined $end3_new) and (defined $start4_new)) {
								if (abs($start2-$end4) < 3000) {
									if (($start3 < $end2) and ($start4 < $end3)) {
										if ((($mark1==0) or ((defined $intron1_aa) and ($intron1_aa=~ /\*/))) and (($mark2==0) or ((defined $intron2_aa) and ($intron2_aa=~ /\*/))) and (($intron1_length > 30) and ($intron2_length > 30))) {# (two introns) spliced exon,not 3x or have stop codon
											print $out_annotation "     "."gene"."            "."complement(".$end4."..".$start2.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."complement(join(".$end4."..".$start4_new.",".$end3_new."..".$start3_new.",".$end2_new."..".$start2."))"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
											$gene_number_seq{$name}++;
										}elsif ((($mark1==1) or ((defined $intron1_aa) and ($intron1_aa!~ /\*/))) and (($mark2==1) or ((defined $intron2_aa) and ($intron2_aa!~ /\*/))) and (($intron1_length < 30) and ($intron2_length < 30))) {# (no intron) joined exon,no intron or no stop codon
											print $out_annotation "     "."gene"."            "."complement(".$end4."..".$start2.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."complement(".$end4."..".$start2.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (negative two-intron PCG) lost two introns!\n";
										}elsif ((($mark1==0) or ((defined $intron1_aa) and ($intron1_aa=~ /\*/))) and (($mark2==1) or ((defined $intron2_aa) and ($intron2_aa!~ /\*/))) and (($intron1_length > 30) and ($intron2_length < 30))) {# (loss intron2)
											print $out_annotation "     "."gene"."            "."complement(".$end4."..".$start2.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."complement(join(".$end4."..".$start3_new.",".$end2_new."..".$start2."))"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (negative two-intron PCG) lost intron 2!\n";
										}elsif ((($mark1==1) or ((defined $intron1_aa) and ($intron1_aa!~ /\*/))) and (($mark2==0) or ((defined $intron2_aa) and ($intron2_aa=~ /\*/))) and (($intron1_length < 30) and ($intron2_length > 30))) {# (loss intron1)
											print $out_annotation "     "."gene"."            "."complement(".$end4."..".$start2.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."complement(join(".$end4."..".$start4_new.",".$end3_new."..".$start2."))"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (negative two-intron PCG) lost intron 1!\n";
										}
									}elsif (($start3 > $end2) or ($start4 > $end3)) {
										print $logfile "Warning: $name (negative two-intron PCG) has not been annotated!\n";
									}
								#}elsif ((!defined $end2_new) or (!defined $start3_new) or (!defined $end3_new) or (!defined $start4_new)) {
								}elsif (abs($start2-$end4) > 3000) {
									print $out_annotation "     "."gene"."            "."complement(".$end3_new."..".$start3_new.")\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."CDS"."             "."complement(".$end3_new."..".$start3_new.")"."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/codon_start=1"."\n";
									print $out_annotation "                     "."/transl_table=11"."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
									$gene_number_seq{$name}++;
									print $logfile "Warning: $name (negative two-intron PCG) must be checked!\n";
								}
							}elsif((grep {$_=~ /\*/} @aa_exon) or (($reverse_start_codon ne "ATG") and ($reverse_start_codon ne "GTG")) or (($reverse_stop_codon ne "TAA") and ($reverse_stop_codon ne "TAG") and ($reverse_stop_codon ne "TGA"))){# non-standard start or stop codon for _gene
								my $start_left=$start2-60;
								my $start_right;
								if (($length_cp-$start2)>=9000) {
									$start_right=$start2+9000;
								}elsif (($length_cp-$start2)<9000){
									$start_right=$start2+int(($length_cp-$start2)/3)*3;
								}
								my $length_exon1=($start2-$end2+1);
								my $length_exon2=($start3-$end3+1);
								my $length_exon12=$length_exon1+$length_exon2;
								my $end_left;
								if ($start4>=9002) {
									if ($length_exon12 % 3==0){
										$end_left=$start4-9000;
									}
									if ($length_exon12 % 3==1){
										$end_left=$start4-9002;
									}
									if ($length_exon12 % 3==2){
										$end_left=$start4-9001;
									}
								}elsif ($start4<9002){
									if ($length_exon12 % 3==0){
										$end_left=$start4-(int($start4/3)*3);
									}
									if ($length_exon12 % 3==1){
										$end_left=$start4-(int($start4/3)*3+2);
									}
									if ($length_exon12 % 3==2){
										$end_left=$start4-(int($start4/3)*3+1);
									}
								}
								my $str1=substr($sequence,($start_left-1),($start_right-$start_left+1));
								my $str2=substr($sequence,($end_left-1),($start4-$end_left+1));
								my $length1=length $str1;
								my $length2=length $str2;
								my $rev_com1=reverse $str1;
								$rev_com1=~ tr/ACGTacgt/TGCAtgca/;
								my $rev_com2=reverse $str2;
								$rev_com2=~ tr/ACGTacgt/TGCAtgca/;
								my $aa1;#left range for start codon
								for (my $i=0;$i<($length1-3);$i+=3){# delete stop codon
									my $codon=substr ($rev_com1,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa1.=$hash_codon{$codon};
									}else{
										$aa1.="X";
										my $m=$i+1;
										#print "Bad codon $codon in position $m of gene $name in species $contig!\n";
									}
								}
								my @aa1=split //,$aa1;
								my (@star1,@start1);
								foreach (0..$#aa1){
									if (($aa1[$_]=~ /\*/) and ($_ < 3000)){
										push @star1,$_;
									}
									if (($aa1[$_]=~ /M/) and ($_ <= 3060)){
										push @start1,$_;
									}
								}
								my $right_star_position=pop @star1;
								until ($right_star_position<(3000+$length_exon1/3)){
									$right_star_position=pop @star1;
								}
								my @start2;
								foreach (@start1){
									push @start2,$_ if ((defined $right_star_position) and ($_ > $right_star_position))
								}
								my $right_start_position=$start2[0];
								my $aa2;#right range for stop codon
								my $j;
								if ($length_exon1 % 3==0){
									$j=0;
								}
								if ($length_exon1 % 3==1){
									$j=2;
								}
								if ($length_exon1 % 3==2){
									$j=1;
								}
								for (my $k=$j;$k<($length2-3);$k+=3){# delete stop codon
									my $codon=substr ($rev_com2,$k,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa2.=$hash_codon{$codon};
									}else{
										$aa2.="X";
										my $n=$k+1;
										#print "Bad codon $codon in position $n of gene $name in species $contig!\n";
									}
								}
								my @aa2=split //,$aa2;
								my @star2;
								foreach (0..$#aa2){
									if ($aa2[$_]=~ /\*/){
										push @star2,$_;
									}
								}
								my $left_star_position=shift @star2;

								my $seq_string1=substr($sequence,($end2-61),120);
								my $seq_string2=substr($sequence,($start3-60),120);
								my $seq_string3=substr($sequence,($end3-61),120);
								my $seq_string4=substr($sequence,($start4-60),120);
								my $rev_seq_string1=reverse $seq_string1;
								$rev_seq_string1=~ tr/ACGTacgt/TGCAtgca/;
								my $rev_seq_string2=reverse $seq_string2;
								$rev_seq_string2=~ tr/ACGTacgt/TGCAtgca/;
								my $rev_seq_string3=reverse $seq_string3;
								$rev_seq_string3=~ tr/ACGTacgt/TGCAtgca/;
								my $rev_seq_string4=reverse $seq_string4;
								$rev_seq_string4=~ tr/ACGTacgt/TGCAtgca/;
								my $aa_seq_string1;
								my $aa_seq_string2;
								my $aa_seq_string3;
								my $aa_seq_string4;
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($rev_seq_string1,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string1.=$hash_codon{$codon};
									}else{
										$aa_seq_string1.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($rev_seq_string2,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string2.=$hash_codon{$codon};
									}else{
										$aa_seq_string2.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($rev_seq_string3,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string3.=$hash_codon{$codon};
									}else{
										$aa_seq_string3.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								for (my $i=0;$i<120;$i+=3){# delete stop codon
									my $codon=substr ($rev_seq_string4,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa_seq_string4.=$hash_codon{$codon};
									}else{
										$aa_seq_string4.="X";
										my $j=$i+1;
										#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}

								if ((defined $right_start_position) and (defined $left_star_position)){
									my $start=$start_right-$right_start_position * 3;
									my $end;
									if ($length_exon12 % 3==0){
										$end=($start4+1)-($left_star_position * 3+3);
									}
									if ($length_exon12 % 3==1){
										$end=($start4-1)-($left_star_position * 3+3);
									}
									if ($length_exon12 % 3==2){
										$end=($start4)-($left_star_position * 3+3);
									}

									my ($str1,$str2,$str3,$str,$rev_com,$length,$end2_new,$start3_new,$end3_new,$start4_new);
									my ($string_intron1,$string_intron2,$intron1_length,$intron2_length,$string_intron1_rev_com,$string_intron2_rev_com);
									my $mark1=0;
									my $mark2=0;
									if ($length_ref_exon1 % 3==0) {
										for (my $i=$match;$i<$a;$i+=1) {
											my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string1,$j,$match);
												if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
													$end2_new=$end2+60-($j+$match)*3-($i-$match)*3;
													$marks++;
												}
											}
											if (defined $end2_new){
												last;
											}
										}
										if (!defined $end2_new) {
											$end2_new=$end2;
										}
										for (my $i=0;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string2,-$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$start3_new=$start3-60+$j*3+$i*3;
													$marks++;
												}
											}
											if (defined $start3_new){
												last;
											}
										}
										if (!defined $start3_new) {
											$start3_new=$start3;
										}
										$intron1_length=$end2_new-$start3_new-1;
										if ($intron1_length==0) {
											$mark1=1;
										}
										if (($intron1_length!=0) and ($intron1_length % 3==0)) {
											$string_intron1=substr ($sequence,$start3_new,$intron1_length);
											$intron1_length=length $string_intron1;
											$string_intron1_rev_com=reverse $string_intron1;
											$string_intron1_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark1=2;
										}

										if ($length_ref_exon2 % 3==0) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3+60-($j+$match)*3-($i-$match)*3;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4-60+$j*3+$i*3;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4;
											}
											$intron2_length=$end3_new-$start4_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,$start4_new,$intron2_length);
												$intron2_length=length $string_intron2;
												$string_intron2_rev_com=reverse $string_intron2;
												$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark2=2;
											}
										}elsif ($length_ref_exon2 % 3==1) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3+60-($j+$match)*3-($i-$match)*3-1;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3-1;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4-60+$j*3+$i*3+2;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4+2;
											}
											$intron2_length=$end3_new-$start4_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,($start4_new+1),$intron2_length);
												$intron2_length=length $string_intron2;
												$string_intron2_rev_com=reverse $string_intron2;
												$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark2=2;
											}
										}elsif ($length_ref_exon2 % 3==2) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3+60-($j+$match)*3-($i-$match)*3-2;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3-2;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4-60+$j*3+$i*3+1;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4+1;
											}
											$intron2_length=$end3_new-$start4_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,($start4_new-1),$intron2_length);
												$intron2_length=length $string_intron2;
												$string_intron2_rev_com=reverse $string_intron2;
												$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark2=2;
											}
										}
									}elsif ($length_ref_exon1 % 3==1) {
										for (my $i=$match;$i<$a;$i+=1) {
											my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string1,$j,$match);
												if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
													$end2_new=$end2+60-($j+$match)*3-($i-$match)*3-1;
													$marks++;
												}
											}
											if (defined $end2_new){
												last;
											}
										}
										if (!defined $end2_new) {
											$end2_new=$end2-1;
										}
										for (my $i=0;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string2,-$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$start3_new=$start3-60+$j*3+$i*3+2;
													$marks++;
												}
											}
											if (defined $start3_new){
												last;
											}
										}
										if (!defined $start3_new) {
											$start3_new=$start3+2;
										}
										$intron1_length=$end2_new-$start3_new-1;
										if ($intron1_length==0) {
											$mark1=1;
										}
										if (($intron1_length!=0) and ($intron1_length % 3==0)) {
											$string_intron1=substr ($sequence,($start3_new+1),$intron1_length);
											$intron1_length=length $string_intron1;
											$string_intron1_rev_com=reverse $string_intron1;
											$string_intron1_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark1=2;
										}

										if ($length_ref_exon2 % 3==0) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3+60-($j+$match)*3-($i-$match)*3-1;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3-1;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4-60+$j*3+$i*3+2;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4+2;
											}
											$intron2_length=$end3_new-$start4_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,($start4_new+1),$intron2_length);
												$intron2_length=length $string_intron2;
												$string_intron2_rev_com=reverse $string_intron2;
												$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark2=2;
											}
										}elsif ($length_ref_exon2 % 3==1) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3+60-($j+$match)*3-($i-$match)*3-2;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3-2;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4-60+$j*3+$i*3+1;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4+1;
											}
											$intron2_length=$end3_new-$start4_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,($start4_new-1),$intron2_length);
												$intron2_length=length $string_intron2;
												$string_intron2_rev_com=reverse $string_intron2;
												$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark2=2;
											}
										}elsif ($length_ref_exon2 % 3==2) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3+60-($j+$match)*3-($i-$match)*3;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4-60+$j*3+$i*3;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4;
											}
											$intron2_length=$end3_new-$start4_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,$start4_new,$intron2_length);
												$intron2_length=length $string_intron2;
												$string_intron2_rev_com=reverse $string_intron2;
												$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark2=2;
											}
										}
									}elsif ($length_ref_exon1 % 3==2) {
										for (my $i=$match;$i<$a;$i+=1) {
											my $match_aa_ref_exon1=substr ($aa_ref_exon1,-$i,$match);
											my $marks=0;
											for (my $j=0;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string1,$j,$match);
												if (($repeat eq $match_aa_ref_exon1) and ($marks == 0)) {
													$end2_new=$end2+60-($j+$match)*3-($i-$match)*3-2;
													$marks++;
												}
											}
											if (defined $end2_new){
												last;
											}
										}
										if (!defined $end2_new) {
											$end2_new=$end2-2;
										}
										for (my $i=0;$i<$b;$i+=1) {
											my $match_aa_ref_exon2=substr ($aa_ref_exon2,$i,$match);
											my $marks=0;
											for (my $j=$match;$j<40;$j+=1) {
												my $repeat=substr ($aa_seq_string2,-$j,$match);
												if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
													$start3_new=$start3-60+$j*3+$i*3+1;
													$marks++;
												}
											}
											if (defined $start3_new){
												last;
											}
										}
										if (!defined $start3_new) {
											$start3_new=$start3+1;
										}
										$intron1_length=$end2_new-$start3_new-1;
										if ($intron1_length==0) {
											$mark1=1;
										}elsif (($intron1_length!=0) and ($intron1_length % 3==0)) {
											$string_intron1=substr ($sequence,($start3_new-1),$intron1_length);
											$intron1_length=length $string_intron1;
											$string_intron1_rev_com=reverse $string_intron1;
											$string_intron1_rev_com=~ tr/ACGTacgt/TGCAtgca/;
											$mark1=2;
										}

										if ($length_ref_exon2 % 3==0) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3+60-($j+$match)*3-($i-$match)*3-2;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3-2;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4-60+$j*3+$i*3+1;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4+1;
											}
											$intron2_length=$end3_new-$start4_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,($start4_new-1),$intron2_length);
												$intron2_length=length $string_intron2;
												$string_intron2_rev_com=reverse $string_intron2;
												$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark2=2;
											}
										}elsif ($length_ref_exon2 % 3==1) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3+60-($j+$match)*3-($i-$match)*3;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4-60+$j*3+$i*3;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4;
											}
											$intron2_length=$end3_new-$start4_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,$start4_new,$intron2_length);
												$intron2_length=length $string_intron2;
												$string_intron2_rev_com=reverse $string_intron2;
												$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark2=2;
											}
										}elsif ($length_ref_exon2 % 3==2) {
											for (my $i=$match;$i<$b;$i+=1) {
												my $match_aa_ref_exon2=substr ($aa_ref_exon2,-$i,$match);
												my $marks=0;
												for (my $j=0;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string3,$j,$match);
													if (($repeat eq $match_aa_ref_exon2) and ($marks == 0)) {
														$end3_new=$end3+60-($j+$match)*3-($i-$match)*3-1;
														$marks++;
													}
												}
												if (defined $end3_new){
													last;
												}
											}
											if (!defined $end3_new) {
												$end3_new=$end3-1;
											}
											for (my $i=0;$i<$c;$i+=1) {
												my $match_aa_ref_exon3=substr ($aa_ref_exon3,$i,$match);
												my $marks=0;
												for (my $j=$match;$j<40;$j+=1) {
													my $repeat=substr ($aa_seq_string4,-$j,$match);
													if (($repeat eq $match_aa_ref_exon3) and ($marks == 0)) {
														$start4_new=$start4-60+$j*3+$i*3+2;
														$marks++;
													}
												}
												if (defined $start4_new){
													last;
												}
											}
											if (!defined $start4_new) {
												$start4_new=$start4+2;
											}
											$intron2_length=$end3_new-$start4_new-1;
											if ($intron2_length==0) {
												$mark2=1;
											}
											if (($intron2_length!=0) and ($intron2_length % 3==0)) {
												$string_intron2=substr ($sequence,($start4_new+1),$intron2_length);
												$intron2_length=length $string_intron2;
												$string_intron2_rev_com=reverse $string_intron2;
												$string_intron2_rev_com=~ tr/ACGTacgt/TGCAtgca/;
												$mark2=2;
											}
										}
									}

									$str1=substr($sequence,($end2_new-1),($start-($end2_new-1)));
									$str2=substr($sequence,($end3_new-1),($start3_new-($end3_new-1)));
									$str3=substr($sequence,($end-1),($start4_new-($end-1)));
									$str=$str1.$str2.$str3;
									$length=length $str;
									$rev_com=reverse $str;
									$rev_com=~ tr/ACGTacgt/TGCAtgca/;

									my $aa;
									for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
										my $codon=substr ($rev_com,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa.=$hash_codon{$codon};
										}else{
											$aa.="X";
											my $j=$i+1;
											#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}

									my $intron1_aa;
									if ($mark1==2) {
										for (my $i=0;$i<$intron1_length;$i+=3){
											my $codon=substr ($string_intron1_rev_com,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$intron1_aa.=$hash_codon{$codon};
											}else{
												$intron1_aa.="X";
												my $j=$i+1;
												#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}
									}
									my $intron2_aa;
									if ($mark2==2) {
										for (my $i=0;$i<$intron2_length;$i+=3){
											my $codon=substr ($string_intron2_rev_com,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$intron2_aa.=$hash_codon{$codon};
											}else{
												$intron2_aa.="X";
												my $j=$i+1;
												#print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}
									}

									#if ((defined $end2_new) and (defined $start3_new) and (defined $end3_new) and (defined $start4_new)) {
									if (abs($start2-$end4) < 3000) {
										if (($start3 < $end2) and ($start4 < $end3)) {
											if ((($mark1==0) or ((defined $intron1_aa) and ($intron1_aa=~ /\*/))) and (($mark2==0) or ((defined $intron2_aa) and ($intron2_aa=~ /\*/))) and (($intron1_length > 30) and ($intron2_length > 30))) {# (two introns) spliced exon,not 3x or have stop codon
												print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."complement(join(".$end."..".$start4_new.",".$end3_new."..".$start3_new.",".$end2_new."..".$start."))"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
												$gene_number_seq{$name}++;
											}elsif ((($mark1==1) or ((defined $intron1_aa) and ($intron1_aa!~ /\*/))) and (($mark2==1) or ((defined $intron2_aa) and ($intron2_aa!~ /\*/))) and (($intron1_length < 30) and ($intron2_length < 30))) {# (no intron) joined exon,no intron or no stop codon
												print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."complement(".$end."..".$start.")\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
												$gene_number_seq{$name}++;
												print $logfile "Warning: $name (negative two-intron PCG) lost two introns!\n";
											}elsif ((($mark1==0) or ((defined $intron1_aa) and ($intron1_aa=~ /\*/))) and (($mark2==1) or ((defined $intron2_aa) and ($intron2_aa!~ /\*/))) and (($intron1_length > 30) and ($intron2_length < 30))) {# (loss intron2)
												print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."complement(join(".$end."..".$start3_new.",".$end2_new."..".$start."))"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
												$gene_number_seq{$name}++;
												print $logfile "Warning: $name (negative two-intron PCG) lost intron 2!\n";
											}elsif ((($mark1==1) or ((defined $intron1_aa) and ($intron1_aa!~ /\*/))) and (($mark2==0) or ((defined $intron2_aa) and ($intron2_aa=~ /\*/))) and (($intron1_length < 30) and ($intron2_length > 30))) {# (loss intron1)
												print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "     "."CDS"."             "."complement(join(".$end."..".$start4_new.",".$end3_new."..".$start."))"."\n";
												print $out_annotation "                     "."/gene=\"$name\""."\n";
												print $out_annotation "                     "."/codon_start=1"."\n";
												print $out_annotation "                     "."/transl_table=11"."\n";
												print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
												#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
												$gene_number_seq{$name}++;
												print $logfile "Warning: $name (negative two-intron PCG) lost intron 1!\n";
											}
										}elsif (($start3 > $end2) or ($start4 > $end3)) {
											print $logfile "Warning: $name (negative two-intron PCG) has not been annotated!\n";
										}
									#}elsif ((!defined $end2_new) or (!defined $start3_new) or (!defined $end3_new) or (!defined $start4_new)) {
									}elsif (abs($start2-$end4) > 3000) {
										print $out_annotation "     "."gene"."            "."complement(".$end3_new."..".$start3_new.")\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."complement(".$end3_new."..".$start3_new.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
										$gene_number_seq{$name}++;
										print $logfile "Warning: $name (negative two-intron PCG) must be checked!\n";
									}
								}elsif ((!defined $right_start_position) and (defined $left_star_position)){
									my $end;
									if ($length_exon12 % 3==0){
										$end=($start4+1)-($left_star_position * 3+3);
									}
									if ($length_exon12 % 3==1){
										$end=($start4-1)-($left_star_position * 3+3);
									}
									if ($length_exon12 % 3==2){
										$end=($start4)-($left_star_position * 3+3);
									}
									if (abs($start2-$end4) < 3000) {
										if (($start3 < $end2) and ($start4 < $end3)) {
											print $out_annotation "     "."gene"."            "."complement(".$end."..".$start2.")\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "     "."CDS"."             "."complement(join(".$end."..".$start4.",".$end3."..".$start3.",".$end2."..".$start2."))"."\n";
											print $out_annotation "                     "."/gene=\"$name\""."\n";
											print $out_annotation "                     "."/codon_start=1"."\n";
											print $out_annotation "                     "."/transl_table=11"."\n";
											print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
											#print $out_annotation "                     "."/translation=\"$aa\""."\n";
											$gene_number_seq{$name}++;
											print $logfile "Warning: $name (negative two-intron PCG) has non-canonical start codon!\n";
										}elsif (($start3 > $end2) or ($start4 > $end3)) {
											print $logfile "Warning: $name (negative two-intron PCG) has not been annotated!\n";
										}
									}elsif (abs($start2-$end4) > 3000) {
										print $out_annotation "     "."gene"."            "."complement(".$end3."..".$start3.")\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."complement(".$end3."..".$start3.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
										$gene_number_seq{$name}++;
										print $logfile "Warning: $name (negative two-intron PCG) must be checked!\n";
									}
								}
							}else{
								print "The exon length of $name (negative two-intron PCG) is not an intergral multiple of 3!\n";
							}
						}
					}
				}
			}
		}
	}
	print $out_annotation "ORIGIN      \n";
	$sequence=lc $sequence;
	for (my $x=0;$x<($length_cp);$x+=60) {
		my $y=$x+1;
		my $number=length $y;
		if ($number==1) {
			print $out_annotation "        $y";
		}elsif ($number==2) {
			print $out_annotation "       $y";
		}elsif ($number==3) {
			print $out_annotation "      $y";
		}elsif ($number==4) {
			print $out_annotation "     $y";
		}elsif ($number==5) {
			print $out_annotation "    $y";
		}elsif ($number==6) {
			print $out_annotation "   $y";
		}
		my $row1=substr($sequence,$x,60);
		my $length=length $row1;
		my $row2;
		for (my $z=0;$z<$length;$z+=10) {
			$row2=substr($row1,$z,10);
			print $out_annotation " $row2";
		}
		print $out_annotation "\n";
	}
	print $out_annotation "//\n";
	close $out_annotation;


	{
	if ($link eq "N") {
		open (my $input_temp,"<","$output_directory/$temp.gb");
		open (my $output_temp,">","$output_directory/$filename.gb");
		while (<$input_temp>) {
			$_=~ s/\r|\n//g;
			$_=~ s/rps12\+1/rps12/g;
			$_=~ s/rps12\+2/rps12/g;
			print $output_temp "$_\n";
		}
		close $input_temp;
		close $output_temp;
		unlink("$output_directory/$temp.gb");
	}elsif ($link eq "Y") {
		my (@row_array,$species_name,$length,$element,@genearray,@output1);
		my $mark=0;
		open (my $input_temp,"<","$output_directory/$temp.gb");
		while (<$input_temp>){
			$_=~ s/\r|\n//g;
			@row_array=split /\s+/,$_;
			if (/^LOCUS/i){
				$species_name=$row_array[1];
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
				}
				$element=$row_array[2];
				$mark=1;
			}elsif(/^\s+\/gene="(rps12.*)"/ and $mark == 1){
				$element=$species_name.":".$1.":".$element;
				push @genearray,$element;
				$element=();
				$mark=0;
			}elsif(/ {5}CDS {13}/){
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
				}
				$element=$row_array[2];
				$mark=2;
			}elsif(/^\s+\/gene="(rps12.*)"/ and $mark == 2){
				$element=$species_name.":".$1.":".$element;
				push @genearray,$element;
				$element=();
				$mark=0;
			}
		}
		close $input_temp;
		#print Dumper \@genearray;


		foreach (@genearray){
			my @array=split /:/,$_;
			push @output1,"$array[0]\t$array[1]\t$array[2]\n";
		}
		@row_array=();
		@genearray=();
		#print Dumper \@output1;


		#sort position information of rps12 gene
		my $col=3;
		my (%sort,@output2);
		foreach (@output1){
			my @row=split /\t/,$_;
			$sort{$_}=$row[$col];
		}
		foreach (sort {$sort{$a} <=> $sort{$b}} keys %sort){
			push @output2,"$_\n";
		}
		@output1=();
		#print Dumper \@output2;


		#separately output gene and CDS features
		my (%SP01,%GENE01,%STRAND01,%START01,%END01,%TYPE01,%STRAND02,%START02,%END02,%TYPE02,@output3,@output4);
		my $cnt1=0;
		foreach (@output2) {
			$_=~ s/\r|\n//g;
			$cnt1++;
			($SP01{$cnt1},$GENE01{$cnt1},$STRAND01{$cnt1},$START01{$cnt1},$END01{$cnt1},$TYPE01{$cnt1},$STRAND02{$cnt1},$START02{$cnt1},$END02{$cnt1},$TYPE02{$cnt1})=(split /\s+/,$_)[0,1,2,3,4,5,6,7,8,9];
			if ($TYPE01{$cnt1} eq "gene") {
				push @output3,$_;
			}elsif ($TYPE01{$cnt1} eq "CDS") {
				push @output4,$_;
			}
		}
		#print Dumper \@output2;
		@output2=();
		#print Dumper \%STRAND02;

		my @rps12;


		#output gene feature of rps12
		my (%SP1,%GENE1,%STRAND1,%START1,%END1,%TYPE1);
		my $cnt2=0;
		foreach (@output3) {
			$_=~ s/\r|\n//g;
			$cnt2++;
			($SP1{$cnt2},$GENE1{$cnt2},$STRAND1{$cnt2},$START1{$cnt2},$END1{$cnt2},$TYPE1{$cnt2})=(split /\s+/,$_)[0,1,2,3,4,5];
		}
		#print Dumper \%STRAND1;

		if ($cnt2 == 3) {
			foreach (1..$cnt2) {
				my $last1=$_;
				my $last2=$_+1;
				my $last3=$_+2;
				#print $last1."\t".$last2."\t".$last3."\n";
				my $length_rps12_1=abs($END1{$last1}-$START1{$last1}+1);
				my $length_rps12_2=abs($END1{$last2}-$START1{$last2}+1);
				my $length_rps12_3=abs($END1{$last3}-$START1{$last3}+1);
				#print $length_rps12_1."\t".$length_rps12_2."\t".$length_rps12_3."\n";
				#print $STRAND1{$last1}."\t".$STRAND1{$last2}."\t".$STRAND1{$last3}."\n";
				if (($STRAND1{$last1} eq "-") and ($STRAND1{$last2} eq "-") and ($STRAND1{$last3} eq "+")) {
					if ((($length_rps12_1 > 80) and ($length_rps12_1 < 200)) and ($length_rps12_2 > 200) and ($length_rps12_3 > 200)) {
						push @rps12,"     gene            complement(join(".$START1{$last2}."..".$END1{$last2}.",".$START1{$last1}."..".$END1{$last1}."))"."\n";
						push @rps12,"                     /gene=\"rps12\""."\n";
						push @rps12,"     gene            join(complement(".$START1{$last1}."..".$END1{$last1}."),".$START1{$last3}."..".$END1{$last3}.")"."\n";
						push @rps12,"                     /gene=\"rps12\""."\n";
					}
				}elsif (($STRAND1{$last1} eq "+") and ($STRAND1{$last2} eq "+") and ($STRAND1{$last3} eq "-")) {
					if ((($length_rps12_1 > 80) and ($length_rps12_1 < 200)) and ($length_rps12_2 > 200) and ($length_rps12_3 > 200)) {
						push @rps12,"     gene            join(".$START1{$last1}."..".$END1{$last1}.",".$START1{$last2}."..".$END1{$last2}.")"."\n";
						push @rps12,"                     /gene=\"rps12\""."\n";
						push @rps12,"     gene            join(".$START1{$last1}."..".$END1{$last1}.",complement(".$START1{$last3}."..".$END1{$last3}."))"."\n";
						push @rps12,"                     /gene=\"rps12\""."\n";
					}
				}elsif (($STRAND1{$last1} eq "-") and ($STRAND1{$last2} eq "+") and ($STRAND1{$last3} eq "-")) {
					if ((($length_rps12_1 > 80) and ($length_rps12_1 < 200)) and ($length_rps12_2 > 200) and ($length_rps12_3 > 200)) {
						push @rps12,"     gene            join(complement(".$START1{$last1}."..".$END1{$last1}."),".$START1{$last2}."..".$END1{$last2}.")"."\n";
						push @rps12,"                     /gene=\"rps12\""."\n";
						push @rps12,"     gene            complement(join(".$START1{$last3}."..".$END1{$last3}.",".$START1{$last1}."..".$END1{$last1}."))"."\n";
						push @rps12,"                     /gene=\"rps12\""."\n";
					}
				}elsif (($STRAND1{$last1} eq "+") and ($STRAND1{$last2} eq "-") and ($STRAND1{$last3} eq "+")) {
					if ((($length_rps12_1 > 80) and ($length_rps12_1 < 200)) and ($length_rps12_2 > 200) and ($length_rps12_3 > 200)) {
						push @rps12,"     gene            join(".$START1{$last1}."..".$END1{$last1}.",complement(".$START1{$last2}."..".$END1{$last2}."))"."\n";
						push @rps12,"                     /gene=\"rps12\""."\n";
						push @rps12,"     gene            join(".$START1{$last1}."..".$END1{$last1}.",".$START1{$last3}."..".$END1{$last3}.")"."\n";
						push @rps12,"                     /gene=\"rps12\""."\n";
					}
				}
				last if ($last3 == 3);
			}
		}elsif ($cnt2 == 2) {
			foreach (1..$cnt2) {
				my $last1=$_;
				my $last2=$_+1;
				#print $last1."\t".$last2."\n";
				my $length_rps12_1=abs($END1{$last1}-$START1{$last1}+1);
				my $length_rps12_2=abs($END1{$last2}-$START1{$last2}+1);
				#print $length_rps12_1."\t".$length_rps12_2."\n";
				#print $STRAND1{$last1}."\t".$STRAND1{$last2}."\n";
				if (($STRAND1{$last1} eq "-") and ($STRAND1{$last2} eq "-")) {
					if ((($length_rps12_1 > 80) and ($length_rps12_1 < 200)) and ($length_rps12_2 > 200)) {
						push @rps12,"     gene            complement(join(".$START1{$last2}."..".$END1{$last2}.",".$START1{$last1}."..".$END1{$last1}."))"."\n";
						push @rps12,"                     /gene=\"rps12\""."\n";
					}
				}elsif (($STRAND1{$last1} eq "+") and ($STRAND1{$last2} eq "+")) {
					if ((($length_rps12_1 > 80) and ($length_rps12_1 < 200)) and ($length_rps12_2 > 200)) {
						push @rps12,"     gene            join(".$START1{$last1}."..".$END1{$last1}.",".$START1{$last2}."..".$END1{$last2}.")"."\n";
						push @rps12,"                     /gene=\"rps12\""."\n";
					}
				}elsif (($STRAND1{$last1} eq "-") and ($STRAND1{$last2} eq "+")) {
					if ((($length_rps12_1 > 80) and ($length_rps12_1 < 200)) and ($length_rps12_2 > 200)) {
						push @rps12,"     gene            join(complement(".$START1{$last1}."..".$END1{$last1}."),".$START1{$last2}."..".$END1{$last2}.")"."\n";
						push @rps12,"                     /gene=\"rps12\""."\n";
					}elsif (($length_rps12_1 > 200) and ($length_rps12_2 > 200)) {
						push @rps12,"     gene            complement(".$START1{$last1}."..".$END1{$last1}.")"."\n";
						push @rps12,"                     /gene=\"rps12\""."\n";
						push @rps12,"     gene            ".$START1{$last2}."..".$END1{$last2}."\n";
						push @rps12,"                     /gene=\"rps12\""."\n";
					}
				}elsif (($STRAND1{$last1} eq "+") and ($STRAND1{$last2} eq "-")) {
					if ((($length_rps12_1 > 80) and ($length_rps12_1 < 200)) and ($length_rps12_2 > 200)) {
						push @rps12,"     gene            join(".$START1{$last1}."..".$END1{$last1}.",complement(".$START1{$last2}."..".$END1{$last2}."))"."\n";
						push @rps12,"                     /gene=\"rps12\""."\n";
					}elsif (($length_rps12_1 > 200) and ($length_rps12_2 > 200)) {
						push @rps12,"     gene            ".$START1{$last1}."..".$END1{$last1}."\n";
						push @rps12,"                     /gene=\"rps12\""."\n";
						push @rps12,"     gene            complement(".$START1{$last2}."..".$END1{$last2}.")"."\n";
						push @rps12,"                     /gene=\"rps12\""."\n";
					}
				}
				last if ($last2 == 2);
			}
		}elsif ($cnt2 == 1) {
			foreach (1..$cnt2) {
				my $last1=$_;
				#print $last1."\n";
				my $length_rps12_1=abs($END1{$last1}-$START1{$last1}+1);
				#print $length_rps12_1."\n";
				#print $STRAND1{$last1}."\n";
				if (($STRAND1{$last1} eq "-")) {
					if ((($length_rps12_1 > 80) and ($length_rps12_1 < 200))) {
						push @rps12,"     gene            complement(".$START1{$last1}."..".$END1{$last1}.")"."\n";
						push @rps12,"                     /gene=\"rps12\""."\n";
					}elsif ($length_rps12_1 > 200) {
						push @rps12,"     gene            complement(".$START1{$last1}."..".$END1{$last1}.")"."\n";
						push @rps12,"                     /gene=\"rps12\""."\n";
					}
				}elsif (($STRAND1{$last1} eq "+")) {
					if ((($length_rps12_1 > 80) and ($length_rps12_1 < 200)) ) {
						push @rps12,"     gene            ".$START1{$last1}."..".$END1{$last1}."\n";
						push @rps12,"                     /gene=\"rps12\""."\n";
					}elsif ($length_rps12_1 > 200) {
						push @rps12,"     gene            ".$START1{$last1}."..".$END1{$last1}."\n";
						push @rps12,"                     /gene=\"rps12\""."\n";
					}
				}
				last if ($last1 == 1);
			}
		}




		#output CDS feature of rps12
		my (%SP2,%GENE2,%STRAND2,%START2,%END2,%TYPE2,%STRAND3,%START3,%END3,%TYPE3);
		my $cnt3=0;
		foreach (@output4) {
			$_=~ s/\r|\n//g;
			$cnt3++;
			($SP2{$cnt3},$GENE2{$cnt3},$STRAND2{$cnt3},$START2{$cnt3},$END2{$cnt3},$TYPE2{$cnt3},$STRAND3{$cnt3},$START3{$cnt3},$END3{$cnt3},$TYPE3{$cnt3})=(split /\s+/,$_)[0,1,2,3,4,5,6,7,8,9];
		}
		#print Dumper \%STRAND3;

		if ($cnt3 == 3) {
			foreach (1..$cnt3) {
				my $last1=$_;
				my $last2=$_+1;
				my $last3=$_+2;
				#print $last1."\t".$last2."\t".$last3."\n";
				my $length_rps12_1=abs($END2{$last1}-$START2{$last1}+1);
				my $length_rps12_2=abs($END2{$last2}-$START2{$last2}+1);
				my $length_rps12_3=abs($END2{$last3}-$START2{$last3}+1);
				#print $length_rps12_1."\t".$length_rps12_2."\t".$length_rps12_3."\n";
				#print $STRAND2{$last1}."\t".$STRAND2{$last2}."\t".$STRAND2{$last3}."\n";
				my ($length_rps12_4,$length_rps12_5);
				if (defined $STRAND3{$last2}) {
					$length_rps12_4=abs($END3{$last2}-$START3{$last2}+1);
				}
				if (defined $STRAND3{$last3}) {
					$length_rps12_5=abs($END3{$last3}-$START3{$last3}+1);
				}
				if (($STRAND2{$last1} eq "-") and ($STRAND2{$last2} eq "-") and ($STRAND2{$last3} eq "+")) {
					if ((defined $STRAND3{$last2}) and (defined $STRAND3{$last3})) {
						if ((($length_rps12_1 > 80) and ($length_rps12_1 < 200)) and ($length_rps12_2 < 80) and ($length_rps12_3 > 200) and ($length_rps12_4 > 200) and ($length_rps12_5 < 80)) {
							push @rps12,"     CDS             complement(join(".$START2{$last2}."..".$END2{$last2}.",".$START3{$last2}."..".$END3{$last2}.",".$START2{$last1}."..".$END2{$last1}."))"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
							push @rps12,"     CDS             join(complement(".$START2{$last1}."..".$END2{$last1}."),".$START2{$last3}."..".$END2{$last3}.",".$START3{$last3}."..".$END3{$last3}.")"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}elsif ((!defined $STRAND3{$last2}) and (!defined $STRAND3{$last3})) {
						if ((($length_rps12_1 > 80) and ($length_rps12_1 < 200)) and ($length_rps12_2 > 200) and ($length_rps12_3 > 200)) {
							push @rps12,"     CDS             complement(join(".$START2{$last2}."..".$END2{$last2}.",".$START2{$last1}."..".$END2{$last1}."))"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
							push @rps12,"     CDS             join(complement(".$START2{$last1}."..".$END2{$last1}."),".$START2{$last3}."..".$END2{$last3}.")"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}elsif ((defined $STRAND3{$last2}) and (!defined $STRAND3{$last3})) {
						if ((($length_rps12_1 > 80) and ($length_rps12_1 < 200)) and ($length_rps12_2 < 80) and ($length_rps12_3 > 200) and ($length_rps12_4 > 200)) {
							push @rps12,"     CDS             complement(join(".$START2{$last2}."..".$END2{$last2}.",".$START3{$last2}."..".$END3{$last2}.",".$START2{$last1}."..".$END2{$last1}."))"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
							push @rps12,"     CDS             join(complement(".$START2{$last1}."..".$END2{$last1}."),".$START2{$last3}."..".$END2{$last3}.")"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}elsif ((!defined $STRAND3{$last2}) and (defined $STRAND3{$last3})) {
						if ((($length_rps12_1 > 80) and ($length_rps12_1 < 200)) and ($length_rps12_2 < 80) and ($length_rps12_3 > 200) and ($length_rps12_5 < 80)) {
							push @rps12,"     CDS             complement(join(".$START2{$last2}."..".$END2{$last2}.",".$START2{$last1}."..".$END2{$last1}."))"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
							push @rps12,"     CDS             join(complement(".$START2{$last1}."..".$END2{$last1}."),".$START2{$last3}."..".$END2{$last3}.",".$START3{$last3}."..".$END3{$last3}.")"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}
				}elsif (($STRAND2{$last1} eq "+") and ($STRAND2{$last2} eq "+") and ($STRAND2{$last3} eq "-")) {
					if ((defined $STRAND3{$last2}) and (defined $STRAND3{$last3})) {
						if ((($length_rps12_1 > 80)and ($length_rps12_1 < 200)) and ($length_rps12_2 > 200) and ($length_rps12_3 < 80) and ($length_rps12_4 < 80) and ($length_rps12_5 > 200)) {
							push @rps12,"     CDS             join(".$START2{$last1}."..".$END2{$last1}.",".$START2{$last2}."..".$END2{$last2}.",".$START3{$last2}."..".$END3{$last2}.")"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
							push @rps12,"     CDS             join(".$START2{$last1}."..".$END2{$last1}.",complement(".$START2{$last3}."..".$END2{$last3}.",".$START3{$last3}."..".$END3{$last3}."))"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}elsif ((!defined $STRAND3{$last2}) and (!defined $STRAND3{$last3})) {
						if ((($length_rps12_1 > 80)and ($length_rps12_1 < 200)) and ($length_rps12_2 > 200) and ($length_rps12_3 > 200)) {
							push @rps12,"     CDS             join(".$START2{$last1}."..".$END2{$last1}.",".$START2{$last2}."..".$END2{$last2}.")"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
							push @rps12,"     CDS             join(".$START2{$last1}."..".$END2{$last1}.",complement(".$START2{$last3}."..".$END2{$last3}."))"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}elsif ((defined $STRAND3{$last2}) and (!defined $STRAND3{$last3})) {
						if ((($length_rps12_1 > 80)and ($length_rps12_1 < 200)) and ($length_rps12_2 > 200) and ($length_rps12_3 < 80) and ($length_rps12_4 < 80)) {
							push @rps12,"     CDS             join(".$START2{$last1}."..".$END2{$last1}.",".$START2{$last2}."..".$END2{$last2}.",".$START3{$last2}."..".$END3{$last2}.")"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
							push @rps12,"     CDS             join(".$START2{$last1}."..".$END2{$last1}.",complement(".$START2{$last3}."..".$END2{$last3}."))"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}elsif ((!defined $STRAND3{$last2}) and (defined $STRAND3{$last3})) {
						if ((($length_rps12_1 > 80)and ($length_rps12_1 < 200)) and ($length_rps12_2 > 200) and ($length_rps12_3 < 80) and ($length_rps12_5 > 200)) {
							push @rps12,"     CDS             join(".$START2{$last1}."..".$END2{$last1}.",".$START2{$last2}."..".$END2{$last2}.")"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
							push @rps12,"     CDS             join(".$START2{$last1}."..".$END2{$last1}.",complement(".$START2{$last3}."..".$END2{$last3}.",".$START3{$last3}."..".$END3{$last3}."))"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}
				}elsif (($STRAND2{$last1} eq "-") and ($STRAND2{$last2} eq "+") and ($STRAND2{$last3} eq "-")) {
					if ((defined $STRAND3{$last2}) and (defined $STRAND3{$last3})) {
						if ((($length_rps12_1 > 80) and ($length_rps12_1 < 200)) and ($length_rps12_2 > 200) and ($length_rps12_3 < 80) and ($length_rps12_4 < 80) and ($length_rps12_5 > 200)) {
							push @rps12,"     CDS             join(complement(".$START2{$last1}."..".$END2{$last1}."),".$START2{$last2}."..".$END2{$last2}.",".$START3{$last2}."..".$END3{$last2}.")"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
							push @rps12,"     CDS             complement(join(".$START2{$last3}."..".$END2{$last3}.",".$START3{$last3}."..".$END3{$last3}.",".$START2{$last1}."..".$END2{$last1}."))"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}elsif ((!defined $STRAND3{$last2}) and (!defined $STRAND3{$last3})) {
						if ((($length_rps12_1 > 80) and ($length_rps12_1 < 200)) and ($length_rps12_2 > 200) and ($length_rps12_3 > 200)) {
							push @rps12,"     CDS             join(complement(".$START2{$last1}."..".$END2{$last1}."),".$START2{$last2}."..".$END2{$last2}.")"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
							push @rps12,"     CDS             complement(join(".$START2{$last3}."..".$END2{$last3}.",".$START2{$last1}."..".$END2{$last1}."))"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}elsif ((defined $STRAND3{$last2}) and (!defined $STRAND3{$last3})) {
						if ((($length_rps12_1 > 80) and ($length_rps12_1 < 200)) and ($length_rps12_2 > 200) and ($length_rps12_3 < 80) and ($length_rps12_4 < 80)) {
							push @rps12,"     CDS             join(complement(".$START2{$last1}."..".$END2{$last1}."),".$START2{$last2}."..".$END2{$last2}.",".$START3{$last2}."..".$END3{$last2}.")"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
							push @rps12,"     CDS             complement(join(".$START2{$last3}."..".$END2{$last3}.",".$START2{$last1}."..".$END2{$last1}."))"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}elsif ((!defined $STRAND3{$last2}) and (defined $STRAND3{$last3})) {
						if ((($length_rps12_1 > 80) and ($length_rps12_1 < 200)) and ($length_rps12_2 > 200) and ($length_rps12_3 < 80) and ($length_rps12_5 > 200)) {
							push @rps12,"     CDS             join(complement(".$START2{$last1}."..".$END2{$last1}."),".$START2{$last2}."..".$END2{$last2}.")"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
							push @rps12,"     CDS             complement(join(".$START2{$last3}."..".$END2{$last3}.",".$START3{$last3}."..".$END3{$last3}.",".$START2{$last1}."..".$END2{$last1}."))"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}
				}elsif (($STRAND2{$last1} eq "+") and ($STRAND2{$last2} eq "-") and ($STRAND2{$last3} eq "+")) {
					if ((defined $STRAND3{$last2}) and (defined $STRAND3{$last3})) {
						if ((($length_rps12_1 > 80) and ($length_rps12_1 < 200)) and ($length_rps12_2 < 80) and ($length_rps12_3 > 200) and ($length_rps12_4 > 200) and ($length_rps12_5 < 80)) {
							push @rps12,"     CDS             join(".$START2{$last1}."..".$END2{$last1}.",complement(".$START2{$last2}."..".$END2{$last2}.",".$START3{$last2}."..".$END3{$last2}."))"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
							push @rps12,"     CDS             join(".$START2{$last1}."..".$END2{$last1}.",".$START2{$last3}."..".$END2{$last3}.",".$START3{$last3}."..".$END3{$last3}.")"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}elsif ((!defined $STRAND3{$last2}) and (!defined $STRAND3{$last3})) {
						if ((($length_rps12_1 > 80) and ($length_rps12_1 < 200)) and ($length_rps12_2 > 200) and ($length_rps12_3 > 200)) {
							push @rps12,"     CDS             join(".$START2{$last1}."..".$END2{$last1}.",complement(".$START2{$last2}."..".$END2{$last2}."))"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
							push @rps12,"     CDS             join(".$START2{$last1}."..".$END2{$last1}.",".$START2{$last3}."..".$END2{$last3}.")"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}elsif ((defined $STRAND3{$last2}) and (!defined $STRAND3{$last3})) {
						if ((($length_rps12_1 > 80) and ($length_rps12_1 < 200)) and ($length_rps12_2 < 80) and ($length_rps12_3 > 200) and ($length_rps12_4 > 200)) {
							push @rps12,"     CDS             join(".$START2{$last1}."..".$END2{$last1}.",complement(".$START2{$last2}."..".$END2{$last2}.",".$START3{$last2}."..".$END3{$last2}."))"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
							push @rps12,"     CDS             join(".$START2{$last1}."..".$END2{$last1}.",".$START2{$last3}."..".$END2{$last3}.")"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}elsif ((!defined $STRAND3{$last2}) and (defined $STRAND3{$last3})) {
						if ((($length_rps12_1 > 80) and ($length_rps12_1 < 200)) and ($length_rps12_2 < 80) and ($length_rps12_3 > 200) and ($length_rps12_5 < 80)) {
							push @rps12,"     CDS             join(".$START2{$last1}."..".$END2{$last1}.",complement(".$START2{$last2}."..".$END2{$last2}."))"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
							push @rps12,"     CDS             join(".$START2{$last1}."..".$END2{$last1}.",".$START2{$last3}."..".$END2{$last3}.",".$START3{$last3}."..".$END3{$last3}.")"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}
				}
				last if ($last3 == 3);
			}
		}elsif ($cnt3 == 2) {
			foreach (1..$cnt3) {
				my $last1=$_;
				my $last2=$_+1;
				#print $last1."\t".$last2."\n";
				my $length_rps12_1=abs($END2{$last1}-$START2{$last1}+1);
				my $length_rps12_2=abs($END2{$last2}-$START2{$last2}+1);
				#print $length_rps12_1."\t".$length_rps12_2."\n";
				#print $STRAND2{$last1}."\t".$STRAND2{$last2}."\n";
				my ($length_rps12_4,$length_rps12_5);
				if (defined $STRAND3{$last1}) {
					$length_rps12_4=abs($END3{$last1}-$START3{$last1}+1);
				}
				if (defined $STRAND3{$last2}) {
					$length_rps12_5=abs($END3{$last2}-$START3{$last2}+1);
				}
				if (($STRAND2{$last1} eq "-") and ($STRAND2{$last2} eq "-")) {
					if ((!defined $STRAND3{$last1}) and (defined $STRAND3{$last2})) {
						if ((($length_rps12_1 > 80)and ($length_rps12_1 < 200)) and ($length_rps12_2 < 80) and ($length_rps12_5 > 200)) {
							push @rps12,"     CDS             complement(join(".$START2{$last2}."..".$END2{$last2}.",".$START3{$last2}."..".$END3{$last2}.",".$START2{$last1}."..".$END2{$last1}."))"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}elsif ((!defined $STRAND3{$last1}) and (!defined $STRAND3{$last2})) {
						if ((($length_rps12_1 > 80)and ($length_rps12_1 < 200)) and ($length_rps12_2 > 200)) {
							push @rps12,"     CDS             complement(join(".$START2{$last2}."..".$END2{$last2}.",".$START2{$last1}."..".$END2{$last1}."))"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}
				}elsif (($STRAND2{$last1} eq "-") and ($STRAND2{$last2} eq "+")) {
					if ((!defined $STRAND3{$last1}) and (defined $STRAND3{$last2})) {
						if ((($length_rps12_1 > 80)and ($length_rps12_1 < 200)) and ($length_rps12_2 > 200) and ($length_rps12_5 < 80)) {
							push @rps12,"     CDS             join(complement(".$START2{$last1}."..".$END2{$last1}."),".$START2{$last2}."..".$END2{$last2}.",".$START3{$last2}."..".$END3{$last2}.")"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}elsif ((defined $STRAND3{$last1}) and (defined $STRAND3{$last2})) {
						if (($length_rps12_1 < 80) and ($length_rps12_2 > 200) and ($length_rps12_4 > 200) and ($length_rps12_5 < 80)) {
							push @rps12,"     CDS             complement(join(".$START2{$last1}."..".$END2{$last1}.",".$START3{$last1}."..".$END3{$last1}."))"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
							push @rps12,"     CDS             join(".$START2{$last2}."..".$END2{$last2}.",".$START3{$last2}."..".$END3{$last2}.")"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}elsif ((!defined $STRAND3{$last1}) and (!defined $STRAND3{$last2})) {
						if ((($length_rps12_1 > 80) and ($length_rps12_1 < 200)) and ($length_rps12_2 > 200)) {
							push @rps12,"     CDS             join(complement(".$START2{$last1}."..".$END2{$last1}."),".$START2{$last2}."..".$END2{$last2}.")"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}elsif (($length_rps12_1 > 200) and ($length_rps12_2 > 200)) {
							push @rps12,"     CDS             complement(".$START2{$last1}."..".$END2{$last1}.")"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
							push @rps12,"     CDS             ".$START2{$last2}."..".$END2{$last2}."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}
				}elsif (($STRAND2{$last1} eq "+") and ($STRAND2{$last2} eq "+")) {
					if ((!defined $STRAND3{$last1}) and (defined $STRAND3{$last2})) {
						if ((($length_rps12_1 > 80)and ($length_rps12_1 < 200)) and ($length_rps12_2 > 200) and ($length_rps12_5 < 80)) {
							push @rps12,"     CDS             join(".$START2{$last1}."..".$END2{$last1}.",".$START2{$last2}."..".$END2{$last2}.",".$START3{$last2}."..".$END3{$last2}.")"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}elsif ((!defined $STRAND3{$last1}) and (!defined $STRAND3{$last2})) {
						if ((($length_rps12_1 > 80)and ($length_rps12_1 < 200)) and ($length_rps12_2 > 200)) {
							push @rps12,"     CDS             join(".$START2{$last1}."..".$END2{$last1}.",".$START2{$last2}."..".$END2{$last2}.")"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}
				}elsif (($STRAND2{$last1} eq "+") and ($STRAND2{$last2} eq "-")) {
					if ((!defined $STRAND3{$last1}) and (defined $STRAND3{$last2})) {
						if ((($length_rps12_1 > 80) and ($length_rps12_1 < 200)) and ($length_rps12_2 < 80) and ($length_rps12_5 > 200)) {
							push @rps12,"     CDS             join(".$START2{$last1}."..".$END2{$last1}.",complement(".$START2{$last2}."..".$END2{$last2}.",".$START3{$last2}."..".$END3{$last2}."))"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}elsif ((defined $STRAND3{$last1}) and (defined $STRAND3{$last2})) {
						if (($length_rps12_1 > 200) and ($length_rps12_2 < 80) and ($length_rps12_4 < 80) and ($length_rps12_5 > 200)) {
							push @rps12,"     CDS             join(".$START2{$last1}."..".$END2{$last1}.",".$START3{$last1}."..".$END3{$last1}.")"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
							push @rps12,"     CDS             complement(join(".$START2{$last2}."..".$END2{$last2}.",".$START3{$last2}."..".$END3{$last2}."))"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}elsif ((!defined $STRAND3{$last1}) and (!defined $STRAND3{$last2})) {
						if ((($length_rps12_1 > 80) and ($length_rps12_1 < 200)) and ($length_rps12_2 > 200)) {
							push @rps12,"     CDS             join(".$START2{$last1}."..".$END2{$last1}.",complement(".$START2{$last2}."..".$END2{$last2}."))"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}elsif (($length_rps12_1 > 200) and ($length_rps12_2 > 200)) {
							push @rps12,"     CDS             ".$START2{$last1}."..".$END2{$last1}."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
							push @rps12,"     CDS             complement(".$START2{$last2}."..".$END2{$last2}.")"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}
				}
				last if ($last2 == 2);
			}
		}elsif ($cnt3 == 1) {
			foreach (1..$cnt3) {
				my $last1=$_;
				#print $last1."\n";
				my $length_rps12_1=abs($END2{$last1}-$START2{$last1}+1);
				#print $length_rps12_1."\n";
				#print $STRAND2{$last1}."\n";
				my $length_rps12_4;
				if (defined $STRAND3{$last1}) {
					$length_rps12_4=abs($END3{$last1}-$START3{$last1}+1);
				}
				if ($STRAND2{$last1} eq "-") {
					if (!defined $STRAND3{$last1}) {
						if (($length_rps12_1 > 80) and ($length_rps12_1 < 200)) {
							push @rps12,"     CDS             complement(".$START2{$last1}."..".$END2{$last1}.")"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}elsif ($length_rps12_1 > 200) {
							push @rps12,"     CDS             complement(".$START2{$last1}."..".$END2{$last1}.")"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}elsif (defined $STRAND3{$last1}) {
						if (($length_rps12_1 < 80) and ($length_rps12_4 > 200)) {
							push @rps12,"     CDS             complement(join(".$START2{$last1}."..".$END2{$last1}.",".$START3{$last1}."..".$END3{$last1}."))"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}
				}elsif ($STRAND2{$last1} eq "+") {
					if (!defined $STRAND3{$last1}) {
						if ($length_rps12_1 > 200) {
							push @rps12,"     CDS             ".$START2{$last1}."..".$END2{$last1}."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}elsif (($length_rps12_1 > 80) and ($length_rps12_1 < 200)) {
							push @rps12,"     CDS             ".$START2{$last1}."..".$END2{$last1}."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}elsif (defined $STRAND3{$last1}) {
						if (($length_rps12_1 > 200) and ($length_rps12_4 < 80)) {
							push @rps12,"     CDS             join(".$START2{$last1}."..".$END2{$last1}.",".$START3{$last1}."..".$END3{$last1}.")"."\n";
							push @rps12,"                     /gene=\"rps12\""."\n";
							push @rps12,"                     /codon_start=1"."\n";
							push @rps12,"                     /transl_table=11"."\n";
							push @rps12,"                     /product=\"ribosomal protein S12\""."\n";
						}
					}
				}
				last if ($last1 == 1);
			}
		}


		open (my $input_temp2,"<","$output_directory/$temp.gb");
		my $marks=0;
		my $flag=0;
		my (@all1,@all2);
		while (<$input_temp2>){
			#$_=~ s/\r|\n//g;
			if ($_=~ /ORIGIN/){
				$flag=1;
			}
			if (/ {5}gene {12}/){
				$marks=1;
				push @all1,$_;
			}elsif(/^\s+\/gene="(rps12.*)"/ and $marks == 1){
				$marks=0;
				delete $all1[-1];
			}elsif(/ {5}CDS {13}/){
				$marks=2;
				push @all1,$_;
			}elsif(/^\s+\/gene="(rps12.*)"/ and $marks == 2){
				$marks=3;
				delete $all1[-1];
			}elsif (/^\s+\/codon_start=1/ and $marks == 3) {
				$marks=4;
			}elsif (/^\s+\/transl_table=11/ and $marks == 4) {
				$marks=5;
			}elsif (/^\s+\/product="ribosomal protein S12"/ and $marks == 5) {
				$marks=0;
			}elsif ($flag == 0) {
				push @all1,$_;
			}
			if ($flag==1){
				push @all2,$_;
			}
		}
		close $input_temp2;

		open (my $output_temp,">","$output_directory/$filename.gb");
		print $output_temp @all1;
		print $output_temp @rps12;
		print $output_temp @all2;
		close $output_temp;
		unlink("$output_directory/$temp.gb");
	}
	}




	%hash=();
	%hash_exon=();
	%hash_exon_length=();
	%hash_exon_sequence=();

	my ($length_ref,$length_seq,%gene_number_seqs);
	$length_ref=keys %gene_number_ref;
	foreach my $key1 (keys %gene_number_seq) {
		if ($key1=~ /rps12/) {
			$key1=~ s/rps12\+1/rps12/g;
			$key1=~ s/rps12\+2/rps12/g;
			$gene_number_seqs{$key1}++;
		}else{
			$gene_number_seqs{$key1}=$gene_number_seq{$key1};
		}
	}
	$length_seq=keys %gene_number_seqs;
	print $logfile "Total number of genes in the reference plastome(s): $length_ref.\n";
	print $logfile "Total number of genes annotated in the target plastome: $length_seq.\n";
	print $logfile "All gene names from the reference plastome(s) that were not annotated in the target plastome:\n";
	foreach my $key2 (keys %gene_number_ref) {
		if (!exists $gene_number_seqs{$key2}) {
			print $logfile "$key2\t";
		}
	}
	print $logfile "\n\n";
	close $logfile;
	my $now6=&gettime;
	print "$now6 || Finish annotating the $j sequence: $filename";
	print "\n";
}

#unlink("all2.fasta");
unlink("$reference1_random");
unlink("$reference2_random");
unlink("$reference3_random");
unlink("$reference4_random");
%hash_codon=();
%hash_product=();

my $now7=&gettime;
print "$now7 || Finish annotating all sequences and total elapsed time is: ",time-$now1," seconds!";
print "\n";


########################################
##subroutines
########################################
sub getdate {
	my ($sec,$min,$hour,$day,$mon,$year,$weekday,$yeardate,$savinglightday)=(localtime(time));
	my %hash_month=("01"=>"JAN","02"=>"FEB","03"=>"MAR","04"=>"APR","05"=>"MAY","06"=>"JUN","07"=>"JUL","08"=>"AUG","09"=>"SEP","10"=>"OCT","11"=>"NOV","12"=>"DEC");
	$year+=1900;
	$mon=($mon<9)?"0".($mon+1):($mon+1);
	$day=($day<10)?"0$day":$day;
	$hour=($hour<10)?"0$hour":$hour;
	$min=($min<10)?"0$min":$min;
	$sec=($sec<10)?"0$sec":$sec;

	my $now="$day-$hash_month{$mon}-$year";
}

sub gettime {
	my ($sec,$min,$hour,$day,$mon,$year,$weekday,$yeardate,$savinglightday)=(localtime(time));
	$year+=1900;
	$mon=($mon<9)?"0".($mon+1):($mon+1);
	$day=($day<10)?"0$day":$day;
	$hour=($hour<10)?"0$hour":$hour;
	$min=($min<10)?"0$min":$min;
	$sec=($sec<10)?"0$sec":$sec;

	my $now="$year.$mon.$day $hour:$min:$sec";
}

sub argument{
	my @options=("help|h","reference|r:s","target|t:s","ir|i:i","pidentity|p:i","link|l:s","redundancy|d:s","qcoverage|q:s","out|o:s","form|f:s","warning|w:s");
	my %options;
	GetOptions(\%options,@options);
	exec ("pod2usage $0") if ((keys %options)==0 or $options{'h'} or $options{'help'});
	if(!exists $options{'reference'}){
		print "***ERROR: No reference directory is assigned!!!\n";
		exec ("pod2usage $0");
	}
	if(!exists $options{'target'}){
		print "***ERROR: No target directory is assigned!!!\n";
		exec ("pod2usage $0");
	}
	if (-f $options{'reference'}) {
		print "***ERROR: Input directory name for reference, not file name!!!\n";
		exec ("pod2usage $0");
	}
	if (-f $options{'target'}) {
		print "***ERROR: Input directory name for target, not file name!!!\n";
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

    PGA_v2.pl Plastid Genome Annotator version 2.0

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

    Plastid Genome Annotator version 2.0

=head1 SYNOPSIS

    PGA_v2.pl -r -t [-i -p -l -d -q -o -f -w]
    Plastid Genome Annotator version 2.0
    Copyright (C) 2023 Xiao-Jian Qu
    Please contact <quxiaojian@sdnu.edu.cn>, if you have any bugs or questions.

    [-h -help]         help information.
    [-r -reference]    required: (default: reference) input directory name containing GenBank-
                       formatted file(s) that from the same or close families.
    [-t -target]       required: (default: target) input directory name containing FASTA-
                       formatted file(s) that will be annotated.
    [-i -ir]           optional: (default: 1000) minimum allowed inverted-repeat (IR) length.
    [-p -pidentity]    optional: (default: 40) any PCGs with a TBLASTN percent identity less
                       than this value will be listed in the log file and will not be annotated.
    [-l -link]         optional: (default: Y) if Y, the exons of the trans-splicing gene rps12
                       will be linked; if N, the exons of the rps12 gene will not be linked.
    [-d -redundancy]   optional: (default: N) if N, any PCGs with a query coverage per annotated
                       PCG less or greater than each of the two values (<1,>1) in the following
                       parameter will not be annotated; if Y, any PCGs with a query coverage per
                       annotated PCG less or greater than each of the two values (<1,>1) in the
                       following parameter will be annotated, which can be used to trace the
                       evolutionary remains (pseudo-genes) of functional PCGs.
    [-q -qcoverage]    optional: (default: 0.5,2.0) any PCGs with a query coverage per annotated
                       PCG less or greater than each of these two values (<1,>1) will be listed
                       in the log file.
    [-o -out]          optional: (default: gb) output directory name.
    [-f -form]         optional: (default: circular) circular or linear form for FASTA-formatted
                       file.
    [-w -warning]      optional: (default: warning) log file name containing warning information
                       for annotated GenBank-formatted file(s).

=cut
