#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Find;
use Data::Dumper;
$|=1;

my $global_options=&argument();
my $input_directory=&default("input","input");
my $inverted_repeat=&default("1000","minir");
my $similarity=&default("90","similarity");
my $min_length=&default("15","length");
my $evalue=&default("10","evalue");
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


####extract genes and their positions
####extract_outergenes_innergenes_label_position
my $pattern1=".gb";
my $pattern2=".gbk";
my $pattern3=".gbf";
my @filenames1;
find(\&target1,$input_directory);
sub target1{
	if (/$pattern1/ or /$pattern2/ or /$pattern3/){
		push @filenames1,"$File::Find::name";
	}
	return;
}
my @filenames=@filenames1;

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
	my (@row_array,$species_name,$length,$element,@genearray,@output1);
	my $mark=0;
	#my $tick=0;
	open (my $in_gb3,"<",$target_name_temp_random2);
	while (<$in_gb3>){
		$_=~ s/\r|\n//g;
		@row_array=split /\s+/,$_;
		if (/^LOCUS/i){
			#$species_name=$row_array[1];
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
				#$tick=2;
			}elsif($row_array[2]=~ /^(\d+..>\d+)$/){# positive split no-intron gene
				$row_array[2]="\+\t$1\t$row_array[1]";
				$row_array[2]=~ s/\..>/\t/g;
				#$tick=1;
			}elsif($row_array[2]=~/^complement\(<(\d+..\d+)\)$/){# negative split no-intron gene
				$row_array[2]="-\t$1\t$row_array[1]";
				$row_array[2]=~ s/\../\t/g;
				#$tick=1;
			}elsif($row_array[2]=~/^complement\((\d+..>\d+)\)$/){# negative split no-intron gene
				$row_array[2]="-\t$1\t$row_array[1]";
				$row_array[2]=~ s/\..>/\t/g;
				#$tick=2;
			}
			$element=$row_array[2];
			$mark=1;
		}elsif(/^\s+\/gene="(.*)"/ and $mark == 1){
			$element=$species_name.":".$1.":".$element;
			push @genearray,$element;
			$element=();
			$mark=0;
		}
	}
	close $in_gb3;

	foreach (@genearray){
		my @array=split /:/,$_;
		push @output1,"$array[0]\t$array[1]\t$array[2]\n";
	}
	@row_array=();
	@genearray=();


	#generate data_species_length.txt file with position of species and length
	#Rosa_roxburghii 16000 21000 156,749bp
	#Rosa_roxburghii 23000 28000 Rosa_roxburghii
	#16000=plastome length/9.8
	#21000=plastome length/7.5
	#23000=plastome length/6.8
	#28000=plastome length/5.6
	open(my $out_species_length,">","$output_directory/$filename\_species_length.txt");
	my $position_species_length1=int($length/9.8);
	my $position_species_length2=int($length/7.5);
	my $position_species_length3=int($length/6.8);
	my $position_species_length4=int($length/5.6);

	my $length_temp=$length;
	my $lengthspecies=substr($length_temp,0,index($length_temp,substr($length_temp,-3,3))).",".substr($length_temp,-3,3);

	#my $length_temp1=$length;
	#my $length_temp2=$length;
	#substr($length_temp1,-3,3,",");
	#my $lengthspecies=$length_temp1.substr($length_temp2,-3,3);

	print $out_species_length $filename."\t".$position_species_length1."\t".$position_species_length2."\t".$lengthspecies."bp"."\n";
	print $out_species_length $filename."\t".$position_species_length3."\t".$position_species_length4."\t".$filename."\n";
	close $out_species_length;

	#karyotype.txt
	#chr - Rosa_roxburghii Rosa_roxburghii 0 156749 black
	open(my $out_karyotype,">","$output_directory/$filename\_karyotype.txt");
	print $out_karyotype "chr - ".$filename." ".$filename." "."0"." ".$length." "."black"."\n";
	close $out_karyotype;


	#put_fasta_sequence_in_array and generate fasta sequence file from gb file
    my $flag=0;
    my @sequence;
	my (@fas1,@fas2);
	open(my $in_gb4,"<",$target_name_temp_random2);
	open(my $out_gb4,">","$input_directory/$filename.fasta");
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
		$_=~ s/\s*//g;
		$_=~ s/\d+//g;
		push @fas1,$_;
	}
    my $fas1=join "",@fas1;
	print $out_gb4 ">".$filename."\n".$fas1."\n";
	close $out_gb4;

	my (@fasta1,@fasta2);
    push @fasta1,$filename,$fas1;
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
			push @output2,($SP1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_});
		}elsif (defined $STRAND2{$_} ne "") {
			if (($GENE1{$_} ne "rps12")) {
				if (($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START1{$_} < $START2{$_})){
					push @output2,($SP1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_});
					push @output2,($SP1{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_});
				}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START1{$_} > $START2{$_})){
					push @output2,($SP1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_});
					push @output2,($SP1{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_});
				}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START1{$_} < $START2{$_})){
					push @output2,($SP1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_});
					push @output2,($SP1{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_});
				}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START1{$_} > $START2{$_})){
					push @output2,($SP1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_});
					push @output2,($SP1{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_});
				}
			}elsif (($GENE1{$_} eq "rps12")) {
				if (($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START1{$_} < $START2{$_})){
					push @output2,($SP1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_});
					push @output2,($SP1{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_});
				}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START1{$_} > $START2{$_})){
					push @output2,($SP1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_});
					push @output2,($SP1{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_});
				}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START1{$_} < $START2{$_})){
					push @output2,($SP1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_});
					push @output2,($SP1{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_});
				}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START1{$_} > $START2{$_})){
					push @output2,($SP1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_});
					push @output2,($SP1{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_});

				}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "+") and ($START1{$_} < $START2{$_})){
					push @output2,($SP1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_});
					push @output2,($SP1{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_});
				}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "+") and ($START1{$_} > $START2{$_})){
					push @output2,($SP1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_});
					push @output2,($SP1{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_});
				}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "-") and ($START1{$_} < $START2{$_})){
					push @output2,($SP1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_});
					push @output2,($SP1{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_});
				}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "-") and ($START1{$_} > $START2{$_})){
					push @output2,($SP1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_});
					push @output2,($SP1{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$GENE1{$_}."\t".$STRAND2{$_});
				}
			}
		}
	}
	@output1=();


	#sort_bed_file
	my $col1=1;
	my $col2=4;
	my $col3=3;
	my %sort;
	foreach (@output2){
		my @row=split /\t/,$_;
		$sort{$_}=$row[$col1];
	}

	#output forward genes and reverse genes
	open (my $out_bed1,">","$output_directory/$filename\_outergenes.label.txt");
	open (my $out_bed2,">","$output_directory/$filename\_innergenes.label.txt");
	open (my $out_bed3,">","$output_directory/$filename\_outergenes.position.txt");
	open (my $out_bed4,">","$output_directory/$filename\_innergenes.position.txt");
	foreach (sort {$sort{$a} <=> $sort{$b}} keys %sort){
		my @row=split /\t/,$_;
		if ($row[$col2] eq "+") {
			if ($row[1] < $row[2]) {
				print $out_bed1 $row[0]."\t".$row[1]."\t".$row[2]."\t".$row[3]."\n";
				if ($row[$col3]=~ /accD/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=chr19"."\n";
				}elsif ($row[$col3]=~ /^atp/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=lgreen"."\n";
				}elsif ($row[$col3]=~ /ccsA/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=chr19"."\n";
				}elsif ($row[$col3]=~ /cemA/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=chr19"."\n";
				}elsif ($row[$col3]=~ /^chl/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=chr19"."\n";
				}elsif ($row[$col3]=~ /clpP/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=orange"."\n";
				}elsif ($row[$col3]=~ /infA/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=chr19"."\n";
				}elsif ($row[$col3]=~ /matK/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=orange"."\n";
				}elsif ($row[$col3]=~ /^ndh/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=chr10"."\n";
				}elsif ($row[$col3]=~ /^pet/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=chrun"."\n";
				}elsif ($row[$col3]=~ /^psa/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=chr13"."\n";
				}elsif ($row[$col3]=~ /^psb/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=dgreen"."\n";
				}elsif ($row[$col3]=~ /rbcL/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=green"."\n";
				}elsif ($row[$col3]=~ /^rpl/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=chr1"."\n";
				}elsif ($row[$col3]=~ /^rpo/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=chr4"."\n";
				}elsif ($row[$col3]=~ /^rps/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=dorange"."\n";
				}elsif ($row[$col3]=~ /^rrn/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=chr5"."\n";
				}elsif ($row[$col3]=~ /^trn/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=chr14"."\n";
				}elsif ($row[$col3]=~ /^ycf/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=vlyellow"."\n";
				}
			}elsif ($row[1] > $row[2]) {
				print $out_bed1 $row[0]."\t".$row[1]."\t".$length_temp."\t".$row[3]."\n";
				print $out_bed1 $row[0]."\t"."1"."\t".$row[2]."\t".$row[3]."\n";
				if ($row[$col3]=~ /accD/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=chr19"."\n";
					print $out_bed3 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=chr19"."\n";
				}elsif ($row[$col3]=~ /^atp/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=lgreen"."\n";
					print $out_bed3 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=lgreen"."\n";
				}elsif ($row[$col3]=~ /ccsA/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=chr19"."\n";
					print $out_bed3 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=chr19"."\n";
				}elsif ($row[$col3]=~ /cemA/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=chr19"."\n";
					print $out_bed3 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=chr19"."\n";
				}elsif ($row[$col3]=~ /^chl/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=chr19"."\n";
					print $out_bed3 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=chr19"."\n";
				}elsif ($row[$col3]=~ /clpP/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=orange"."\n";
					print $out_bed3 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=orange"."\n";
				}elsif ($row[$col3]=~ /infA/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=chr19"."\n";
					print $out_bed3 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=chr19"."\n";
				}elsif ($row[$col3]=~ /matK/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=orange"."\n";
					print $out_bed3 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=orange"."\n";
				}elsif ($row[$col3]=~ /^ndh/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=chr10"."\n";
					print $out_bed3 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=chr10"."\n";
				}elsif ($row[$col3]=~ /^pet/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=chrun"."\n";
					print $out_bed3 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=chrun"."\n";
				}elsif ($row[$col3]=~ /^psa/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=chr13"."\n";
					print $out_bed3 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=chr13"."\n";
				}elsif ($row[$col3]=~ /^psb/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=dgreen"."\n";
					print $out_bed3 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=dgreen"."\n";
				}elsif ($row[$col3]=~ /rbcL/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=green"."\n";
					print $out_bed3 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=green"."\n";
				}elsif ($row[$col3]=~ /^rpl/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=chr1"."\n";
					print $out_bed3 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=chr1"."\n";
				}elsif ($row[$col3]=~ /^rpo/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=chr4"."\n";
					print $out_bed3 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=chr4"."\n";
				}elsif ($row[$col3]=~ /^rps/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=dorange"."\n";
					print $out_bed3 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=dorange"."\n";
				}elsif ($row[$col3]=~ /^rrn/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=chr5"."\n";
					print $out_bed3 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=chr5"."\n";
				}elsif ($row[$col3]=~ /^trn/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=chr14"."\n";
					print $out_bed3 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=chr14"."\n";
				}elsif ($row[$col3]=~ /^ycf/) {
					print $out_bed3 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=vlyellow"."\n";
					print $out_bed3 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=vlyellow"."\n";
				}
			}
		}elsif ($row[$col2] eq "-") {
			if ($row[1] < $row[2]) {
				print $out_bed2 $row[0]."\t".$row[1]."\t".$row[2]."\t".$row[3]."\n";
				if ($row[$col3]=~ /accD/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=chr19"."\n";
				}elsif ($row[$col3]=~ /^atp/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=lgreen"."\n";
				}elsif ($row[$col3]=~ /ccsA/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=chr19"."\n";
				}elsif ($row[$col3]=~ /cemA/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=chr19"."\n";
				}elsif ($row[$col3]=~ /^chl/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=chr19"."\n";
				}elsif ($row[$col3]=~ /clpP/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=orange"."\n";
				}elsif ($row[$col3]=~ /infA/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=chr19"."\n";
				}elsif ($row[$col3]=~ /matK/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=orange"."\n";
				}elsif ($row[$col3]=~ /^ndh/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=chr10"."\n";
				}elsif ($row[$col3]=~ /^pet/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=chrun"."\n";
				}elsif ($row[$col3]=~ /^psa/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=chr13"."\n";
				}elsif ($row[$col3]=~ /^psb/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=dgreen"."\n";
				}elsif ($row[$col3]=~ /rbcL/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=green"."\n";
				}elsif ($row[$col3]=~ /^rpl/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=chr1"."\n";
				}elsif ($row[$col3]=~ /^rpo/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=chr4"."\n";
				}elsif ($row[$col3]=~ /^rps/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=dorange"."\n";
				}elsif ($row[$col3]=~ /^rrn/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=chr5"."\n";
				}elsif ($row[$col3]=~ /^trn/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=chr14"."\n";
				}elsif ($row[$col3]=~ /^ycf/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$row[2]."\t"."fill_color=vlyellow"."\n";
				}
			}elsif ($row[1] > $row[2]) {
				print $out_bed2 $row[0]."\t".$row[1]."\t".$length_temp."\t".$row[3]."\n";
				print $out_bed2 $row[0]."\t"."1"."\t".$row[2]."\t".$row[3]."\n";
				if ($row[$col3]=~ /accD/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=chr19"."\n";
					print $out_bed4 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=chr19"."\n";
				}elsif ($row[$col3]=~ /^atp/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=lgreen"."\n";
					print $out_bed4 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=lgreen"."\n";
				}elsif ($row[$col3]=~ /ccsA/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=chr19"."\n";
					print $out_bed4 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=chr19"."\n";
				}elsif ($row[$col3]=~ /cemA/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=chr19"."\n";
					print $out_bed4 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=chr19"."\n";
				}elsif ($row[$col3]=~ /^chl/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=chr19"."\n";
					print $out_bed4 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=chr19"."\n";
				}elsif ($row[$col3]=~ /clpP/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=orange"."\n";
					print $out_bed4 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=orange"."\n";
				}elsif ($row[$col3]=~ /infA/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=chr19"."\n";
					print $out_bed4 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=chr19"."\n";
				}elsif ($row[$col3]=~ /matK/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=orange"."\n";
					print $out_bed4 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=orange"."\n";
				}elsif ($row[$col3]=~ /^ndh/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=chr10"."\n";
					print $out_bed4 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=chr10"."\n";
				}elsif ($row[$col3]=~ /^pet/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=chrun"."\n";
					print $out_bed4 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=chrun"."\n";
				}elsif ($row[$col3]=~ /^psa/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=chr13"."\n";
					print $out_bed4 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=chr13"."\n";
				}elsif ($row[$col3]=~ /^psb/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=dgreen"."\n";
					print $out_bed4 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=dgreen"."\n";
				}elsif ($row[$col3]=~ /rbcL/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=green"."\n";
					print $out_bed4 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=green"."\n";
				}elsif ($row[$col3]=~ /^rpl/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=chr1"."\n";
					print $out_bed4 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=chr1"."\n";
				}elsif ($row[$col3]=~ /^rpo/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=chr4"."\n";
					print $out_bed4 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=chr4"."\n";
				}elsif ($row[$col3]=~ /^rps/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=dorange"."\n";
					print $out_bed4 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=dorange"."\n";
				}elsif ($row[$col3]=~ /^rrn/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=chr5"."\n";
					print $out_bed4 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=chr5"."\n";
				}elsif ($row[$col3]=~ /^trn/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=chr14"."\n";
					print $out_bed4 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=chr14"."\n";
				}elsif ($row[$col3]=~ /^ycf/) {
					print $out_bed4 $row[0]."\t".$row[1]."\t".$length_temp."\t"."fill_color=vlyellow"."\n";
					print $out_bed4 $row[0]."\t"."1"."\t".$row[2]."\t"."fill_color=vlyellow"."\n";
				}
			}
		}
	}
	close $out_bed1;
	close $out_bed2;
	close $out_bed3;
	close $out_bed4;
	@output2=();




	####generate ticks.conf file
	open(my $out_ticks,">","$output_directory/$filename\_ticks.conf");
	print $out_ticks "show_ticks          = yes"."\n";
	print $out_ticks "show_tick_labels    = yes"."\n";
	print $out_ticks "show_grid           = yes"."\n\n";

	print $out_ticks "<ticks>"."\n";
	print $out_ticks "skip_first_label     = no"."\n";
	print $out_ticks "skip_last_label      = no"."\n";
	print $out_ticks "radius               = dims(ideogram,radius_outer)"."\n";
	print $out_ticks "tick_separation      = 2p"."\n";
	print $out_ticks "min_label_distance_to_edge = 10p"."\n";
	print $out_ticks "label_separation = 5p"."\n";
	print $out_ticks "label_offset     = 5p"."\n";
	print $out_ticks "multiplier       = 1e-3"."\n";
	print $out_ticks "color            = black"."\n\n";

	print $out_ticks "<tick>"."\n";
	print $out_ticks "spacing        = 1u"."\n";
	print $out_ticks "size           = 15p"."\n";
	print $out_ticks "thickness      = 2p"."\n";
	print $out_ticks "orientation    = in"."\n";
	print $out_ticks "</tick>"."\n\n";

	print $out_ticks "<tick>"."\n";
	print $out_ticks "spacing        = 5u"."\n";
	print $out_ticks "size           = 20p"."\n";
	print $out_ticks "thickness      = 4p"."\n";
	print $out_ticks "show_label     = yes"."\n";
	print $out_ticks "label_size     = 17p"."\n";
	print $out_ticks "label_offset   = 60p"."\n";
	print $out_ticks "format         = %d"."\n\n";

	print $out_ticks "grid           = yes"."\n";
	print $out_ticks "grid_start     = 0.78r"."\n";
	print $out_ticks "grid_end       = 0.98r"."\n";
	print $out_ticks "grid_color     = grey"."\n";
	print $out_ticks "grid_thickness = 2p"."\n";
	print $out_ticks "suffix         = kb"."\n";
	print $out_ticks "orientation    = in"."\n";
	print $out_ticks "</tick>"."\n\n";

	print $out_ticks "</ticks>"."\n";
	close $out_ticks;




	####generate ideogram.conf file
	open(my $out_ideogram,">","$output_directory/$filename\_ideogram.conf");
	print $out_ideogram "<ideogram>"."\n\n";

	print $out_ideogram "<spacing>"."\n";
	print $out_ideogram "default = 0r"."\n";
	print $out_ideogram "break   = 0"."\n";
	print $out_ideogram "</spacing>"."\n\n\n";

	print $out_ideogram "radius           = 0.8r"."\n";
	print $out_ideogram "thickness        = 2p"."\n";
	print $out_ideogram "fill             = yes"."\n\n";

	print $out_ideogram "fill_color       = gpos100"."\n";
	print $out_ideogram "stroke_thickness = 3"."\n";
	print $out_ideogram "stroke_color     = black"."\n\n";

	print $out_ideogram "</ideogram>"."\n";
	close $out_ideogram;




	####generate visualization.conf file
	open(my $out_circos,">","$output_directory/$filename\_visualization.conf");
	print $out_circos "<<include etc/colors_fonts_patterns.conf>>"."\n";
	print $out_circos "<<include etc/housekeeping.conf>>"."\n\n";

	print $out_circos "<<include $filename\_ideogram.conf>>"."\n";
	print $out_circos "<<include $filename\_ticks.conf>>"."\n\n";

	print $out_circos "<image>"."\n";
	print $out_circos "<<include etc/image.conf>>"."\n";
	print $out_circos "#angle_orientation* = clockwise"."\n";
	print $out_circos "angle_orientation* = counterclockwise"."\n";
	print $out_circos "#background* = transparent"."\n";
	print $out_circos "#radius* = 1500p"."\n";
	print $out_circos "angle_offset* = 8.6"."\n";
	print $out_circos "</image>"."\n\n";

	print $out_circos "karyotype   = $output_directory/$filename\_karyotype.txt"."\n";
	print $out_circos "chromosomes_units           = 1000"."\n";
	print $out_circos "chromosomes_display_default = yes"."\n\n\n\n\n";

	print $out_circos "<plots>"."\n\n";

	print $out_circos "<plot>"."\n";
	print $out_circos "type  = text"."\n";
	print $out_circos "file  = $output_directory/$filename\_species_length.txt"."\n";
	print $out_circos "r0    = 1.25r"."\n";
	print $out_circos "r1    = 1.85r"."\n";
	print $out_circos "label_parallel = yes"."\n";
	print $out_circos "color = black"."\n";
	print $out_circos "label_font = default"."\n";
	print $out_circos "label_size = 35"."\n";
	print $out_circos "</plot>"."\n\n\n";

	print $out_circos "<plot>"."\n";
	print $out_circos "type  = text"."\n";
	print $out_circos "file  = $output_directory/$filename\_outergenes.label.txt"."\n\n";

	print $out_circos "r0    = 1.07r"."\n";
	print $out_circos "r1    = 1.22r"."\n\n";

	print $out_circos "label_font = light"."\n";
	print $out_circos "label_size = 20p"."\n\n";

	print $out_circos "rpadding   = 0p"."\n\n";

	print $out_circos "label_snuggle         = yes"."\n";
	print $out_circos "max_snuggle_distance  = 1r"."\n";
	print $out_circos "snuggle_tolerance     = 0.3r"."\n";
	print $out_circos "snuggle_sampling      = 3"."\n";
	print $out_circos "snuggle_refine        = yes"."\n";
	print $out_circos "</plot>"."\n\n\n";

	print $out_circos "<plot>"."\n";
	print $out_circos "type  = text"."\n";
	print $out_circos "file  = $output_directory/$filename\_innergenes.label.txt"."\n\n";

	print $out_circos "r0    = 1.01r"."\n";
	print $out_circos "r1    = 1.16r"."\n\n";

	print $out_circos "label_font = light"."\n";
	print $out_circos "label_size = 20p"."\n\n";

	print $out_circos "rpadding   = 0p"."\n\n";

	print $out_circos "label_snuggle         = yes"."\n";
	print $out_circos "max_snuggle_distance  = 1r"."\n";
	print $out_circos "snuggle_tolerance     = 0.3r"."\n";
	print $out_circos "snuggle_sampling      = 3"."\n";
	print $out_circos "snuggle_refine        = yes"."\n";
	print $out_circos "</plot>"."\n\n\n";

	print $out_circos "<plot>"."\n";
	print $out_circos "type = highlight"."\n";
	print $out_circos "file = $output_directory/$filename\_quadripartite.IRb.IRa.txt"."\n";
	print $out_circos "r0 = 0.996r"."\n";
	print $out_circos "r1 = 1.004r"."\n";
	print $out_circos "</plot>"."\n\n\n";

	print $out_circos "<plot>"."\n";
	print $out_circos "type = histogram"."\n";
	print $out_circos "file = $output_directory/$filename\_GCcontent_background.hist.txt"."\n";
	print $out_circos "r0   = 0.78r"."\n";
	print $out_circos "r1   = 0.88r"."\n";
	print $out_circos "extend_bin = no"."\n";
	print $out_circos "stroke_type = outline"."\n";
	print $out_circos "thickness   = 0.01r"."\n";
	print $out_circos "color      = vlgrey"."\n";
	print $out_circos "fill_color = vlgrey"."\n";
	print $out_circos "</plot>"."\n\n\n";

	print $out_circos "<plot>"."\n";
	print $out_circos "type = histogram"."\n";
	print $out_circos "file = $output_directory/$filename\_GCcontent_foreground.hist.txt"."\n";
	print $out_circos "r0   = 0.78r"."\n";
	print $out_circos "r1   = 0.88r"."\n";
	print $out_circos "extend_bin = no"."\n";
	print $out_circos "stroke_type = outline"."\n";
	print $out_circos "thickness   = 0.01r"."\n";
	print $out_circos "color      = grey"."\n";
	print $out_circos "fill_color = grey"."\n";
	print $out_circos "</plot>"."\n\n\n";

	print $out_circos "<plot>"."\n";
	print $out_circos "type  = text"."\n";
	print $out_circos "file  = $output_directory/$filename\_quadripartite.label.txt"."\n";
	print $out_circos "r0    = 0.740r"."\n";
	print $out_circos "r1    = 0.800r"."\n";
	print $out_circos "label_parallel = yes"."\n";
	print $out_circos "color = black"."\n";
	print $out_circos "label_font = light"."\n";
	print $out_circos "label_size = 29"."\n";
	print $out_circos "</plot>"."\n\n\n";

	print $out_circos "<plot>"."\n";
	print $out_circos "type = histogram"."\n";
	print $out_circos "file = $output_directory/$filename\_quadripartite.junction.txt"."\n";
	print $out_circos "r0   = 0.725r"."\n";
	print $out_circos "r1   = 0.765r"."\n";
	print $out_circos "extend_bin = no"."\n";
	print $out_circos "stroke_type = outline"."\n";
	print $out_circos "thickness   = 5r"."\n";
	print $out_circos "color      = black"."\n";
	print $out_circos "fill_color = black"."\n";
	print $out_circos "</plot>"."\n\n\n";

	print $out_circos "<plot>"."\n";
	print $out_circos "type = histogram"."\n";
	print $out_circos "file = $output_directory/$filename\_repeatsFR.hist.txt"."\n";
	print $out_circos "r0   = 0.67r"."\n";
	print $out_circos "r1   = 0.72r"."\n";
	print $out_circos "extend_bin = no"."\n";
	print $out_circos "stroke_type = outline"."\n";
	print $out_circos "thickness   = 0.5r"."\n";
	print $out_circos "color      = vvdgrey"."\n";
	print $out_circos "fill_color = vdgrey"."\n";
	print $out_circos "</plot>"."\n\n\n";

	print $out_circos "<plot>"."\n";
	print $out_circos "type = histogram"."\n";
	print $out_circos "file = $output_directory/$filename\_repeatsIR.hist.txt"."\n";
	print $out_circos "r0   = 0.67r"."\n";
	print $out_circos "r1   = 0.72r"."\n";
	print $out_circos "extend_bin = no"."\n";
	print $out_circos "stroke_type = outline"."\n";
	print $out_circos "thickness   = 0.5r"."\n";
	print $out_circos "color      = chr4"."\n";
	print $out_circos "fill_color = chr5"."\n";
	print $out_circos "</plot>"."\n\n\n";

	print $out_circos "</plots>"."\n\n\n\n\n";

	print $out_circos "<highlights>"."\n";
	print $out_circos "z = 0"."\n";
	print $out_circos "fill_color = white"."\n\n";

	print $out_circos "<highlight>"."\n";
	print $out_circos "file       = $output_directory/$filename\_outergenes.position.txt"."\n";
	print $out_circos "r0         = 1.00r"."\n";
	print $out_circos "r1         = 1.06r"."\n";
	print $out_circos "stroke_color = black"."\n";
	print $out_circos "stroke_thickness = 1"."\n";
	print $out_circos "</highlight>"."\n\n";

	print $out_circos "<highlight>"."\n";
	print $out_circos "file       = $output_directory/$filename\_innergenes.position.txt"."\n";
	print $out_circos "r0         = 0.94r"."\n";
	print $out_circos "r1         = 1.00r"."\n";
	print $out_circos "stroke_color = black"."\n";
	print $out_circos "stroke_thickness = 1"."\n";
	print $out_circos "</highlight>"."\n\n";

	print $out_circos "<highlight>"."\n";
	print $out_circos "file = $output_directory/$filename\_quadripartite.highlight.txt"."\n";
	print $out_circos "r0 = 0.725r"."\n";
	print $out_circos "r1 = 0.730r"."\n";
	print $out_circos "show_label = yes"."\n";
	print $out_circos "stroke_color = black"."\n";
	print $out_circos "stroke_thickness = 1"."\n";
	print $out_circos "</highlight>"."\n\n";

	print $out_circos "</highlights>"."\n\n\n\n\n";

	print $out_circos "<links>"."\n\n";

	print $out_circos "<link>"."\n";
	print $out_circos "file          = $output_directory/$filename\_repeatsFR.link.txt"."\n";
	print $out_circos "radius        = 0.66r"."\n";
	print $out_circos "bezier_radius = 0r"."\n";
	print $out_circos "color         = black"."\n";
	print $out_circos "thickness     = 2"."\n";
	print $out_circos "</link>"."\n\n";

	print $out_circos "<link>"."\n";
	print $out_circos "file          = $output_directory/$filename\_repeatsIR.link.txt"."\n";
	print $out_circos "radius        = 0.66r"."\n";
	print $out_circos "bezier_radius = 0r"."\n";
	print $out_circos "color         = red"."\n";
	print $out_circos "thickness     = 2"."\n";
	print $out_circos "</link>"."\n\n";

	print $out_circos "</links>"."\n";

	close $out_circos;
}




####count GC content
my $pattern4=".fasta";
my $pattern5=".fas";
my $pattern6=".fa";
my @filenames2;
find(\&target2,$input_directory);
sub target2{
	if (/$pattern4/ or /$pattern5/ or /$pattern6/){
		push @filenames2,"$File::Find::name";
	}
	return;
}
my @filenames3=@filenames2;

while (@filenames2) {
	my $filename_fasta=shift @filenames2;
	my $target_name=substr($filename_fasta,0,rindex($filename_fasta,"\."));
	$target_name=~ s/\s/_/g;
	my $filename=substr($target_name,rindex($target_name,"\/")+1);

	open(my $in_fasta,"<",$filename_fasta);
	open(my $out_fasta1,">","$output_directory/$filename\_GCcontent_background.hist.txt");
	open(my $out_fasta2,">","$output_directory/$filename\_GCcontent_foreground.hist.txt");

	my ($header,$sequence);
	my @array;
	while (defined ($header=<$in_fasta>) && defined ($sequence=<$in_fasta>)){
		$header=~ s/\r|\n//g;
		$sequence=~ s/\r|\n//g;
		my $length=length $sequence;
		for (my $i=0;$i<($length-100+1);$i+=100) {
			my $subsequence=substr($sequence,$i,100);
			my @nt=split //,$subsequence;
			my $gc=0;
			my $j=0;
			foreach my $key (@nt) {
				if (($key eq "G") or ($key eq "C") or ($key eq "g") or ($key eq "c")) {
					$j++;
					$gc=$j/100;
				}
			}
			push @array,$gc;
		}
	}

	print $out_fasta1 $filename."\t"."1"."\t"."2"."\t"."0"."\n";
	print $out_fasta2 $filename."\t"."1"."\t"."2"."\t"."1"."\n";
	my $length=$#array+1;
	my $j=0;
	my $cnt=0;
	for (my $i=1;$i<($length*100);$i+=100) {
		$j=$i+99;
		print $out_fasta1 $filename."\t"."$i"."\t"."$j"."\t"."1"."\n";
		print $out_fasta2 $filename."\t"."$i"."\t"."$j"."\t".$array[$cnt]."\n";
		$cnt++;
	}

	close $in_fasta;
	close $out_fasta1;
	close $out_fasta2;
}




####blast for repeats
####detection_of_longest_IR_and_repeats_for_FR_and_IR
while (@filenames3) {
	my $filename_fasta=shift @filenames3;
	my $target_name=substr($filename_fasta,0,rindex($filename_fasta,"\."));
	$target_name=~ s/\s/_/g;
	my $filename=substr($target_name,rindex($target_name,"\/")+1);
	#my $target_path=$target_name;
	#$target_path=~ s/$filename//g;
	#my $target_name_random = $target_path."/".$random."_".$filename;
	#my $IR_temp_random = $target_path."/".$random."_IR_temp";
	#my $blast_IR_temp_random = $target_path."/".$random."_IR_temp_blast";
	#my $blast_repeats_temp_random = $target_path."/".$random."_repeats_temp_blast";
	my $target_name_random = $input_directory."/".$random."_".$filename;
	my $IR_temp_random = $input_directory."/".$random."_IR_temp";
	my $blast_IR_temp_random = $input_directory."/".$random."_IR_temp_blast";
	my $blast_repeats_temp_random = $input_directory."/".$random."_repeats_temp_blast";

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
		$header=~ s/\`/_/g;
		$header=~ s/\~/_/g;
		$header=~ s/\!/_/g;
		$header=~ s/\@/_/g;
		$header=~ s/\#/_/g;
		$header=~ s/\$/_/g;
		$header=~ s/\%/_/g;
		$header=~ s/\^/_/g;
		$header=~ s/\&/_/g;
		$header=~ s/\*/_/g;
		$header=~ s/\(/_/g;
		$header=~ s/\)/_/g;
		$header=~ s/\-/_/g;
		$header=~ s/\+/_/g;
		$header=~ s/\=/_/g;
		$header=~ s/\{/_/g;
		$header=~ s/\}/_/g;
		$header=~ s/\[/_/g;
		$header=~ s/\]/_/g;
		$header=~ s/\|/_/g;
		$header=~ s/\\/_/g;
		$header=~ s/\:/_/g;
		$header=~ s/\;/_/g;
		$header=~ s/\"/_/g;
		$header=~ s/\'/_/g;
		$header=~ s/\</_/g;
		$header=~ s/\>/_/g;
		$header=~ s/\,/_/g;
		$header=~ s/\?/_/g;
		$header=~ s/\//_/g;
		$sequence=shift @fasta;
		$sequence= uc $sequence;
		$length_cp=length $sequence;
	}
	#my $rev_coms=reverse $sequence;
	#$rev_coms=~ tr/ACGTacgt/TGCAtgca/;


	#IRb and IRa
	my $cmd1="makeblastdb.exe -in $IR_temp_random -hash_index -dbtype nucl -blastdb_version 4";
	my $cmd2="blastn.exe -task blastn -query $IR_temp_random -db $IR_temp_random -outfmt 6 -perc_identity 99 -out $blast_IR_temp_random";
	my $cmd3="makeblastdb -in $IR_temp_random -hash_index -dbtype nucl -blastdb_version 4";
	my $cmd4="blastn -task blastn -query $IR_temp_random -db $IR_temp_random -outfmt 6 -perc_identity 99 -out $blast_IR_temp_random";

	#repeats including FR(forward repeats) and IR(inverted repeats)
	my $cmd5="makeblastdb.exe -in $target_name_random -hash_index -dbtype nucl -blastdb_version 4";
	my $cmd6="blastn.exe -task blastn -query $target_name_random -db $target_name_random -outfmt 6 -perc_identity $similarity -evalue $evalue -out $blast_repeats_temp_random";
	my $cmd7="makeblastdb -in $target_name_random -hash_index -dbtype nucl -blastdb_version 4";
	my $cmd8="blastn -task blastn -query $target_name_random -db $target_name_random -outfmt 6 -perc_identity $similarity -evalue $evalue -out $blast_repeats_temp_random";


	if ($osname eq "MSWin32") {
		system("$cmd1");#IRb and IRa
		system("$cmd2");#IRb and IRa
		system("$cmd5");#repeats FR and IR
		system("$cmd6");#repeats FR and IR
	}elsif ($osname eq "cygwin") {
		system("$cmd3");#IRb and IRa
		system("$cmd4");#IRb and IRa
		system("$cmd7");#repeats FR and IR
		system("$cmd8");#repeats FR and IR
	}elsif ($osname eq "linux") {
		system("$cmd3");#IRb and IRa
		system("$cmd4");#IRb and IRa
		system("$cmd7");#repeats FR and IR
		system("$cmd8");#repeats FR and IR
	}elsif ($osname eq "darwin") {
		system("$cmd3");#IRb and IRa
		system("$cmd4");#IRb and IRa
		system("$cmd7");#repeats FR and IR
		system("$cmd8");#repeats FR and IR
	}

	unlink ("$target_name_random");
	unlink ("$target_name_random.nhd");
	unlink ("$target_name_random.nhi");
	unlink ("$target_name_random.nhr");
	unlink ("$target_name_random.nin");
	unlink ("$target_name_random.nog");
	unlink ("$target_name_random.nsd");
	unlink ("$target_name_random.nsi");
	unlink ("$target_name_random.nsq");
	unlink ("$IR_temp_random");
	unlink ("$IR_temp_random.nhd");
	unlink ("$IR_temp_random.nhi");
	unlink ("$IR_temp_random.nhr");
	unlink ("$IR_temp_random.nin");
	unlink ("$IR_temp_random.nog");
	unlink ("$IR_temp_random.nsd");
	unlink ("$IR_temp_random.nsi");
	unlink ("$IR_temp_random.nsq");


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
	my ($JLb,$JSb,$JSa);
	$JLb=$JLB-1;
	$JSb=$JSB+1;
	$JSa=$JSA-1;

	my ($junction1,$junction2,$junction3,$junction4);
	$junction1=$JLb-1;
	$junction2=$JSB-1;
	$junction3=$JSa-1;
	$junction4=$JLA-1;

	open (my $out_quadripartite1,">","$output_directory/$filename\_quadripartite.label.txt");
	open (my $out_quadripartite2,">","$output_directory/$filename\_quadripartite.highlight.txt");
	open (my $out_quadripartite3,">","$output_directory/$filename\_quadripartite.junction.txt");
	open (my $out_quadripartite4,">","$output_directory/$filename\_quadripartite.IRb.IRa.txt");
	if (defined $JLB) {
		print $out_quadripartite1 $filename."\t"."1"."\t".$JLb."\t"."LSC"."\n";
		print $out_quadripartite1 $filename."\t".$JLB."\t".$JSB."\t"."IRb"."\n";
		print $out_quadripartite1 $filename."\t".$JSb."\t".$JSa."\t"."SSC"."\n";
		print $out_quadripartite1 $filename."\t".$JSA."\t".$JLA."\t"."IRa"."\n";

		print $out_quadripartite2 $filename."\t"."1"."\t".$JLb."\t"."fill_color=black,z=0"."\n";
		print $out_quadripartite2 $filename."\t".$JLB."\t".$JSB."\t"."fill_color=black,z=0"."\n";
		print $out_quadripartite2 $filename."\t".$JSb."\t".$JSa."\t"."fill_color=black,z=0"."\n";
		print $out_quadripartite2 $filename."\t".$JSA."\t".$JLA."\t"."fill_color=black,z=0"."\n";

		print $out_quadripartite3 $filename."\t"."1"."\t"."2"."\t"."0"."\n";
		print $out_quadripartite3 $filename."\t".$junction1."\t".$JLb."\t"."10"."\n";
		print $out_quadripartite3 $filename."\t".$junction2."\t".$JSB."\t"."10"."\n";
		print $out_quadripartite3 $filename."\t".$junction3."\t".$JSa."\t"."10"."\n";
		print $out_quadripartite3 $filename."\t".$junction4."\t".$JLA."\t"."10"."\n";

		print $out_quadripartite4 $filename."\t".$JLB."\t".$JSB."\t"."fill_color=black,z=0"."\n";
		print $out_quadripartite4 $filename."\t".$JSA."\t".$JLA."\t"."fill_color=black,z=0"."\n";

	}
	close $out_quadripartite1;
	close $out_quadripartite2;
	close $out_quadripartite3;
	close $out_quadripartite4;


	##remove duplicates of repeats (FR and IR)
	open(my $input_repeats,"<",$blast_repeats_temp_random);
	open(my $output_repeats1,">","$output_directory/$filename\_repeatsFR.link.txt");
	open(my $output_repeats2,">","$output_directory/$filename\_repeatsFR.hist.txt");
	open(my $output_repeats3,">","$output_directory/$filename\_repeatsIR.link.txt");
	open(my $output_repeats4,">","$output_directory/$filename\_repeatsIR.hist.txt");
	my $length_max=2000;#set maxminum repeat length to 2000
	my (%repeats1,%repeats2,%repeats3);
	my @repeats;
	while (<$input_repeats>) {
		$_=~ s/\r|\n//g;
		my ($qseqid,$sseqid,$pident,$length,$mismatch,$gapopen,$qstart,$qend,$sstart,$send,$evalue,$bitscore)=split(/\t/,$_);
		if (($length >= $min_length) and ($length <= $length_max)) {#length limitation
			if (($qstart.$qend ne $sstart.$send) and ($qstart.$qend ne $send.$sstart)) {
				$repeats1{$qstart."\t".$qend."\t".$sstart."\t".$send}++;#all repeats
				$repeats2{$sstart."\t".$send."\t".$qstart."\t".$qend}++;#one type of duplicate
				$repeats3{$send."\t".$sstart."\t".$qend."\t".$qstart}++;#another type of duplicate
			}
		}
		#$qseqid."\t".$qstart."\t".$qend."\t".$sseqid."\t".$sstart."\t".$send
		#push @blast_hit,$qseqid."\t".$sseqid."\t".$sseq if ($length >= $min_length);
	}
	foreach my $key (keys %repeats1) {
		if (exists $repeats2{$key}) {
			my ($qstart,$qend,$sstart,$send)=split /\t/,$key;
			if ($qstart < $sstart) {#remove one type of duplicate
				push @repeats,$key;
			}
		}elsif (exists $repeats3{$key}) {
			my ($qstart,$qend,$sstart,$send)=split /\t/,$key;
			if ($qstart < $sstart) {#remove another type of duplicate
				push @repeats,$key;
			}
		}else{
			push @repeats,$key;#no duplicate
		}
	}
	print $output_repeats2 $filename."\t"."1"."\t"."2"."\t"."0"."\n";
	print $output_repeats4 $filename."\t"."1"."\t"."2"."\t"."0"."\n";
	foreach my $key (@repeats) {
		my ($qstart,$qend,$sstart,$send)=split /\t/,$key;
		if ($sstart < $send) {#repeats for FR
			print $output_repeats1 $filename."\t".$qstart."\t".$qend."\t".$filename."\t".$sstart."\t".$send."\n";
			print $output_repeats2 $filename."\t".$qstart."\t".$qend."\t"."10"."\n".$filename."\t".$sstart."\t".$send."\t"."10"."\n";
		}elsif ($sstart > $send) {#repeats for IR
			print $output_repeats3 $filename."\t".$qstart."\t".$qend."\t".$filename."\t".$sstart."\t".$send."\n";
			print $output_repeats4 $filename."\t".$qstart."\t".$qend."\t"."10"."\n".$filename."\t".$sstart."\t".$send."\t"."10"."\n";
		}
	}
	close $input_repeats;
	close $output_repeats1;
	close $output_repeats2;
	close $output_repeats3;
	close $output_repeats4;
	unlink ("$blast_repeats_temp_random");
}




####batch run circos
#perl bin\circos -conf test\input\circos.conf -outputdir test\output -outputfile circos -svg -png

while (@filenames) {
	my $filename_gb=shift @filenames;
	my $target_name=substr($filename_gb,0,rindex($filename_gb,"\."));
	$target_name=~ s/\s/_/g;
	my $filename=substr($target_name,rindex($target_name,"\/")+1);

	my $cmd1="perl bin\\circos -conf $output_directory\\$filename\_visualization.conf -outputdir $output_directory -outputfile $filename\_visualization -svg -png";
	#my $cmd2="perl bin\\circos -conf $output_directory\\$filename\_visualization.conf -outputdir $output_directory -outputfile $filename\_visualization -debug_group textplace";
	system("$cmd1");
	#system("$cmd2");
}




####functions
sub argument{
	my @options=("help|h","input|i:s","minir|m:i","similarity|s:i","length|l:i","evalue|e:s","output|o:s");
	my %options;
	GetOptions(\%options,@options);
	exec ("pod2usage $0") if ((keys %options)==0 or $options{'h'} or $options{'help'});
	if(!exists $options{'input'}){
		print "***ERROR: No input directory is assigned!!!\n";
		exec ("pod2usage $0");
	}
	if (-f $options{'input'}) {
		print "***ERROR: Input directory name, not file name!!!\n";
		exec ("pod2usage $0");
	}
	if(!exists $options{'output'}){
		print "***ERROR: No output directory is assigned!!!\n";
		exec ("pod2usage $0");
	}
	if (-f $options{'output'}) {
		print "***ERROR: Output directory name, not file name!!!\n";
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

    visualization_v1.pl Generate circular plastome map with GC content, quadripartite structure and repeats version 1.0

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

    Generate circular plastome map with GC content, quadripartite structure and repeats

=head1 SYNOPSIS

    visualization_v1.pl -i [-m -s -l -e] -o
    Copyright (C) 2025 Xiao-Jian Qu
    Please contact <quxiaojian@sdnu.edu.cn>, if you have any bugs or questions.

    [-h -help]         help information.
    [-i -input]        required: (default: input) input directory name containing GenBank-formatted file(s).
    [-m -minir]        optional: (default: 1000) minimum allowed length for the longest inverted-repeat (IR).
    [-s -similarity]   optional: (default: 90) minimum allowed percent identity for blast search of dispersed repeats.
    [-l -length]       optional: (default: 15) minimum allowed length for dispersed repeats.
    [-e -evalue]       optional: (default: 10) evalue for blast search of dispersed repeats.
    [-o -output]       required: (default: output) output directory name.

=cut
