#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Find;
use Data::Dumper;
$|=1;

my $global_options=&argument();
my $input_reference=&default("reference.gb","reference");
my $input_target=&default("target.gb","target");
my $inverted_repeat=&default("1000","minir");
my $angle=&default("-90","angle");
my $pattern=&default("1to1","pattern");
my $ribbon=&default("no","ribbon");
my $crest=&default("0.50","crest");
my $bezier=&default("0.25","bezier");
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
my $left_filename_path=substr($input_reference,0,rindex($input_reference,"\."));
$left_filename_path=~ s/\s/_/g;
my $left_filename=substr($left_filename_path,rindex($left_filename_path,"\\")+1);
my $left_path=$left_filename_path;
$left_path=~ s/$left_filename//g;

my $right_filename_path=substr($input_target,0,rindex($input_target,"\."));
$right_filename_path=~ s/\s/_/g;
my $right_filename=substr($right_filename_path,rindex($right_filename_path,"\\")+1);
my $right_path=$right_filename_path;
$right_path=~ s/$right_filename//g;


my @filenames1=($input_reference,$input_target);
#my @filenames2=("$left_path/$left_filename.fasta","$right_path/$right_filename.fasta");
my @filenames2;
while (@filenames1) {
	my $filename_gb=shift @filenames1;
	my $target_name=substr($filename_gb,0,rindex($filename_gb,"\."));
	$target_name=~ s/\s/_/g;
	my $filename=substr($target_name,rindex($target_name,"\\")+1);
	my $target_path=$target_name;
	#$target_path=~ s/$filename//g;
	$target_path=substr($target_name,0,rindex($target_name,"\\"));
	my $target_name_temp_random1 = $target_path."/".$random."_".$filename."_temp1";
	my $target_name_temp_random2 = $target_path."/".$random."_".$filename."_temp2";

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




	#generate species_length.txt file with position of species and length
	#Amborella_trichopoda	16600	21691	162,686bp
	#Amborella_trichopoda	23924	29051	Amborella_trichopoda
	#Rosa_roxburghii	15994	20899	156,749bp
	#Rosa_roxburghii	23051	27990	Rosa_roxburghii

	#16000=plastome length/9.8
	#21000=plastome length/7.5
	#23000=plastome length/6.8
	#28000=plastome length/5.6
	open(my $out_species_length,">>","$output_directory/species_length.txt");
	my $position_species_length1=int($length/8.8);
	my $position_species_length2=int($length/5.5);
	my $position_species_length3=int($length/5.8);
	my $position_species_length4=int($length/3.6);

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
	open(my $out_karyotype,">>","$output_directory/karyotype.txt");
	print $out_karyotype "chr - ".$filename." ".$filename." "."0"." ".$length." "."black"."\n";
	close $out_karyotype;




	#put_fasta_sequence_in_array and generate fasta sequence file from gb file
	my $flag=0;
	my @sequence;
	my (@fas1,@fas2);
	open(my $in_gb4,"<",$target_name_temp_random2);
	open(my $out_gb4,">","$target_path/$filename.fasta");
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
	push @filenames2,"$target_path\\$filename.fasta";

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
	my %order;
	foreach (@output2){
		my @row=split /\t/,$_;
		$order{$_}=$row[$col1];
	}

	#output forward genes and reverse genes
	open (my $out_bed1,">","$output_directory/$filename\_outergenes.label.txt");
	open (my $out_bed2,">","$output_directory/$filename\_innergenes.label.txt");
	open (my $out_bed3,">","$output_directory/$filename\_outergenes.position.txt");
	open (my $out_bed4,">","$output_directory/$filename\_innergenes.position.txt");
	foreach (sort {$order{$a} <=> $order{$b}} keys %order){
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
}




####blast for repeats
####detection_of_longest_IR_and_repeats_for_FR_and_IR
while (@filenames2) {
	my $filename_fasta=shift @filenames2;
	my $target_name=substr($filename_fasta,0,rindex($filename_fasta,"\."));
	$target_name=~ s/\s/_/g;
	my $filename=substr($target_name,rindex($target_name,"\\")+1);
	my $target_path=$target_name;
	#$target_path=~ s/$filename//g;
	$target_path=substr($target_name,0,rindex($target_name,"\\"));
	my $target_name_random = $target_path."/".$random."_".$filename;
	my $IR_temp_random = $target_path."/".$random."_IR_temp";
	my $blast_IR_temp_random = $target_path."/".$random."_IR_temp_blast";

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


	if ($osname eq "MSWin32") {
		system("$cmd1");#IRb and IRa
		system("$cmd2");#IRb and IRa
	}elsif ($osname eq "cygwin") {
		system("$cmd3");#IRb and IRa
		system("$cmd4");#IRb and IRa
	}elsif ($osname eq "linux") {
		system("$cmd3");#IRb and IRa
		system("$cmd4");#IRb and IRa
	}elsif ($osname eq "darwin") {
		system("$cmd3");#IRb and IRa
		system("$cmd4");#IRb and IRa
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

	my ($JLB,$JSB,$JLA,$JSA);
	if (%IR != 0) {
		my (@IR_length,@boundary);
		foreach my $key (sort {$b <=> $a} keys %IR) {
			push @IR_length,$key;
			push @boundary,$IR{$key};
		}
		my ($AA,$BB,$CC,$DD)=split /\t/,$boundary[0];
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
	}

	#open (my $out_quadripartite1,">","$output_directory/$filename\_quadripartite.label.txt");
	#open (my $out_quadripartite2,">","$output_directory/$filename\_quadripartite.highlight.txt");
	#open (my $out_quadripartite3,">","$output_directory/$filename\_quadripartite.junction.txt");
	open (my $out_quadripartite4,">","$output_directory/$filename\_quadripartite.IRb.IRa.txt");
	if (defined $JLB) {
		#print $out_quadripartite1 $filename."\t"."1"."\t".$JLb."\t"."LSC"."\n";
		#print $out_quadripartite1 $filename."\t".$JLB."\t".$JSB."\t"."IRb"."\n";
		#print $out_quadripartite1 $filename."\t".$JSb."\t".$JSa."\t"."SSC"."\n";
		#print $out_quadripartite1 $filename."\t".$JSA."\t".$JLA."\t"."IRa"."\n";

		#print $out_quadripartite2 $filename."\t"."1"."\t".$JLb."\t"."fill_color=black,z=0"."\n";
		#print $out_quadripartite2 $filename."\t".$JLB."\t".$JSB."\t"."fill_color=black,z=0"."\n";
		#print $out_quadripartite2 $filename."\t".$JSb."\t".$JSa."\t"."fill_color=black,z=0"."\n";
		#print $out_quadripartite2 $filename."\t".$JSA."\t".$JLA."\t"."fill_color=black,z=0"."\n";

		#print $out_quadripartite3 $filename."\t"."1"."\t"."2"."\t"."0"."\n";
		#print $out_quadripartite3 $filename."\t".$junction1."\t".$JLb."\t"."10"."\n";
		#print $out_quadripartite3 $filename."\t".$junction2."\t".$JSB."\t"."10"."\n";
		#print $out_quadripartite3 $filename."\t".$junction3."\t".$JSa."\t"."10"."\n";
		#print $out_quadripartite3 $filename."\t".$junction4."\t".$JLA."\t"."10"."\n";

		print $out_quadripartite4 $filename."\t".$JLB."\t".$JSB."\t"."fill_color=black,z=0"."\n";
		print $out_quadripartite4 $filename."\t".$JSA."\t".$JLA."\t"."fill_color=black,z=0"."\n";
	}
	#close $out_quadripartite1;
	#close $out_quadripartite2;
	#close $out_quadripartite3;
	close $out_quadripartite4;
}




####generate genes.link.txt, genes.unique.reference.txt, and genes.unique.target.txt files
open(my $in_innergenes_reference,"<","$output_directory/$left_filename\_innergenes.label.txt");
open(my $in_outergenes_reference,"<","$output_directory/$left_filename\_outergenes.label.txt");
open(my $in_innergenes_target,"<","$output_directory/$right_filename\_innergenes.label.txt");
open(my $in_outergenes_target,"<","$output_directory/$right_filename\_outergenes.label.txt");

my $column1=1;
my $column2=3;
my (%hash_reference1,%hash_target1,%hash_reference2,%hash_target2);
while (<$in_innergenes_reference>) {
	$_=~ s/\r|\n//g;
	my @row=split /\t/,$_;
	$hash_reference1{$_}=$row[$column1];
	$hash_reference2{$row[$column2]}++;
}
while (<$in_outergenes_reference>) {
	$_=~ s/\r|\n//g;
	my @row=split /\t/,$_;
	$hash_reference1{$_}=$row[$column1];
	$hash_reference2{$row[$column2]}++;
}
while (<$in_innergenes_target>) {
	$_=~ s/\r|\n//g;
	my @row=split /\t/,$_;
	$hash_target1{$_}=$row[$column1];
	$hash_target2{$row[$column2]}++;
}
while (<$in_outergenes_target>) {
	$_=~ s/\r|\n//g;
	my @row=split /\t/,$_;
	$hash_target1{$_}=$row[$column1];
	$hash_target2{$row[$column2]}++;
}

close $in_innergenes_reference;
close $in_outergenes_reference;
close $in_innergenes_target;
close $in_outergenes_target;




open(my $in_quadripartite_reference,"<","$output_directory/$left_filename\_quadripartite.IRb.IRa.txt") if (-s "$output_directory/$left_filename\_quadripartite.IRb.IRa.txt");
open(my $in_quadripartite_target,"<","$output_directory/$right_filename\_quadripartite.IRb.IRa.txt") if (-s "$output_directory/$right_filename\_quadripartite.IRb.IRa.txt");

my ($junctionLB_reference,$junctionSB_reference,$junctionSA_reference,$junctionLA_reference);
my ($junctionLB_target,$junctionSB_target,$junctionSA_target,$junctionLA_target);
if (-s "$output_directory/$left_filename\_quadripartite.IRb.IRa.txt") {
	my ($first_row,$second_row);
	while (defined($first_row=<$in_quadripartite_reference>) and defined($second_row=<$in_quadripartite_reference>)) {
		$first_row=~ s/\r|\n//g;
		$second_row=~ s/\r|\n//g;
		my @row1=split /\t/,$first_row;
		my @row2=split /\t/,$second_row;
		$junctionLB_reference=$row1[1];
		$junctionSB_reference=$row1[2];
		$junctionSA_reference=$row2[1];
		$junctionLA_reference=$row2[2];
	}
}
if (-s "$output_directory/$right_filename\_quadripartite.IRb.IRa.txt") {
	my ($first_row,$second_row);
	while (defined($first_row=<$in_quadripartite_target>) and defined($second_row=<$in_quadripartite_target>)) {
		$first_row=~ s/\r|\n//g;
		$second_row=~ s/\r|\n//g;
		my @row1=split /\t/,$first_row;
		my @row2=split /\t/,$second_row;
		$junctionLB_target=$row1[1];
		$junctionSB_target=$row1[2];
		$junctionSA_target=$row2[1];
		$junctionLA_target=$row2[2];
	}
}

close $in_quadripartite_reference if (-s "$output_directory/$left_filename\_quadripartite.IRb.IRa.txt");
close $in_quadripartite_target if (-s "$output_directory/$right_filename\_quadripartite.IRb.IRa.txt");




#my $pattern="1tomore";
#my $pattern="1to1";
if ($pattern eq "1tomore") {
	open(my $out_genes_link,">","$output_directory/genes.link.txt");
	foreach my $key1 (sort {$hash_reference1{$a} <=> $hash_reference1{$b}} keys %hash_reference1){
		my @row1=split /\t/,$key1;
		foreach my $key2 (sort {$hash_target1{$a} <=> $hash_target1{$b}} keys %hash_target1) {
			my @row2=split /\t/,$key2;
			if (($row1[3] ne "rps12") and ($row2[3] ne "rps12")) {#non-rps12 genes
				if ($row2[3] eq $row1[3]) {
					print $out_genes_link $row1[0]."\t".$row1[1]."\t".$row1[2]."\t".$row2[0]."\t".$row2[1]."\t".$row2[2]."\n";
				}
			}elsif (($row1[3] eq "rps12") and ($row2[3] eq "rps12")) {#rps12 gene
				my $length_reference=abs($row1[2]-$row1[1]);
				my $length_target=abs($row2[2]-$row2[1]);
				if (($length_reference < 200) and ($length_target < 200)) {
					print $out_genes_link $row1[0]."\t".$row1[1]."\t".$row1[2]."\t".$row2[0]."\t".$row2[1]."\t".$row2[2]."\n";
				}elsif (($length_reference > 200) and ($length_target > 200)) {
					print $out_genes_link $row1[0]."\t".$row1[1]."\t".$row1[2]."\t".$row2[0]."\t".$row2[1]."\t".$row2[2]."\n";
				}
			}
		}
	}
	close $out_genes_link;

	open(my $out_genes_unique_reference,">","$output_directory/genes.unique.reference.txt");
	open(my $out_genes_unique_target,">","$output_directory/genes.unique.target.txt");
	foreach my $key1 (sort {$hash_reference1{$a} <=> $hash_reference1{$b}} keys %hash_reference1) {
		my @row1=split /\t/,$key1;
		if (!exists $hash_target2{$row1[3]}) {
			print $out_genes_unique_reference $row1[0]."\t".$row1[1]."\t".$row1[2]."\t"."10"."\n";
		}
	}
	foreach my $key2 (sort {$hash_target1{$a} <=> $hash_target1{$b}} keys %hash_target1) {
		my @row2=split /\t/,$key2;
		if (!exists $hash_reference2{$row2[3]}) {
			print $out_genes_unique_target $row2[0]."\t".$row2[1]."\t".$row2[2]."\t"."10"."\n";
		}
	}
	close $out_genes_unique_reference;
	close $out_genes_unique_target;

}elsif ($pattern eq "1to1") {
	open(my $out_genes_link,">","$output_directory/genes.link.txt");
	my $cnt1=-1;
	my (%SP1,%START1,%END1,%GENE1);
	foreach my $key (sort {$hash_reference1{$a} <=> $hash_reference1{$b}} keys %hash_reference1){
		$key=~ s/\r|\n//g;
		$cnt1++;
		($SP1{$cnt1},$START1{$cnt1},$END1{$cnt1},$GENE1{$cnt1})=(split /\t/,$key)[0,1,2,3];
	}

	my $cnt2=-1;
	my (%SP2,%START2,%END2,%GENE2);
	foreach my $key (sort {$hash_target1{$a} <=> $hash_target1{$b}} keys %hash_target1){
		$key=~ s/\r|\n//g;
		$cnt2++;
		($SP2{$cnt2},$START2{$cnt2},$END2{$cnt2},$GENE2{$cnt2})=(split /\t/,$key)[0,1,2,3];
	}

	my ($last1,$last2,$last3,$last4);
	if (($cnt1 != -1) and ($cnt2 != -1)) {
		foreach my $key1 (0..$cnt1) {#0-14
			my $last1_temp=$key1-1;
			if ($last1_temp == -1) {
				$last1=$cnt1;
			}elsif ($last1_temp != -1) {
				$last1=$key1-1;
			}
			my $last2_temp=$key1+1;
			if ($last2_temp == ($cnt1+1)) {
				$last2=0;
			}elsif ($last2_temp != ($cnt1+1)) {
				$last2=$key1+1;
			}

			my $cnt=0;
			foreach my $key2 (0..$cnt2) {#0-14
				my $last3_temp=$key2-1;
				if ($last3_temp == -1) {
					$last3=$cnt2;
				}elsif ($last1_temp != -1) {
					$last3=$key2-1;
				}
				my $last4_temp=$key2+1;
				if ($last4_temp == ($cnt2+1)) {
					$last4=0;
				}elsif ($last4_temp != ($cnt2+1)) {
					$last4=$key2+1;
				}

				if (($GENE1{$key1} ne "rps12") and ($GENE2{$key2} ne "rps12")) {#non-rps12 genes
					##if matched, then change the patterns of 1 to 1, or 1 to 2, or 2 to 1, or 2 to 2
					##to the patterns of 1 to 1, or 1 to 1, or 1 to 1, or 1 to 1, respectively.
					if (($GENE2{$key2} eq $GENE1{$key1}) and ($hash_reference2{$GENE1{$key1}} == 1) and ($hash_target2{$GENE2{$key2}} == 1)) {
						#occurrence 1 time in reference, occurrence 1 time in target
						#gene1 gene2 gene3 gene10 gene11 gene12 gene13 gene14
						#print $GENE1{$key1}."\t".$GENE2{$key2}."\n";
						print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
					}elsif (($GENE2{$key2} eq $GENE1{$key1}) and ($hash_reference2{$GENE1{$key1}} == 1) and ($hash_target2{$GENE2{$key2}} == 2)) {
						#occurrence 1 time in reference, occurrence 2 times in target
						#change gene4 gene4 gene5 gene5 to gene4 gene5
						#print $GENE1{$key1}."\t".$GENE2{$key2}."\n";
						my $i=$cnt;
						if ($i == 0) {
							if (($GENE2{$last3} eq $GENE1{$last1}) and ($GENE2{$last4} eq $GENE1{$last2})) {#positive order
								$i++;
								print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
							}elsif (($GENE2{$last3} eq $GENE1{$last1}) and ($GENE2{$last4} ne $GENE1{$last2})) {#positive order
								$i++;
								print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
							}elsif (($GENE2{$last3} ne $GENE1{$last1}) and ($GENE2{$last4} eq $GENE1{$last2})) {#positive order
								$i++;
								print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
							}
							if (($GENE2{$last3} ne $GENE1{$last1}) and ($GENE2{$last4} ne $GENE1{$last2})) {#positive order
								if (defined $junctionLB_target) {
									if ((($END2{$key2} >= $junctionLB_target) and ($START2{$key2} <= $junctionSB_target))) {#genes in reference vs genes in IRb of target
										#print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
									}elsif ((($END2{$key2} >= $junctionSA_target) and ($START2{$key2} <= $junctionLA_target))) {#genes in reference vs genes in IRa of target
										#print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
									}
								}
							}
						}elsif (($i != 1) or ($i != 2)) {
							#print $GENE1{$key1}."\t".$GENE2{$key2}."\n";
							#if (($GENE2{$last3} eq $GENE1{$last2}) and ($GENE2{$last4} eq $GENE1{$last1})) {#negative order
							#	print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
							#}elsif (($GENE2{$last3} eq $GENE1{$last2}) and ($GENE2{$last4} ne $GENE1{$last1})) {#negative order
							#	print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
							#}elsif (($GENE2{$last3} ne $GENE1{$last2}) and ($GENE2{$last4} eq $GENE1{$last1})) {#negative order
							#	print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
							#}
						}
					}elsif (($GENE2{$key2} eq $GENE1{$key1}) and ($hash_reference2{$GENE1{$key1}} == 2) and ($hash_target2{$GENE2{$key2}} == 1)) {
						#occurrence 2 times in reference, occurrence 1 time in target
						#change gene8 gene9 gene9 gene8 to gene8 gene9
						#print $GENE1{$key1}."\t".$GENE2{$key2}."\n";
						my $i=$cnt;
						if ($i == 0) {
							if (($GENE2{$last3} eq $GENE1{$last1}) and ($GENE2{$last4} eq $GENE1{$last2})) {#positive order
								$i++;
								print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
							}elsif (($GENE2{$last3} eq $GENE1{$last1}) and ($GENE2{$last4} ne $GENE1{$last2})) {#positive order
								$i++;
								print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
							}elsif (($GENE2{$last3} ne $GENE1{$last1}) and ($GENE2{$last4} eq $GENE1{$last2})) {#positive order
								$i++;
								print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
							}
							if (($GENE2{$last3} ne $GENE1{$last1}) and ($GENE2{$last4} ne $GENE1{$last2})) {#positive order
								if (defined $junctionLB_reference) {
									if ((($END1{$key1} >= $junctionLB_reference) and ($START1{$key1} <= $junctionSB_reference))) {#genes in IRb of reference vs genes in target
										#print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
									}elsif ((($END1{$key1} >= $junctionSA_reference) and ($START1{$key1} <= $junctionLA_reference))) {#genes in IRa of reference vs genes in target
										#print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
									}
								}
							}
						}elsif (($i != 1) or ($i != 2)) {
							#print $GENE1{$key1}."\t".$GENE2{$key2}."\n";
							#if (($GENE2{$last3} eq $GENE1{$last2}) and ($GENE2{$last4} eq $GENE1{$last1})) {#negative
							#	print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
							#}elsif (($GENE2{$last3} eq $GENE1{$last2}) and ($GENE2{$last4} ne $GENE1{$last1})) {#negative
							#	print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
							#}elsif (($GENE2{$last3} ne $GENE1{$last2}) and ($GENE2{$last4} eq $GENE1{$last1})) {#negative
							#	print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
							#}
						}
					}elsif (($GENE2{$key2} eq $GENE1{$key1}) and ($hash_reference2{$GENE1{$key1}} == 2) and ($hash_target2{$GENE2{$key2}} == 2)) {
						#occurrence 2 times in reference, occurrence 2 times in target
						#change gene6 gene6 gene7 gene7 gene7 gene7 gene6 gene6 to gene6 gene7 gene7 gene6
						#print $GENE1{$key1}."\t".$GENE2{$key2}."\n";
						my $i=$cnt;
						if ($i == 0) {
							if (($GENE1{$key1} ne "trnH-GUG") and ($GENE2{$key2} ne "trnH-GUG")) {
								if (($GENE2{$last3} eq $GENE1{$last1}) and ($GENE2{$last4} eq $GENE1{$last2})) {#positive order
									$i++;
									print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
								}elsif (($GENE2{$last3} eq $GENE1{$last1}) and ($GENE2{$last4} ne $GENE1{$last2})) {#positive order
									$i++;
									print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
								}elsif (($GENE2{$last3} ne $GENE1{$last1}) and ($GENE2{$last4} eq $GENE1{$last2})) {#positive order
									$i++;
									print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
								}
								if (($GENE2{$last3} ne $GENE1{$last1}) and ($GENE2{$last4} ne $GENE1{$last2})) {#positive order
									if ((defined $junctionLB_reference) and (defined $junctionLB_target)) {
										if ((($END1{$key1} >= $junctionLB_reference) and ($START1{$key1} <= $junctionSB_reference)) and (($END2{$key2} >= $junctionLB_target) and ($START2{$key2} <= $junctionSB_target))) {#genes in IRb of reference vs genes in IRb of target
											print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
										}elsif ((($END1{$key1} >= $junctionSA_reference) and ($START1{$key1} <= $junctionLA_reference)) and (($END2{$key2} >= $junctionSA_target) and ($START2{$key2} <= $junctionLA_target))) {#genes in IRa of reference vs genes in IRa of target
											print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
										}
									}
								}
							}elsif (($GENE1{$key1} eq "trnH-GUG") and ($GENE2{$key2} eq "trnH-GUG")) {
								print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
							}

						}elsif (($i != 1) or ($i != 2)) {
							#print $GENE1{$key1}."\t".$GENE2{$key2}."\n";
							if (($GENE2{$last3} eq $GENE1{$last2}) and ($GENE2{$last4} eq $GENE1{$last1})) {#negative
								print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
							}elsif (($GENE2{$last3} eq $GENE1{$last2}) and ($GENE2{$last4} ne $GENE1{$last1})) {#negative
								print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
							}elsif (($GENE2{$last3} ne $GENE1{$last2}) and ($GENE2{$last4} eq $GENE1{$last1})) {#negative
								print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
							}
						}
					}
				}elsif (($GENE1{$key1} eq "rps12") and ($GENE2{$key2} eq "rps12")) {#rps12 gene
					my $length_reference=abs($END1{$key1}-$START1{$key1});
					my $length_target=abs($END2{$key2}-$START2{$key2});
					if (($length_reference < 200) and ($length_target < 200)) {
						print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
					}elsif (($length_reference > 200) and ($length_target > 200)) {
						if (($GENE2{$last3} eq $GENE1{$last1}) and ($GENE2{$last4} eq $GENE1{$last2})) {#positive order
							print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
						}elsif (($GENE2{$last3} eq $GENE1{$last1}) and ($GENE2{$last4} ne $GENE1{$last2})) {#positive order
							print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
						}elsif (($GENE2{$last3} ne $GENE1{$last1}) and ($GENE2{$last4} eq $GENE1{$last2})) {#positive order
							print $out_genes_link $SP1{$key1}."\t".$START1{$key1}."\t".$END1{$key1}."\t".$SP2{$key2}."\t".$START2{$key2}."\t".$END2{$key2}."\n";
						}
					}
				}
			}
		}
	}
	close $out_genes_link;

	open(my $out_genes_unique_reference,">","$output_directory/genes.unique.reference.txt");
	open(my $out_genes_unique_target,">","$output_directory/genes.unique.target.txt");
	foreach my $key1 (sort {$hash_reference1{$a} <=> $hash_reference1{$b}} keys %hash_reference1) {
		my @row1=split /\t/,$key1;
		if (!exists $hash_target2{$row1[3]}) {
			print $out_genes_unique_reference $row1[0]."\t".$row1[1]."\t".$row1[2]."\t"."10"."\n";
		}
	}
	foreach my $key2 (sort {$hash_target1{$a} <=> $hash_target1{$b}} keys %hash_target1) {
		my @row2=split /\t/,$key2;
		if (!exists $hash_reference2{$row2[3]}) {
			print $out_genes_unique_target $row2[0]."\t".$row2[1]."\t".$row2[2]."\t"."10"."\n";
		}
	}
	close $out_genes_unique_reference;
	close $out_genes_unique_target;

}




####generate ticks.conf file
open(my $out_ticks,">","$output_directory/ticks.conf");
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
print $out_ticks "spacing        = 1u"."\n";#small tick
print $out_ticks "size           = 15p"."\n";
print $out_ticks "thickness      = 2p"."\n";
print $out_ticks "orientation    = in"."\n";
print $out_ticks "</tick>"."\n\n";

print $out_ticks "<tick>"."\n";
print $out_ticks "spacing        = 10u"."\n";#big tick
print $out_ticks "size           = 20p"."\n";
print $out_ticks "thickness      = 4p"."\n";
print $out_ticks "orientation    = in"."\n";

print $out_ticks "show_label     = yes"."\n";#big tick label
print $out_ticks "label_size     = 15p"."\n";#big tick label length
print $out_ticks "label_offset   = 0p"."\n";#big tick label radius
print $out_ticks "format         = %d"."\n\n";
print $out_ticks "suffix         = kb"."\n";#big tick label suffix
print $out_ticks "color          = black"."\n\n";

print $out_ticks "grid           = no"."\n";#big tick grid line
print $out_ticks "grid_start     = 0.88r"."\n";
print $out_ticks "grid_end       = 0.98r"."\n";
print $out_ticks "grid_color     = grey"."\n";
print $out_ticks "grid_thickness = 2p"."\n";
print $out_ticks "</tick>"."\n\n";

print $out_ticks "</ticks>"."\n";
close $out_ticks;




####generate ideogram.conf file
open(my $out_ideogram,">","$output_directory/ideogram.conf");
print $out_ideogram "<ideogram>"."\n\n";

print $out_ideogram "<spacing>"."\n";
print $out_ideogram "default = 0.01r"."\n";
print $out_ideogram "break   = 0"."\n";
print $out_ideogram "</spacing>"."\n\n\n";

print $out_ideogram "radius           = 0.8r"."\n";
print $out_ideogram "thickness        = 2p"."\n";
print $out_ideogram "fill             = yes"."\n\n";

print $out_ideogram "fill_color       = gpos100"."\n";
print $out_ideogram "stroke_thickness = 2"."\n";
print $out_ideogram "stroke_color     = black"."\n\n";

print $out_ideogram "</ideogram>"."\n";
close $out_ideogram;




####generate visualization.conf file
open(my $out_circos,">","$output_directory/homology.conf");
print $out_circos "<<include etc/colors_fonts_patterns.conf>>"."\n";
print $out_circos "<<include etc/housekeeping.conf>>"."\n\n";

print $out_circos "<<include ideogram.conf>>"."\n";
print $out_circos "<<include ticks.conf>>"."\n\n";

print $out_circos "<image>"."\n";
print $out_circos "<<include etc/image.conf>>"."\n";
print $out_circos "angle_orientation* = counterclockwise"."\n";
print $out_circos "#background* = transparent"."\n";
print $out_circos "#radius* = 1500p"."\n";
print $out_circos "angle_offset* = $angle"."\n";
print $out_circos "</image>"."\n\n";

print $out_circos "karyotype   = $output_directory/karyotype.txt"."\n";
print $out_circos "chromosomes_units           = 1000"."\n";
print $out_circos "chromosomes_display_default = yes"."\n";
print $out_circos "chromosomes_reverse = /$right_filename/"."\n\n\n\n\n";

print $out_circos "<plots>"."\n\n";

print $out_circos "<plot>"."\n";
print $out_circos "type  = text"."\n";
print $out_circos "file  = $output_directory/species_length.txt"."\n";
print $out_circos "r0    = 1.20r"."\n";
print $out_circos "r1    = 1.70r"."\n";
print $out_circos "label_parallel = yes"."\n";
print $out_circos "color = black"."\n";
print $out_circos "label_font = default"."\n";
print $out_circos "label_size = 30"."\n";
print $out_circos "</plot>"."\n\n\n";

print $out_circos "<plot>"."\n";
print $out_circos "type  = text"."\n";
print $out_circos "file  = $output_directory/$left_filename\_outergenes.label.txt"."\n\n";

print $out_circos "r0    = 1.07r"."\n";
print $out_circos "r1    = 1.22r"."\n\n";

print $out_circos "label_font = light"."\n";
print $out_circos "label_size = 15p"."\n\n";

print $out_circos "rpadding   = 0p"."\n\n";

print $out_circos "label_snuggle         = yes"."\n";
print $out_circos "max_snuggle_distance  = 1r"."\n";
print $out_circos "snuggle_tolerance     = 0.3r"."\n";
print $out_circos "snuggle_sampling      = 3"."\n";
print $out_circos "snuggle_refine        = yes"."\n";
print $out_circos "</plot>"."\n\n\n";

print $out_circos "<plot>"."\n";
print $out_circos "type  = text"."\n";
print $out_circos "file  = $output_directory/$left_filename\_innergenes.label.txt"."\n\n";

print $out_circos "r0    = 1.01r"."\n";
print $out_circos "r1    = 1.16r"."\n\n";

print $out_circos "label_font = light"."\n";
print $out_circos "label_size = 15p"."\n\n";

print $out_circos "rpadding   = 0p"."\n\n";

print $out_circos "label_snuggle         = yes"."\n";
print $out_circos "max_snuggle_distance  = 1r"."\n";
print $out_circos "snuggle_tolerance     = 0.3r"."\n";
print $out_circos "snuggle_sampling      = 3"."\n";
print $out_circos "snuggle_refine        = yes"."\n";
print $out_circos "</plot>"."\n\n\n";

print $out_circos "<plot>"."\n";
print $out_circos "type  = text"."\n";
print $out_circos "file  = $output_directory/$right_filename\_outergenes.label.txt"."\n\n";

print $out_circos "r0    = 1.07r"."\n";
print $out_circos "r1    = 1.22r"."\n\n";

print $out_circos "label_font = light"."\n";
print $out_circos "label_size = 15p"."\n\n";

print $out_circos "rpadding   = 0p"."\n\n";

print $out_circos "label_snuggle         = yes"."\n";
print $out_circos "max_snuggle_distance  = 1r"."\n";
print $out_circos "snuggle_tolerance     = 0.3r"."\n";
print $out_circos "snuggle_sampling      = 3"."\n";
print $out_circos "snuggle_refine        = yes"."\n";
print $out_circos "</plot>"."\n\n\n";

print $out_circos "<plot>"."\n";
print $out_circos "type  = text"."\n";
print $out_circos "file  = $output_directory/$right_filename\_innergenes.label.txt"."\n\n";

print $out_circos "r0    = 1.01r"."\n";
print $out_circos "r1    = 1.16r"."\n\n";

print $out_circos "label_font = light"."\n";
print $out_circos "label_size = 15p"."\n\n";

print $out_circos "rpadding   = 0p"."\n\n";

print $out_circos "label_snuggle         = yes"."\n";
print $out_circos "max_snuggle_distance  = 1r"."\n";
print $out_circos "snuggle_tolerance     = 0.3r"."\n";
print $out_circos "snuggle_sampling      = 3"."\n";
print $out_circos "snuggle_refine        = yes"."\n";
print $out_circos "</plot>"."\n\n\n";

my $genes_unique_reference="$output_directory/genes.unique.reference.txt";
my $genes_unique_target="$output_directory/genes.unique.target.txt";
if (-s $genes_unique_reference) {
	print $out_circos "<plot>"."\n";
	print $out_circos "type = scatter"."\n";#line scatter histogram heatmap text tile connector
	print $out_circos "file = $output_directory/genes.unique.reference.txt"."\n";
	print $out_circos "r0 = 0.900r"."\n";
	print $out_circos "r1 = 0.931r"."\n";
	print $out_circos "min = 10"."\n";
	print $out_circos "max = 10"."\n";
	print $out_circos "glyph = triangle"."\n";#circle triangle rectangle
	print $out_circos "glyph_size = 18p"."\n";
	print $out_circos "color = chr5"."\n";
	#print $out_circos "stroke_color = black"."\n";
	#print $out_circos "stroke_thickness = 1p"."\n";
	print $out_circos "</plot>"."\n\n\n";
}

if (-s $genes_unique_target) {
	print $out_circos "<plot>"."\n";
	print $out_circos "type = scatter"."\n";#line scatter histogram heatmap text tile connector
	print $out_circos "file = $output_directory/genes.unique.target.txt"."\n";
	print $out_circos "r0 = 0.900r"."\n";
	print $out_circos "r1 = 0.931r"."\n";
	print $out_circos "min = 10"."\n";
	print $out_circos "max = 10"."\n";
	print $out_circos "glyph = triangle"."\n";#circle triangle rectangle
	print $out_circos "glyph_size = 18p"."\n";
	print $out_circos "color = chr5"."\n";
	#print $out_circos "stroke_color = black"."\n";
	#print $out_circos "stroke_thickness = 1p"."\n";
	print $out_circos "</plot>"."\n\n\n";
}

print $out_circos "<plot>"."\n";
print $out_circos "type = highlight"."\n";
print $out_circos "file = $output_directory/$left_filename\_quadripartite.IRb.IRa.txt"."\n";
print $out_circos "r0 = 0.996r"."\n";
print $out_circos "r1 = 1.004r"."\n";
print $out_circos "</plot>"."\n\n\n";

print $out_circos "<plot>"."\n";
print $out_circos "type = highlight"."\n";
print $out_circos "file = $output_directory/$right_filename\_quadripartite.IRb.IRa.txt"."\n";
print $out_circos "r0 = 0.996r"."\n";
print $out_circos "r1 = 1.004r"."\n";
print $out_circos "</plot>"."\n\n\n";

print $out_circos "</plots>"."\n\n\n\n\n";


print $out_circos "<highlights>"."\n";
print $out_circos "z = 0"."\n";
print $out_circos "fill_color = white"."\n\n";

print $out_circos "<highlight>"."\n";
print $out_circos "file       = $output_directory/$left_filename\_outergenes.position.txt"."\n";
print $out_circos "r0         = 1.00r"."\n";
print $out_circos "r1         = 1.06r"."\n";
print $out_circos "stroke_color = black"."\n";
print $out_circos "stroke_thickness = 1"."\n";
print $out_circos "</highlight>"."\n\n";

print $out_circos "<highlight>"."\n";
print $out_circos "file       = $output_directory/$left_filename\_innergenes.position.txt"."\n";
print $out_circos "r0         = 0.94r"."\n";
print $out_circos "r1         = 1.00r"."\n";
print $out_circos "stroke_color = black"."\n";
print $out_circos "stroke_thickness = 1"."\n";
print $out_circos "</highlight>"."\n\n";

print $out_circos "<highlight>"."\n";
print $out_circos "file       = $output_directory/$right_filename\_outergenes.position.txt"."\n";
print $out_circos "r0         = 1.00r"."\n";
print $out_circos "r1         = 1.06r"."\n";
print $out_circos "stroke_color = black"."\n";
print $out_circos "stroke_thickness = 1"."\n";
print $out_circos "</highlight>"."\n\n";

print $out_circos "<highlight>"."\n";
print $out_circos "file       = $output_directory/$right_filename\_innergenes.position.txt"."\n";
print $out_circos "r0         = 0.94r"."\n";
print $out_circos "r1         = 1.00r"."\n";
print $out_circos "stroke_color = black"."\n";
print $out_circos "stroke_thickness = 1"."\n";
print $out_circos "</highlight>"."\n\n";

print $out_circos "</highlights>"."\n\n\n\n\n";


print $out_circos "<links>"."\n\n";

print $out_circos "<link>"."\n";
print $out_circos "file          = $output_directory/genes.link.txt"."\n";
print $out_circos "radius        = 0.922r"."\n";
print $out_circos "bezier_radius = 0r"."\n";
print $out_circos "color         = green"."\n";
print $out_circos "thickness     = 2r"."\n";

if (($ribbon eq "no") or ($ribbon eq "n") or ($ribbon eq "NO") or ($ribbon eq "No") or ($ribbon eq "N")) {
	print $out_circos "ribbon        = no"."\n";#no
}elsif (($ribbon eq "yes") or ($ribbon eq "y") or ($ribbon eq "YES") or ($ribbon eq "Yes") or ($ribbon eq "Y")) {
	print $out_circos "ribbon        = yes"."\n";#yes
}

if ($crest == 0.01) {
	print $out_circos "crest         = 0.01"."\n";#0.01
}elsif ($crest == 0.25) {
	print $out_circos "crest         = 0.25"."\n";#0.25
}elsif ($crest == 0.50) {
	print $out_circos "crest         = 0.50"."\n";#0.50, display best
}elsif ($crest == 0.75) {
	print $out_circos "crest         = 0.75"."\n";#0.75
}elsif ($crest == 1.00) {
	print $out_circos "crest         = 1.00"."\n";#1.00, circos default
}

if ($bezier == 0.01) {
	print $out_circos "bezier_radius_purity = 0.01"."\n";#0.01
}elsif ($bezier == 0.25) {
	print $out_circos "bezier_radius_purity = 0.25"."\n";#0.25, display best
}elsif ($bezier == 0.50) {
	print $out_circos "bezier_radius_purity = 0.50"."\n";#0.50
}elsif ($bezier == 0.75) {
	print $out_circos "bezier_radius_purity = 0.75"."\n";#0.75, default
}elsif ($bezier == 1.00) {
	print $out_circos "bezier_radius_purity = 1.00"."\n";#1.00
}

print $out_circos "</link>"."\n\n";

print $out_circos "</links>"."\n";

close $out_circos;




####batch run circos
#perl bin\circos -conf test\input\circos.conf -outputdir test\output -outputfile circos -svg -png

my $command1="perl bin\\circos -conf $output_directory\\homology.conf -outputdir $output_directory -outputfile homology -svg -png";
#my $command2="perl bin\\circos -conf $output_directory\\homology.conf -outputdir $output_directory -outputfile homology -debug_group textplace";
system("$command1");
#system("$command2");




####functions
sub argument{
	my @options=("help|h","reference|r:s","target|t:s","minir|m:i","angle|a:i","pattern|p:s","ribbon|i:s","crest|c:s","bezier|b:s","output|o:s");
	my %options;
	GetOptions(\%options,@options);
	exec ("pod2usage $0") if ((keys %options)==0 or $options{'h'} or $options{'help'});
	if(!exists $options{'reference'}){
		print "***ERROR: No reference filename is assigned!!!\n";
		exec ("pod2usage $0");
	}
	if(!exists $options{'target'}){
		print "***ERROR: No target filename is assigned!!!\n";
		exec ("pod2usage $0");
	}
	#if(!exists $options{'minir'}){
	#	print "***ERROR: No minimum allowed length for the longest inverted-repeat (IR) is assigned!!!\n";
	#	exec ("pod2usage $0");
	#}
	#if(!exists $options{'angle'}){
	#	print "***ERROR: No rotation angle of the homologous gene map is assigned!!!\n";
	#	exec ("pod2usage $0");
	#}
	#if(!exists $options{'pattern'}){
	#	print "***ERROR: No pattern for linking syntenic genes is assigned!!!\n";
	#	exec ("pod2usage $0");
	#}
	#if(!exists $options{'ribbon'}){
	#	print "***ERROR: No ribbon (yes) or line (no) for linking syntenic genes is assigned!!!\n";
	#	exec ("pod2usage $0");
	#}
	#if(!exists $options{'crest'}){
	#	print "***ERROR: No crest peak distance of linking lines or ribbons is assigned!!!\n";
	#	exec ("pod2usage $0");
	#}
	#if(!exists $options{'bezier'}){
	#	print "***ERROR: No bezier curve degree of linking lines or ribbons is assigned!!!\n";
	#	exec ("pod2usage $0");
	#}
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

    gene_homology_v1.pl Connect homologous genes between reference and target plastomes version 1.0

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

    Connect homologous genes between reference and target plastomes

=head1 SYNOPSIS

    gene_homology_v1.pl -r -t -m -a -p -i -c -b -o
    Copyright (C) 2025 Xiao-Jian Qu
    Please contact <quxiaojian@sdnu.edu.cn>, if you have any bugs or questions.

    [-h -help]         help information.
    [-r -reference]    required: (default: reference.gb) input reference file of GenBank-formatted file.
    [-t -target]       required: (default: target.gb) input target file of GenBank-formatted file.
    [-m -minir]        required: (default: 1000) minimum allowed length for the longest inverted-repeat (IR).
    [-a -angle]        required: (default: -90) any rotation angle of the homologous gene map, -90=270.
    [-p -pattern]      required: (default: 1to1) patterns (1to1/1tomore) for linking homologous genes.
    [-i -ribbon]       required: (default: no) ribbon (yes/y) or line (no/n) for linking homologous genes.
    [-c -crest]        required: (default: 0.50) crest peak distance of linking lines or ribbons, and the
                       alternative values are [0.01/0.25/0.50/0.75/1.00].
    [-b -bezier]       required: (default: 0.25) bezier curve degree of linking lines or ribbons, and the
                       alternative values are [0.01/0.25/0.50/0.75/1.00].
    [-o -output]       required: (default: output) output directory name.

=cut
