#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Find;
use Data::Dumper;
$|=1;

my $global_options=&argument();
my $input_directory=&default("input","input");
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
my @filenames;
find(\&target,$input_directory);
sub target{
	if (/$pattern1/ or /$pattern2/ or /$pattern3/){
		push @filenames,"$File::Find::name";
	}
	return;
}


while (@filenames) {
	my $filename_gb=shift @filenames;
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
	open (my $in_gb3,"<",$target_name_temp_random2);
	while (<$in_gb3>){
		$_=~ s/\r|\n//g;
		@row_array=split /\s+/,$_;
		if (/^LOCUS/i){
			$species_name=$filename;
			$length=$row_array[2];
		}elsif(/ {5}CDS {13}/ or / {5}tRNA {12}/ or / {5}rRNA {12}/){
			if ($row_array[2]=~ /^\d+..\d+$/){#
				$row_array[2]="\+\t$row_array[2]\t$row_array[1]";
				$row_array[2]=~ s/\../\t/g;
			}elsif($row_array[2]=~/^complement\((\d+..\d+)\)$/){#
				$row_array[2]="-\t$1\t$row_array[1]";
				$row_array[2]=~ s/\../\t/g;
			}elsif($row_array[2]=~ /^join\((\d+..\d+),(\d+..\d+)\)$/) {#
				$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
				$row_array[2]=~ s/\../\t/g;
			}elsif($row_array[2]=~ /^join\((\d+..\d+),(\d+..\d+),(\d+..\d+)\)$/) {#
				$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
				$row_array[2]=~ s/\../\t/g;
			}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){#
				$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]";
				$row_array[2]=~ s/\../\t/g;
			}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),(\d+..\d+)\)$/){#
				$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
				$row_array[2]=~ s/\../\t/g;
			}elsif($row_array[2]=~ /^join\((\d+..\d+),complement\((\d+..\d+)\)\)$/){#
				$row_array[2]="+\t$1\t$row_array[1]\t-\t$2\t$row_array[1]";
				$row_array[2]=~ s/\../\t/g;
			}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){#
				$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
				$row_array[2]=~ s/\../\t/g;
			}elsif($row_array[2]=~ /^complement\(join\((\d+..\d+),(\d+..\d+)\)\)$/){#
				$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]";
				$row_array[2]=~ s/\../\t/g;
			}elsif($row_array[2]=~ /^complement\(join\((\d+..\d+),(\d+..\d+),(\d+..\d+)\)\)$/){#
				$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
				$row_array[2]=~ s/\../\t/g;
			}elsif($row_array[2]=~ /^order\((\d+..\d+),(\d+..\d+)\)$/){#
				$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
				$row_array[2]=~ s/\../\t/g;
			}elsif($row_array[2]=~ /^order\((\d+..\d+),(\d+..\d+),(\d+..\d+)\)$/){#
				$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
				$row_array[2]=~ s/\../\t/g;
			}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){#
				$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]";
				$row_array[2]=~ s/\../\t/g;
			}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){#
				$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
				$row_array[2]=~ s/\../\t/g;
			}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),(\d+..\d+)\)$/){#
				$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
				$row_array[2]=~ s/\../\t/g;
			}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),(\d+..\d+),(\d+..\d+)\)$/){#
				$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
				$row_array[2]=~ s/\../\t/g;
			}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),(\d+..\d+),(\d+..\d+)\)$/){#
				$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
				$row_array[2]=~ s/\../\t/g;
			}elsif($row_array[2]=~ /^join\((\d+..\d+),(\d+..\d+),complement\((\d+..\d+)\)\)$/) {#
				$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
				$row_array[2]=~ s/\../\t/g;
			}elsif($row_array[2]=~ /^<(\d+..\d+)$/){# positive split no-intron gene
				$row_array[2]="\+\t$1\t$row_array[1]";
				$row_array[2]=~ s/\../\t/g;
			}elsif($row_array[2]=~ /^(\d+..>\d+)$/){# positive split no-intron gene
				$row_array[2]="\+\t$1\t$row_array[1]";
				$row_array[2]=~ s/\..>/\t/g;
			}elsif($row_array[2]=~/^complement\(<(\d+..\d+)\)$/){# negative split no-intron gene
				$row_array[2]="-\t$1\t$row_array[1]";
				$row_array[2]=~ s/\../\t/g;
			}elsif($row_array[2]=~/^complement\((\d+..>\d+)\)$/){# negative split no-intron gene
				$row_array[2]="-\t$1\t$row_array[1]";
				$row_array[2]=~ s/\..>/\t/g;
			}
			$element=$row_array[2];
			$mark=1;
		}elsif(/^\s+\/gene="(.*)"/ and $mark == 1){
			$element=$1.":".$element;
			push @genearray,$element;
			$element=();
			$mark=0;
		}
	}
	close $in_gb3;

	foreach (@genearray){
		my @array=split /:/,$_;
		push @output1,"$array[0]\t$array[1]\n";
	}
	unlink ("$target_name_temp_random1");
	unlink ("$target_name_temp_random2");


	#sort_bed_file
	my $col=2;
	my (%sort,@output_temp);
	foreach (@output1){
		my @row=split /\t/,$_;
		$sort{$_}=$row[$col];
	}
	foreach (sort {$sort{$a} <=> $sort{$b}} keys %sort){
		push @output_temp,"$_\n";
	}
	@output1=();


	#edit_bed_file
	my (%GENE1,%STRAND1,%START1,%END1,%TYPE1,%STRAND2,%START2,%END2,%TYPE2,%STRAND3,%START3,%END3,%TYPE3,@output2);
	my $cnt1=0;
	foreach (@output_temp) {
		$_=~ s/\r|\n//g;
		$cnt1++;
		($GENE1{$cnt1},$STRAND1{$cnt1},$START1{$cnt1},$END1{$cnt1},$TYPE1{$cnt1},$STRAND2{$cnt1},$START2{$cnt1},$END2{$cnt1},$TYPE2{$cnt1},$STRAND3{$cnt1},$START3{$cnt1},$END3{$cnt1},$TYPE3{$cnt1})=(split /\s+/,$_)[0,1,2,3,4,5,6,7,8,9,10,11,12];
	}
	@output_temp=();

	foreach (1..$cnt1) {
		if (defined $STRAND2{$_} eq "") {
			if ($TYPE1{$_} eq "CDS") {
				if ($STRAND1{$_} eq "-") {
					push @output2,("<"."\t".$START1{$_}."\t".$END1{$_}."\t".$GENE1{$_}."\n");
					push @output2,($START1{$_}."\t".$END1{$_}."\t"."exon\n");
				}elsif ($STRAND1{$_} eq "+") {
					push @output2,(">"."\t".$START1{$_}."\t".$END1{$_}."\t".$GENE1{$_}."\n");
					push @output2,($START1{$_}."\t".$END1{$_}."\t"."exon\n");
				}
			}elsif ($TYPE1{$_} eq "tRNA") {
				if ($STRAND1{$_} eq "-") {
					push @output2,("<"."\t".$START1{$_}."\t".$END1{$_}."\t".$GENE1{$_}."\n");
					push @output2,($START1{$_}."\t".$END1{$_}."\t"."utr\n");
				}elsif ($STRAND1{$_} eq "+") {
					push @output2,(">"."\t".$START1{$_}."\t".$END1{$_}."\t".$GENE1{$_}."\n");
					push @output2,($START1{$_}."\t".$END1{$_}."\t"."utr\n");
				}
			}elsif ($TYPE1{$_} eq "rRNA") {
				if ($STRAND1{$_} eq "-") {
					push @output2,("<"."\t".$START1{$_}."\t".$END1{$_}."\t".$GENE1{$_}."\n");
					push @output2,($START1{$_}."\t".$END1{$_}."\t"."utr\n");
				}elsif ($STRAND1{$_} eq "+") {
					push @output2,(">"."\t".$START1{$_}."\t".$END1{$_}."\t".$GENE1{$_}."\n");
					push @output2,($START1{$_}."\t".$END1{$_}."\t"."utr\n");
				}
			}
		}elsif ((defined $STRAND2{$_} ne "") and (defined $STRAND3{$_} eq "")) {
			if ($TYPE1{$_} eq "CDS") {
				if (($STRAND1{$_} eq "-") and ($START1{$_} < $START2{$_})){
					push @output2,("<"."\t".$START1{$_}."\t".$END2{$_}."\t".$GENE1{$_}."\n");
					push @output2,($START1{$_}."\t".$END1{$_}."\t"."exon\n");
					push @output2,($START2{$_}."\t".$END2{$_}."\t"."exon\n");
				}elsif(($STRAND1{$_} eq "-") and ($START1{$_} > $START2{$_})){
					push @output2,("<"."\t".$START2{$_}."\t".$END1{$_}."\t".$GENE1{$_}."\n");
					push @output2,($START2{$_}."\t".$END2{$_}."\t"."exon\n");
					push @output2,($START1{$_}."\t".$END1{$_}."\t"."exon\n");
				}elsif(($STRAND1{$_} eq "+") and ($START1{$_} < $START2{$_})){
					push @output2,(">"."\t".$START1{$_}."\t".$END2{$_}."\t".$GENE1{$_}."\n");
					push @output2,($START1{$_}."\t".$END1{$_}."\t"."exon\n");
					push @output2,($START2{$_}."\t".$END2{$_}."\t"."exon\n");
				}elsif(($STRAND1{$_} eq "+") and ($START1{$_} > $START2{$_})){
					push @output2,(">"."\t".$START2{$_}."\t".$END1{$_}."\t".$GENE1{$_}."\n");
					push @output2,($START2{$_}."\t".$END2{$_}."\t"."exon\n");
					push @output2,($START1{$_}."\t".$END1{$_}."\t"."exon\n");
				}
			}elsif ($TYPE1{$_} eq "tRNA") {
				if (($STRAND1{$_} eq "-") and ($START1{$_} < $START2{$_})){
					push @output2,("<"."\t".$START1{$_}."\t".$END2{$_}."\t".$GENE1{$_}."\n");
					push @output2,($START1{$_}."\t".$END1{$_}."\t"."utr\n");
					push @output2,($START2{$_}."\t".$END2{$_}."\t"."utr\n");
				}elsif(($STRAND1{$_} eq "-") and ($START1{$_} > $START2{$_})){
					push @output2,("<"."\t".$START2{$_}."\t".$END1{$_}."\t".$GENE1{$_}."\n");
					push @output2,($START2{$_}."\t".$END2{$_}."\t"."utr\n");
					push @output2,($START1{$_}."\t".$END1{$_}."\t"."utr\n");
				}elsif(($STRAND1{$_} eq "+") and ($START1{$_} < $START2{$_})){
					push @output2,(">"."\t".$START1{$_}."\t".$END2{$_}."\t".$GENE1{$_}."\n");
					push @output2,($START1{$_}."\t".$END1{$_}."\t"."utr\n");
					push @output2,($START2{$_}."\t".$END2{$_}."\t"."utr\n");
				}elsif(($STRAND1{$_} eq "+") and ($START1{$_} > $START2{$_})){
					push @output2,(">"."\t".$START2{$_}."\t".$END1{$_}."\t".$GENE1{$_}."\n");
					push @output2,($START2{$_}."\t".$END2{$_}."\t"."utr\n");
					push @output2,($START1{$_}."\t".$END1{$_}."\t"."utr\n");
				}
			}elsif ($TYPE1{$_} eq "rRNA") {
				if (($STRAND1{$_} eq "-") and ($START1{$_} < $START2{$_})){
					push @output2,("<"."\t".$START1{$_}."\t".$END2{$_}."\t".$GENE1{$_}."\n");
					push @output2,($START1{$_}."\t".$END1{$_}."\t"."utr\n");
					push @output2,($START2{$_}."\t".$END2{$_}."\t"."utr\n");
				}elsif(($STRAND1{$_} eq "-") and ($START1{$_} > $START2{$_})){
					push @output2,("<"."\t".$START2{$_}."\t".$END1{$_}."\t".$GENE1{$_}."\n");
					push @output2,($START2{$_}."\t".$END2{$_}."\t"."utr\n");
					push @output2,($START1{$_}."\t".$END1{$_}."\t"."utr\n");
				}elsif(($STRAND1{$_} eq "+") and ($START1{$_} < $START2{$_})){
					push @output2,(">"."\t".$START1{$_}."\t".$END2{$_}."\t".$GENE1{$_}."\n");
					push @output2,($START1{$_}."\t".$END1{$_}."\t"."utr\n");
					push @output2,($START2{$_}."\t".$END2{$_}."\t"."utr\n");
				}elsif(($STRAND1{$_} eq "+") and ($START1{$_} > $START2{$_})){
					push @output2,(">"."\t".$START2{$_}."\t".$END1{$_}."\t".$GENE1{$_}."\n");
					push @output2,($START2{$_}."\t".$END2{$_}."\t"."utr\n");
					push @output2,($START1{$_}."\t".$END1{$_}."\t"."utr\n");
				}
			}
		}elsif ((defined $STRAND2{$_} ne "") and (defined $STRAND3{$_} ne "")) {
			if (($STRAND1{$_} eq "-") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_})){
				push @output2,("<"."\t".$START1{$_}."\t".$END3{$_}."\t".$GENE1{$_}."\n");
				push @output2,($START1{$_}."\t".$END1{$_}."\t"."exon\n");
				push @output2,($START2{$_}."\t".$END2{$_}."\t"."exon\n");
				push @output2,($START3{$_}."\t".$END3{$_}."\t"."exon\n");
			}elsif(($STRAND1{$_} eq "-") and ($START1{$_} > $START2{$_}) and ($START2{$_} > $START3{$_})){
				push @output2,("<"."\t".$START3{$_}."\t".$END1{$_}."\t".$GENE1{$_}."\n");
				push @output2,($START3{$_}."\t".$END3{$_}."\t"."exon\n");
				push @output2,($START2{$_}."\t".$END2{$_}."\t"."exon\n");
				push @output2,($START1{$_}."\t".$END1{$_}."\t"."exon\n");
			}elsif(($STRAND1{$_} eq "-") and ($START1{$_} > $START3{$_}) and ($START3{$_} > $START2{$_})){
				push @output2,("<"."\t".$START2{$_}."\t".$END1{$_}."\t".$GENE1{$_}."\n");
				push @output2,($START2{$_}."\t".$END2{$_}."\t"."exon\n");
				push @output2,($START3{$_}."\t".$END3{$_}."\t"."exon\n");
				push @output2,($START1{$_}."\t".$END1{$_}."\t"."exon\n");
			}elsif(($STRAND1{$_} eq "+") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_})){
				push @output2,(">"."\t".$START1{$_}."\t".$END3{$_}."\t".$GENE1{$_}."\n");
				push @output2,($START1{$_}."\t".$END1{$_}."\t"."exon\n");
				push @output2,($START2{$_}."\t".$END2{$_}."\t"."exon\n");
				push @output2,($START3{$_}."\t".$END3{$_}."\t"."exon\n");
			}elsif(($STRAND1{$_} eq "+") and ($START1{$_} > $START2{$_}) and ($START2{$_} > $START3{$_})){
				push @output2,(">"."\t".$START3{$_}."\t".$END1{$_}."\t".$GENE1{$_}."\n");
				push @output2,($START3{$_}."\t".$END3{$_}."\t"."exon\n");
				push @output2,($START2{$_}."\t".$END2{$_}."\t"."exon\n");
				push @output2,($START1{$_}."\t".$END1{$_}."\t"."exon\n");
			}elsif(($STRAND1{$_} eq "+") and ($START1{$_} > $START3{$_}) and ($START3{$_} > $START2{$_})){
				push @output2,(">"."\t".$START2{$_}."\t".$END1{$_}."\t".$GENE1{$_}."\n");
				push @output2,($START2{$_}."\t".$END2{$_}."\t"."exon\n");
				push @output2,($START3{$_}."\t".$END3{$_}."\t"."exon\n");
				push @output2,($START1{$_}."\t".$END1{$_}."\t"."exon\n");
			}
		}
	}

	#output_bed_file
	open (my $out_bed,">","$output_directory/$filename\_mVISTA.txt");
	foreach (@output2){
		print $out_bed $_;
	}
	close $out_bed;
	@output2=();
}




##function
sub argument{
	my @options=("help|h","input|i:s","output|o:s");
	my %options;
	GetOptions(\%options,@options);
	exec ("pod2usage $0") if ((keys %options)==0 or $options{'h'} or $options{'help'});
	if(!exists $options{'input'}){
		print "***ERROR: No input directory is assigned!!!\n";
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

    mVISTA_v1.pl version 1.0

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

    generate input file for mVISTA alignment

=head1 SYNOPSIS

    mVISTA_v1.pl -i -o
    example: perl generate_input_file_for_mVISTA_alignment.pl -i input -o output
    Copyright (C) 2025 Xiao-Jian Qu
    Please contact <quxiaojian@sdnu.edu.cn>, if you have any bugs or questions.

    [-h -help]         help information.
    [-i -input]        required: (default: input) input directory name containing GenBank flatfiles.
    [-o -output]       required: (default: output) output directory name.

=cut
