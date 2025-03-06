#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Find;
use Data::Dumper;
$|=1;

my $global_options=&argument();
my $input_directory=&default("input","input");
my $site=&default("1","site");
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


my %hash_codon=("TAA"=>"*","TAG"=>"*","TGA"=>"*","TCA"=>"S","TCC"=>"S","TCG"=>"S","TCT"=>"S","TTC"=>"F","TTT"=>"F","TTA"=>"L","TTG"=>"L","TAC"=>"Y","TAT"=>"Y","TGC"=>"C","TGT"=>"C","TGG"=>"W","CTA"=>"L","CTC"=>"L","CTG"=>"L","CTT"=>"L","CCA"=>"P","CCC"=>"P","CCG"=>"P","CCT"=>"P","CAC"=>"H","CAT"=>"H","CAA"=>"Q","CAG"=>"Q","CGA"=>"R","CGC"=>"R","CGG"=>"R","CGT"=>"R","ATA"=>"I","ATC"=>"I","ATT"=>"I","ATG"=>"M","ACA"=>"T","ACC"=>"T","ACG"=>"T","ACT"=>"T","AAC"=>"N","AAT"=>"N","AAA"=>"K","AAG"=>"K","AGC"=>"S","AGT"=>"S","AGA"=>"R","AGG"=>"R","GTA"=>"V","GTC"=>"V","GTG"=>"V","GTT"=>"V","GCA"=>"A","GCC"=>"A","GCG"=>"A","GCT"=>"A","GAC"=>"D","GAT"=>"D","GAA"=>"E","GAG"=>"E","GGA"=>"G","GGC"=>"G","GGG"=>"G","GGT"=>"G");
#my %hash_codon=("TCA"=>"S","TCC"=>"S","TCG"=>"S","TCT"=>"S","TTC"=>"F","TTT"=>"F","TTA"=>"L","TTG"=>"L","TAC"=>"Y","TAT"=>"Y","TGC"=>"C","TGT"=>"C","TGG"=>"W","CTA"=>"L","CTC"=>"L","CTG"=>"L","CTT"=>"L","CCA"=>"P","CCC"=>"P","CCG"=>"P","CCT"=>"P","CAC"=>"H","CAT"=>"H","CAA"=>"Q","CAG"=>"Q","CGA"=>"R","CGC"=>"R","CGG"=>"R","CGT"=>"R","ATA"=>"I","ATC"=>"I","ATT"=>"I","ATG"=>"M","ACA"=>"T","ACC"=>"T","ACG"=>"T","ACT"=>"T","AAC"=>"N","AAT"=>"N","AAA"=>"K","AAG"=>"K","AGC"=>"S","AGT"=>"S","AGA"=>"R","AGG"=>"R","GTA"=>"V","GTC"=>"V","GTG"=>"V","GTT"=>"V","GCA"=>"A","GCC"=>"A","GCG"=>"A","GCT"=>"A","GAC"=>"D","GAT"=>"D","GAA"=>"E","GAG"=>"E","GGA"=>"G","GGC"=>"G","GGG"=>"G","GGT"=>"G");

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
	my ($header,$sequence);
	while (defined ($header=<$input>) and defined ($sequence=<$input>)){
		$header=~ s/\r|\n//g;
		$sequence=~ s/\r|\n//g;

		my $length=length $sequence;
		my $remainder=$length % 3;
		if ($remainder != 0) {
			open(my $warning,">>","$output_directory/bad_codons.txt");
			print $warning "The length of gene $filename in species $header is not a multiple of 3!\n";
			close $warning;
		}elsif ($remainder == 0) {
			if ($site == 1) {#1=site1 + site2 + site3
				my (@site1,@site2,@site3);
				for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length($sequence)-1) or (length($sequence)-2) is OK
				#for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
					my $site1=substr ($sequence,$i,1);
					my $site2=substr ($sequence,$i+1,1);
					my $site3=substr ($sequence,$i+2,1);
					push @site1,$site1;
					push @site2,$site2;
					push @site3,$site3;

					my $codon=substr ($sequence,$i,3);
					$codon=uc $codon;
					if (exists $hash_codon{$codon}){
						if (($hash_codon{$codon} eq "*") and ($i != ($length-3))) {
							my $j=$i+1;
							open(my $warning,">>","$output_directory/bad_codons.txt");
							print $warning "Bad codon $codon in position $j of $filename in species $header!\n";
							close $warning;
						}
					}else{
						my $j=$i+1;
						open(my $warning,">>","$output_directory/bad_codons.txt");
						print $warning "Bad codon $codon in position $j of $filename in species $header!\n";
						close $warning;
					}
				}

				open(my $output1,">>","$output_directory/$filename\_site1.fasta");
				open(my $output2,">>","$output_directory/$filename\_site2.fasta");
				open(my $output3,">>","$output_directory/$filename\_site3.fasta");

				print $output1 $header."\n";
				my $site1=join "",@site1;
				print $output1 $site1."\n";

				print $output2 $header."\n";
				my $site2=join "",@site2;
				print $output2 $site2."\n";

				print $output3 $header."\n";
				my $site3=join "",@site3;
				print $output3 $site3."\n";

				close $output1;
				close $output2;
				close $output3;
			}elsif ($site == 2) {#2=site12 + site3
				my (@site1,@site2);
				for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length($sequence)-1) or (length($sequence)-2) is OK
				#for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
					my $site1=substr ($sequence,$i,2);
					my $site2=substr ($sequence,$i+2,1);
					push @site1,$site1;
					push @site2,$site2;

					my $codon=substr ($sequence,$i,3);
					$codon=uc $codon;
					if (exists $hash_codon{$codon}){
						if (($hash_codon{$codon} eq "*") and ($i != ($length-3))) {
							my $j=$i+1;
							open(my $warning,">>","$output_directory/bad_codons.txt");
							print $warning "Bad codon $codon in position $j of $filename in species $header!\n";
							close $warning;
						}
					}else{
						my $j=$i+1;
						open(my $warning,">>","$output_directory/bad_codons.txt");
						print $warning "Bad codon $codon in position $j of $filename in species $header!\n";
						close $warning;
					}
				}

				open(my $output1,">>","$output_directory/$filename\_site12.fasta");
				open(my $output2,">>","$output_directory/$filename\_site3.fasta");

				print $output1 $header."\n";
				my $site1=join "",@site1;
				print $output1 $site1."\n";

				print $output2 $header."\n";
				my $site2=join "",@site2;
				print $output2 $site2."\n";

				close $output1;
				close $output2;
			}elsif ($site == 3) {#3=site1 + site23.
				my (@site1,@site2);
				for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length($sequence)-1) or (length($sequence)-2) is OK
				#for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
					my $site1=substr ($sequence,$i,1);
					my $site2=substr ($sequence,$i+1,2);
					push @site1,$site1;
					push @site2,$site2;

					my $codon=substr ($sequence,$i,3);
					$codon=uc $codon;
					if (exists $hash_codon{$codon}){
						if (($hash_codon{$codon} eq "*") and ($i != ($length-3))) {
							my $j=$i+1;
							open(my $warning,">>","$output_directory/bad_codons.txt");
							print $warning "Bad codon $codon in position $j of $filename in species $header!\n";
							close $warning;
						}
					}else{
						my $j=$i+1;
						open(my $warning,">>","$output_directory/bad_codons.txt");
						print $warning "Bad codon $codon in position $j of $filename in species $header!\n";
						close $warning;
					}
				}

				open(my $output1,">>","$output_directory/$filename\_site1.fasta");
				open(my $output2,">>","$output_directory/$filename\_site23.fasta");

				print $output1 $header."\n";
				my $site1=join "",@site1;
				print $output1 $site1."\n";

				print $output2 $header."\n";
				my $site2=join "",@site2;
				print $output2 $site2."\n";

				close $output1;
				close $output2;
			}
		}
	}
	close $input;
	unlink("$output_filename_temp");
}




##function
sub argument{
	my @options=("help|h","input|i:s","site|s:i","output|o:s");
	my %options;
	GetOptions(\%options,@options);
	exec ("pod2usage $0") if ((keys %options)==0 or $options{'h'} or $options{'help'});
	if(!exists $options{'input'}){
		print "***ERROR: No input directory is assigned!!!\n";
		exec ("pod2usage $0");
	}elsif(!exists $options{'site'}){
		print "***ERROR: No site position is assigned!!!\n";
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

    extract_sites_v1.pl version 1.0

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

    extract codon sites of PCGs, such as first/second/third codon site, first two codon sites, last two codon sites

=head1 SYNOPSIS

    extract_sites_v1.pl -i -s -o
    Copyright (C) 2025 Xiao-Jian Qu
    Please contact <quxiaojian@sdnu.edu.cn>, if you have any bugs or questions.

    [-h -help]   help information.
    [-i -input]  required: (default: input) input directory.
    [-s -site]   required: (default: 1) 1/2/3, extraction pattern of codon sites,
	             1=site1 + site2 + site3, 2=site12 + site3, 3=site1 + site23.
    [-o -output] required: (default: output) output directory.

=cut
