#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Find;
use FindBin qw($Bin);
use Data::Dumper;
$|=1;

my $global_options=&argument();
my $input_directory=&default("input","input");
my $runs=&default("100","runs");
my $model=&default("GTRGAMMA","model");
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


while (@filenames) {
	my $filename_fa=shift @filenames;
	my $filename_fasta=$filename_fa;
	$filename_fasta=~ s/\s/_/g;
	rename($filename_fa,$filename_fasta);
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


	open(my $fasta,"<",$output_filename_temp);
	my ($header,$sequence,@species,$length,%species);
	while (defined ($header=<$fasta>) and defined ($sequence=<$fasta>)) {
		$header=~ s/\r|\n//g;
		$sequence=~ s/\r|\n//g;
		$length=length $sequence;
		if ($header=~ />(.+)/) {
			push @species,$1;
			$species{$1}++;
		}
	}
	close $fasta;
	unlink("$output_filename_temp");

	my $warning="$output_directory/warning.txt";
	foreach my $key (keys %species) {
		if ($species{$key} > 1) {
			open(my $output_warning,">>",$warning);
			print $output_warning "Duplicated species name $key for the gene $filename!\n";
			close $output_warning;
		}
	}


	my $species=@species;
	if (($species <= 50) and ($length <= 10000)) {
		if (($runs == 100) or ($runs == 200) or ($runs == 500) or ($runs == 1000)) {
			system("$Bin/scripts/raxmlHPC.exe -f a -x 12345 -p 12345 -k -N $runs -m $model -s $Bin/$filename_fasta -n $filename.result -w $Bin/$output_directory");
			#$Bin/scripts/raxmlHPC-PTHREADS-AVX.exe -f a -x 12345 -p 12345 -k -T 10 -# 100 -m GTRGAMMA -s input/atpA.fasta -n atpA.result -w output
		}
		#if ("$filename_fasta.reduced") {
		#	unlink("$filename_fasta.reduced");
		#}
	}elsif ($species > 50) {
		open(my $output_warning,">>",$warning);
		print $output_warning "Species number larger than 50 for the gene $filename!\n";
		close $output_warning;
	}elsif ($length > 10000) {
		open(my $output_warning,">>",$warning);
		print $output_warning "Sequence length longer than 10 kb for the gene $filename!\n";
		close $output_warning;
	}
}




##function##
sub argument{
	my @options=("help|h","input|i:s","runs|r:i","model|m:s","output|o:s");
	my %options;
	GetOptions(\%options,@options);
	exec ("pod2usage $0") if ((keys %options)==0 or $options{'h'} or $options{'help'});
	if(!exists $options{'input'}){
		print "***ERROR: No input directory is assigned!!!\n";
		exec ("pod2usage $0");
	}elsif(!exists $options{'runs'}){
		print "***ERROR: No number of runs is assigned!!!\n";
		exec ("pod2usage $0");
	}elsif(!exists $options{'model'}){
		print "***ERROR: No evolutionary model is assigned!!!\n";
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

    batch_RAxML_v1.pl version 1.0

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

    batch RAxML phylogenetic reconstruction

=head1 SYNOPSIS

    batch_RAxML_v1.pl -i -r -m -o
    Copyright (C) 2025 Xiao-Jian Qu
    Please contact <quxiaojian@sdnu.edu.cn>, if you have any bugs or questions.

    [-h -help]   help information.
    [-i -input]  required: (default: input) input directory containing .fasta sequence matrices.
    [-r -runs]   required: (default: 100) number of runs, i.e., 100/200/500/1000.
    [-m -model]  required: (default: GTRGAMMA) evolutionary models, i.e.,
	             GTRGAMMA/GTRGAMMAI/GTRCAT/GTRCATI, PROTGAMMAGTR/PROTCATGTR.
    [-o -output] required: (default: output) output directory.

=cut
