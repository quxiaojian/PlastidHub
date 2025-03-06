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


my $pattern1=".fasta";
my $pattern2=".fas";
my $pattern3=".fa";
my $pattern4=".fsa";
my $pattern5=".phylip";
my $pattern6=".phy";
my @filenames;
find(\&target,$input_directory);
sub target{
	if (/$pattern1/ or /$pattern2/ or /$pattern3/ or /$pattern4/ or /$pattern5/ or /$pattern6/){
		push @filenames,"$File::Find::name";
	}
	return;
}


my @array;
while (@filenames) {
	my $filename_fasta_phy=shift @filenames;
	my $target_name=substr($filename_fasta_phy,0,rindex($filename_fasta_phy,"\."));
	$target_name=~ s/\s/_/g;
	my $filename=substr($target_name,rindex($target_name,"\/")+1);
	my $output_filename_temp="$input_directory/$filename\_temp";
	my $genename=$filename;

	if ($filename_fasta_phy=~ /\.phy/) {
		open (my $input,"<",$filename_fasta_phy);
		while (<$input>) {
			$_=~ s/\r|\n//g;
			if ($_=~ /^\s(\d+)\s(\d+)/) {
				push @array,"$genename\t$2";
			}
		}
		close $input;
	}elsif (($filename_fasta_phy=~ /\.fa/) or ($filename_fasta_phy=~ /\.fsa/)) {
		open(my $input_temp,"<",$filename_fasta_phy);
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


		open (my $input,"<",$output_filename_temp);
		my $j=0;
		my ($header,$sequence,$length);
		while (defined ($header=<$input>) and defined ($sequence=<$input>)) {
			$header=~ s/\r|\n//g;
			$sequence=~ s/\r|\n//g;
			$length=length $sequence;
			push @array,"$genename\t$length";
			$j++;
			last if ($j > 0);
		}
		close $input;
		unlink("$output_filename_temp");
	}
}


my $sum=0;
open(my $output,">","$output_directory/PartitionFinder.txt");
foreach (@array) {
	my ($gene,$length)=split(/\s+/,$_);
	$sum=$sum+$length;
	my $sumlength=$sum+1-$length;
	print $output "$gene = $sumlength-$sum;\n";
}
close $output;




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
		print "***ERROR: No output filename is assigned!!!\n";
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

    generate_PartitionFinder_v1.pl version 1.0

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

    generate PartitionFinder format from alignment matrices

=head1 SYNOPSIS

    generate_PartitionFinder_v1.pl -i -o
    Copyright (C) 2025 Xiao-Jian Qu
    Please contact <quxiaojian@sdnu.edu.cn>, if you have any bugs or questions.

    [-h -help]         help information.
    [-i -input]        required: (default: input) input directory name containing fasta or phylip alignment matrices.
    [-o -output]       required: (default: output) output file name.

=cut
