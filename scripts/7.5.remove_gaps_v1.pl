#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Find;
use Data::Dumper;
$|=1;
#array of hash plus transposed matrix(two-dimensional array)

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
	open(my $output1,">","$output_directory/$filename.fasta");
	open(my $output2,">","$output_directory/gapsites_$filename.txt");

	my ($header,$sequence,$length,@id,@array);
	while (defined ($header=<$input>) && defined ($sequence=<$input>)) {
		$header=~ s/\r|\n//g;
		$sequence=~ s/\r|\n//g;
		$sequence=uc $sequence;
		$length=length $sequence;
		push @id,$header;

		for (0..$length-1){
			push @{$array[$_]},substr($sequence,$_,1);
		}
	}
	#print Dumper \@array;
	close $input;


	my @position;
	my @nongap_array;#nongap sites,first array_ref is fisrt column
	for my $column (0..$#array){#1 of 36 columns
		for my $nucleotide ($array[$column]){#1 of 4 nt in each column
			if (((grep{$_ eq "A"} @$nucleotide) or (grep{$_ eq "T"} @$nucleotide) or (grep{$_ eq "G"} @$nucleotide) or (grep{$_ eq "C"} @$nucleotide)) or ((grep{$_ eq "S"} @$nucleotide) or (grep{$_ eq "F"} @$nucleotide) or (grep{$_ eq "L"} @$nucleotide) or (grep{$_ eq "Y"} @$nucleotide) or (grep{$_ eq "C"} @$nucleotide) or (grep{$_ eq "W"} @$nucleotide) or (grep{$_ eq "P"} @$nucleotide) or (grep{$_ eq "H"} @$nucleotide) or (grep{$_ eq "Q"} @$nucleotide) or (grep{$_ eq "I"} @$nucleotide) or (grep{$_ eq "M"} @$nucleotide) or (grep{$_ eq "T"} @$nucleotide) or (grep{$_ eq "K"} @$nucleotide) or (grep{$_ eq "R"} @$nucleotide) or (grep{$_ eq "V"} @$nucleotide) or (grep{$_ eq "A"} @$nucleotide) or (grep{$_ eq "D"} @$nucleotide) or (grep{$_ eq "E"} @$nucleotide) or (grep{$_ eq "G"} @$nucleotide) or (grep{$_ eq "N"} @$nucleotide))) {
				push @nongap_array,[@$nucleotide];
				delete $array[$column];#gap sites, undefined
			}else {
				push @position,$column;
			}
		}
	}
	#print Dumper \@nongap_array;
	@array=grep {defined $_} @array;#gap sites,remove undefined elements


	#nongap sites,first array_ref is first row
	my @transposed_nongap_array=map {my $x=$_;[map {$nongap_array[$_][$x]} (0..$#nongap_array)]} (0..$#{$nongap_array[0]});
	for (0..(@id-1)){
		print $output1 $id[$_],"\n",@{$transposed_nongap_array[$_]},"\n";
	}
	close $output1;


	#gap position in alignment matrix
	print $output2 "Gap position(s) in alignment matrix:\n";
	if (@position != 0) {
		foreach (@position) {
			my $position=$_+1;
			print $output2 $position." ";
		}
		print $output2 "\n\n";
		print $output2 "Gap site(s) in alignment matrix (i.e. -):\n";
	}

	#gap sites,first array_ref is first row
	if (@array != 0) {
		my @transposed_gap_array=map {my $x=$_;[map {$array[$_][$x]} (0..$#array)]} (0..$#{$array[0]});
		for (0..(@id-1)){
			print $output2 $id[$_],"\n",@{$transposed_gap_array[$_]},"\n";
		}
		close $output2;
	}else {
		close $output2;
		unlink("$output_directory/gapsites_$filename.txt");
	}
	unlink("$output_filename_temp");
}




##function
sub argument{
	my @options=("help|h","input|i:s","output|o:s");
	my %options;
	GetOptions(\%options,@options);
	exec ("pod2usage $0") if ((keys %options)==0 or $options{'h'} or $options{'help'});
	if(!exists $options{'input'}){
		print "***ERROR: No input file is assigned!!!\n";
		exec ("pod2usage $0");
	}elsif(!exists $options{'output'}){
		print "***ERROR: No output file is assigned!!!\n";
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

    remove_gap_sites_from_alignment_matrices_v1.pl version 1.0

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

    remove gap sites from alignment matrices

=head1 SYNOPSIS

    remove_gap_sites_from_alignment_matrices_v1.pl -i -o
    Copyright (C) 2025 Xiao-Jian Qu
    Please contact <quxiaojian@sdnu.edu.cn>, if you have any bugs or questions.

    [-h -help]         help information.
    [-i -input]        required: (default: input) input directory containing .fasta format alignment matrices.
    [-o -output]       required: (default: output) output directory.

=cut
