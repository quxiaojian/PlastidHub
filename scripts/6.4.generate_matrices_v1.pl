#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Find;
use Data::Dumper;
$|=1;

my $global_options=&argument();
my $fastafilename=&default("input","input");
my $genefilename=&default("name","name");
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


my $genename;
my @names;
if ($genefilename ne "name") {
	open ($genename,"<",$genefilename);
	while (<$genename>) {
		$_=~ s/\r|\n//g;
		push @names,$_;
	}
	close $genename;
}


my $output_filename_temp="$fastafilename\_temp";
open(my $input_temp,"<",$fastafilename);
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


open (my $temp,"<",$output_filename_temp);
my %genes;
while (<$temp>) {
	$_=~ s/\r|\n//g;
	if ($_=~ />(.+?)_(.+)/) {
		$genes{$1}++;
	}
}
close $temp;


open (my $file,"<",$output_filename_temp);
my ($header,$sequence);
while (defined ($header=<$file>) and defined ($sequence=<$file>)) {
	$header=~ s/\r|\n//g;
	$sequence=~ s/\r|\n//g;
	if (@names != 0) {
		foreach my $name (@names) {#gene names in the genename.txt
			if ($header=~ /^>(.+?)_/) {
				if ($1 eq $name) {
					$header=~ s/>(.+?)_/>/g;
					open (my $output,">>","$output_directory/$name.fasta");
					print $output "$header\n$sequence\n";
					close $output;
				}
			}
			#if ($header=~ m/$name/) {
			#	open (my $output,">>","$output_directory/$name.fasta");
			#	print $output "$header\n$sequence\n";
			#	close $output;
			#}
		}
	}elsif (@names == 0) {#gene names in the .fasta file
		foreach my $gene (keys %genes) {
			if ($header=~ /^>(.+?)_/) {
				if ($1 eq $gene) {
					$header=~ s/>(.+?)_/>/g;
					open (my $output,">>","$output_directory/$gene.fasta");
					print $output "$header\n$sequence\n";
					close $output;
				}
			}
			#if ($header=~ m/$gene/) {
			#	open (my $output,">>","$output_directory/$gene.fasta");
			#	print $output "$header\n$sequence\n";
			#	close $output;
			#}
		}
	}
}
close $file;
unlink("$output_filename_temp");




##function##
sub argument{
	my @options=("help|h","input|i:s","name|n:s","output|o:s");
	my %options;
	GetOptions(\%options,@options);
	exec ("pod2usage $0") if ((keys %options)==0 or $options{'h'} or $options{'help'});
	if(!exists $options{'input'}){
		print "***ERROR: No input filename is assigned!!!\n";
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

    generate_matrices_v1.pl version 1.0

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

    generate sequence matrices from one fasta file

=head1 SYNOPSIS

    generate_matrices_v1.pl -i [-n] -o
    Copyright (C) 2025 Xiao-Jian Qu
    Please contact <quxiaojian@sdnu.edu.cn>, if you have any bugs or questions.

    [-h -help]         help information.
    [-i -input]        required: (default: input) input filename containing sequences from multiple species.
    [-n -name]         optional: (default: codingname.txt/noncodingname.txt) input filename containing sequence names or not assign this file.
    [-o -output]       required: (default: output) output directory.

=cut
