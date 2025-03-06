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
	open(my $output,">","$output_directory/$filename.fasta");
	my ($header, $sequence);
	while (defined ($header=<$input>) && defined ($sequence=<$input>)) {
		$header=~ s/\r|\n//g;
		$sequence=~ s/\r|\n//g;
		my $rev_com=reverse $sequence;
		$rev_com=~ tr/ACGTacgt/TGCAtgca/;
		my $count_dash = $header =~ tr/-/-/;
		my $temp;
		if ($count_dash == 5) {#5>trnK-UUU-1-trnK-UUU-1_
			my ($part1, $part2, $part3, $part4, $part5, $part6, $part7, $part8, $part9);
			if ($header =~ /^>(trn(f?[A-Z])-([A-Z]{3})-(\d))-(trn(f?[A-Z])-([A-Z]{3})-(\d))_(.+)/) {
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$part4 = $4;
				$part5 = $5;
				$part6 = $6;
				$part7 = $7;
				$part8 = $8;
				$part9 = $9;
				if ($part2 lt $part6) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part2 gt $part6) {
					$temp = ">".$part5."-".$part1."_".$part9."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part2 eq $part6) {
					if ($part3 lt $part7) {
						print $output "$header\n";
						print $output "$sequence\n";
					}elsif ($part3 gt $part7) {
						$temp = ">".$part5."-".$part1."_".$part9."\n";
						print $output $temp;
						print $output "$rev_com\n";
					}elsif ($part3 eq $part7) {
						if ($part4 <= $part8) {
							print $output "$header\n";
							print $output "$sequence\n";
						}elsif ($part4 > $part8) {
							$temp = ">".$part5."-".$part1."_".$part9."\n";
							print $output $temp;
							print $output "$rev_com\n";
						}
					}
				}
			}
		}elsif ($count_dash == 4) {
			my ($part1, $part2, $part3, $part4, $part5, $part6, $part7, $part8, $part9);
			if ($header =~ /^>(trn(f?[A-Z])-([A-Z]{3}))-(trn(f?[A-Z])-([A-Z]{3})-(\d))_(.+)/) {#4>trnK-UUU-trnK-UUU-1_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$part4 = $4;
				$part5 = $5;
				$part6 = $6;
				$part7 = $7;
				$part8 = $8;
				if ($part2 lt $part5) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part2 gt $part5) {
					$temp = ">".$part4."-".$part1."_".$part8."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part2 eq $part5) {
					if ($part3 lt $part6) {
						print $output "$header\n";
						print $output "$sequence\n";
					}elsif ($part3 gt $part6) {
						$temp = ">".$part4."-".$part1."_".$part8."\n";
						print $output $temp;
						print $output "$rev_com\n";
					}elsif ($part3 eq $part6) {
						print $output "$header\n";
						print $output "$sequence\n";
					}
				}
			}elsif ($header =~ /^>(trn(f?[A-Z])-([A-Z]{3})-(\d))-(trn(f?[A-Z])-([A-Z]{3}))_(.+)/) {#4>trnK-UUU-1-trnK-UUU_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$part4 = $4;
				$part5 = $5;
				$part6 = $6;
				$part7 = $7;
				$part8 = $8;
				if ($part2 lt $part6) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part2 gt $part6) {
					$temp = ">".$part5."-".$part1."_".$part8."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part2 eq $part6) {
					if ($part3 lt $part7) {
						print $output "$header\n";
						print $output "$sequence\n";
					}elsif ($part3 gt $part7) {
						$temp = ">".$part5."-".$part1."_".$part8."\n";
						print $output $temp;
						print $output "$rev_com\n";
					}elsif ($part3 eq $part7) {
						$temp = ">".$part5."-".$part1."_".$part8."\n";
						print $output $temp;
						print $output "$rev_com\n";
					}
				}
			}elsif ($header =~ /^>((t)rnf?[A-Z]-([A-Z]{3})-(\d))-(([a-su-z])[a-zA-Z0-9]+-(\d))_(.+)/) {#4>trnK-UUU-1-matK(rpl32/rrn16)-1_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$part4 = $4;
				$part5 = $5;
				$part6 = $6;
				$part7 = $7;
				$part8 = $8;
				if ($part2 lt $part6) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part2 gt $part6) {
					$temp = ">".$part5."-".$part1."_".$part8."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part2 eq $part6) {
					print "ERROR!";
				}
			}elsif ($header =~ /^>(([a-su-z])[a-zA-Z0-9]+-(\d))-((t)rnf?[A-Z]-([A-Z]{3})-(\d))_(.+)/) {#4>matK(rpl32/rrn16)-1-trnK-UUU-1_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$part4 = $4;
				$part5 = $5;
				$part6 = $6;
				$part7 = $7;
				$part8 = $8;
				if ($part2 lt $part5) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part2 gt $part5) {
					$temp = ">".$part4."-".$part1."_".$part8."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part2 eq $part5) {
					print "ERROR!";
				}
			}
		}elsif ($count_dash == 3) {
			my ($part1, $part2, $part3, $part4, $part5, $part6, $part7, $part8, $part9);
			if ($header =~ /^>(trn(f?[A-Z])-([A-Z]{3}))-(trn(f?[A-Z])-([A-Z]{3}))_(.+)/) {#3>trnK-UUU-trnK-UUU_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$part4 = $4;
				$part5 = $5;
				$part6 = $6;
				$part7 = $7;
				if ($part2 lt $part5) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part2 gt $part5) {
					$temp = ">".$part4."-".$part1."_".$part7."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part2 eq $part5) {
					if ($part3 lt $part6) {
						print $output "$header\n";
						print $output "$sequence\n";
					}elsif ($part3 gt $part6) {
						$temp = ">".$part4."-".$part1."_".$part7."\n";
						print $output $temp;
						print $output "$rev_com\n";
					}elsif ($part3 eq $part6) {
						print $output "$header\n";
						print $output "$sequence\n";
					}
				}
			}elsif ($header =~ /^>((t)rnf?[A-Z]-([A-Z]{3})-(\d))-(([a-su-z])[a-zA-Z0-9]+)_(.+)/) {#3>trnK-UUU-1-matK(rpl32/rrn16)_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$part4 = $4;
				$part5 = $5;
				$part6 = $6;
				$part7 = $7;
				if ($part2 lt $part6) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part2 gt $part6) {
					$temp = ">".$part5."-".$part1."_".$part7."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part2 eq $part6) {
					print "ERROR!";
				}
			}elsif ($header =~ /^>(([a-su-z])[a-zA-Z0-9]+)-((t)rnf?[A-Z]-([A-Z]{3})-(\d))_(.+)/) {#3>matK(rpl32/rrn16)-trnK-UUU-1_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$part4 = $4;
				$part5 = $5;
				$part6 = $6;
				$part7 = $7;
				if ($part2 lt $part4) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part2 gt $part4) {
					$temp = ">".$part3."-".$part1."_".$part7."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part2 eq $part4) {
					print "ERROR!";
				}
			}elsif ($header =~ /^>((t)rnf?[A-Z]-([A-Z]{3}))-(([a-su-z])[a-zA-Z0-9]+-(\d))_(.+)/) {#3>trnK-UUU-matK(rpl32/rrn16)-1_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$part4 = $4;
				$part5 = $5;
				$part6 = $6;
				$part7 = $7;
				if ($part2 lt $part5) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part2 gt $part5) {
					$temp = ">".$part4."-".$part1."_".$part7."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part2 eq $part5) {
					print "ERROR!";
				}
			}elsif ($header =~ /^>(([a-su-z])[a-zA-Z0-9]+-(\d))-((t)rnf?[A-Z]-([A-Z]{3}))_(.+)/) {#3>matK(rpl32/rrn16)-1-trnK-UUU_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$part4 = $4;
				$part5 = $5;
				$part6 = $6;
				$part7 = $7;
				if ($part2 lt $part5) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part2 gt $part5) {
					$temp = ">".$part4."-".$part1."_".$part7."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part2 eq $part5) {
					print "ERROR!";
				}
			}elsif ($header =~ /^>(([a-su-z][a-zA-Z0-9]+)-(\d+))-(([a-su-z][a-zA-Z0-9]+)-(\d))_(.+)/) {#3>matK(rpl32/rrn16)-1-matK(rpl32/rrn16)-1_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$part4 = $4;
				$part5 = $5;
				$part6 = $6;
				$part7 = $7;
				if ($part2 lt $part5) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part2 gt $part5) {
					$temp = ">".$part4."-".$part1."_".$part7."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part2 eq $part5) {
					if ($part3 <= $part6) {
						print $output "$header\n";
						print $output "$sequence\n";
					}elsif ($part3 > $part6) {
						$temp = ">".$part4."-".$part1."_".$part7."\n";
						print $output $temp;
						print $output "$rev_com\n";
					}
				}
			}elsif ($header =~ /^>(trn[A-Z]-[A-Z]{3}-\d)-(rrn4.5)_(.+)/) {#3>trnK-UUU-1-rrn4.5_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$temp = ">".$part2."-".$part1."_".$part3."\n";
				print $output $temp;
				print $output "$rev_com\n";
			}elsif ($header =~ /^>(rrn4.5)-(trn[A-Z]-[A-Z]{3}-\d)_(.+)/) {#3>rrn4.5-trnK-UUU-1_
				print $output $header;
				print $output "$sequence\n";
			}
		}elsif ($count_dash == 2) {
			my ($part1, $part2, $part3, $part4, $part5, $part6, $part7, $part8, $part9);
			if ($header =~ /^>((t)rnf?[A-Z]-([A-Z]{3}))-(([a-su-z])[a-zA-Z0-9]+)_(.+)/) {#2>trnK-UUU-matK(rpl32/rrn16)_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$part4 = $4;
				$part5 = $5;
				$part6 = $6;
				if ($part2 lt $part5) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part2 gt $part5) {
					$temp = ">".$part4."-".$part1."_".$part6."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part2 eq $part5) {
					print "ERROR!";
				}
			}elsif ($header =~ /^>(([a-su-z])[a-zA-Z0-9]+)-((t)rnf?[A-Z]-([A-Z]{3}))_(.+)/) {#2>matK(rpl32/rrn16)-trnK-UUU_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$part4 = $4;
				$part5 = $5;
				$part6 = $6;
				if ($part2 lt $part4) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part2 gt $part4) {
					$temp = ">".$part3."-".$part1."_".$part6."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part2 eq $part4) {
					print "ERROR!";
				}
			}elsif ($header =~ /^>(([a-su-z][a-zA-Z]+)-(\d))-([a-su-z][a-zA-Z]+)_(.+)/) {#2>matK-1-matK_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$part4 = $4;
				$part5 = $5;
				if ($part2 lt $part4) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part2 gt $part4) {
					$temp = ">".$part4."-".$part1."_".$part5."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part2 eq $part4) {
					$temp = ">".$part4."-".$part1."_".$part5."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}
			}elsif ($header =~ /^>(([a-su-z][a-zA-Z]+)-(\d))-([a-su-z][a-zA-Z]+[0-9]+)_(.+)/) {#2>matK-1-rpl32(rrn16)_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$part4 = $4;
				$part5 = $5;
				if ($part2 lt $part4) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part2 gt $part4) {
					$temp = ">".$part4."-".$part1."_".$part5."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part2 eq $part4) {
					print "ERROR!";
				}
			}elsif ($header =~ /^>(([a-su-z][a-zA-Z]+[0-9]+)-(\d))-([a-su-z][a-zA-Z]+)_(.+)/) {#2>rpl32(rrn16)-1-matK_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$part4 = $4;
				$part5 = $5;
				if ($part2 lt $part4) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part2 gt $part4) {
					$temp = ">".$part4."-".$part1."_".$part5."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part2 eq $part4) {
					print "ERROR!";
				}
			}elsif ($header =~ /^>((([a-su-z][a-zA-Z]+)([0-9]+))-(\d))-(([a-su-z][a-zA-Z]+)([0-9]+))_(.+)/) {#2>rpl32(rrn16)-1-rpl32(rrn16)_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$part4 = $4;
				$part5 = $5;
				$part6 = $6;
				$part7 = $7;
				$part8 = $8;
				$part9 = $9;
				if ($part3 lt $part7) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part3 gt $part7) {
					$temp = ">".$part6."-".$part1."_".$part9."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part3 eq $part7) {
					if ($part4 < $part8) {
						print $output "$header\n";
						print $output "$sequence\n";
					}elsif ($part4 >= $part8) {
						$temp = ">".$part6."-".$part1."_".$part9."\n";
						print $output $temp;
						print $output "$rev_com\n";
					}
				}
			}elsif ($header =~ /^>([a-su-z][a-zA-Z]+)-(([a-su-z][a-zA-Z]+)-(\d))_(.+)/) {#2>matK-matK-1_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$part4 = $4;
				$part5 = $5;
				if ($part1 lt $part3) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part1 gt $part3) {
					$temp = ">".$part2."-".$part1."_".$part5."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part1 eq $part3) {
					print $output "$header\n";
					print $output "$sequence\n";
				}
			}elsif ($header =~ /^>([a-su-z][a-zA-Z]+[0-9]+)-(([a-su-z][a-zA-Z]+)-(\d))_(.+)/) {#2>rpl32(rrn16)-matK-1_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$part4 = $4;
				$part5 = $5;
				if ($part1 lt $part3) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part1 gt $part3) {
					$temp = ">".$part2."-".$part1."_".$part5."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part1 eq $part3) {
					print "ERROR!";
				}
			}elsif ($header =~ /^>([a-su-z][a-zA-Z]+)-(([a-su-z][a-zA-Z]+[0-9]+)-(\d))_(.+)/) {#2>matK-rpl32(rrn16)-1_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$part4 = $4;
				$part5 = $5;
				if ($part1 lt $part3) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part1 gt $part3) {
					$temp = ">".$part2."-".$part1."_".$part5."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part1 eq $part3) {
					print "ERROR!";
				}
			}elsif ($header =~ /^>(([a-su-z][a-zA-Z]+)([0-9]+))-((([a-su-z][a-zA-Z]+)([0-9]+))-(\d))_(.+)/) {#2>rpl32(rrn16)-rpl32-1_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$part4 = $4;
				$part5 = $5;
				$part6 = $6;
				$part7 = $7;
				$part8 = $8;
				$part9 = $9;
				if ($part2 lt $part6) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part2 gt $part6) {
					$temp = ">".$part4."-".$part1."_".$part9."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part2 eq $part6) {
					if ($part3 <= $part7) {
						print $output "$header\n";
						print $output "$sequence\n";
					}elsif ($part3 > $part7) {
						$temp = ">".$part4."-".$part1."_".$part9."\n";
						print $output $temp;
						print $output "$rev_com\n";
					}
				}
			}elsif ($header =~ /^>(trn[A-Z]-[A-Z]{3})-(rrn4.5)_(.+)/) {#2>trnK-UUU-rrn4.5_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$temp = ">".$part2."-".$part1."_".$part3."\n";
				print $output $temp;
				print $output "$rev_com\n";
			}elsif ($header =~ /^>(rrn4.5)-(trn[A-Z]-[A-Z]{3})_(.+)/) {#2>rrn4.5-trnK-UUU_
				print $output $header;
				print $output "$sequence\n";
			}elsif ($header =~ /^>(([a-su-z][a-zA-Z]+)-(\d))-(rrn4.5)_(.+)/) {#2>matK-1-rrn4.5_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$part4 = $4;
				$part5 = $5;
				if ($part2 lt $part4) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part2 gt $part4) {
					$temp = ">".$part4."-".$part1."_".$part5."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part2 eq $part4) {
					print "ERROR!";
				}
			}elsif ($header =~ /^>(rrn4.5)-(([a-su-z][a-zA-Z]+)-(\d))_(.+)/) {#2>rrn4.5-matK-1_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$part4 = $4;
				$part5 = $5;
				if ($part1 lt $part3) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part1 gt $part3) {
					$temp = ">".$part2."-".$part1."_".$part5."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part1 eq $part3) {
					print "ERROR!";
				}
			}elsif ($header =~ /^>((([a-su-z][a-zA-Z]+)[0-9]+)-(\d))-((rrn)(4.5))_(.+)/) {#2>rpl32(rrn16)-1-rrn4.5_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$part4 = $4;
				$part5 = $5;
				$part6 = $6;
				$part7 = $7;
				$part8 = $8;
				if ($part3 lt $part6) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part3 gt $part6) {
					$temp = ">".$part5."-".$part1."_".$part8."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part3 eq $part6) {
					if ($part4 < $part7) {
						print $output "$header\n";
						print $output "$sequence\n";
					}elsif ($part4 >= $part7) {
						$temp = ">".$part5."-".$part1."_".$part8."\n";
						print $output $temp;
						print $output "$rev_com\n";
					}
				}
			}elsif ($header =~ /^>((rrn)(4.5))-((([a-su-z][a-zA-Z]+)[0-9]+)-(\d))_(.+)/) {#2>rrn4.5-rpl32(rrn16)-1_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$part4 = $4;
				$part5 = $5;
				$part6 = $6;
				$part7 = $7;
				$part8 = $8;
				if ($part2 lt $part6) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part2 gt $part6) {
					$temp = ">".$part4."-".$part1."_".$part8."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part2 eq $part6) {
					if ($part3 <= $part7) {
						print $output "$header\n";
						print $output "$sequence\n";
					}elsif ($part3 > $part7) {
						$temp = ">".$part4."-".$part1."_".$part8."\n";
						print $output $temp;
						print $output "$rev_com\n";
					}
				}
			}
		}elsif ($count_dash == 1) {
			my ($part1, $part2, $part3, $part4, $part5, $part6, $part7, $part8, $part9);
			if ($header =~ /^>([a-su-z][a-zA-Z0-9]+)-([a-su-z][a-zA-Z0-9]+)_(.+)/) {#1>matK(rpl32/rrn16)-matK(rpl32/rrn16)_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				if ($part1 lt $part2) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part1 gt $part2) {
					$temp = ">".$part2."-".$part1."_".$part3."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part1 eq $part2) {
					print $output "$header\n";
					print $output "$sequence\n";
				}
			}elsif ($header =~ /^>([a-su-z][a-zA-Z]+)-((rrn)4.5)_(.+)/) {#1>matK-rrn4.5_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$part4 = $4;
				if ($part1 lt $part3) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part1 gt $part3) {
					$temp = ">".$part2."-".$part1."_".$part4."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part1 eq $part3) {
					print "ERROR!";
				}
			}elsif ($header =~ /^>((rrn)4.5)-([a-su-z][a-zA-Z]+)_(.+)/) {#1>rrn4.5-matK_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$part4 = $4;
				if ($part2 lt $part3) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part2 gt $part3) {
					$temp = ">".$part3."-".$part1."_".$part4."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part2 eq $part3) {
					print "ERROR!";
				}
			}elsif ($header =~ /^>(([a-su-z][a-zA-Z]+)([0-9]+))-((rrn)(4.5))_(.+)/) {#1>rpl32(rrn16)-rrn4.5_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$part4 = $4;
				$part5 = $5;
				$part6 = $6;
				$part7 = $7;
				if ($part2 lt $part5) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part2 gt $part5) {
					$temp = ">".$part4."-".$part1."_".$part7."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part2 eq $part5) {
					if ($part3 < $part6) {
						print $output "$header\n";
						print $output "$sequence\n";
					}elsif ($part3 > $part6) {
						$temp = ">".$part4."-".$part1."_".$part7."\n";
						print $output $temp;
						print $output "$rev_com\n";
					}elsif ($part3 == $part6) {
						print "ERROR!";
					}
				}
			}elsif ($header =~ /^>((rrn)(4.5))-(([a-su-z][a-zA-Z]+)([0-9]+))_(.+)/) {#1>rrn4.5-rpl32(rrn16)_
				$part1 = $1;
				$part2 = $2;
				$part3 = $3;
				$part4 = $4;
				$part5 = $5;
				$part6 = $6;
				$part7 = $7;
				if ($part2 lt $part5) {
					print $output "$header\n";
					print $output "$sequence\n";
				}elsif ($part2 gt $part5) {
					$temp = ">".$part4."-".$part1."_".$part7."\n";
					print $output $temp;
					print $output "$rev_com\n";
				}elsif ($part2 eq $part5) {
					if ($part3 < $part6) {
						print $output "$header\n";
						print $output "$sequence\n";
					}elsif ($part3 > $part6) {
						$temp = ">".$part4."-".$part1."_".$part7."\n";
						print $output $temp;
						print $output "$rev_com\n";
					}elsif ($part3 == $part6) {
						print "ERROR!";
					}
				}
			}elsif ($header =~ /^>(rrn4.5)-(rrn4.5)_(.+)/) {#1>rrn4.5-rrn4.5_
				print $output "$header\n";
			}
		}
	}
	close $input;
	close $output;
	unlink("$output_filename_temp");
}




####
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

    sort_genenames_v1.pl version 1.0

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

    sort flanking genenames and sequences for intergenic regions

=head1 SYNOPSIS

    sort_genenames_v1.pl -i -o
    Copyright (C) 2025 Xiao-Jian Qu
    Please contact <quxiaojian@sdnu.edu.cn>, if you have any bugs or questions.

    [-h -help]         help information.
    [-i -input]        required: (default: input) input directory containing fasta files.
    [-o -output]       required: (default: output) output directory.

=cut
