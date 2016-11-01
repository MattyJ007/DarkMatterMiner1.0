#!/usr/bin/perl -w

# script was written by Xuewen Wang
# gmat: Genome_wide MicroSatelllite Analysis Tool 
# originally designed for mining repeat motif like microsatllite or SSR 
# and do statistic analyzing on chromosme and whole genome level, especially faster in huge genome 
# 
# Usage: $0 -i [DNA sequence file] -o [motif searching results] -r [minimum repeated times] 
#        -m [minimum motif length] -x [maximum motif length] -s [highlight motif in sequence 0 or 1]
#run e.g. :
#          perl gmat.pl -r 5 -m 2 -x 10 -s 0 -i testseq.fasta
#result files started with input file name with suffix .fms, .ssr, .sat1,.sat2
&usageinfo();

# get command options
my %commandin = @ARGV;
if ((scalar @ARGV)%2 != 0){print "arguments must be in pair";}
my $refin=$commandin{"-i"}; #information of input sequence file
my $fmseqfile=$refin.".fms"; # formated sequence file
my $ssrfile=$refin.".ssr"; # default output file name if no output file name is given
my $ssrfiletmp=$refin.".ssrtmp";
my $motifmin=$commandin{"-m"}||2;# minimum length of any repeat unit, default value is 2
my $motifmax=$commandin{"-x"}||10; # maximum length of any repeat unit,, default value is 10
my $motif_times_min=$commandin{"-r"}||5; # minimum repeated times of a repeat unit, default value is 5
my $seqhight=$commandin{"-s"}||0; #highlight the repeated unit, also called motif, in each sequence which the motif is located


open (OUTLOG,">info.log")|| die (" LOG file writing failed, please check");
&runtime(OUTLOG);
print OUTLOG "In this run, parameters setting are: -r $motif_times_min -m $motifmin -x $motifmax -s $seqhight -i $refin\n";

#formating for huge genome sequence
print "Please wait during formating...\n\n"; #command line run information
system ("perl formatchunk.pl  -i $refin -f 1 -len 0 -o $fmseqfile");
my $fileinfo1= "Sequences were successfully formatted. 
file: $fmseqfile.
file: $refin.sat1\n\n";
print OUTLOG $fileinfo1;
print $fileinfo1;
#echo formated completed

&runtime(OUTLOG);
print "Mining microsatellite...\n";
system ("perl gssr.pl -r $motif_times_min -m $motifmin -x $motifmax -s $seqhight -i $fmseqfile -o $ssrfiletmp");
system ("perl gssrtrim.pl -i $ssrfiletmp -o $ssrfile");
my $fileinfo2="Microsatellite data were provided.
file:  $ssrfile\n\n";
print OUTLOG $fileinfo2;
print $fileinfo2;
#echo ssr analysis completed

print "Statistical Analysing...\n";
&runtime(OUTLOG);
system ("perl gsts.pl -i $ssrfile");
my $fileinfo3="statistic results were provided.
file:  $refin.sat2\n\n";
print OUTLOG $fileinfo3;
print $fileinfo3;
print "Done.\n";
unlink $ssrfiletmp;
&runtime(OUTLOG);

# remove temperary files
#unlink $fmseqfile;
 
sub usageinfo
	{# print program name and usage for help if no input aguments available in command line
	my @usage=(); # showing content on how to use the programme
	$usage[0]="Usage: searching the repeated motif such as SSR in a genome sequence file.\n";
	$usage[1]=" perl gmat.pl -r 5 -m 2 -x 10 -s 0 -i testseq.fasta \n";
	$usage[2]=" for    help: perl $0 ; \n";
	$usage[3]=" for running: perl $0 -i [formated sequence file]  -r [minimum repeated times] 
	-m [minimum motif length] -x [maximum motif length] -s [highlight motif in sequence 0 or 1] \n";
	$usage[4]="Author: Xuewen Wang\n";
	$usage[5]="year 2012\n";
	
	unless(@ARGV){print @usage; exit;} 
 }

sub runtime(){
	my $OUTfile=shift @_;
	my $local_time = gmtime();
	print {$OUTfile} "$0 was run. Current time is $local_time\n";
}

exit;# exit whole programe
