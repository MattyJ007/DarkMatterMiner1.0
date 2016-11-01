#!/usr/bin/perl -w
use strict;

# script was written by Xuewen Wang, originally designed for parsing SSR output from gssr.pl via
#        removing repated loci information 
#         and merge loci in overlap region 

# this version run very fast and less memory is required. no bug was found.
#Usage: $0 -i [motif information] -o [trimmed motif searching results] 
# print "Software $0 written by Xuewen Wang from 2012 Jan.  
#current version: 2012 Dec
# The function is to post parsing SSR loci information
#run e.g. :
#perl gssr.pl -i data\testseq.fasta.ssrtemp -o data\testseq.fasta.ssr



&usageinfo();

# get command options
my %commandin = @ARGV;
if ((scalar @ARGV)%2 != 0){print "arguments must be in pair";}
my $refin=$commandin{"-i"}; #information of input sequence file
my $outfiledefault=$refin.".ssr"; # default output file name if no output file name is given
my $outputfile=$commandin{"-o"}||$outfiledefault; #information of input sequence file

#file status check and prepare
 unless(-e $refin){print "input file does not exist in the working directoty, or the path is not given\n";}
 open (DNAhand, "<$refin")|| die (" $0 filed open failed, pleasse check");
 open(OUT, ">$outputfile.tmp")|| die (" $0 Write file failed, pleasse check");
 
 print OUT "Name\tSeq_Len\tStartPos\tEndPos\tRepetitions\tMotif\n"; #prepare head
 
 my $myvaluein="";
 my $mykeyid="";
 my %allseqinformation=();
 my %allseqinfor=();

while ( <DNAhand>){	# one line at at time
	 chomp ;
	 if ($_=~/^>.*/){
			(my $Name,my $Seq_chunk,my $StartPos,my $EndPos,my $Repetitions,my $tMotif,my $Seq_len)=split("\t",);
			 $myvaluein = join("\t", $Name,$Seq_len,$StartPos,$EndPos,$Repetitions,$tMotif);
			 $mykeyid = join("\$\$", $Name,$StartPos);
			 
			 if(! exists $allseqinformation{$mykeyid}){
				  $allseqinformation{$mykeyid}=$EndPos; 
			 }elsif($allseqinformation{$mykeyid}<$EndPos){ 
			 # to find far end in the right side of overlap region
			      $allseqinformation{$mykeyid}=$EndPos;
			 }
			 			 
			 # hashing all information based on $mykeyid
                 if(! exists $allseqinfor{$mykeyid}){
				     $allseqinfor{$mykeyid}=$myvaluein; 
				 }
	 } #end if
	 
} #ending while readin

#output temporary results: 
     foreach my $keyout(sort keys %allseqinformation){
	     print OUT $allseqinfor{$keyout}, "\n";	
     }
	 %allseqinformation=(); %allseqinfor=();
close DNAhand;
close OUT;



#################################################
# resovle the left far end of the overlaped region
     open (DNAhand, "<$outputfile.tmp")|| die (" $0 filed open failed, pleasse check");
     open(OUT, ">$outputfile")|| die (" $0 Write file failed, pleasse check"); 
     print OUT "Name\tSeq_Len\tStartPos\tEndPos\tRepetitions\tMotif\n"; #prepare head

     while ( <DNAhand>){	# one line at at time
		chomp ;
		if ($_=~/^>.*/){
			(my $Name,my $Seq_len,my $StartPos,my $EndPos,my $Repetitions,my $tMotif)=split("\t",);
			 $myvaluein = $_;
			 $mykeyid = join("\$\$", $Name,$EndPos);			 
			 
			 if(! exists $allseqinformation{$mykeyid}){
				  $allseqinformation{$mykeyid}=$StartPos; 
			 }elsif($allseqinformation{$mykeyid}>$StartPos){ 
			 # to find far end in the right side of overlap region
			      $allseqinformation{$mykeyid}=$StartPos;
			 }
			 			 
			 # hashing all information based on $mykeyid
                 if(! exists $allseqinfor{$mykeyid}){
				     $allseqinfor{$mykeyid}=$myvaluein; 
				 }
	 } #end if
	 
} #ending while readin

#output temporary results: 
     foreach my $keyout(sort keys %allseqinformation){
	     print OUT $allseqinfor{$keyout}, "\n";	
     }

close DNAhand;
close OUT;
unlink "$outputfile.tmp";

sub usageinfo
	{# print program name and usage for help if no input aguments available in command line
	my @usage=(); # showing content on how to use the programme
	$usage[0]="Usage: posr parsing SSR loci file output from gssr.pl\n";
	$usage[1]=" for    help: perl $0 \n";
	$usage[2]=" for running: perl $0 -i [ssr output file for gssr.pl] -o [final SSR file name]  \n";
	$usage[3]=" e.g. perl gssr.pl -i data\testseq.fasta.ssrtemp -o data\testseq.fasta.ssr\n";
	$usage[4]="Author: Xuewen Wang\n";
	$usage[5]="year 2012\n";
	unless(@ARGV){print @usage; exit;} 
 }


exit;# exit whole programe
