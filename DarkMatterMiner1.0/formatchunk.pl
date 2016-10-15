#!/usr/bin/perl -w
use strict;
&usageinfo();
# script written by Xuewen Wang, first version of this script was written in 2009, 
#chunking sequence to the expected size i.e 5000 000 each line for genome sequence
#tested by Xuewen Wang for several years in handreds of cases and working very nicely
#usage: perl formatchunk.pl -i [input DNA sequence file name]  -o [output file name]
# usage example: perl formatchunk.pl -i ssg_groupSeqfasta.fasta  -o formated.ssg_groupSeqfasta.fasta
#programmed by Xuewen Wang, previous version, 2009, 2010 Mar, 2011 Feb,2012 Aug; 
# Current version:  2012 Dec;

# get command options
my %commandin = @ARGV;
if ((scalar @ARGV)%2 != 0){print "arguments must in pair";}
my $refin=$commandin{"-i"}; 
#input sequence file

my $outputfile=$commandin{"-o"}; 
#output sequence file

#file check and prepare
 unless(-e $refin){print "Data file containg data for extracting does not 
 exist in the working directoty\n";}
 open(OUT, ">$outputfile")|| die (" formated file writing failed, pleasse check");
 my $satfile=$outputfile.".sat1";
 open(OUTsat, ">$satfile")|| die (" sat1 file writing failed, pleasse check");


#my $filename = 'ssg_groupSeqfasta.fasta'; 
 my $filename=$refin;
open(FILE, $filename) || die("Couldn't read file $filename\n");  

local $/ = "\n>";  # read by FASTA tag >
my $seqformat="false";
my $seqno=0;
my $seqlen=0;
my $eachseqlen=0;
my $startpos =0;
my $chunkvalue=3000000;
my $chunksize=$chunkvalue;
my $overlap=0;
my $ii=0;
my $seq="";
my $testlength=0;

while (my $seqinfo = <FILE>) {
chomp $seqinfo;
     
    $seqinfo =~ s/^>*(.+)\n//;  #  match and remove FASTA head info
     my $seqID=">".$1; #   head line
	 $seqno ++ if($1);
	 $seqformat="TRUE" if($1);
	 $seqinfo=~ s/[0-9\n\s]//g;  # remove newlines, numbers and white space
	 $eachseqlen=length $seqinfo;
	 $seqlen=$seqlen+$eachseqlen;
	 #print $seqID, "\n",$seqinfo,"\n";
	 print OUT  $seqID,"\n";	
	 #print "seqLen", $eachseqlen,"chunksize:$chunksize\n";
     
	 
	 # to get the chunks
	    if( $eachseqlen <= $chunksize){
			print OUT  $seqinfo, "\n";
			#print length $seqinfo, "shortseq\n";
		}else{		#($eachseqlen > $chunksize)
			while($startpos <= $eachseqlen-1){              
				$seq = substr($seqinfo, $startpos, $chunksize);
				$ii =$ii+1;				
			    print OUT  $seq, "\n";
				$startpos = $startpos+$chunksize-$overlap;
				if ($startpos+$chunksize >= $eachseqlen-1){$chunksize=$eachseqlen-$startpos+1;} 
				# final end part seq
				$testlength=$testlength+length $seq;
				#print $testlength, "bp longseq in substri \n";
			} # ending while		
			
		}
	 #reset for next
		$startpos =0; 
		$testlength=0;
		$chunksize=$chunkvalue;
     
} #end while

&runtime(\*OUTsat); #for run log information
if($seqformat eq "TRUE"){
	# output formated sequenced once it is a qualified sequence
	print OUTsat "statistic summary of input sequence\(s\)",  "\n";
	print OUTsat "total input sequences #: ", $seqno, "\n"; 
	print OUTsat "total length (bp) after formatted: ", $seqlen, "\n";
	print OUTsat "total chunking sites after formatted: ", $ii, "\n";
}else{ 
	# check the input sequence is qualified or not
	print OUTsat "input sequence format is not qualified. check the presence of sign '>' \n";
	exit;
}


close FILE;
close OUT;
close OUTsat;

sub usageinfo
 {# print program name and usage for help if no input aguments available in command line
 my @usage=(); # showing content on how to use the programme
 $usage[0]="Usage: perl formatchunk.pl -i ssg_groupSeqfasta.fasta  -o formated.ssg_groupSeqfasta.fasta\n";
 $usage[1]=" for    help: perl $0 ; \n";
 $usage[2]=" for running: perl $0  -i [input DNA sequence file name]  -o [output file name] \n";
 $usage[3]="Author: Xuewen Wang\n";
 $usage[4]="year 2009-2012\n";
 unless(@ARGV){print @usage; exit;} 
 }

sub runtime() {
my $OUTfile=shift @_;
my $local_time = gmtime();
print {$OUTfile} "$0 was run and results were yielded at $local_time\n";
}
exit;

