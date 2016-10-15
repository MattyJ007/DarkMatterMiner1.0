#!/usr/bin/perl -w
use strict;

# script was written by Xuewen Wang, originally designed for 
#        looking repeated motif e.g. SSR in fasta file
# input file should in fasta formatted:   
#    >ID in oneline 
#    seq in next lines
# this version run very fast and less memory is required. no bug was found.
#Usage: $0 -i [formated sequence file] -o [motif searching results] -r [minimum repeated times] 
#       -m [minimum motif length] -x [maximum motif length] -s [highlight motif in sequence 0 or 1]
# print "Software $0 written by Xuewen Wang from 2012 Jan.  
#       The function is to find motif and return motif type, length, position from a fasta file. 
#       such as SSR\n" ;
#run e.g. :
#       perl gssr.pl -r 5 -m 2 -x 10 -s 0 -i testseq.fasta -o testseq.fasta.ssr
#faster speed for long string DNA seq , for example 50Mb or longer in genome DNA
#treat short seq fragment at a time and then nest segment, with short end seq overlapped
&usageinfo();

# get command options
my %commandin = @ARGV;
if ((scalar @ARGV)%2 != 0){print "arguments must be in pair";}
my $refin=$commandin{"-i"}; #information of input sequence file
my $outfiledefault=$refin.".ssr"; # default output file name if no output file name is given
my $outputfile=$commandin{"-o"}||$outfiledefault; #information of output sequence file
my $motif_min=$commandin{"-m"}||2;# minimum length of any repeat unit, default value is 2
my $motif_max=$commandin{"-x"}||10; # maximum length of any repeat unit,, default value is 10
my $motif_times_min=$commandin{"-r"}||5; # minimum repeated times of a repeat unit, default value is 5
my $seqhight=$commandin{"-s"}||0; #highlight all the repeated sequence consisting of motif in each source sequence 

#file status check and prepare
 unless(-e $refin){print "input file does not exist in the working directoty, or the path is not given\n";}
 open (DNAhand, "<$refin")|| die (" $0 filed open failed, pleasse check");
 open(OUT, ">$outputfile")|| die (" $0 Write file failed, pleasse check");
 my $runlog=$0."\.log";
 open (OUTstat,">>$runlog")|| die (" LOG file writing failed, please check");
 print OUT "Name\tSeq_Len\tStartPos\tEndPos\tRepetitions\tMotif\n"; #prepare head
 
#initial
 my $seqID="";
 my $lengthseq=0;
 my $Sequences="";
 my $Sequences_sub="";
 my $positionm=0;
 my $motiflen=0;
 my $readincount=0;
 my $newID="";
 my $oldID="";
 my $leftm=0;
 my $rightm=0;
 my @SSRloci="";
 my $SSRinfo="";
 my $totallength=0;
			
 $motif_times_min=$motif_times_min-1; # recaculat the motif minimum repeated times
#start to read in sequence data in while loop


while ( $Sequences=<DNAhand>){	# one line at at time
	chomp $Sequences;
	
	#get seq ID and seq	    				
		if($Sequences=~m/^>.+/){ 
		    
			$seqID=$Sequences;
            $positionm=0;
            $Sequences_sub="";
			$Sequences="";
            $leftm=0;
            $rightm=0;
            $readincount=0;
			
			
             #output SSR loci for current sequence
	         foreach $SSRinfo(@SSRloci){
			   (print OUT $SSRinfo, "\t",$totallength,"\n" )if(length $SSRinfo != 0);
	         }
			 @SSRloci="";	#reset
			 $totallength=0;
			
		}else{
			$Sequences=~ s/[0-9\n\s]//g;		# remove newlines, numbers and white space				
			$totallength +=length $Sequences;
			$Sequences=$Sequences_sub.$Sequences;
			$readincount =$readincount+1;
			$lengthseq=length $Sequences;
			
		#} #else to do while search:	looking for motif	
		
	
	
	#looking for motif
	# get the motif return in ()
	while($Sequences=~ /([ACGT]{$motif_min,$motif_max}?)\1{$motif_times_min,}/gix){ #search perfect 
		#\1 to double the unit and repeated times, real repeated times should plus one from first match in ()
		#backref to get the same unit stored in $1
		#print OUT $motif_min, $motif_max, $motif_times_min, "\n";	
		my $motif=uc $1;
		#print $1,"*1st match*\t";
        $motiflen=length $motif;		
		my $removelen=int($motiflen/2)+1 ;
		$removelen =1 if ($removelen <=1);
		$leftm=$-[0]+1; #  left position, start position in current segment
		$rightm= $+[0]; #  right position in current segment
	    #print $leftm, "\t", $rightm," 1st matchposition\t";
		
		#to remove
		unless ($motif=~/^([ACGT]{1,$removelen})\1+$/ix){#motif while 
			#print $1,"*2*\t";
		#filter out the motif which will be duplicated in later longer motif , 
		#eg. A is motif, AA should not be motif if the setting $mitif_min is changed to 2.
		#e.g. AC is motif, ACAC should not be motif if the setting $mitif_min is set to four
			#my $microsate= $&; # detailed sequence of microsatellite		
			#my $repeat_times = length($microsate)/(length $motif);			
						
			my $repeat_times=($rightm-$leftm+1)/$motiflen;
			$leftm=$leftm+$positionm;
			$rightm=$rightm+$positionm;
			 #print "Positionm: $positionm\t";
			 #print $leftm, "\t", $rightm," after 2nd matchposition\n";
			
			#temperary store SSR loci for each sequence
			 push (@SSRloci,join("\t",$seqID,$readincount, $leftm,$rightm,$repeat_times,$motif));			
			
			#output results, tab delimited
			#print  "Name\tSeq_Len\tStartPos\tEndPos\trepetitions\tSub_unit\n"; #prepare head
			#print  $seqID,"\t",$lengthseq,"\t";
			#print  $leftm,"\t"; # print left position, start position
			#print  $rightm,"\t"; #print right position
			#print  $repeat_times, "\t"; #"repeated times: ",
			#print  $motif, "\n"; #"motif: ",
			
			
			#output highlighted motif in sequence which the motif is coming from
			#output original sequence if asked
			if ($seqhight == 1){
				print OUT lc $` ; # print seq before left position, start position
				#$microsate= $&;
				print OUT uc $& ; #"Microsatellites: ",
				print OUT lc $',"\n" ; #print sequence after right position
			} #end if
			
		 } #ending unless 
	    } #ending while for getting motif 
		  
	     
		 
		 # Reset the position for next line, set reconsider length=20
		 my $reconsiderlen=20;
         $positionm = $positionm + $lengthseq - $reconsiderlen + 1;
		 #if is new seqID, set to 0;

        # Discard the data in the buffer, except for a portion at the end
        # so patterns that appear across line breaks are not missed
        $Sequences_sub = substr($Sequences, $lengthseq - $reconsiderlen + 1, $reconsiderlen - 1);
	    #if is new seqID, set $Sequences_sub to 0;
	
	} # end  else to do while search:	looking for motif
	
	if (eof(DNAhand)){ #print data for last seqID
			         foreach $SSRinfo(@SSRloci){
						(print OUT $SSRinfo, "\t" )if(length $SSRinfo != 0);
						 (print OUT $totallength, "\n")if(length $SSRinfo != 0);
					}
					
	          
    }
} #ending while readin

&runtime(\*OUTstat); #for run log information
close DNAhand;
close OUT;
close OUTstat;
 
sub usageinfo
	{# print program name and usage for help if no input aguments available in command line
	my @usage=(); # showing content on how to use the programme \n";
	$usage[1]=" for    help: perl $0 ; \n";
	$usage[2]=" for running: perl $0 -i [formated sequence file] -o [motif searching results] 
	-r [minimum repeated times] -m [minimum motif length] -x [maximum motif length] 
	-s [highlight motif in sequence 0 or 1] \n";
	$usage[3]="Author: Xuewen Wang\n";
	$usage[4]="year 2012\n";
	unless(@ARGV){print @usage; exit;} 
 }

sub runtime(){
	my $OUTfile=shift @_;
	my $local_time = gmtime();
	print {$OUTfile} "$0 was run and results were yielded at $local_time\n";
}

exit;# exit whole programe
