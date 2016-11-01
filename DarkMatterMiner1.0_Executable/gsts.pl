#!/usr/bin/perl -w
use strict;
# statement
	# script written by Xuewen Wang July, Dec, 2012, copy right reserved
	# original designed for statistic analysis for result output of motif , 
	# e.g. SSR, from gssrtrim.pl mining  script written by Xuewen Wang.
	# the following statistic results or tables were provided after running this script.
	# all results in tabular format and can be open in any spreadsheet such as MS Excell software;
	# table 1 provides total occurence of each mer based on length
	# table 2 provides total occurence of each motif which does NOT count the reverse complement motif based on sequence of each motif
	# table 3 provides total occurence of each grouped motif which counts the reverse complement motif based on sequence of each motif
	# table 4 provides total motifs occurence in each sequence or chromosome 
	# 	and the distribution of motif frequency in each sequence or chromosome
	# 	and the total sequence or chromosomes, total motifs, average frequence
	# this version run very fast and less memory is required. no bug was found.
	#current version: 2012 Dec.

# Usage: $0 -i [ssr result file]
# e.g. perl gsts.pl -i testfasta.fasta.ssr
	&usageinfo();

# get command options
	my %commandin = @ARGV;
	if ((scalar @ARGV)%2 != 0){print "arguments must be in pair";}
	my $refin=$commandin{"-i"}; #information of input file
	my $outputfile=$refin.".sat2";
	#my $refin="testseq.fasta.ssr"; # test file

# file status check and preparation
	open (DNAhand,$refin)|| die (" $0 file opening failed, please check");
	open (OUTstat,">$outputfile")|| die (" $0 file writing failed, please check");
# definition
	my $word="";
	my $freqref="";
	my $freqref_name ="";
	my $freqref_subunit="";
	my $freqref_unitlen="";
	my $freqref_subunit_motif="";
	my @name=(); 
	my @subunit=(); 
	my @unitlen=();
	my %seq_len=();

# read in data
	while ($word=<DNAhand>){ 	
		chomp $word;
		if($word =~ /^\>.+/){ # only read in the ssr information line
			my @cellinfor=split("\t",$word);
			# get information of >ID (name), Sub_unit, and Sub_unit length 
			# and then do statistic analysis
				push (@name, $cellinfor[0]); #name of ID, ok
				if (!exists $seq_len{$cellinfor[0]}){
					$seq_len{$cellinfor[0]}=$cellinfor[1];
				}# get sequence length
				push (@subunit, $cellinfor[5]); # Sub_unit type, ok
				push (@unitlen, length $cellinfor[5]); # Sub_unit length: -mers, e.g. dimer, ok
			
		}else{
		next; #ignore the highlighted lines if -s 1 was used in motif mining
		}
		
	} #end while

# do statistic analysis
	$freqref_unitlen = &fre_ct_batch(\@unitlen); #count mers occurence
	$freqref_subunit = &fre_ct_batch(\@subunit); #count motif ocurrence
	$freqref_subunit_motif = &fre_ct_batch_motif($freqref_subunit); #count motif ocurrence, pair grouped
	$freqref_name = &fre_ct_batch(\@name); #count ocurrence of ID of sequence 

# output results
	&runtime(\*OUTstat);# for running log
	
	#for table 1
	print OUTstat "Table 1\n";
	print OUTstat "Motif(-mer)\t total\n";
	&val_ascend($freqref_unitlen,\*OUTstat); #sort and print
	print OUTstat "\n";
	
	#for table 2
	print OUTstat "Table 2\n";
	print OUTstat "Motif\t total\n";
	&val_ascend($freqref_subunit, \*OUTstat); #sort and print
	print OUTstat "\n";	
	
	#for table 3
	print OUTstat "Table 3 paired\n";
	print OUTstat "Grouped_Motif\t total\n";
	&val_ascend($freqref_subunit_motif,\*OUTstat); #sort and print
	print OUTstat "\n";
	
	#for table 4
	print OUTstat "Table 4\n";
	print OUTstat "SeqID\tTotal_Motifs\tSeqSize\tFrequency\(Motifs\/Mb\)" ,"\n";
	&val_ascend_mb($freqref_name,\*OUTstat,\%seq_len);
	
#close and exit
	close 	DNAhand;
	exit;

#sub routines
############################
	sub val_ascend(){ #sort
	#to print key, value out of %{$freqref};
	#print "final output sorted by values (ocurrence)in ascending order:\n";
	#print "words (key)\t","ocurrence\n";
	
		my ($fref,$outfile)=@_;
		my %refhash=%{$fref};
		my $value_total=0;
		my $key_ct=0;
		foreach  my $key(sort {$refhash{$b} <=> $refhash{$a}} keys %refhash){
			print {$outfile} $key ,"\t", $refhash{$key},"\n";
			$value_total =$refhash{$key}+$value_total;#count the sum
			$key_ct =$key_ct+1; #coutn the total occurence of key
		}
		
		#print out results
		print {$outfile} "total_above\ttotal_above\n";
		print {$outfile} $key_ct, "\t",$value_total,"\n";
	
	} #end sub

#sort and do frequence caculation for each sequence
############################
	sub val_ascend_mb(){ #sort
	#to print key, value out of %{$freqref};
	#print "final output sorted by values (ocurrence) of sequence in ascending order:\n";
	#print "sequence\t","ocurrence in each sequence\t","sequence length\t","frequence per mbp\n";
	
	#definition
	my ($fref,$outfile,$seq_len_ref)=@_;
	my %refhash=%{$fref};
	my $value_total=0;
	my $key_ct=0;
	my $len_ct=0;
	my %seq_lensub=%{$seq_len_ref};
	
	#sort and caculate the statistic results
	foreach  my $key(sort {$refhash{$b} <=> $refhash{$a}} keys %refhash){ 
	# key is ID name, value is occurrence	
	# see previous in main :  $seq_len{$cellinfor[0]};
			#$seq_lensub{$key) is the lenght of each sequence
			print {$outfile} $key ,"\t", $refhash{$key}, "\t", $seq_lensub{$key}, "\t", $refhash{$key}/$seq_lensub{$key}*1000000,"\n";
			$value_total =$refhash{$key}+$value_total;#count the sum of occurence (value)
			$key_ct =$key_ct+1; #count the total occurence of key
			$len_ct=$len_ct+$seq_lensub{$key}; #count total length of all sequences
				
	} # end foreach
	
	# print out summary of each of first three collumns, and the average frequence of all sequence
	print {$outfile} "total_above", "\t","total_above","\t","total_above","\t", "average_frequency","\n";
	if($len_ct!=0){
		print {$outfile} $key_ct, "\t",$value_total,"\t",$len_ct,"\t", $value_total/$len_ct*1000000,"\n";}
} #end sub

#count occurrence
############################
	sub fre_ct_batch(){ 
		#definition
		my @cells=@{shift @_}; #ok
		my %allseqinformation= ();
		my $sword="";
	
		#count occurence based on the key value
		foreach $sword (@cells){
			if(! exists $allseqinformation{$sword}){$allseqinformation{$sword}=1;}
			else{$allseqinformation{$sword}++;}
		#print $sword ,"\t"; # word content
		#print $allseqinformation{$sword},"\n"; #word count
		} #end each
		return \%allseqinformation;
	} # end sub

#sub count grouped motif distribution
############################
	sub fre_ct_batch_motif(){ # pass a hash reference to here, grouping
		#definition
		my %refhash=%{shift @_}; #ok, only motif passed here
		my %refhashc=(); #paired value
		
		#find the grouped value and then do statistic
		foreach  my $key( keys %refhash){
			#get complemented index in $key 
			my $keytmp=$key;
			$keytmp=~tr/AGTC/TCAG/;
			my $key_cr=reverse $keytmp;
			my $pairkey=$key."\/".$key_cr;
	
			if($key ne $key_cr){ # keep the same value , no addition when equal: $key ne $key_cr
				if(exists $refhash{$key_cr}){			
					( $refhashc{$pairkey}=$refhash{$key_cr}+$refhash{$key})if($refhash{$key}!= -1); #get the paired value
					$refhash{$key_cr}=-1; # marked complemented key value			
				}else{ # not exist $refhash{$key_cr}
					( $refhashc{$pairkey}=$refhash{$key})if($refhash{$key}!= -1); #get the paired value: only it self
				}
			}else{ # equal keys means there is no complemented key value here
				($refhashc{$pairkey}=$refhash{$key})if($refhash{$key}!= -1); #the paired value is only itself
			}
	
			#print "pairvalue  ", $pairkey, "\t", $refhashc{$pairkey},"\n"; # only one of  duplicated reverse pairs eg.ag/ct, ct/ag have value because more pairs key than caculated values
		} # end each
		return \%refhashc;
	} #end sub

#sub 
############################
	sub runtime() {
		my $OUTfile=shift @_;
		my $local_time = gmtime();
		print {$OUTfile} "$0 was run and results were yielded at $local_time\n";
	} # end sub

#sub 
############################
	sub usageinfo{
	# print program name and usage for help if no input aguments available in command line
		my @usage=(); # showing content on how to use the programme
		$usage[0]="Usage: statistic analysis for result output from motif script gssrtrim.pl.\n";
		$usage[1]=" for    help: perl $0 ; \n";
		$usage[2]=" for running: perl $0 -i [ssr result file name from gssrtrim.pl] \n";
		$usage[3]="Author: Xuewen Wang\n";
		$usage[4]="year 2012\n";
		unless(@ARGV){print @usage; exit;} 
 } #end sub
 exit;