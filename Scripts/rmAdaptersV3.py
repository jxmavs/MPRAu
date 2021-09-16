## Dustin Griesemer
##
## 02/15/15
##
## rmAdapters.py
##
## usage: python rmAdapters.py fastq_file_path -o output_file_path -d discard_file_path [-Lbarcode Lbarcode] [-Rbarcode Rbarcode] [-Radapter Radapter] [-seqxN] [-barxN]
##
## example: python rmAdapters.py /idi/sabeti-data/dustin/sequencing/AALE9_2/fastq/AALE9.AACTTGAC.extendedFrags.fastq -o trimmed.expos.fastq -d discarded.expos.fastq 
##
## Removes adapters from an MPRA fastq file, cataloging barcodes in the process. Outputs an updated fastq file, with the adapters trimmed and the barcode in the sequence identifier. If the -seqxN option is used, sequences with N's are removed. If the -barxN option is used, sequences with N's in the barcode position are removed.

import sys
import os
import argparse
import fuzzysearch as fs
from Bio import SeqIO
from Bio.Seq import Seq
import Levenshtein as lv
import gzip

def Ladapter_trim(seq,Lbar,Rbar,Nbar=6,Lbar_dist=3,Rbar_dist=3,Bar_clamp=2):
#Takes in sequence in string form, outputs index at which trimming should occur and barcode (or error) sequence
#Avoids using fuzzysearch for efficiency
#Searches for an exact match to 10bp on either the left or right side of the barcode.
#Finds Lbar & Rbar based on position and performs Levenshtein distance check.
#Nbar: expected length of barcode
#Finds Barcode, requiring an exact match of Bar_clamp bases on either side.
	if Lbar[-10:] in seq:
		Lbar_pos=seq.find(Lbar[-10:])-len(Lbar)+10
		Bar_pos=Lbar_pos+len(Lbar)
		Rbar_pos=Bar_pos+Nbar
	elif Rbar[0:10] in seq:
		Rbar_pos=seq.find(Rbar[0:10])
		Bar_pos=Rbar_pos-Nbar
		Lbar_pos=Bar_pos-len(Lbar)
	else:
		return 'Lbar/Rbar_not_found',0 #Couldn't find an exact match on either left or right of barcode
	Lbar_seq=seq[Lbar_pos:(Lbar_pos+len(Lbar))]
	Bar_seq=seq[Bar_pos:(Bar_pos+Nbar)]
	Rbar_seq=seq[Rbar_pos:(Rbar_pos+len(Rbar))]
	if lv.distance(Lbar,Lbar_seq)>Lbar_dist:
		return 'Lbar_incorrect',0 #Lbar sequence was too different from expected
	if lv.distance(Rbar,Rbar_seq)>Rbar_dist:
		return 'Rbar_incorrect',0 #Rbar sequence was too different from expected
	if seq[(Bar_pos-Bar_clamp):Bar_pos]!=Lbar[-Bar_clamp:]:
		return 'Lclamp_mismatch',0 #clamp left of barcode mismatch
	if seq[Rbar_pos:(Rbar_pos+Bar_clamp)]!=Rbar[0:Bar_clamp]:
		return 'Rclamp_mismatch',0 #clamp right of barcode mismatch
#	return Bar_seq,(Rbar_pos+len(Rbar))
	return Bar_seq,Rbar_pos #don't trim Rbar


def Ladapter_trim_fuzzy(seq,Lbar,Rbar,Nbar=6):
#takes in sequence in string form, outputs index at which trimming should occur and barcode sequence
	matches_l=fs.find_near_matches(Lbar,seq,max_l_dist=3)
	if matches_l==[]:
		return 'Lbar not found',0
	else:
		dist_l=[i.dist for i in matches_l]
		min_index_l=dist_l.index(min(dist_l))
		match_l=matches_l[min_index_l]

	matches_r=fs.find_near_matches(Rbar,seq,max_l_dist=3)
	if matches_r==[]:
		return 'Rbar not found',0
	else:
		dist_r=[i.dist for i in matches_r]
		min_index_r=dist_r.index(min(dist_r))
		match_r=matches_r[min_index_r]

	N=match_r.start-match_l.end
	if N!=Nbar:
		return 'Barcode incorrect length',0
	else:
		bar=seq[match_l.end:match_r.start]
		trim_index=match_r.end
			
	return bar,int(trim_index)


def Radapter_trim(seq,Radp):
#takes in sequence in string form, outputs index at which trimming should occur
#avoids using fuzzysearch as much as possible for efficiency
	trim_index=seq.find(Radp) #requires an exact match
	if trim_index==-1:
		trim_index=len(seq) #if an exact match is not found, don't alter the sequence
	else:
		trim_index=trim_index+len(Radp) #trim AFTER Radapter
	return int(trim_index)


def Radapter_trim_fuzzy(seq,Radp):
#takes in sequence in string form, outputs index at which trimming should occur
	matches=fs.find_near_matches(Radp,seq,max_l_dist=3)
	if matches==[]:
		return len(seq)
	else:
		dist=[i.dist for i in matches]
		min_index=dist.index(min(dist))
		match=matches[min_index]
		trim_index=match.start
		
	return int(trim_index)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq_file_path',type=str)
    parser.add_argument('-o',type=str,help='output file path for trimmed sequences')
    parser.add_argument('-d',type=str,help='output file path for discarded sequences')
    parser.add_argument('-readType',type=str,help='type of read sequence (1) read 1, (2) read 2, (3) paired end')
    #parser.add_argument('-Lbar',default='CGAGCTGTACAAGTAATTCTAGTTG',type=str,help='sequence left of barcode') #PGKGiga
    parser.add_argument('-Lbar',default='GTTTAAAGCCCAACGCTAGTC',type=str,help='sequence left of barcode') #Bmt
    parser.add_argument('-Rbar',default='CGAGCTCGCTAGCCT',type=str,help='sequence right of barcode')
    parser.add_argument('-Radp',default='AGATCGGAAGAGCGTCG',type=str,help='sequence of right adapter')
    parser.add_argument('-barxN',action='store_true')
    parser.add_argument('-seqxN',action='store_true')
    parser.add_argument('-s',type=str,help='output file path for stats')
    parsed = parser.parse_args()

    fastq_file = gzip.open(parsed.fastq_file_path,'rb')
    trimmed_file = open(parsed.o,'w')
    trimmed=[]
    discard_file = open(parsed.d,'w')
    discarded=[]

    num_total_seq=0
    num_total_discard=0
    num_total_discard_non_N=0
    #	it=0
    #	Nline=0
    for record in SeqIO.parse(fastq_file,'fastq'):
    #		it=it+1
    #		if it==100000:
    #			Nline=Nline+10000
    #			print Nline
    #			it=0			
    #			break
        num_total_seq+=1
        if parsed.readType=="1":
            #reverse comeplmeent
            #changed 1/29/19
            #prev_annotations=record.letter_annotations
            #record.letter_annotations = {}
            #record.seq=Seq(str(record.seq.reverse_complement()))
            #record.letter_annotations=prev_annotations
            #seq=str(record.seq)
            #do this to also reverse complement the letter annotations
            description=record.description
            record=record.reverse_complement()
            record.description=description
            seq=str(record.seq)
        else:
            #seq=str(record.seq)
            seq=str(record.seq)

        [barcode,l_index] = Ladapter_trim(seq,parsed.Lbar,parsed.Rbar)
        record.description = record.description + ' bar:%s' %(barcode) #add barcode to description
        record.description='_'.join(record.description.split(' ')) #convert ' ' to '_' so samtools doesn't cut them off
        record.id=record.description #so separate fields aren't added
        record.name=record.description #so separate fields aren't added
        
        if parsed.seqxN: #if the seqxN option is used, skip sequences with N's anywhere.
            if 'N' in seq:
                num_total_discard+=1
                discarded=[record]
                SeqIO.write(discarded,discard_file,'fastq')
                continue

        #[barcode,l_index] = Ladapter_trim(seq,parsed.Lbar,parsed.Rbar)
        if parsed.barxN: #if the barxN option is used, skip sequences with N in the barcode position.
            if 'N' in barcode:
                num_total_discard+=1
                discarded=[record]
                SeqIO.write(discarded,discard_file,'fastq')
                continue

        #record.description = record.description + ' bar:%s' %(barcode) #add barcode to description
        #record.description='_'.join(record.description.split(' ')) #convert ' ' to '_' so samtools doesn't cut them off
        #record.id=record.description #so separate fields aren't added
        #record.name=record.description #so separate fields aren't added
        if l_index==0: #failed to identify the barcode sequence
            num_total_discard+=1
            num_total_discard_non_N+=1 
            discarded=[record]
            SeqIO.write(discarded,discard_file,'fastq')
            continue
       
        #don't trim if just read 1, because the sequence in the beginning is immediately the oligo variable sequence
        if parsed.readType=="1":
            r_index=len(seq)
        else:
            r_index = Radapter_trim(seq,parsed.Radp)
        record=record[l_index:r_index]
    #		record=record[l_index:] #don't trim Radapter

        trimmed=[record]
        SeqIO.write(trimmed,trimmed_file,'fastq')
    #		trimmed.append(record)

    #	SeqIO.write(discarded,discard_file,'fastq')
    #	SeqIO.write(trimmed,trimmed_file,'fastq')	

    fastq_file.close()
    trimmed_file.close()
    discard_file.close()
    
    with open(parsed.s, "w") as out_stats_file:
        out_stats_file.write("Num_Total_Seq\tNum_Total_Discard\tNum_Total_Discard_Non_N\n")
        out_stats_file.write("{0}\t{1}\t{2}\n".format(num_total_seq, num_total_discard, num_total_discard_non_N))

def one_record():
#	Lbar='CGAGCTGTACAAGTAATTCTAGTTG' #PGKGiga
	Lbar='GTTTAAAGCCCAACGCTAGTC' #BmtGiga
	Rbar='CGAGCTCGCTAGCCT'
	Radp='AGATCGGAAGAGCGTCG'
	
	fastq_file=open('/idi/sabeti-data/dustin/sequencing/AALE9_2/fastq/AALE9.AACTTGAC.extendedFrags.fastq','r')
	record=SeqIO.parse(fastq_file,'fastq').next()

	seq=str(record.seq)
		
	[barcode,l_index] = Ladapter_trim(seq,Lbar,Rbar)
	if l_index==0: #failed to identify the barcode sequence
		print 'DISCARDED SEQUENCE'
		print record.format('fastq')
		fastq_file.close()
		return
	record.description = record.description + ' bar:%s' %(barcode)
	record.description='_'.join(record.description.split(' '))
	record.id=record.description
	record.name=record.description
	
	r_index = Radapter_trim(seq,Radp)
	record=record[l_index:r_index]
#	record=record[l_index:] #don't trim Radapter

	print 'TRIMMED SEQUENCE'
	print record.format('fastq')
	
	fastq_file.close()


if __name__ == '__main__':
	if len(sys.argv)==1:
		one_record()
	else:
		main()
