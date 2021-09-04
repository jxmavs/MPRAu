## 06/14/16
##
## Dustin Griesemer
##
## filterSAM.py -inputF SAM_input_path -passF pass_reads_path -failF fail_reads_path

## Example: python /home/unix/dusting/src/MPRA/seq/filterSAM.py -inputF /idi/sabeti-data/dustin/sequencing/HLBCXX/BWA/Bmt1_CTGCGGAT.merged.trimmed_seqxN_M_qF_md.sam -passF /idi/sabeti-data/dustin/sequencing/HLBCXX/BWA/Bmt1_CTGCGGAT.merged.trimmed_seqxN_M_qF_md_filtered.sam -failF /idi/sabeti-data/dustin/sequencing/HLBCXX/BWA/Bmt1_CTGCGGAT.merged.trimmed_seqxN_M_qF_md_failed.sam
##
## Filter SAM file for read homology
##

import sys
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('-inputF')
parser.add_argument('-passF')
parser.add_argument('-failF')
parser.add_argument('-infoF')
parser.add_argument('-s')
parsed = parser.parse_args()

input_file=open(parsed.inputF,'r')
info_file=open(parsed.infoF,'r')

pass_file=open(parsed.passF,'w')
fail_file=open(parsed.failF,'w')

length={}

for line in info_file:
	if line[0]=='#':
		continue
	vals=line.strip().split('\t')
	seqID=vals[0]
	seq_length=vals[4]
	length[seqID]=int(seq_length)
	
num_pass_homology_filter=0
num_fail_homology_filter=0

for line in input_file:
    if line[0]=='@':
        pass_file.write(line)
        fail_file.write(line)
        continue
    vals=line.strip().split('\t')
    seqID=vals[2]
    seq_length=length[seqID]
    MD_tag=vals[12].split(':')[-1]
    MD_vals=re.findall('[0-9]+|\^[A-Z]+|[A-Z]+',MD_tag)
    match_count=0
    mismatch_count=0
    for entry in MD_vals:
        if entry=='':
            continue
        try:
            match_count=match_count+int(entry)
        except ValueError:
            if entry[0]=='^':
                entry=entry[1:]
            mismatch_count=mismatch_count+len(entry)
    MD_length=match_count+mismatch_count
    MD_max=min(seq_length,123) #only 123bp can be read at most
    if MD_length>MD_max:
        percent_match=(float(MD_max)-float(mismatch_count))/float(MD_max)
    else:
        percent_match=float(match_count)/float(MD_max)
    if percent_match>0.95:
        num_pass_homology_filter+=1
        pass_file.write(line)
    else:
        num_fail_homology_filter+=1
        fail_file.write(line)

with open(parsed.s, "w") as out_stats_file:
    out_stats_file.write("Num_Pass_Filter\tNum_Fail_Filter\n")
    out_stats_file.write("{0}\t{1}\n".format(num_pass_homology_filter, num_fail_homology_filter))
