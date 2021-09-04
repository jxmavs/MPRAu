## 01/16/19
##
## James Xue
##
## compileReads.py
##
## usage: python compileReads.py -sam sam_file_path -fasta fasta_file_path [-output_dir output_file_directory] [-output_pre output_file_prefix]
##
## Compiles table of oligo & barcode counts from filtered sam file. Uses fasta file for oligo sequence names.
## Incorporates tag for type of oligo (ref, alt, CXCL2, random)
##
## Header: oligo_name, tag, indices (ordered)

import sys
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-sam',help='filtered sam file')
parser.add_argument('-fasta',help='fasta alignment file')
parser.add_argument('-output_dir',help='output file directory')
parser.add_argument('-output_pre',help='output file prefix')
#parser.add_argument('-nobars',help='collapse barcodes',action='store_true')
parsed = parser.parse_args()

if parsed.output_dir is None:
	parsed.output_dir='/'.join(parsed.sam.split('/')[:-1])
if parsed.output_pre is None:
	parsed.output_pre='.'.join(parsed.sam.split('/')[-1].split('.')[:-1])

sam_file=open(parsed.sam,'r') #output of SAMtools idxstats function
fasta_file=open(parsed.fasta,'r') #fasta alignment file to retrieve oligo names

#if parsed.nobars:
#	output_file_path=parsed.output_dir+'/'+parsed.output_pre+'.counts.nobars.table'
#else:
#	output_file_path=parsed.output_dir+'/'+parsed.output_pre+'.counts.table'
output_file_path=parsed.output_dir+'/'+parsed.output_pre+'.counts.table'

output_file=open(output_file_path,'w')

bases=['A','C','G','T']
indices=[a+b+c+d+e+f for a in bases for b in bases for c in bases for d in bases for e in bases for f in bases]

refDict={}
#n=0
for line in fasta_file:
	if line[0]=='>':
#		n=n+1
#		print n
		Rname=line.strip()[1:]
		refDict[Rname]={}
		for index in indices:
			refDict[Rname][index]=0
	else:
		continue

#lineN=0
for line in sam_file:
	if line[0]=='@':
		continue
#	else:
#		lineN=lineN+1
#	if lineN>100000:
#		break
	vals=line.strip().split('\t')
	Qname=vals[0]
	Rname=vals[2]
    #barcode is last 6 bp of read name
	index=Qname[-6:]
	refDict[Rname][index]=refDict[Rname][index]+1

#print 'Finished refDict'

#write out header
output_file.write('oligo_name\t')
for ind in range(len(indices)-1):
    index=indices[ind]
    output_txt='%s\t' %(index)
    output_file.write(output_txt)

index=indices[len(indices)-1]
output_txt='%s\n' %(index)
output_file.write(output_txt)

for oligo_name in sorted(refDict.keys()):
    output_txt='%s\t' %(oligo_name)
    output_file.write(output_txt)
    for ind in range(len(indices)-1):
        index=indices[ind]
        output_txt='%s\t' %(refDict[oligo_name][index])
        output_file.write(output_txt)
    index=indices[len(indices)-1]
    output_txt='%s\n' %(refDict[oligo_name][index])
    output_file.write(output_txt)

output_file.close()
