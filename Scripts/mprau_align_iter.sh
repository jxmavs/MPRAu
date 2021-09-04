file_name=$1
file_path=$2
temp_path=$3
out_path=$4
bwa_ind=$5
alignment_fasta_file=$6
info_file=$7
numCores=$8

#use a specific version of python 2.7 due to Seq not being installed on python_2.7.9
python_dir=/broad/software/free/Linux/redhat_6_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin
#python -c "import sys; print sys.executable"

#see /Users/jxue/Documents/Dustin_Project/users_dustin_sequencing_HLBCXX_BWA_execute_alignment_merge_trim_M_qF_UGER for an example script of the pipeline I follow

#trim the fastq files

#cut off first bp
script_dir=/idi/sabeti-scratch/jxue/dustin_project/CMS-GWAS/scripts

#gunzip -c ${file_path}/${file_name}.fastq.gz > ${temp_path}/${file_name}.fastq

#num_initial_reads=$( gunzip -c ${temp_path}/${file_name}.fastq.gz | wc -l | awk '{print $1/4}' ) 

################## TRIM #############################
#also adds barcode info to read name
#V3 adds on stats info and reads in gzipped file
#uncomment
#${python_dir}/python ${script_dir}/rmAdaptersV3.py ${file_path}/${file_name}.fastq.gz -o ${temp_path}/${file_name}.merged.trimmed.fastq -d ${temp_path}/${file_name}.merged.discarded.fastq -readType 1 -barxN -s ${temp_path}/${file_name}.rmAdapters.stats

num_initial_reads=$(awk 'NR>1{print $1}' ${temp_path}/${file_name}.rmAdapters.stats) 

num_trimmed_discarded=$(awk 'NR>1{print $2}' ${temp_path}/${file_name}.rmAdapters.stats) 

num_trimmed_discarded_non_n=$(awk 'NR>1{print $3}' ${temp_path}/${file_name}.rmAdapters.stats)

num_trimmed_kept=$((num_initial_reads-num_trimmed_discarded))

#num_trimmed_kept=$( wc -l ${temp_path}/${file_name}.merged.trimmed.fastq | awk '{print $1/4}' )
#num_trimmed_discarded=$( wc -l ${temp_path}/${file_name}.merged.discarded.fastq | awk '{print $1/4}' )

#note that sequences in discarded.fastq don't contain sequences with N's in barcode
#so num_trimmed_discarded+num_trimmed_kept is less than num_initial_reads

#uncomment
#gzip ${temp_path}/${file_name}.merged.discarded.fastq 
#gzip ${temp_path}/${file_name}.merged.trimmed.fastq

#rm -f ${temp_path}/${file_name}.fastq

########### process stuff here #########


bwa mem -M ${bwa_ind} -t ${numCores} ${temp_path}/${file_name}.merged.trimmed.fastq.gz > ${temp_path}/${file_name}.merged.trimmed.sam

#rm -f ${temp_path}/${file_name}.merged.trimmed.fastq

############# post alignment processing ##########################

########## filter SAM #################

#256 removes duplicates

samtools view -h -S -q 1 -F 256 ${temp_path}/${file_name}.merged.trimmed.sam -o ${temp_path}/${file_name}.merged.trimmed.qF.sam -U ${temp_path}/${file_name}.merged.trimmed.qF.unaligned.sam

#don't need to uniq since there are no sequences with the same read name in the file
num_aligned=$( grep -v @ ${temp_path}/${file_name}.merged.trimmed.qF.sam |  wc -l | awk '{print $1}' )

#uniq due to reads multimapping
#there is overlap of ids in the unaligned file with the ids in the aligned file due to multimapping, in a multimap, one read is kept, the other random read is sent to unaligned
#https://www.biostars.org/p/304614/
#num_unaligned=$( grep -v @ ${temp_path}/${file_name}.merged.trimmed.qF.unaligned.sam | awk '{print $1}' | sort | uniq | wc -l | awk '{print $1}' )

#keep a record of the bam file after alignment, as bwa randomly splits reads that can align to multiple places 

samtools view -S -b ${temp_path}/${file_name}.merged.trimmed.sam > ${temp_path}/${file_name}.merged.trimmed.bam

rm -f ${temp_path}/${file_name}.merged.trimmed.sam

samtools view -S -b ${temp_path}/${file_name}.merged.trimmed.qF.unaligned.sam > ${temp_path}/${file_name}.merged.trimmed.qF.unaligned.bam

rm -f ${temp_path}/${file_name}.merged.trimmed.qF.unaligned.sam 

################### construct the calmd flag #################

#remove previous indeces
#rm -f ${alignment_fasta_file}.fai
samtools calmd ${temp_path}/${file_name}.merged.trimmed.qF.sam ${alignment_fasta_file} > ${temp_path}/${file_name}.merged.trimmed.qF.md.sam

rm -f ${temp_path}/${file_name}.merged.trimmed.qF.sam

#### NOTE:  I can potentially optimize below by reading from a bam file and writing to one
################### filter for homology#################################

#make sure to make seq parameters table first if needed
#difference is the info_file being read as input instead of being hard coded
#V3 adds on stats info
${python_dir}/python ${script_dir}/filterSAMV3.py -infoF ${info_file} -inputF ${temp_path}/${file_name}.merged.trimmed.qF.md.sam -passF ${temp_path}/${file_name}.merged.trimmed.qF.md_filtered.sam -failF ${temp_path}/${file_name}.merged.trimmed.qF.md_failed.sam -s ${temp_path}/${file_name}.homologyFilter.stats

#num_pass_homology_filter=$( grep -v @ ${temp_path}/${file_name}.merged.trimmed.qF.md_filtered.sam | wc -l | awk '{print $1}' )

#num_fail_homology_filter=$( grep -v @ ${temp_path}/${file_name}.merged.trimmed.qF.md_failed.sam | wc -l | awk '{print $1}' )

num_pass_homology_filter=$(awk 'NR>1{print $1}' ${temp_path}/${file_name}.homologyFilter.stats)
num_fail_homology_filter=$(awk 'NR>1{print $2}' ${temp_path}/${file_name}.homologyFilter.stats)


samtools view -S -b ${temp_path}/${file_name}.merged.trimmed.qF.md.sam  > ${temp_path}/${file_name}.merged.trimmed.qF.md.bam

rm -f ${temp_path}/${file_name}.merged.trimmed.qF.md.sam


samtools view -S -b ${temp_path}/${file_name}.merged.trimmed.qF.md_failed.sam > ${temp_path}/${file_name}.merged.trimmed.qF.md_failed.bam

rm -f ${temp_path}/${file_name}.merged.trimmed.qF.md_failed.sam


######## get the tables, with barcodes and without #####################
#simply add up the counts
${python_dir}/python ${script_dir}/compileReads_no_dup_sep.py -sam ${temp_path}/${file_name}.merged.trimmed.qF.md_filtered.sam -fasta ${alignment_fasta_file}

samtools view -S -b ${temp_path}/${file_name}.merged.trimmed.qF.md_filtered.sam > ${temp_path}/${file_name}.merged.trimmed.qF.md_filtered.bam

rm -f ${temp_path}/${file_name}.merged.trimmed.qF.md_filtered.sam

################# create counts per oligo table ###################################################
#sum up all the barcode counts per row across all columns
#no header in file, just name and count
awk 'NR>1{sum=0; for(i=2; i<=NF; i++) {sum += $i}; print $1"\t"sum }' ${temp_path}/${file_name}.merged.trimmed.qF.md_filtered.counts.table > ${temp_path}/${file_name}.merged.trimmed.qF.md_filtered.counts.nobars.table.parsed

############ create a barcodes only count table  #########################
#do i-1 to ignore the oligo column
awk 'NR==1{for (i=2;i<=NF;i++) names[i-1]=$i} NR>1{for (i=2;i<=NF;i++) count[i-1]+=$i} END{for (i=1;i<=NF-1;i++){print names[i]"\t"count[i]}}' ${temp_path}/${file_name}.merged.trimmed.qF.md_filtered.counts.table > ${temp_path}/${file_name}.merged.trimmed.qF.md_filtered.counts.barcodes_only.table.parsed

#rm ${temp_path}/${file_name}.merged.trimmed.qF.md_filtered.counts.table

########### output summary statistics into a file ##########################
#num_aligned should equal num_pass_homology_filter plus num_fail_homology_filter

echo -e "${num_initial_reads}\t${num_trimmed_kept}\t${num_trimmed_discarded}\t${num_trimmed_discarded_non_n}\t${num_aligned}\t${num_pass_homology_filter}\t${num_fail_homology_filter}" > ${temp_path}/${file_name}_alignment_stats.txt
