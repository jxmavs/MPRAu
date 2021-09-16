file_name=$1
temp_path=$2
bwa_ind=$3
alignment_fasta_file=$4
info_file=$5
numCores=$6


################## trim reads #############################

#also adds barcode info to read name
python ${script_dir}/rmAdaptersV3.py ${temp_path}/${file_name}.fastq.gz -o ${temp_path}/${file_name}.merged.trimmed.fastq -d ${temp_path}/${file_name}.merged.discarded.fastq -readType 1 -barxN -s ${temp_path}/${file_name}.rmAdapters.stats

num_initial_reads=$(awk 'NR>1{print $1}' ${temp_path}/${file_name}.rmAdapters.stats) 

num_trimmed_discarded=$(awk 'NR>1{print $2}' ${temp_path}/${file_name}.rmAdapters.stats) 

num_trimmed_discarded_non_n=$(awk 'NR>1{print $3}' ${temp_path}/${file_name}.rmAdapters.stats)

num_trimmed_kept=$((num_initial_reads-num_trimmed_discarded))

#note that sequences in discarded.fastq don't contain sequences with N's in barcode
#so num_trimmed_discarded+num_trimmed_kept is less than num_initial_reads

gzip ${temp_path}/${file_name}.merged.discarded.fastq 
gzip ${temp_path}/${file_name}.merged.trimmed.fastq

########### process stuff here #########

bwa mem -M ${bwa_ind} -t ${numCores} ${temp_path}/${file_name}.merged.trimmed.fastq.gz > ${temp_path}/${file_name}.merged.trimmed.sam

############# post alignment processing ##########################

########## filter SAM #################

#256 removes duplicates

samtools view -h -S -q 1 -F 256 ${temp_path}/${file_name}.merged.trimmed.sam -o ${temp_path}/${file_name}.merged.trimmed.qF.sam -U ${temp_path}/${file_name}.merged.trimmed.qF.unaligned.sam

num_aligned=$( grep -v @ ${temp_path}/${file_name}.merged.trimmed.qF.sam |  wc -l | awk '{print $1}' )

#keep a record of the bam file after alignment, as bwa randomly splits reads that can align to multiple places 

samtools view -S -b ${temp_path}/${file_name}.merged.trimmed.sam > ${temp_path}/${file_name}.merged.trimmed.bam

rm -f ${temp_path}/${file_name}.merged.trimmed.sam

samtools view -S -b ${temp_path}/${file_name}.merged.trimmed.qF.unaligned.sam > ${temp_path}/${file_name}.merged.trimmed.qF.unaligned.bam

rm -f ${temp_path}/${file_name}.merged.trimmed.qF.unaligned.sam 

################### construct the calmd flag #################

#remove previous indeces
rm -f ${alignment_fasta_file}.fai

samtools calmd ${temp_path}/${file_name}.merged.trimmed.qF.sam ${alignment_fasta_file} > ${temp_path}/${file_name}.merged.trimmed.qF.md.sam

rm -f ${temp_path}/${file_name}.merged.trimmed.qF.sam

################### filter for homology#################################

python ${script_dir}/filterSAMV3.py -infoF ${info_file} -inputF ${temp_path}/${file_name}.merged.trimmed.qF.md.sam -passF ${temp_path}/${file_name}.merged.trimmed.qF.md_filtered.sam -failF ${temp_path}/${file_name}.merged.trimmed.qF.md_failed.sam -s ${temp_path}/${file_name}.homologyFilter.stats

num_pass_homology_filter=$(awk 'NR>1{print $1}' ${temp_path}/${file_name}.homologyFilter.stats)
num_fail_homology_filter=$(awk 'NR>1{print $2}' ${temp_path}/${file_name}.homologyFilter.stats)

samtools view -S -b ${temp_path}/${file_name}.merged.trimmed.qF.md.sam  > ${temp_path}/${file_name}.merged.trimmed.qF.md.bam

rm -f ${temp_path}/${file_name}.merged.trimmed.qF.md.sam


samtools view -S -b ${temp_path}/${file_name}.merged.trimmed.qF.md_failed.sam > ${temp_path}/${file_name}.merged.trimmed.qF.md_failed.bam

rm -f ${temp_path}/${file_name}.merged.trimmed.qF.md_failed.sam


######## get the tables, with barcodes and without #####################

#simply add up the counts
python ${script_dir}/compileReads_no_dup_sep.py -sam ${temp_path}/${file_name}.merged.trimmed.qF.md_filtered.sam -fasta ${alignment_fasta_file}

samtools view -S -b ${temp_path}/${file_name}.merged.trimmed.qF.md_filtered.sam > ${temp_path}/${file_name}.merged.trimmed.qF.md_filtered.bam

rm -f ${temp_path}/${file_name}.merged.trimmed.qF.md_filtered.sam

################# create counts per oligo table ###################################################

#sum up all the barcode counts per row across all columns
awk 'NR>1{sum=0; for(i=2; i<=NF; i++) {sum += $i}; print $1"\t"sum }' ${temp_path}/${file_name}.merged.trimmed.qF.md_filtered.counts.table > ${temp_path}/${file_name}.merged.trimmed.qF.md_filtered.counts.nobars.table.parsed

############ create a barcodes only count table  #########################

#do i-1 to ignore the oligo column
awk 'NR==1{for (i=2;i<=NF;i++) names[i-1]=$i} NR>1{for (i=2;i<=NF;i++) count[i-1]+=$i} END{for (i=1;i<=NF-1;i++){print names[i]"\t"count[i]}}' ${temp_path}/${file_name}.merged.trimmed.qF.md_filtered.counts.table > ${temp_path}/${file_name}.merged.trimmed.qF.md_filtered.counts.barcodes_only.table.parsed

rm -f ${temp_path}/${file_name}.merged.trimmed.qF.md_filtered.counts.table

########### output summary statistics into a file ##########################

echo -e "${num_initial_reads}\t${num_trimmed_kept}\t${num_trimmed_discarded}\t${num_trimmed_discarded_non_n}\t${num_aligned}\t${num_pass_homology_filter}\t${num_fail_homology_filter}" > ${temp_path}/${file_name}_alignment_stats.txt
