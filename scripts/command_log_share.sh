######################### set pre-variables & download ##########################

#change directory to downloaded git repository 
initial_data_dir=/Users/jxue/Documents/Dustin_Project/CMS-GWAS/Scripts/github/MPRAu
header_name=mprau_novaseq_12_15_18
#numCores to specify for alignment
numCores=1

alignment_fasta_file=${initial_data_dir}/data/GWASrewritepos_CMS_alignment_file.fasta
info_file=${initial_data_dir}/data/GWASrewritepos_CMS_alignment_file_align.table
bwa_ind=${initial_data_dir}/data/GWASrewritepos_CMS_alignment_file_align

#directory where the final output files will be stored
out_dir=${initial_data_dir}/analysis

#temp_dir are intermediate files created that can be deleted later, useful for debugging purposes
temp_dir=${initial_data_dir}/temp

script_dir=${initial_data_dir}/scripts

#qsub template file, modify if using other submission system
qsub_template_file=${script_dir}/runTasksParallelShare.sh

#make index file
bwa index ${alignment_fasta_file} -p ${bwa_ind}

#raw file
#first column is file_path, second column is lane number, third column is read number (read 1 or read 2), third column is p7 index, fourth column is p5 index
initial_file_info=${initial_data_dir}/data/mprau_novaseq_12_15_18_read_meta_file.txt

mkdir -p ${temp_dir}
mkdir -p ${out_dir}

#grab MPRAu files only from the master file list, and keep read 1 only (read 2 is more noisy, look at the fastqc reports)

#create reference file for creating the DESeq2 count/conditions files later
awk '{if($3=="HEK293"){print $NF"\t"$3"_"$5} else{print $NF"\t"$3}}' ${initial_file_info} > ${temp_dir}/${header_name}_cell_type_info.txt

initial_conditions_file=${temp_dir}/${header_name}_cell_type_info.txt

#download read 1 files, second field is read num
awk -v temp_dir="${temp_dir}" -v header_name="${header_name}" '$2==1{print temp_dir"/"header_name"_"$NF".fastq.gz\t"$1}' ${initial_file_info} > ${temp_dir}/${header_name}_download_info.txt
bash ${script_dir}/download_general_submit.sh ${temp_dir}/${header_name}_download_info.txt ${temp_dir} ${qsub_template_file} ${script_dir}

############# alignment ##################################

awk -v header_name="${header_name}" ${initial_file_info} '{print header_name"_"$NF}' > ${temp_dir}/${header_name}_file_name_list_fastq_ids.txt 

bash ${script_dir}/mprau_align_submit.sh ${temp_dir}/${header_name}_file_name_list_fastq_ids.txt ${bwa_ind} ${temp_dir} ${header_name} ${alignment_fasta_file} ${info_file} ${numCores} ${qsub_template_file} ${script_dir}

############# make the count table, to be used later for DESeq2 ##################

out_header_name=${header_name}_nobar

deseq_count_file_names=${temp_dir}/deseq_${out_header_name}_file_names.txt

#first column is id for file, second column is link to actual count file, third column contains the deseq condition (i.e. cell type)

file_suffix="merged.trimmed.qF.md_filtered.counts.nobars.table.parsed"
awk -v temp_dir="${temp_dir}" -v header_name="${header_name}" -v file_suffix="${file_suffix}" '{print $1"\t"temp_dir"/"header_name"_"$1"."file_suffix"\t"$2}' ${initial_conditions_file} > ${deseq_count_file_names}

bash ${script_dir}/deseq_merge_count_data_mprau.sh ${deseq_count_file_names} ${out_header_name} ${out_dir} ${temp_dir}

######## (optional) make a table for the hexamer barcode counts ########

out_header_name=${header_name}_bar

deseq_count_file_names=${temp_dir}/deseq_${out_header_name}_file_names.txt

#first column is id for file, second column is link to actual count file, third column contains the deseq condition (i.e. cell type)

file_suffix="merged.trimmed.qF.md_filtered.counts.barcodes_only.table.parsed"
awk -v temp_dir="${temp_dir}" -v header_name="${header_name}" -v file_suffix="${file_suffix}" '{print $1"\t"temp_dir"/"header_name"_"$1"."file_suffix"\t"$2}' ${initial_conditions_file} > ${deseq_count_file_names}

bash ${script_dir}/deseq_merge_count_data_mprau.sh ${deseq_count_file_names} ${out_header_name} ${out_dir} ${temp_dir}

### combine the alignment statistics files ##########

rm -f ${out_dir}/${header_name}_all_alignment_stats.txt
#put in header
echo -e "Sample\tNum_Initial_Reads\tNum_Trimmed_Kept\tNum_Trimmed_Discarded\tNum_Trimmed_Discarded_Non_N\tNum_Aligned\tNum_Pass_Homology_Filter\tNum_Fail_Homology_Filter" > ${out_dir}/${header_name}_all_alignment_stats.txt
while read -a array
do
    file_name="${array[0]}"
    awk -v file_name=${file_name} '{print file_name"\t"$0}' ${temp_dir}/${file_name}_alignment_stats.txt >> ${out_dir}/${header_name}_all_alignment_stats.txt

done<${out_dir}/${header_name}_file_name_list_fastq_ids.txt

