######################### alignment ##########################

data_dir=/idi/sabeti-scratch/jxue/dustin_project/CMS-GWAS
alignment_fasta_file=${data_dir}/GWASrewritepos_CMS_alignment_file.fasta


out_path=${data_dir}/Novaseq_12-15-18/analysis
bwa_ind=${data_dir}/GWASrewritepos_CMS_alignment_file_align
file_path=${data_dir}/Novaseq_12-15-18/data
temp_dir=/broad/hptmp/jxue/nova_12-15-18/mprau

script_dir=${data_dir}/scripts
header_name=mprau_novaseq_12_15_18

info_file=${data_dir}/GWASrewritepos_CMS_alignment_file_align.table
numCores=1

#make index file
bwa index ${alignment_fasta_file} -p ${data_dir}/GWASrewritepos_CMS_alignment_file_align

#raw file
#first column is file_path, second column is lane number, third column is read number (read 1 or read 2), third column is p7 index, fourth column is p5 index
initial_file_info=/idi/sabeti-scratch/jxue/novaseq_mprau_mpradel/files/novaseq_12_15_18/novaseq_12_15_18_file_info_local_parsed.txt

mkdir -p ${temp_dir}
mkdir -p ${out_path}
mkdir -p ${file_path}

#first column of ${file_path}/${header_name}_file_name_list.txt is the file path
#second column is the output file name
#grab MPRAu files only from the master file list, and keep read 1 only (read 2 is more noisy, look at the fastqc reports)
grep MPRAu ${initial_file_info} | awk '$4==1{print}'  > ${file_path}/${header_name}_file_name_list.txt

#make the master condition file/header file names referencing fastq

awk '{print $2}' ${file_path}/${header_name}_file_name_list.txt | sort | uniq > ${file_path}/${header_name}_file_name_list_unique_ids.txt

#add read info here
awk -v header_name=${header_name} '{print $2"_1_"header_name}' ${file_path}/${header_name}_file_name_list.txt | sort | uniq > ${file_path}/${header_name}_file_name_list_fastq_ids.txt

cat ${file_path}/${header_name}_file_name_list_unique_ids.txt | awk -F"_" '{if(NF==3){print $0"\t"$2} else if(NF==4 || NF==5){print $0"\t"$2"_"$3} else{print $0"\t"$2"_"$3"_"$4"_"$5"_"$6}}' > ${file_path}/${header_name}_cell_type_info.txt

master_condition_file_sorted=${file_path}/${header_name}_cell_type_info.txt

#gunzip, combine all lanes into one file (if files are split across lanes), then gzip the file

bash ${script_dir}/combine_gzip_mult_seq_files_submit.sh ${file_path}/${header_name}_file_name_list.txt ${header_name} ${temp_dir} ${file_path}


#ls -l ${file_path}/*.fastq | awk '{print $NF}' | awk -F"/" '{print $NF}' | awk -F"." '{print $1}' >  ${file_names_file}
#note that the output is to temp_path
bash ${script_dir}/mprau_align_submit.sh ${file_path}/${header_name}_file_name_list_fastq_ids.txt ${out_path} ${bwa_ind} ${file_path} ${temp_dir} ${header_name} ${alignment_fasta_file} ${info_file} ${numCores}


############# make the deseq count table ##################

header_name=mprau_novaseq_12_15_18
out_header_name=mprau_novaseq_12_15_18_nobar

deseq_count_file_names=${temp_dir}/deseq_${out_header_name}_file_names.txt
readNum=1

#add on count file names

rm -f ${deseq_count_file_names}_bar_temp_2
while read -a array
do

    file_ID="${array[0]}"
    echo -e "${temp_dir}/${file_ID}_${readNum}_${header_name}.merged.trimmed.qF.md_filtered.counts.nobars.table.parsed" >> ${deseq_count_file_names}_bar_temp_2
done<${master_condition_file_sorted}

#first column is id for file, second column is link to actual count file, thrid column contains the deseq condition (i.e. cell type)
paste <(awk '{print $1}' ${master_condition_file_sorted}) <(cat ${deseq_count_file_names}_bar_temp_2 ) > ${deseq_count_file_names}_bar_temp_3
paste <(cat ${deseq_count_file_names}_bar_temp_3) <( awk '{print $2}' ${master_condition_file_sorted} ) > ${deseq_count_file_names}

rm ${deseq_count_file_names}_bar_temp_2 ${deseq_count_file_names}_bar_temp_3



bash ${script_dir}/deseq_merge_count_data_mprau.sh ${deseq_count_file_names} ${out_header_name} ${file_path} ${temp_dir}


######## (optional) make a table for the hexamer barcode counts ########

header_name=mprau_novaseq_12_15_18
out_header_name=mprau_novaseq_12_15_18_bar

deseq_count_file_names=${temp_dir}/deseq_${out_header_name}_file_names.txt
readNum=1

rm -f ${deseq_count_file_names}_bar_temp_2
while read -a array
do

    file_ID="${array[0]}"
    echo -e "${temp_dir}/${file_ID}_${readNum}_${header_name}.merged.trimmed.qF.md_filtered.counts.barcodes_only.table.parsed" >> ${deseq_count_file_names}_bar_temp_2
done<${master_condition_file_sorted}

#first column is id for file, second column is link to actual count file, thrid column contains the deseq condition (i.e. cell type)
paste <(awk '{print $1}' ${master_condition_file_sorted}) <(cat ${deseq_count_file_names}_bar_temp_2 ) > ${deseq_count_file_names}_bar_temp_3
paste <(cat ${deseq_count_file_names}_bar_temp_3) <( awk '{print $2}' ${master_condition_file_sorted} ) > ${deseq_count_file_names}

rm ${deseq_count_file_names}_bar_temp_2 ${deseq_count_file_names}_bar_temp_3

bash ${script_dir}/deseq_merge_count_data_mprau.sh ${deseq_count_file_names} ${out_header_name} ${file_path} ${temp_dir}


### combine the alignment statistics files ##########

header_name=mprau_novaseq_12_15_18
rm -f ${temp_dir}/${header_name}_all_alignment_stats.txt
#put in header
echo -e "Sample\tNum_Initial_Reads\tNum_Trimmed_Kept\tNum_Trimmed_Discarded\tNum_Trimmed_Discarded_Non_N\tNum_Aligned\tNum_Pass_Homology_Filter\tNum_Fail_Homology_Filter" > ${temp_dir}/${header_name}_all_alignment_stats.txt
while read -a array
do
    file_name="${array[0]}"
    awk -v file_name=${file_name} '{print file_name"\t"$0}' ${temp_dir}/${file_name}_alignment_stats.txt >> ${temp_dir}/${header_name}_all_alignment_stats.txt

done<${file_path}/${header_name}_file_name_list_fastq_ids.txt

