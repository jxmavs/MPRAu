#file names list, first column is the sequence file location, second column is the file name to combine all the files into
file_names_file=$1
header=$2
temp_dir=$3
file_path=$4
#############################################################################################################

mkdir -p ${temp_dir}

script_dir=/idi/sabeti-scratch/jxue/dustin_project/CMS-GWAS/scripts

#get list of unique file names
awk '{print $2"\t"$4}' ${file_names_file} | sort | uniq > ${file_names_file}.unique.ids

rm -f ${temp_dir}/all_runs_${header}.txt
#submit jobs
while read -a array
do

    unique_file_id="${array[0]}"
    read_num="${array[1]}"
    
    awk -v unique_id=${unique_file_id} -v read_num="${read_num}" '($2==unique_id) && ($4==read_num){print}' ${file_names_file} > ${file_names_file}.subset.${unique_file_id}.${read_num}
    
    subset_file_name=${file_names_file}.subset.${unique_file_id}.${read_num}

    echo "bash ${script_dir}/combine_gzip_mult_seq_files_iter.sh ${subset_file_name} ${header} ${file_path} ${unique_file_id}_${read_num}" >> ${temp_dir}/all_runs_${header}.txt

done<${file_names_file}.unique.ids

#submit all jobs
del_dir=/idi/sabeti-scratch/jxue/ape_project/deletions_project
totalMem=5G
totalNumTasks=$(wc -l ${temp_dir}/all_runs_${header}.txt  | awk '{print $1}')
totalTime=50:00:00
numTasksAtOnceLimit=10
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g; s/#$ -tc numTasksAtOnceLimit/#$ -tc ${numTasksAtOnceLimit}/g" ${del_dir}/runTasksEvalLimit.sh > ${temp_dir}/combine_gzip_mult_seq_files_submit_${header}.sh
rm -f ${temp_dir}/combine_gzip_mult_seq_files_submit_${header}.log
qsub -o ${temp_dir}/combine_gzip_mult_seq_files_submit_${header}.log  ${temp_dir}/combine_gzip_mult_seq_files_submit_${header}.sh ${temp_dir}/all_runs_${header}.txt
