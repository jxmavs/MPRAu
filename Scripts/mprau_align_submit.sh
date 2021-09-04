file_names_file=$1
out_path=$2
bowtie_ind=$3
file_path=$4
temp_dir=$5
header_name=$6
alignment_fasta_file=$7
info_file=$8
numCores=$9
#############################################################################################################

mkdir -p ${temp_dir}

#how much do we trim off beginning/end, remove barcode and adapter sequences
#cutAMTEnd removes CGTCAAGCGGCCAGTT (16 bp) adapter sequence from end
#cutAMTBegin removes barcode (20bp barcode sequence) + adapter TCTAGAGGTTCGTCGACGCGATCGCAGGAGCCGCAGTG (38 bp long) 
script_dir=/idi/sabeti-scratch/jxue/dustin_project/CMS-GWAS/scripts

rm -f ${temp_dir}/all_runs_${header_name}.txt
#submit jobs
while read -a array
do

    file_name="${array[0]}"
    echo "bash ${script_dir}/mprau_align_iter.sh ${file_name} ${file_path} ${temp_dir} ${out_path} ${bowtie_ind} ${alignment_fasta_file} ${info_file} ${numCores}" >> ${temp_dir}/all_runs_${header_name}.txt

done<${file_names_file}

#submit all jobs
del_dir=/idi/sabeti-scratch/jxue/ape_project/deletions_project
totalMem=50G
totalNumTasks=$(wc -l ${temp_dir}/all_runs_${header_name}.txt  | awk '{print $1}')
totalTime=50:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g; s/#$ -pe smp numCores/#$ -pe smp ${numCores}/g; s/#$ -binding linear:numCores/#$ -binding linear:${numCores}/g" ${del_dir}/runTasksParallel.sh > ${temp_dir}/mprau_align_submit_${header_name}.sh
rm -f ${temp_dir}/mprau_align_submit_${header_name}.log
qsub -o ${temp_dir}/mprau_align_submit_${header_name}.log  ${temp_dir}/mprau_align_submit_${header_name}.sh ${temp_dir}/all_runs_${header_name}.txt
