file_headers_file=$1
bowtie_ind=$2
temp_dir=$3
header_name=$4
alignment_fasta_file=$5
info_file=$6
numCores=$7
qsub_template_file=$8
script_dir=$9

rm -f ${temp_dir}/${header_name}_mprau_align_submit.txt
#submit jobs
while read -a array
do

    file_header="${array[0]}"
    echo "bash ${script_dir}/mprau_align_iter.sh ${file_header} ${temp_dir} ${bowtie_ind} ${alignment_fasta_file} ${info_file} ${numCores}" >> ${temp_dir}/${header_name}_mprau_align_submit.txt

done<${file_headers_file}

#submit all jobs
totalMem=50G
totalNumTasks=$(wc -l ${temp_dir}/${header_name}_mprau_align_submit.txt  | awk '{print $1}')
totalTime=50:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g; s/#$ -pe smp numCores/#$ -pe smp ${numCores}/g; s/#$ -binding linear:numCores/#$ -binding linear:${numCores}/g" ${qsub_template_file} > ${temp_dir}/${header_name}_mprau_align_submit.sh
rm -f ${temp_dir}/${header_name}_mprau_align_submit.log
qsub -o ${temp_dir}/${header_name}_mprau_align_submit.log  ${temp_dir}/${header_name}_mprau_align_submit.sh ${temp_dir}/${header_name}_mprau_align_submit.txt
