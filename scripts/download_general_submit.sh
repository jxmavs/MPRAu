#this script takes in an input file to download data from a source

header_name=$1
temp_dir=$2
qsub_template_file=$3
script_dir=$4

rm -f ${temp_dir}/${header_name}_download_general.txt

while read -a array
do
    
    file_out_path=${array[0]}
    link=${array[1]}
    
    echo "bash ${script_dir}/download_general_iter.sh ${file_out_path} ${link}" >> ${temp_dir}/${header_name}_download_general.txt 
done<${input_encode_data}

#submit all jobs
totalMem=10G
totalNumTasks=$(wc -l ${temp_dir}/${header_name}_download_general.txt | awk '{print $1}')
totalTime=50:00:00
numCores=1
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g; s/#$ -pe smp numCores/#$ -pe smp ${numCores}/g; s/#$ -binding linear:numCores/#$ -binding linear:${numCores}/g" ${qsub_template_file} > ${temp_dir}/${header_name}_download_general.sh

rm -f ${temp_dir}/${header_name}_download_general.log
qsub -o ${temp_dir}/${header_name}_download_general.log ${temp_dir}/${header_name}_download_general.sh ${temp_dir}/${header_name}_download_general.txt
