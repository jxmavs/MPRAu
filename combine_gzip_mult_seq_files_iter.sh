subset_file_name_list=$1
header=$2
file_path=$3
file_ID=$4

rm -f ${file_path}/${file_ID}_${header}.fastq
while read -a array
do
    orig_file_path="${array[0]}"
    gunzip -c ${orig_file_path} >> ${file_path}/${file_ID}_${header}.fastq
done<${subset_file_name_list}

gzip ${file_path}/${file_ID}_${header}.fastq
rm -f ${subset_file_name_list}
