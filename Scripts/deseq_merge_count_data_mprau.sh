#this script generates the count and condition files for DESEQ

#list files
deseq_count_file_names=$1
header_name=$2
file_path=$3
temp_path=$4
mkdir -p ${temp_path}

iter=1
headerString="Gene\t"

#create the design matrix

rm -f  ${file_path}/${header_name}_deseq_cond_file.txt
#rm -f ${file_path}/simulation_design_mat.txt
 
while read -a array
do
    fileName=${array[1]}
    echo ${fileName}
    #0 is reference, 1 is condition
    deseq_condition=${array[2]}
    headerIter=${array[0]}
    if [[ ${iter} -eq 1 ]];
        then
            #first column contains the column name, second column contains number of unique barcodes, thrid column contains number of non-unique barcodes
            awk '{print $1"\t"$2}' ${fileName} >  ${temp_path}/${header_name}_temp_count_file.txt
            
        else
            paste <(cat ${temp_path}/${header_name}_temp_count_file.txt) <(awk '{print $2}' ${fileName})  > ${temp_path}/${header_name}_temp_count_file_augment.txt
            #echo "paste <(cat ${temp_path}/${header_name}_temp_count_file.txt) <(awk 'NR>1{print $5}' ${fileName}) "
            mv ${temp_path}/${header_name}_temp_count_file_augment.txt ${temp_path}/${header_name}_temp_count_file.txt
    fi
    #build the header string
    headerString+="${headerIter}\t" 

    echo -e "${headerIter}\t${deseq_condition}" >> ${file_path}/${header_name}_deseq_cond_file.txt

    iter=$((${iter}+1)) 

done<${deseq_count_file_names}

awk '{print $1"\t"$3}' ${deseq_count_file_names} > ${file_path}/${header_name}_deseq_cond_file.txt
#add header to count table
echo -e ${headerString} > ${file_path}/${header_name}_deseq.txt
cat ${temp_path}/${header_name}_temp_count_file.txt >> ${file_path}/${header_name}_deseq.txt
rm -f ${temp_path}/${header_name}_temp_count_file.txt
