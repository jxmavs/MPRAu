The master script that runs everything to generate the count table from raw reads is in  Scripts/command_log_share.sh  - the analysis is broken up into components of alignment and generating count tables.  

I also output a final statistics table which gives me info on how many reads are filtered out, etc.

After generating the count data, to generate the processed data (DESeq2 fold changes/p-values), use Scripts/generate_processed_DESeq2_data.R

