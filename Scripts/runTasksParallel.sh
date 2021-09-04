#!/bin/sh
#$ -cwd 
#$ -pe smp numCores
#$ -binding linear:numCores
#$ -l h_vmem=totalMem
#$ -l h_rt=totalTime
#$ -t 1-totalNumTasks
#$ -V
#$ -j y
source /broad/software/scripts/useuse
reuse -q GCC-5.2
reuse -q Python-2.7
reuse -q R-3.4

export LD_LIBRARY_PATH=/broad/software/free/Linux/redhat_6_x86_64/pkgs/r_3.4.0/lib64/R/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/libiconv_1.13.1/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_5.1.0/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_5.1.0/lib64:/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_4.9.0/lib64:/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_4.9.0/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/perl_5.8.9/lib:/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/oracle_instantclient-10.2.0.4.0/instantclient_10_2:/broad/uge/8.4.0/lib/lx-amd64:/broad/software/free/Linux/redhat_6_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/zlib_1.2.6/lib64:/broad/software/free/Linux/redhat_6_x86_64/pkgs/hdf5_1.8.9/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/graphviz_2.28.0/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/db_4.7.25/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/tcltk8.5.9/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/jdk1.8.0_92/jre/lib/amd64/server:/broad/software/free/Linux/redhat_6_x86_64/pkgs/jdk1.8.0_92/lib:/idi/sabeti-scratch/jxue/bin/jellyfish-1.1.11/lib:/idi/sabeti-scratch/jxue/bin/openmpi/lib:/idi/sabeti-scratch/jxue/bin/libncurses5/include:/idi/sabeti-scratch/jxue/bin/packages/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/libiconv_1.13.1/lib:/idi/sabeti-scratch/jxue/bin/lib

#export LD_PRELOAD=/broad/software/free/Linux/redhat_6_x86_64/pkgs/r_3.4.0/lib64/R/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/libiconv_1.13.1/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_5.1.0/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_5.1.0/lib64:/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_4.9.0/lib64:/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_4.9.0/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/perl_5.8.9/lib:/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/oracle_instantclient-10.2.0.4.0/instantclient_10_2:/broad/uge/8.4.0/lib/lx-amd64:/broad/software/free/Linux/redhat_6_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/zlib_1.2.6/lib64:/broad/software/free/Linux/redhat_6_x86_64/pkgs/hdf5_1.8.9/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/graphviz_2.28.0/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/db_4.7.25/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/tcltk8.5.9/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/jdk1.8.0_92/jre/lib/amd64/server:/broad/software/free/Linux/redhat_6_x86_64/pkgs/jdk1.8.0_92/lib:/idi/sabeti-scratch/jxue/bin/jellyfish-1.1.11/lib:/idi/sabeti-scratch/jxue/bin/openmpi/lib:/idi/sabeti-scratch/jxue/bin/libncurses5/include:/idi/sabeti-scratch/jxue/bin/packages/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/libiconv_1.13.1/lib:/idi/sabeti-scratch/jxue/bin/lib
#this script submits a task array and evaluates each of the tasks from a text file
tasklist=$1
runNum=${SGE_TASK_ID}

#run the command from the file
eval $(awk -v runNum=${runNum} 'NR==runNum{print}' ${tasklist})

