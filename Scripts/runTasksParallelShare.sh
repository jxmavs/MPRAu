#!/bin/sh
#$ -cwd 
#$ -pe smp numCores
#$ -binding linear:numCores
#$ -l h_vmem=totalMem
#$ -l h_rt=totalTime
#$ -t 1-totalNumTasks
#$ -V
#$ -j y

#this script submits a task array and evaluates each of the tasks from a text file
tasklist=$1
runNum=${SGE_TASK_ID}

#run the command from the file
eval $(awk -v runNum=${runNum} 'NR==runNum{print}' ${tasklist})

