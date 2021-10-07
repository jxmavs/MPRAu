library(ggplot2)
library(DESeq2) #Differential Expression Analysis

#this script generates the processed DESeq2 log2FoldChange count files

##################### set initial variables ##########################

#change to downloaded directory
primary_file_path="/Users/jxue/Documents/Dustin_Project/CMS-GWAS/Scripts/github/MPRAu"

inputConditionsFileName=paste(primary_file_path, "data", "mprau_novaseq_12_15_18_nobar_deseq_cond_file.txt", sep="/")
inputCountFileName=paste(primary_file_path, "data", "mprau_novaseq_12_15_18_nobar_deseq.txt", sep="/")

#columns are:  treat type, ref type, comparison group, row subset group, name of comparison group 
compareInfoDfFileName=paste(primary_file_path, "data", "novaseq_cell_type_comparisons.txt", sep="/")

out_file_path=paste(primary_file_path, "Processed_Output", sep="/")
if(!file.exists(out_file_path)){
	dir.create(out_file_path)
}
out_prefix="mprau_novaseq_12_15_18"
out_file_full_path=paste(c(out_file_path, out_prefix), collapse="/")

#read in master metafile
masterMetaDfFileName=paste(primary_file_path, "data", "mprau_novaseq_12_15_18_master_metafile_encode.txt", sep="/")

##################### read tables/set initial lists##################################################

countData=read.table(inputCountFileName, header=TRUE, row.names=1, check.names=FALSE)
colData=read.table(inputConditionsFileName, header=FALSE, row.names=1)
colnames(colData)="CellType"
compareInfoDF=read.table(compareInfoDfFileName, header=FALSE)
masterMetaDf=read.table(masterMetaDfFileName, header=TRUE)

rowNames=rownames(countData)
sortedCountDataInd=order(rowNames)

countData=countData[sortedCountDataInd,]

#parse into list of CMS/GWAS/CMS|GWAS ids

CMSGWASIdMapList=list()
CMSGWASIdMapList[["CMS"]]=as.character(masterMetaDf[which(masterMetaDf[,"array"]=="CMS"), "oligo_id"])
CMSGWASIdMapList[["GWAS"]]=as.character(masterMetaDf[which(masterMetaDf[,"array"]=="GWAS"), "oligo_id"])
CMSGWASIdMapList[["CMS|GWAS"]]=as.character(masterMetaDf[which(masterMetaDf[,"array"]=="CMS_GWAS"), "oligo_id"])

########## create new count table by parsing the names ###################
rowNames=rownames(countData)

#remove bashing for now, duplicate column names with a "/"
splitRowNames=strsplit(rowNames, "/")
numRowsDF=sum(sapply(splitRowNames, length))

countDataParsedPre=data.frame(matrix(nrow=numRowsDF, ncol=ncol(countData)))
colnames(countDataParsedPre)=colnames(countData)

newRowNames=rep(NA, numRowsDF)

ind_row=0

for( ind in 1:length(splitRowNames)){

	rowNameIterVec=splitRowNames[[ind]]
	
	for ( ind2 in 1:length(rowNameIterVec)){
	
		rowNameIter=rowNameIterVec[ind2]
		
		#add on info to new data frame
		ind_row=1+ind_row
		countDataParsedPre[ind_row,]=countData[ind,]
		newRowNames[ind_row]=rowNameIter
	}
}

rownames(countDataParsedPre)=newRowNames

countDataParsed=countDataParsedPre[which(!grepl("bash", newRowNames)),]

sortedCountDataInd=order(rowNames)
countDataParsed=countDataParsed[sortedCountDataInd,]

###################################################

#have a seperate countData for each comparison with alignment to a different - this is necessary since the GWAS and CMS array are not pooled

numComparisonSets=max(compareInfoDF[,3])
countDataColNames=colnames(countDataParsed)
countDataRowNames=rownames(countDataParsed)

countDataAllSets=list()
colDataAllSets=list()

for(ind in 1:numComparisonSets){
	
	compareInfoDFSubset=compareInfoDF[which(compareInfoDF[,3]==ind),]
	
	#ID for pulling out specific rows from count data
	rowID=unique(as.character(compareInfoDFSubset[,4]))
	IDs=CMSGWASIdMapList[[rowID]]

	treat_and_nontreat_cols=unique(c(as.character(compareInfoDFSubset[,1]),as.character(compareInfoDFSubset[,2])))
	
	#get col data names
	col_subset_ind=rownames(colData)[which(colData[,1] %in% treat_and_nontreat_cols)]
	row_subset_ind=which(countDataRowNames %in% IDs)
	
	countDataSubset=countDataParsed[row_subset_ind, col_subset_ind]
	colDataSubset=colData[which(colData[,1]%in%treat_and_nontreat_cols),,drop=FALSE]
	
	
	##### remove rows with no ref/alt corresponding pair - this may be due to a GWAS oligo being a CMS oligo in another array #####
	rowNames=rownames(countDataSubset)
	SeqNames=gsub("_alt|_ref", "", rowNames)
	
	#create map of seq names to row names
	
	rowNamesTable=data.frame(SeqNames, rowNames)
	
	seqNamesTableCount=table(SeqNames)
	seqNamesRemove=names(seqNamesTableCount[which(seqNamesTableCount==1)])
	
	rowNamesRemove=as.character(rowNamesTable[which(rowNamesTable[,1]%in%seqNamesRemove),2])

	print(rowNamesRemove)
	
	countDataSubset=countDataSubset[which(!(rowNames%in% rowNamesRemove)),]

	rowNamesSubset=rownames(countDataSubset)
	sortedCountDataSubsetInd=order(rowNamesSubset)
	countDataSubset=countDataSubset[sortedCountDataSubsetInd,]
	
	countDataAllSets[[ind]]=countDataSubset
	colDataAllSets[[ind]]=colDataSubset
}

#*****************************************  using DESEQ2 to calculate FC/Skew ********************************************************#


allDESEQResults=list()

#boolean to make plots or not
plotDESEQ=FALSE
for(ind in 1:numComparisonSets){
	#ind=1
	countDataSubset=countDataAllSets[[ind]]
	compareInfoDFSubset=compareInfoDF[which(compareInfoDF[,3]==ind),]
	#change - to . for DESEQ purposes
	compareInfoDFSubset[,2]=gsub("-",".", compareInfoDFSubset[,2])
	
	colDataSubset=colDataAllSets[[ind]]
	
	ref_ind=which(grepl("ref", rownames(countDataSubset)))
	alt_ind=which(grepl("alt", rownames(countDataSubset)))
	
	#human ref is reference
	Ref_Count_Data=countDataSubset[ref_ind,]
	Alt_Count_Data=countDataSubset[alt_ind,]

	colnames(Ref_Count_Data)=paste(colnames(Ref_Count_Data), "ref", sep="_")
	colnames(Alt_Count_Data)=paste(colnames(Alt_Count_Data), "alt", sep="_")

	rowNames=rownames(Ref_Count_Data)
	#ID=strsplit(as.character(rowNames),  "_")
	SeqNames=gsub("_ref", "", rowNames)

	#check if Alt seq names are the same
	altRowNames=rownames(Alt_Count_Data)
	#altID=strsplit(as.character(altRowNames),  "_")
	altSeqNames=gsub("_alt", "", altRowNames)
	
	print(all(altSeqNames==SeqNames))
	print(altSeqNames[which(!(altSeqNames %in% SeqNames))])
	
	rownames(Ref_Count_Data)=SeqNames
	rownames(Alt_Count_Data)=SeqNames

	#combine columns
	revisedCountData=data.frame(Ref_Count_Data, Alt_Count_Data)

	#revise colData, add Ref/Alt info 
	revisedColData=rbind(colDataSubset, colDataSubset)
	rownames(revisedColData)[1:(nrow(revisedColData)/2)]=colnames(Ref_Count_Data)
	rownames(revisedColData)[(nrow(revisedColData)/2+1):nrow(revisedColData)]=colnames(Alt_Count_Data)
	
	#change - to .
	
	rownames(revisedColData)=gsub("-",".", rownames(revisedColData))
	revisedColData[,"CellType"]=gsub("-",".", revisedColData[,"CellType"])
	
	revisedColData$Ref_Alt=c( rep("Ref", nrow(revisedColData)/2 ), rep("Alt", nrow(revisedColData)/2) )

	revisedColData[,"Ref_Alt"]=factor(revisedColData[,"Ref_Alt"])
	revisedColData[,"CellType"]=factor(revisedColData[,"CellType"])

	#loop and run DESEQ across all cell type comparisons 
	deseqDataTableRevisedAll=data.frame()

	for(i in 1:nrow(compareInfoDFSubset)){
		
		treat_level=as.character(compareInfoDFSubset[i,1])
		base_level=as.character(compareInfoDFSubset[i,2])
	
		revisedColData[,"Ref_Alt"]=relevel(revisedColData[,"Ref_Alt"], "Ref")
		revisedColData[,"CellType"]=relevel(revisedColData[,"CellType"], base_level)
		
		deseqData=DESeqDataSetFromMatrix(countData = revisedCountData, colData = revisedColData, design = ~ Ref_Alt+CellType+Ref_Alt:CellType)

		#betaPrior should automatically be false for interaction terms
		deseqDataResults <- DESeq(deseqData, betaPrior=FALSE, fitType="local")
		
		##############################################################################################################
		
		#I may need to relevel here
		resultsNames(deseqDataResults)
		#"Ref_AltAlt.RNA_DNARNA"
		deseqOutputSkew=results(deseqDataResults, name=paste(c("Ref_AltAlt.CellType", treat_level), collapse=""))
		#contrast=c("RNA_DNA", "RNA", "DNA")
		deseqOutputRef=results(deseqDataResults, contrast=c("CellType", treat_level, base_level))
		#list(c("RNA_DNA_RNA_vs_DNA","Ref_AltAlt.RNA_DNARNA"))
		deseqOutputAlt=results(deseqDataResults, contrast=list(c(paste(c("CellType", treat_level, "vs", base_level), collapse="_"), paste(c("Ref_AltAlt.CellType", treat_level), collapse="") )) )

		#Skew

		deseqOutputDataFrame=data.frame(deseqOutputSkew)
		OriginalName=rownames(deseqOutputDataFrame)
		rownames(deseqOutputDataFrame)=c()

		deseqDataTableRevisedSkew=data.frame(SeqName=OriginalName, deseqOutputDataFrame)
		colnames(deseqDataTableRevisedSkew)[2:ncol(deseqDataTableRevisedSkew)]=paste(colnames(deseqDataTableRevisedSkew)[2:ncol(deseqDataTableRevisedSkew)], "Skew", treat_level, sep="_")

		#Ref 

		deseqOutputDataFrame=data.frame(deseqOutputRef)
		OriginalName=rownames(deseqOutputDataFrame)
		rownames(deseqOutputDataFrame)=c()

		deseqDataTableRevisedRef=data.frame(SeqName=OriginalName, deseqOutputDataFrame)
		colnames(deseqDataTableRevisedRef)[2:ncol(deseqDataTableRevisedRef)]=paste(colnames(deseqDataTableRevisedRef)[2:ncol(deseqDataTableRevisedRef)], "Ref", treat_level, sep="_")

		#Alt

		deseqOutputDataFrame=data.frame(deseqOutputAlt)
		OriginalName=rownames(deseqOutputDataFrame)
		rownames(deseqOutputDataFrame)=c()

		deseqDataTableRevisedAlt=data.frame(SeqName=OriginalName, deseqOutputDataFrame)
		colnames(deseqDataTableRevisedAlt)[2:ncol(deseqDataTableRevisedAlt)]=paste(colnames(deseqDataTableRevisedAlt)[2:ncol(deseqDataTableRevisedAlt)], "Alt", treat_level, sep="_")

		deseqDataTableRevised_1=merge(deseqDataTableRevisedSkew, deseqDataTableRevisedRef, by="SeqName")
		deseqDataTableRevised=merge(deseqDataTableRevised_1, deseqDataTableRevisedAlt, by="SeqName")

		if(nrow(deseqDataTableRevisedAll)==0){
			deseqDataTableRevisedAll=deseqDataTableRevised
		} else{
			deseqDataTableRevisedAll=merge(deseqDataTableRevisedAll, deseqDataTableRevised, by="SeqName")
		}
	
		####################### some prelim DESEQ plots ##########################
	
		if(plotDESEQ){
	
			plotName=paste(c("treat", treat_level, "base", base_level,"Skew"), collapse="_")
			plotTitle=plotName
			plotOutFileName=paste(c(out_file_full_path,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")
			pdf(plotOutFileName)
			par(cex.axis=1, cex.lab=1.5, cex.main=1, cex.sub=1)
			plotMA(deseqOutputSkew, alpha=0.01, ylab="Skew", xlab="Mean of Normalized Counts")
			dev.off()

			plotName=paste(c("treat", treat_level, "base", base_level,"Ref_Activity"), collapse="_")
			plotTitle=plotName
			plotOutFileName=paste(c(out_file_full_path,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")
			pdf(plotOutFileName)
			par(cex.axis=1, cex.lab=1.5, cex.main=1, cex.sub=1)
			plotMA(deseqOutputRef, alpha=0.01, ylab="Ref Activity", xlab="Mean of Normalized Counts" )
			dev.off()

			plotName=paste(c("treat", treat_level, "base", base_level,"Alt_Activity"), collapse="_")
			plotTitle=plotName
			plotOutFileName=paste(c(out_file_full_path,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")
			pdf(plotOutFileName)
			par(cex.axis=1, cex.lab=1.5, cex.main=1, cex.sub=1)
			plotMA(deseqOutputAlt, alpha=0.01, ylab="Alt Activity", xlab="Mean of Normalized Counts" )
			dev.off()

			plotName=paste(c("treat", treat_level, "base", base_level,"disp_fit"), collapse="_")
			plotTitle=plotName
			plotOutFileName=paste(c(out_file_full_path,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")
			pdf(plotOutFileName)
			par(cex.axis=1, cex.lab=1.5, cex.main=1, cex.sub=1)
			plotDispEsts(deseqDataResults)
			dev.off()
			
		}
		
		#####  write out normalized counts & individual deseq2 result data for encode processing #####
		deseqDataResultsTemp=estimateSizeFactors(deseqData)
		deseqNormalizedCounts=counts(deseqDataResultsTemp, normalized=TRUE)
		
		deseqAvgNormalizedCountsOutputDf=data.frame(matrix(nrow=nrow(deseqNormalizedCounts), ncol=5))
		colnames(deseqAvgNormalizedCountsOutputDf)=c("SeqName", paste(treat_level, "RNA_avg_Ref", sep="_"), paste(treat_level, "RNA_avg_Alt", sep="_"), paste(base_level, "avg_Ref", sep="_"), paste(base_level, "avg_Alt", sep="_"))
		
		deseqAvgNormalizedCountsOutputDf[,1]=rownames(deseqNormalizedCounts)
		deseqAvgNormalizedCountsOutputDf[,2]=apply(deseqNormalizedCounts[,which(grepl(paste("MPRAu", treat_level, ".*_ref", sep="_"), colnames(deseqNormalizedCounts)))], 1, mean, na.rm=T)
		deseqAvgNormalizedCountsOutputDf[,3]=apply(deseqNormalizedCounts[,which(grepl(paste("MPRAu", treat_level, ".*_alt", sep="_"), colnames(deseqNormalizedCounts)))], 1, mean, na.rm=T)
		deseqAvgNormalizedCountsOutputDf[,4]=apply(deseqNormalizedCounts[,which(grepl(paste("MPRAu", base_level, ".*_ref", sep="_"), colnames(deseqNormalizedCounts)))], 1, mean, na.rm=T)
		deseqAvgNormalizedCountsOutputDf[,5]=apply(deseqNormalizedCounts[,which(grepl(paste("MPRAu", base_level, ".*_alt", sep="_"), colnames(deseqNormalizedCounts)))], 1, mean, na.rm=T)
		
		out_file_full_path=paste(c(out_file_path, out_prefix), collapse="/")
		textOutFileName=paste(c(out_file_full_path, treat_level, "deseq_results_avg_normalized_counts.txt"), collapse="_")
		write.table(deseqAvgNormalizedCountsOutputDf, quote=F, row.names=F, textOutFileName)
	
		out_file_full_path=paste(c(out_file_path, out_prefix), collapse="/")
		textOutFileName=paste(c(out_file_full_path, treat_level, "deseq_results.txt"), collapse="_")
		write.table(deseqDataTableRevised, quote=F, row.names=F, textOutFileName)
		
		#############################
	}
	
	allDESEQResults[[ind]]=deseqDataTableRevisedAll
	
}


### combine table into one file, plot GWAS/CMS shared oligos concordance ###

#pool CMS and GWAS first
deseq_df_CMS=allDESEQResults[[1]]
colnames_deseq=colnames(deseq_df_CMS)
#keep padj and log2FoldChange
keep_col_ind=which(grepl("padj*|log2FoldChange*|lfcSE*|pvalue*", colnames_deseq))
deseq_df_CMS=deseq_df_CMS[,c(1, keep_col_ind)]

deseq_df_GWAS=allDESEQResults[[2]]
colnames_deseq=colnames(deseq_df_GWAS)
#keep padj and log2FoldChange
keep_col_ind=which(grepl("padj*|log2FoldChange*|lfcSE*|pvalue*", colnames_deseq))
deseq_df_GWAS=deseq_df_GWAS[,c(1, keep_col_ind)]

deseq_df_GWAS_CMS_overlap=merge(deseq_df_GWAS, deseq_df_CMS, by="SeqName")

## create the merged table, just get GWAS hits for any overlaps with CMS ####
deseq_df_CMS_nonoverlap_GWAS=deseq_df_CMS[which(!(deseq_df_CMS[,"SeqName"]%in%deseq_df_GWAS[,"SeqName"])),]

#change colnames of all

colnames(deseq_df_GWAS)=gsub("_GWAS", "" , colnames(deseq_df_GWAS))
colnames(deseq_df_CMS_nonoverlap_GWAS)=gsub("_CMS", "" , colnames(deseq_df_CMS_nonoverlap_GWAS))

deseq_df_HEK293FT=rbind(deseq_df_GWAS, deseq_df_CMS_nonoverlap_GWAS)

deseq_df_all=deseq_df_HEK293FT
#collapse to single data frame
for(ind in 3:length(allDESEQResults)){
	
	deseq_df_iter=allDESEQResults[[ind]]
	colnames_deseq=colnames(deseq_df_iter)
	#keep padj and log2FoldChange
	keep_col_ind=which(grepl("padj*|log2FoldChange*|lfcSE*|pvalue*", colnames_deseq))
	
	deseq_df_all=merge(deseq_df_all, deseq_df_iter[, c(1, keep_col_ind)], by="SeqName")
}

#write out deseq table only
out_file_full_path=paste(c(out_file_path, out_prefix), collapse="/")
textOutFileName=paste(c(out_file_full_path, "_deseq_results_all_table.txt"), collapse="")
write.table(deseq_df_all, quote=F, row.names=F, textOutFileName)

#write out deseq table without deletion bashing
not_del_ind=which(!grepl("del", as.character(deseq_df_all[,1])))
deseq_df_all_del_bash_removed=deseq_df_all[not_del_ind,]
out_file_full_path=paste(c(out_file_path, out_prefix), collapse="/")
textOutFileName=paste(c(out_file_full_path, "_deseq_results_all_table_no_del_bash.txt"), collapse="")
write.table(deseq_df_all_del_bash_removed, quote=F, row.names=F, textOutFileName)


#*********************** generate snv tiling results ***********************#

########### preprocess bashing #############################################

countDataParsed=countDataParsedPre

#sort again
rowNames=rownames(countDataParsed)
sortedCountDataInd=order(rowNames)
countDataParsed=countDataParsed[sortedCountDataInd,]

###################################################

numComparisonSets=max(compareInfoDF[,3])
countDataColNames=colnames(countDataParsed)
countDataRowNames=rownames(countDataParsed)

countDataRevisedBashAllSets=list()
colDataAllSets=list()

for(ind in 1:numComparisonSets){
	
	compareInfoDFSubset=compareInfoDF[which(compareInfoDF[,3]==ind),]
	
	#get info, if CMS, then there are no bashing oligos
	subsetInfo=unique(compareInfoDFSubset[,4])
	
	#ID for pulling out specific rows from count data
	rowID=unique(as.character(compareInfoDFSubset[,4]))
	IDs=CMSGWASIdMapList[[rowID]]

	treat_and_nontreat_cols=unique(c(as.character(compareInfoDFSubset[,1]),as.character(compareInfoDFSubset[,2])))
	
	#get col data names
	col_subset_ind=rownames(colData)[which(colData[,1] %in% treat_and_nontreat_cols)]
	row_subset_ind=which(countDataRowNames %in% IDs)
	
	countDataSubset=countDataParsed[row_subset_ind, col_subset_ind]
	colDataSubset=colData[which(colData[,1]%in%treat_and_nontreat_cols),,drop=FALSE]
	
	countDataRevisedBashAllSets[[ind]]=countDataSubset
	
	colDataAllSets[[ind]]=colDataSubset
}



################ use DESEQ2 to calculate FC/Skew (altered bashing version) ###############

out_file_full_path_snp_bash=paste(c(out_file_path, paste(c(out_prefix, "snv_tiling"), collapse="_")), collapse="/")

allDESEQResultsSnpBash=list()

#boolean to make plots or not
plotDESEQ=FALSE

for(ind in 1:numComparisonSets){

	countDataRevisedBashSubset=countDataRevisedBashAllSets[[ind]]
	compareInfoDFSubset=compareInfoDF[which(compareInfoDF[,3]==ind),]
	#change - to . for DESEQ purposes
	compareInfoDFSubset[,2]=gsub("-",".", compareInfoDFSubset[,2])
	
	#get info, if CMS, then there are no bashing oligos
	subsetInfo=unique(compareInfoDFSubset[,4])
	
	colDataSubset=colDataAllSets[[ind]]
	
	countDataRevisedBashSubsetRowNames=rownames(countDataRevisedBashSubset)
	
	#change - to .
	revisedColData=colDataSubset
	rownames(revisedColData)=gsub("-",".", rownames(revisedColData))
	revisedColData[,"CellType"]=gsub("-",".", revisedColData[,"CellType"])
	
	revisedColData[,"CellType"]=factor(revisedColData[,"CellType"])

	#loop and run DESEQ across all cell type comparisons 
	deseqDataTableRevisedAll=data.frame()

	colnames(countDataRevisedBashSubset)=gsub("-",".", colnames(countDataRevisedBashSubset))
	
	for(i in 1:nrow(compareInfoDFSubset)){

		treat_level=as.character(compareInfoDFSubset[i,1])
		base_level=as.character(compareInfoDFSubset[i,2])
	
		revisedColData[,"CellType"]=relevel(revisedColData[,"CellType"], base_level)

		deseqData=DESeqDataSetFromMatrix(countData = countDataRevisedBashSubset, colData = revisedColData, design = ~ CellType)

		#betaPrior should automatically be false for interaction terms
		deseqDataResults <- DESeq(deseqData, betaPrior=FALSE, fitType="local")
		
		resultsNames(deseqDataResults)

		deseqOutput=results(deseqDataResults, contrast=c("CellType", treat_level, base_level))

		deseqOutputDataFrame=data.frame(deseqOutput)
		OriginalName=rownames(deseqOutputDataFrame)
		rownames(deseqOutputDataFrame)=c()

		deseqDataTableRevised=data.frame(SeqName=OriginalName, deseqOutputDataFrame)

		colnames(deseqDataTableRevised)[2:ncol(deseqDataTableRevised)]=paste(colnames(deseqDataTableRevised)[2:ncol(deseqDataTableRevised)], treat_level, sep="_")

		if(nrow(deseqDataTableRevisedAll)==0){
			deseqDataTableRevisedAll=deseqDataTableRevised
		}
		else{
			deseqDataTableRevisedAll=merge(deseqDataTableRevisedAll, deseqDataTableRevised, by="SeqName")
		}
	
		####################### some prelim DESEQ plots ##########################
	
		if(plotDESEQ){

			plotName=paste(c("treat", treat_level, "base", base_level,"Activity"), collapse="_")
			plotTitle=plotName
			plotOutFileName=paste(c(out_file_full_path_snp_bash,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")
			pdf(plotOutFileName)
			par(cex.axis=1, cex.lab=1.5, cex.main=1, cex.sub=1)
			plotMA(deseqOutput, alpha=0.01, ylab="Activity", xlab="Mean of Normalized Counts" )
			dev.off()

			plotName=paste(c("treat", treat_level, "base", base_level,"disp_fit"), collapse="_")
			plotTitle=plotName
			plotOutFileName=paste(c(out_file_full_path_snp_bash,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")
			pdf(plotOutFileName)
			par(cex.axis=1, cex.lab=1.5, cex.main=1, cex.sub=1)
			plotDispEsts(deseqDataResults)
			dev.off()
			
		}
		
		####### write out table results for encode processing later ########
		textOutFileName=paste(c(out_file_full_path_snp_bash, treat_level, "deseq_results.txt"), collapse="_")
		write.table(deseqDataTableRevised, quote=F, row.names=F, textOutFileName)


	}
	
	allDESEQResultsSnpBash[[ind]]=deseqDataTableRevisedAll

}

### combine table into one file ###

#pool CMS and GWAS first
deseq_df_CMS=allDESEQResultsSnpBash[[1]]
colnames_deseq=colnames(deseq_df_CMS)
#keep padj and log2FoldChange
keep_col_ind=which(grepl("padj*|log2FoldChange*|lfcSE*|pvalue*", colnames_deseq))
deseq_df_CMS=deseq_df_CMS[,c(1, keep_col_ind)]

deseq_df_GWAS=allDESEQResultsSnpBash[[2]]
colnames_deseq=colnames(deseq_df_GWAS)
#keep padj and log2FoldChange
keep_col_ind=which(grepl("padj*|log2FoldChange*|lfcSE*|pvalue*", colnames_deseq))
deseq_df_GWAS=deseq_df_GWAS[,c(1, keep_col_ind)]

deseq_df_GWAS_CMS_overlap=merge(deseq_df_GWAS, deseq_df_CMS, by="SeqName")

## create the merged table, just get GWAS hits for any overlaps with CMS ####
deseq_df_CMS_nonoverlap_GWAS=deseq_df_CMS[which(!(deseq_df_CMS[,"SeqName"]%in%deseq_df_GWAS[,"SeqName"])),]

#change colnames of all

colnames(deseq_df_GWAS)=gsub("_GWAS", "" , colnames(deseq_df_GWAS))
colnames(deseq_df_CMS_nonoverlap_GWAS)=gsub("_CMS", "" , colnames(deseq_df_CMS_nonoverlap_GWAS))

deseq_df_HEK293FT=rbind(deseq_df_GWAS, deseq_df_CMS_nonoverlap_GWAS)

deseq_df_all=deseq_df_HEK293FT
#collapse to single data frame
for(ind in 3:length(allDESEQResultsSnpBash)){
	
	deseq_df_iter=allDESEQResultsSnpBash[[ind]]
	colnames_deseq=colnames(deseq_df_iter)
	#keep padj and log2FoldChange
	keep_col_ind=which(grepl("padj*|log2FoldChange*|lfcSE*|pvalue*", colnames_deseq))
	
	deseq_df_all=merge(deseq_df_all, deseq_df_iter[, c(1, keep_col_ind)], by="SeqName")
}


#write out full analysis
#textOutFileName=paste(c(out_file_full_path_snp_bash, "_deseq_results_all_table_full.txt"), collapse="")
#write.table(deseq_df_all, quote=F, row.names=F, textOutFileName)

#write out snv tiling only oligo results
orig_oligo_names=unique(as.character(allBashMetaDf[,"OrigSeqName"]))

bashed_orig_oligo_ind=which(deseq_df_all[,"SeqName"] %in% orig_oligo_names)

deseq_df_all_bashed_orig=deseq_df_all[bashed_orig_oligo_ind,]
deseq_df_all_bash_only=deseq_df_all[which(grepl("bash", as.character(deseq_df_all[,"SeqName"]))),]

deseq_df_all_bash_all=rbind(deseq_df_all_bashed_orig, deseq_df_all_bash_only)
textOutFileName=paste(c(out_file_full_path_snp_bash, "_deseq_results_all_table.txt"), collapse="")
write.table(deseq_df_all, quote=F, row.names=F, textOutFileName)

#******************* generate encode final processed data  *******************#

######### write out the bed files for encode ##########

cell_names=c("HEK293_CMS", "HEK293_GWAS", "HEPG2", "HMEC", "K562", "GM12878", "SKNSH")

cell_plasmid_map_encode=list()
cell_plasmid_map_encode[["GM12878"]]="PB_CMS.GWAS_4_11_Cycles"
cell_plasmid_map_encode[["K562"]]="PB_CMS.GWAS_4_11_Cycles"
cell_plasmid_map_encode[["HEPG2"]]="PB_CMS.GWAS_4_11_Cycles"
cell_plasmid_map_encode[["SKNSH"]]="PB_CMS.GWAS_4_11_Cycles"
cell_plasmid_map_encode[["HEK293_CMS"]]="PB_CMS_4_11_Cycles"
cell_plasmid_map_encode[["HEK293_GWAS"]]="PB_GWAS_3_11_Cycles"
cell_plasmid_map_encode[["HMEC"]]="PB_CMS.GWAS_4_11_Cycles"
cell_plasmid_map_encode[["HUVECS"]]="PB_CMS.GWAS_4_16.5_Cycles"

masterMetaDfVar=masterMetaDf[which(masterMetaDf[,"project"]=="variant"),]
masterMetaDfVar_parsed=unique(masterMetaDfVar[,c("mpra_variant_id", "chrom", "oligo_starts", "oligo_ends", "strand", "var_start", "var_end")])

##### make block information columns: start, end, length, count, and starts #####

start_info=strsplit(as.character(masterMetaDfVar_parsed[,"oligo_starts"]), ",")
end_info=strsplit(as.character(masterMetaDfVar_parsed[,"oligo_ends"]), ",")

block_additional_info_df=data.frame(matrix(nrow=length(start_info), ncol=5))
colnames(block_additional_info_df)=c("chromStart", "chromEnd", "block_count", "block_sizes", "block_starts")

for(ind in 1:length(start_info)){

	start_info_iter=start_info[[ind]]
	end_info_iter=end_info[[ind]]
	
	block_counts_all=length(start_info_iter)
	block_sizes_all=rep(NA, length(start_info_iter))
	block_starts_all=rep(NA, length(start_info_iter))
	
	for(ind2 in 1:length(start_info_iter)){
			
		block_sizes_all[ind2]=as.numeric(end_info_iter[ind2])-as.numeric(start_info_iter[ind2])+1
		block_starts_all[ind2]=as.numeric(start_info_iter[ind2])-1
		
	}
	
	block_sizes_all_str=paste(block_sizes_all, collapse=",")
	block_starts_all_str=paste(block_starts_all, collapse=",")
	
	block_additional_info_df[ind,1]=as.numeric(block_starts_all[1])
	block_additional_info_df[ind,2]=as.numeric(end_info_iter[length(end_info_iter)])
	block_additional_info_df[ind,3]=block_counts_all
	block_additional_info_df[ind,4]=block_sizes_all_str
	block_additional_info_df[ind,5]=block_starts_all_str

}

masterMetaDfVar_parsed=cbind(masterMetaDfVar_parsed, block_additional_info_df)
print(dim(masterMetaDfVar_parsed))
print(dim(block_additional_info_df))

##### merge with initial Variant Info #####

#add empty score column
masterMetaDfVar_parsed[,"score"]=rep("0", nrow(masterMetaDfVar_parsed))

#add on chr
masterMetaDfVar_parsed[,"chrom"]=paste("chr", masterMetaDfVar_parsed[,"chrom"], sep="")

#loop to write out bed file for each cell type
for(ind in 1:length(cell_names)){

	cell_name_iter=cell_names[ind]
	
	AvgNormalizedCountsDfFileName=paste(c(out_file_full_path, cell_name_iter, "deseq_results_avg_normalized_counts.txt"), collapse="_")
	avg_normalized_counts_df_iter=read.table(AvgNormalizedCountsDfFileName, header=T)
		
	output_df_1=merge(masterMetaDfVar_parsed, avg_normalized_counts_df_iter, by.x="mpra_variant_id", by.y="SeqName")
	colnames(output_df_1)[1]="SeqName"
	
	deseqDfFileName=paste(c(out_file_full_path, cell_name_iter, "deseq_results.txt"), collapse="_")
	deseq_df_iter=read.table(deseqDfFileName, header=T)
	
	output_df=merge(deseq_df_iter, output_df_1, by="SeqName")
	
	print(dim(output_df))
	
	l2fc_ref_col=paste("log2FoldChange_Ref", cell_name_iter, sep="_")
	pvalue_ref_col=paste("pvalue_Ref", cell_name_iter, sep="_")
	padj_ref_col=paste("padj_Ref", cell_name_iter, sep="_")
	
	log10_pvalue_ref_col=paste("log10_pvalue_Ref", cell_name_iter, sep="_")
	log10_padj_ref_col=paste("log10_padj_Ref", cell_name_iter, sep="_")
	
	output_df[,log10_pvalue_ref_col]=-log(output_df[,pvalue_ref_col], base=10)
	output_df[,log10_padj_ref_col]=-log(output_df[,padj_ref_col], base=10)
		
	normalized_ref_plasmid_count_col=paste(cell_plasmid_map_encode[[cell_name_iter]], "avg_Ref", sep="_")
	normalized_ref_rna_count_col=paste(cell_name_iter, "RNA_avg_Ref", sep="_")
	
	#create bed entry
	encode_output_df=output_df[,c("chrom", "chromStart", "chromEnd", "SeqName", "score", "strand", l2fc_ref_col, normalized_ref_plasmid_count_col, normalized_ref_rna_count_col, log10_pvalue_ref_col, log10_padj_ref_col, "block_count", "block_sizes", "block_starts")]
	
	#change colnames
	colnames(encode_output_df)=c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "log2FoldChange", "inputCount", "outputCount", "minusLog10Pvalue", "minusLog10Qvalue", "blockCount", "blockSizes", "blockStarts")

	#filter out _2 and NA in p-value fields 
	encode_output_df_parsed=encode_output_df[which(! (is.na(encode_output_df[,"chromStart"]) | is.na(encode_output_df[,"minusLog10Pvalue"])) ),]
	
	#for the duplicated ids
	check_dup_id=apply(encode_output_df_parsed[,c("chrom", "chromStart", "chromEnd", "strand", "blockStarts", "blockSizes")], 1, function(x){paste(x, collapse="_")})
	check_dup_id_table=table(check_dup_id)
	encode_output_df_parsed[,"check_dup_id"]=check_dup_id
	check_dup_ids_look=names(check_dup_id_table[which(check_dup_id_table>1)])

	encode_output_df_parsed_nonduplicated=encode_output_df_parsed[which(!(encode_output_df_parsed[,"check_dup_id"] %in% check_dup_ids_look)),]

	encode_output_df_parsed_duplicated=encode_output_df_parsed[which(encode_output_df_parsed[,"check_dup_id"] %in% check_dup_ids_look),]
	encode_output_df_parsed_duplicated=encode_output_df_parsed_duplicated[order(encode_output_df_parsed_duplicated[,"check_dup_id"]), ]

	encode_output_df_parsed_duplicated_revised=data.frame(matrix(nrow=length(check_dup_ids_look), ncol=ncol(encode_output_df_parsed)))
	colnames(encode_output_df_parsed_duplicated_revised)=c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "log2FoldChange", "inputCount", "outputCount", "minusLog10Pvalue", "minusLog10Qvalue", "blockCount", "blockSizes", "blockStarts", "check_dup_id")

	encode_output_df_parsed_final=rbind(encode_output_df_parsed_nonduplicated, encode_output_df_parsed_duplicated_revised)
	
	#take the one with the least significant signal to be conservative
	for(ind2 in 1:length(check_dup_ids_look)){
		check_dup_ids_look_iter=check_dup_ids_look[ind2]
		encode_output_df_parsed_dup_iter=encode_output_df_parsed[which(encode_output_df_parsed[,"check_dup_id"] %in% check_dup_ids_look_iter),]
		ind_keep=which.min(encode_output_df_parsed_dup_iter[,"minusLog10Pvalue"])
				
		encode_output_df_parsed_final[nrow(encode_output_df_parsed_nonduplicated)+ind2,]=encode_output_df_parsed_dup_iter[ind_keep,]
	}

	#remove check_dup_id field
	encode_output_df_parsed_final=encode_output_df_parsed_final[, 1:(ncol(encode_output_df_parsed_final)-1)]
	
	encode_output_df_parsed_final=encode_output_df_parsed_final[order(as.character(encode_output_df_parsed_final[,"chrom"]), encode_output_df_parsed_final[,"chromStart"], as.character(encode_output_df_parsed_final[,"strand"])), ]
	print(dim(encode_output_df_parsed_final))
		
	out_file_full_path=paste(c(out_file_path, out_prefix), collapse="/")
	textOutFileName=paste(c(out_file_full_path, cell_name_iter, "encode_output.bed"), collapse="_")
	write.table(encode_output_df_parsed_final, quote=F, row.names=F, col.names=F, sep="\t", file=textOutFileName)

}


########## write out in a combined file for oligo activity for encode ############

cell_names=c("HEK293_CMS", "HEK293_GWAS", "HEPG2", "HMEC", "K562", "GM12878", "SKNSH")

for(i in 1:length(cell_names)){

	cell_name_iter=cell_names[i]

	#read mpra variant & del tiling info
	textOutFileName=paste(c(out_file_full_path, cell_name_iter, "deseq_results.txt"), collapse="_")
	deseq_df_iter=read.table(textOutFileName, header=T)

	deseq_ref_col_names=colnames(deseq_df_iter)[which(grepl("padj_Ref|log2FoldChange_Ref|lfcSE_Ref|SeqName|pvalue_Ref", colnames(deseq_df_iter)))]
	deseq_alt_col_names=colnames(deseq_df_iter)[which(grepl("padj_Alt|log2FoldChange_Alt|lfcSE_Alt|SeqName|pvalue_Alt", colnames(deseq_df_iter)))]

	####### process variant results ############
	deseq_df_iter_var_only=deseq_df_iter[which(!grepl("del", deseq_df_iter[,"SeqName"])),]

	deseq_df_iter_var_only_ref=data.frame(SeqName_tag=paste(deseq_df_iter_var_only[,"SeqName"], "ref", sep="_"), deseq_df_iter_var_only[, deseq_ref_col_names])
	deseq_df_iter_var_only_alt=data.frame(SeqName_tag=paste(deseq_df_iter_var_only[,"SeqName"], "alt", sep="_"), deseq_df_iter_var_only[, deseq_alt_col_names])

	#change col names
	colnames(deseq_df_iter_var_only_ref)=gsub("_Ref", "", colnames(deseq_df_iter_var_only_ref))
	colnames(deseq_df_iter_var_only_alt)=gsub("_Alt", "", colnames(deseq_df_iter_var_only_alt))

	deseq_df_iter_var_only_revised=rbind(deseq_df_iter_var_only_ref, deseq_df_iter_var_only_alt)

	#merge with master file to get oligo_ids for variants
	masterMetaDf_var=masterMetaDf[which(masterMetaDf[,"project"]=="variant"),]
	masterMetaDf_var_with_tag=data.frame(SeqName_tag=paste(masterMetaDf_var[,"mpra_variant_id"], masterMetaDf_var[,"tag"], sep="_"), masterMetaDf_var)

	deseq_df_iter_var_only_revised_with_meta=merge(masterMetaDf_var_with_tag, deseq_df_iter_var_only_revised, by="SeqName_tag")

	#filter to keep only columns of interest
	deseq_cols_interest=colnames(deseq_df_iter_var_only_revised_with_meta)[which(grepl("padj|log2FoldChange|lfcSE|pvalue", colnames(deseq_df_iter_var_only_revised_with_meta)))]
	deseq_df_iter_var_only_final=deseq_df_iter_var_only_revised_with_meta[,c("oligo_id", deseq_cols_interest)]
	
	######## remaining preprocess ###################
	
	#CMS array does not have any bashing oligos
	
	if(cell_name_iter=="HEK293_CMS"){
	
		deseq_df_combined_final=merge(deseq_df_iter_var_only_final, masterMetaDf[,c("oligo_id", "oligo_group_id", "project")], by="oligo_id")
		dim(deseq_df_iter_var_only_final)
		dim(deseq_df_combined_final)
		
	} else{
	
		
		####### process del tiling results ##############
		deseq_df_iter_del_only=deseq_df_iter[which(grepl("del", deseq_df_iter[,"SeqName"])),]

		deseq_df_iter_del_only_ref=deseq_df_iter_del_only[, deseq_ref_col_names]
		deseq_df_iter_del_only_alt=deseq_df_iter_del_only[, deseq_alt_col_names]

		#change col names
		colnames(deseq_df_iter_del_only_ref)=gsub("_Ref", "", colnames(deseq_df_iter_del_only_ref))
		colnames(deseq_df_iter_del_only_alt)=gsub("_Alt", "", colnames(deseq_df_iter_del_only_alt))

		#add back ref/alt info into id
		deseq_df_iter_del_only_ref[,"SeqName"]=gsub("del", "ref_del", as.character(deseq_df_iter_del_only_ref[,"SeqName"]))
		deseq_df_iter_del_only_alt[,"SeqName"]=gsub("del", "alt_del", as.character(deseq_df_iter_del_only_alt[,"SeqName"]))

		deseq_df_iter_del_only_revised=rbind(deseq_df_iter_del_only_ref, deseq_df_iter_del_only_alt)
		colnames(deseq_df_iter_del_only_revised)[1]="oligo_id"

		#####  process snv tiling results  ######
		#read snv tiling analysis
		textOutFileName=paste(c(out_file_full_path, "snp_bash", cell_name_iter, "deseq_results.txt"), collapse="_")
		deseq_df_iter_snp_bash=read.table(textOutFileName, header=T)
		
		deseq_col_names=colnames(deseq_df_iter_snp_bash)[which(grepl("padj|log2FoldChange|lfcSE|SeqName|pvalue", colnames(deseq_df_iter_snp_bash)))]
		deseq_df_iter_snp_bash=deseq_df_iter_snp_bash[,deseq_col_names]
	
		#keep only snp tiling ids
		deseq_df_iter_bash_only=deseq_df_iter_snp_bash[which(grepl("bash", as.character(deseq_df_iter_snp_bash[,"SeqName"]))),]
		colnames(deseq_df_iter_bash_only)[1]="oligo_id"

		##### combine into one activity table ####

		deseq_df_combined=rbind(deseq_df_iter_var_only_final, deseq_df_iter_del_only_revised, deseq_df_iter_bash_only)
	
		#final columns: oligo ID, deseq results, oligo_group_id, project 
		deseq_df_combined_final=merge(deseq_df_combined, masterMetaDf[,c("oligo_id", "oligo_group_id", "project")], by="oligo_id")
	
		dim(deseq_df_combined)
		dim(deseq_df_combined_final)
	}
	
	###### write out table ######
	
	textOutFileName=paste(paste(c(out_file_full_path, "activity_file_encode", cell_name_iter), collapse="_"), ".txt", sep="")
	write.table(deseq_df_combined_final, quote=F, row.names=F, textOutFileName)

}

########### write out variant mprau results seperated by cell type ###########

masterMetaDf_var=masterMetaDf[which(masterMetaDf[,"project"]=="variant"),]
masterMetaDf_var_ids=unique(as.character(masterMetaDf_var[,"mpra_variant_id"]))

cell_names=c("HEK293_CMS", "HEK293_GWAS", "HEPG2", "HMEC", "K562", "GM12878", "SKNSH")

for(ind in 1:length(cell_names)){
	cell_name_iter=cell_names[ind]
	
	#read mpra variant & del tiling info
	textOutFileName=paste(c(out_file_full_path, cell_name_iter, "deseq_results.txt"), collapse="_")
	deseq_df_iter=read.table(textOutFileName, header=T)
	
	length(masterMetaDf_var_ids)
	deseq_df_var_final=deseq_df_iter[which(as.character(deseq_df_iter[,"SeqName"]) %in% masterMetaDf_var_ids), ]
	colnames(deseq_df_var_final)[1]="mpra_variant_id"
	dim(deseq_df_var_final)
	
	###### write out table ######
	deseq_col_names=colnames(deseq_df_var_final)[which(grepl("padj|log2FoldChange|lfcSE|SeqName|pvalue", colnames(deseq_df_var_final)))]
	textOutFileName=paste(paste(c(out_file_full_path, "variant_file_encode", cell_name_iter), collapse="_"), ".txt", sep="")
	write.table(deseq_df_var_final[,deseq_col_names], quote=F, row.names=F, textOutFileName)
}

