library(ggplot2)
library(DESeq2) #Differential Expression Analysis

#this script generates the processed DESeq2 log2FoldChange count files

###############################################

#change to downloaded directory
primary_file_path="/Users/jxue/Documents/Dustin_Project/CMS-GWAS/Scripts/github/MPRAu"

inputConditionsFileName=paste(primary_file_path, "Data", "mprau_novaseq_12_15_18_nobar_deseq_cond_file.txt", sep="/")
inputCountFileName=paste(primary_file_path, "Data", "mprau_novaseq_12_15_18_nobar_deseq.txt", sep="/")

#columns are:  treat type, ref type, comparison group, row subset group, name of comparison group 
compareInfoDFFileName=paste(primary_file_path, "Data", novaseq_cell_type_comparisons.txt", sep="/")

out_file_path=paste(primary_file_path, "Processed_Output", sep="/")
if(!file.exists(out_file_path)){
	dir.create(out_file_path)
}
out_prefix="mprau_novaseq_12_15_18"
out_file_full_path=paste(c(out_file_path, out_prefix), collapse="/")

#read in CMS/GWAS ID map
CMSGWASIdMapFileName=paste(primary_file_path, "Data", "GWASrewritepos_CMS_arrayassign", sep="/")

###################################################################################

countData=read.table(inputCountFileName, header=TRUE, row.names=1, check.names=FALSE)
colData=read.table(inputConditionsFileName, header=FALSE, row.names=1)
colnames(colData)="CellType"
compareInfoDF=read.table(compareInfoDFFileName, header=FALSE)
CMSGWASIdMap=read.table(CMSGWASIdMapFileName, header=TRUE)

colnames(CMSGWASIdMap)[5:6]=c("CMS","GWAS")

sortedIdInd=order(CMSGWASIdMap[,2])
CMSGWASIdMap=CMSGWASIdMap[sortedIdInd,]

rowNames=rownames(countData)
sortedCountDataInd=order(rowNames)

countData=countData[sortedCountDataInd,]

#parse into list of CMS/GWAS/CMS|GWAS ids

########## revise CMSGWASIdMap ###################

idNames=as.character(CMSGWASIdMap[,4])

#remove bashing for now, duplicate column names with a "/"
splitIDNames=strsplit(idNames, "/")
numRowsDF=sum(sapply(splitIDNames, length))

CMSGWASIdMapRevised=data.frame(matrix(nrow=numRowsDF, ncol=3))
colnames(CMSGWASIdMapRevised)=c("SeqNames", "CMS", "GWAS")

ind_row=0

for( ind in 1:length(splitIDNames)){

	rowNameIterVec=splitIDNames[[ind]]
	
	for( ind2 in 1:length(rowNameIterVec)){
	
		rowNameIter=rowNameIterVec[ind2]
		
		#add on info to new data frame
		ind_row=1+ind_row
		CMSGWASIdMapRevised[ind_row, 1]=rowNameIter
		CMSGWASIdMapRevised[ind_row, c(2,3)]=CMSGWASIdMap[ind ,c(5,6)]
	}
}

dim(CMSGWASIdMapRevised)
CMSGWASIdMapRevised=unique(CMSGWASIdMapRevised)
dim(CMSGWASIdMapRevised)

###############################

CMSGWASIdMapList=list()
CMSGWASIdMapList[["CMS"]]=as.character(CMSGWASIdMapRevised[which(CMSGWASIdMapRevised[,"CMS"]==1),1])
CMSGWASIdMapList[["GWAS"]]=as.character(CMSGWASIdMapRevised[which(CMSGWASIdMapRevised[,"GWAS"]==1),1])
CMSGWASIdMapList[["CMS|GWAS"]]=as.character(CMSGWASIdMapRevised[,1])

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

########### remove CXCL2, missing, random, include bashing #################

#remove "16_835850_A" annd "rs145747214" - these completely duplicated in everything, including gene/oligo position
#remove rs2008157_5'_End_2 and rs2008159_5'_End_2 -these were not in the original array but dustin put the sequences to complete ref/alt pairs
#remove esv2672388, rs3203358_3'_End_v2, and rs58290679_3'_End_v2 as initialMetaInfo has nothing in these entries
#rs3203358_3'_End_v2 and rs58290679_3'_End_v2 are essentially duplicate variants with same context as rs3203358_3'_End and rs58290679_3'_End
#esv2672388 is misgenerated, in particular esv2672388_alt has zeros in the count table

ignore_cols=c("bash", "CXCL2", "missing", "random", "16_835850_A", "rs145747214", "rs2008157_5'_End_.*2", "rs2008159_5'_End_.*2", "esv2672388", "rs3203358_3'_End_v2", "rs58290679_3'_End_v2") 


ignore_cols_string=paste(ignore_cols, collapse="|")
countDataParsed=countDataParsedPre[which(!grepl(ignore_cols_string, newRowNames)),]

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
	
	
	##### remove rows with no ref/alt - this may be due to a GWAS oligo being a CMS oligo in another array #####
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

	#boolean to make plots or not
	plotDESEQ=FALSE

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
		
		##### added 9/12/21 - write out normalized counts & individual deseq2 result data for encode #####
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

