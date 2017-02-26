

############################# for LT6535 #####################################################
#make design matrix
LT6535_index_no_myc$lev <- paste(LT6535_index_no_myc$Marker, LT6535_index_no_myc$Time, sep=".")
lev <- unique(LT6535_index_no_myc$lev)
fac <- factor(LT6535_index_no_myc$lev, levels = lev)
design <- model.matrix(~0+fac)
colnames(design) <- lev

####################### calculate weights for DEG analysis ##########################################

setwd("D:/data_analysis/microarrays/FLAG_NmrA microarray/JP_raw_data_27-05-2015/DEG analysis")

#To calculate weights, use background corrected and normalised data BEFORE filtering for expressed genes
#For duplicate correlation, need to remove control probes and reorder so rep probes next to each other

#calculate weights for replicate probes (on the same array)
library("statmod")

#calcuating weights based on individual rep spots on array
#Need to reorder probes so replicate probes are together to calculate weights for replicated probes
data.q.order <- data.q.no_controls[order(data.q.no_controls$genes$SystematicName)]

dupcor<-duplicateCorrelation(data.q.order,design=design, ndup=3,spacing=1,block=NULL)

dupcor$consensus

#calcualte array weights; low weights indicate outliers
#analysis is based on empirical reproducibility of gene expression measures from replicate arrays
#returns array level weights

arrayw <- arrayWeights(data.q.no_controls,design=design) 

#arrange arrays in order of line, timepoint and rep: order is in design matrix for graphs
#TdTom, FLAG/ 8h, 16h, 24h/ rep1, rep2, rep3
list <- c(2,8,14,6,12,18,4,10,16,1,7,13,5,11,17,3,9,15)
matrix_arrayw <- as.matrix(arrayw)
arrayw2 <- matrix_arrayw[list,]
arrayw2

#draw barplot of array weights
par(mar=c(8,4,1.09,0.56)) #change margin size of graph to fit bar labels (resizes bottom, left, top, right margin)
barplot(arrayw2,ylab="Weight", col="White",las=2,main="arrayWeights", ylim=c(0,4), 
        names.arg=c("TdTom_8h_rep1","TdTom_8h_rep2","TdTom_8h_rep3","TdTom_16h_rep1","TdTom_16h_rep2","TdTom_16h_rep3",
                    "TdTom_24h_rep1","TdTom_24h_rep2","TdTom_24h_rep3",
                    "FLAG_8h_rep1","FLAG_8h_rep2","FLAG_8h_rep3","FLAG_16h_rep1","FLAG_16h_rep2","FLAG_16h_rep3",
                    "FLAG_24h_rep1","FLAG_24h_rep2","FLAG_24h_rep3"))
abline(h=1,lwd=1,lty=2)


############################# find DEGS ####################################################

#combine weights from duplicateCorrelation and arrayWeights to normalised data 
#lmFit will average over replicate spots - no need to use averep
fitx <- lmFit(data.q.order, design=design, ndups=3, spacing=1, correlation=dupcor$consensus, weights=arrayw)


# Contrast to find genes that differ over the time course between FLAG and TdTOM
cont.dif <- makeContrasts(
  Dif8h = "FLAG2.8h-TdTom.8h",
  Dif16h = "FLAG2.16h-TdTom.16h",
  Dif24h = "FLAG2.24h-TdTom.24h",
  Dif8to16h=(FLAG2.16h-FLAG2.8h)-(TdTom.16h-TdTom.8h),
  Dif16to24h=(FLAG2.24h-FLAG2.16h)-(TdTom.24h-TdTom.16h),
  Dif8to24h=(FLAG2.24h-FLAG2.8h)-(TdTom.24h-TdTom.8h),
  levels=design)

fit <- contrasts.fit(fitx, cont.dif) #This looks for all the contrasts, to check a specific one use e.g. cont.dif[,"Dif24h"]

#can also filter out 25% of lowest intensity probe sets
#use trend=TRUE for eBayes to accomodate mean-variance trend
keep <- fit$Amean>quantile(fit$Amean,prob=0.25)
summary(keep)
fit2 <- eBayes(fit[keep,],trend=TRUE)
results2 <- decideTests(fit2,p.value=0.005)
topTable(fit2,adjust="BH", p.value=0.005)
summary(results2)

############### Extract DEGS to use for Pearson correlation ############################

# Contrast to find genes that differ over the time course between FLAG and TdTOM
cont.dif2 <- makeContrasts(
  TdTom_Dif8to16h=(TdTom.16h-TdTom.8h),
  TdTom_Dif16to24h=(TdTom.24h-TdTom.16h),
  TdTom_Dif8to24h=(TdTom.24h-TdTom.8h),
  FLAG_Dif8to16h=(FLAG2.16h-FLAG2.8h),
  FLAG_Dif16to24h=(FLAG2.24h-FLAG2.16h),
  FLAG_Dif8to24h=(FLAG2.24h-FLAG2.8h),
  levels=design)

fit2 <- contrasts.fit(fitx, cont.dif2) #This looks for all the contrasts, to check a specific one use e.g. cont.dif[,"Dif24h"]

#can also filter out 25% of lowest intensity probe sets
#use trend=TRUE for eBayes to accomodate mean-variance trend
fit3 <- eBayes(fit2[keep,],trend=TRUE)
results3 <- decideTests(fit3)
topTable(fit3,adjust="BH", p.value=0.005)
summary(results3)

#show top DEGs: n=Inf to get whole list, cut off at p value 0.05 and sort by adj P val
#TdTom and FLAG lines and order from most positive to most negative)

#for DEGS in seperate lines, use fitx=fit2, cont.diff=cont.dif2
#dif = comparison you want to make (e.g., TdTom_Dif8to16h)

DEG_table <- function(cont.dif,dif){
  fit_dif <- contrasts.fit(fitx, cont.dif[,dif])
  fit_dif <- eBayes(fit_dif[keep,], trend=TRUE)
  table <- topTable(fit_dif, n=Inf, adjust="BH", p.value=0.005, sort.by="p", resort.by="p")
  return(table)
}

#after making all contrast separately, make a unique list of genes to use to subset data for pearsons correlation

DEG_List <- function(cont_1,cont_2,cont_3){
  if (missing(cont_3)){
    List_1 <- cont_1$SystematicName
    List_2 <- cont_2$SystematicName
    unique_DEGs <- unique(c(List_1,List_2))
    
  }else{
    List_1 <- cont_1$SystematicName
    List_2 <- cont_2$SystematicName
    List_3 <- cont_3$SystematicName
    unique_DEGs <- unique(c(List_1,List_2,List_3))
  }
  return(unique_DEGs)
}

write.table(FLAG_DEG_List,file="D:/data_analysis/microarrays/FLAG_NmrA microarray/DEG analysis/FLAG-NMRA DEG List.txt",
            quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(TdTom_DEG_List,file="D:/data_analysis/microarrays/FLAG_NmrA microarray/DEG analysis/TdTom DEG List.txt",
            quote=FALSE,row.names=FALSE, col.names=FALSE)


####################### make seperate table for FLAG and TdTom and subset DEGs ######################################

normalised_data <- read.csv("D:/data_analysis/microarrays/FLAG_NmrA microarray/JP_raw_data_27-05-2015/JP_raw_data_27-05-2015/LT6535 normalised data mean expression.csv")

TdTom_normalised <- normalised_data[,1:4]
FLAG_normalised <- normalised_data[,c(1,5:7)]

#subset DEGs

TdTom_normalised_DEGs <- TdTom_normalised[TdTom_normalised$SystematicName %in% TdTom_DEG_List,]
rownames(TdTom_normalised_DEGs) <- NULL

FLAG_normalised_DEGs <- FLAG_normalised[FLAG_normalised$SystematicName %in% FLAG_DEG_List,]
rownames(FLAG_normalised_DEGs) <- NULL

write.csv(FLAG_normalised_DEGs,file="D:/data_analysis/microarrays/FLAG_NmrA microarray/DEG analysis/FLAG-NMRA normalised DEGs.csv",
          row.names=FALSE)

write.csv(TdTom_normalised_DEGs,file="D:/data_analysis/microarrays/FLAG_NmrA microarray/DEG analysis/TdTom normalised DEGs.csv",
          row.names=FALSE)

#Use last two files to generate Pearson correlation table #

#########################look at the contrasts seperately###############################

#8h
fit_Dif8h <- contrasts.fit(fitx,cont.dif[,"Dif8h"])
fit_Dif8h  <- eBayes(fit_Dif8h[keep,],trend=TRUE)
topTable_8h <- topTable(fit_Dif8h, n=Inf, adjust="BH", p.value=0.005, sort.by="logFC", resort.by="logFC")

results_Dif8h <- decideTests(fit_Dif8h,p.value=0.05)
summary(results_Dif8h)



######### write tables to file ##################
setwd("D:/data_analysis/microarrays/FLAG_NmrA microarray/DEG analysis")

write.csv(topTable_8hr,file="Results_6545_Dif8h.csv",)
write.csv(topTable_16hr,file="Results_6545_Dif16h.csv",)
write.csv(topTable_24hr,file="Results_6545_Dif24h.csv",)

write.csv(topTable_8to16hr,file="Results_6545_Dif8to16h.csv",)
write.csv(topTable_16to24hr,file="Results_6545_Dif16to24h.csv",)
write.csv(topTable_8to24hr,file="Results_6545_Dif8to24h.csv",)


##########################################################################################
loopy <- colnames(cont.dif)
for (i in loopy){
  print(i) 
  fit4 <- contrasts.fit(fit, cont.dif[,as.vector(i)])
  FDR.2 <- p.adjust(fit4$table$PValue, method="BH") # Could use qvalue from library(qvalue)
  if (sum(FDR.2 < 0.05) > 0){
    print("writing table")
    print(sum(FDR.2 < 0.05))
    write.table(topTags(fit.2, n=sum(FDR.2 < 0.05)), paste(i,".txt", sep=""), sep="\t", quote = FALSE)
  } else {
    print("no DEGs in this comparison")  
  }
}



?topTable
