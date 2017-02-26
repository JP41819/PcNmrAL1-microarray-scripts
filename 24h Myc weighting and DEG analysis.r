
#make design matrix

design <- cbind(TdTom=c(1,0,1,0,1,0),FLAG=c(0,1,0,1,0,1))


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
list <- c(1,3,5,2,4,6)
matrix_arrayw <- as.matrix(arrayw)
arrayw2 <- matrix_arrayw[list,]
arrayw2

#draw barplot of array weights
par(mar=c(8,4,1.09,0.56)) #change margin size of graph to fit bar labels (resizes bottom, left, top, right margin)
barplot(arrayw2,ylab="Weight", col="White",las=2,main="arrayWeights", ylim=c(0,4), 
        names.arg=c("TdTom_rep1","TdTom_rep3","TdTom_rep4",
                    "FLAG_rep1","FLAG_rep3","FLAG_rep4"))
abline(h=1,lwd=1,lty=2)


############################# find DEGS ####################################################

#combine weights from duplicateCorrelation and arrayWeights to normalised data 
#lmFit will average over replicate spots - no need to use averep
fitx <- lmFit(data.q.order, design=design, ndups=3, spacing=1, correlation=dupcor$consensus, weights=arrayw)


# Contrast to find genes that differ over the time course between FLAG and TdTOM
cont.dif <- makeContrasts(Dif = "FLAG-TdTom",levels=design)

fit <- contrasts.fit(fitx, cont.dif) #This looks for all the contrasts, to check a specific one use e.g. cont.dif[,"Dif24h"]


#can also filter out 25% of lowest intensity probe sets
#use trend=TRUE for eBayes to accomodate mean-variance trend

keep <- fit$Amean>quantile(fit$Amean,prob=0.25)
summary(keep)
fit2 <- eBayes(fit[keep,],trend=TRUE)
results2 <- decideTests(fit2,p.value=0.005)
topTable(fit2,adjust="BH", p.value=0.005)
table <- topTable(fit2, n=Inf, adjust="BH", p.value=0.005, sort.by="p", resort.by="p")
summary(results2)

####################### save DEG list for analysis ######################################

write.csv(table, file="LT6535 24h Myc TdTom Vs FLAG-NMRA DEGs 3rep.csv")

