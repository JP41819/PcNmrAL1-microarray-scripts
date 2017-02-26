source("http://www.bioconductor.org/biocLite.R")
biocLite()
biocLite("limma")
library("limma")

file_index<-read.csv("LT6535_microarray_files_index 3rep.csv")

##LOADING FILES AND BG CORRECTION

data <- read.maimages(file_index$file, source="agilent",columns=list(G = "gMedianSignal", Gb = "gBGMeanSignal"))


##################### Remove unwanted probes ########################################

data_probes<-read.csv("D:/data_analysis/microarrays/FLAG_NmrA microarray/JP_raw_data_27-05-2015/JP_raw_data_27-05-2015/data_probes.csv", header=TRUE)

#read in file with probes to be removed
control_probes <- read.table("D:/data_analysis/microarrays/FLAG_NmrA microarray/JP_raw_data_27-05-2015/JP_raw_data_27-05-2015/ctrl_probes.txt")

#subset from the list of probes on the microarray, the probes you want to keep
probes_to_subset<-data_probes[!data_probes$ProbeName %in% control_probes$V1,]

#subset the probes removed and check numbers add up
probes_removed <- data_probes[data_probes$ProbeName %in% control_probes$V1,]

#make a list of the spot number corresponding to probes you want to keep
probes_to_keep <- probes_to_subset[,1]

#remove unwanted probes from EListRaw, then use for normalisation
data_probes_removed <- data[probes_to_keep,]

write.csv(data_probes_removed,file="data_probes_removed.csv")
write.csv(probes_removed,file="probes_removed.csv")
write.csv(probes_to_subset,file="probes_to_subset.csv")

######################## background correct and normalise data #######################################

data.bg <- backgroundCorrect(data_probes_removed, method="normexp") #use wt.fun to flag control and low expressed probes?
data.q <- normalizeBetweenArrays(data.bg,method="quantile")

#keeps all non-control probes
data.q.no_controls<- data.q[data.q$genes$ControlType==0,] # use for weights calculations and DEG analysis

####################### for other analyses ###########################################################

#average out rep spots
data.q.no_controls.avg <- avereps(data.q.no_controls,ID=data.q.no_controls$genes$ProbeName)

#TURN DATA INTO DATAFRAME FOR FOLLOW-UP ANALYSIS
data.q.no_controls.avg.frame <- as.data.frame.EList(data.q.no_controls.avg, row.names = NULL, )

# Calculate means
#rename column names to sample names:
sample_names  <- read.csv("LT6535_microarray_files_index_no_file_extension 3rep.csv")
names(data.q.no_controls.avg.frame)<-sample_names[match(names(data.q.no_controls.avg.frame),sample_names[,"file"]), "sample"]

library(plyr)
names(data.q.no_controls.avg.frame) [5] <- "SystematicName"
rownames(data.q.no_controls.avg.frame) <- data.q.no_controls.avg.frame$SystematicName
data.q.no_controls.avg.frame <- data.q.no_controls.avg.frame[,-(1:5)] #remove unwanted columns


#calculate means over reps samples
#Use file names to subset normalised data 

data.q.no_controls.avg.frame$TdTom <- rowMeans(subset(data.q.no_controls.avg.frame, select = c("TdTom_rep1", "TdTom_rep3","TdTom_rep4")), na.rm = TRUE)
data.q.no_controls.avg.frame$FLAG <- rowMeans(subset(data.q.no_controls.avg.frame, select = c("FLAG_rep1", "FLAG_rep3", "FLAG_rep4")), na.rm = TRUE)

#reorder columns so reps are together
data.q.no_controls.avg.frame <- data.q.no_controls.avg.frame[,c(1,3,5,2,4,6,7,8)]
write.csv(data.q.no_controls.avg.frame,file="LT6535 24h myc normalised data log2 expression 3rep.csv")

# SUBSET MEANS TABLE, ORIGINAL TABLE
data.corrected.means <- subset(data.q.no_controls.avg.frame[,7:8])

write.csv(data.corrected.means,file="LT6535 24h myc normalised data mean log2 expression 3rep.csv")

######################### Plot graphs to see if normalisation is ok ######################################
#PCA
raw_data<-data$E
pr_6535<-prcomp(t(raw_data)) #transpose data so doesn't plot all genes
plot(pr_6535)
summary(pr_6535)

#plot PCA

plot(pr_6535$x[,1],pr_6535$x[,2],pch=18, col=file_index$colour1,main="PCA by line") #1=black=FLAG, 2=red=TdTom
legend(1.75e+06,4e+05, c("TdTom", "FLAG"), cex=0.8, 
       col=c("red","black"), pch=18,); # make legend

#following code labels points, but not sure how to resize font
#text(pr_6535,text(pr_6535$x[,1],pr_6535$x[,2],labels=Samples_6535$Code),pos=3,offset=0.25,size=0.5)

#2=red=TdTom_8h, 3=green=TdTom_16h, 4=blue=TdTom_24h
#6=magenta=FLAG_8h, 1=black=FLAG_16h, 5=cyan=FLAG_24h







