source("http://www.bioconductor.org/biocLite.R")
biocLite()
biocLite("limma")
library("limma")

setwd("D:/data_analysis/microarrays/FLAG_NmrA microarray/JP_raw_data_27-05-2015/JP_raw_data_27-05-2015")

file_index<-read.csv("JP_data_files_extract_file_names.csv")

##LOADING FILES AND BG CORRECTION
#Create index of all files releating to 6535 samples
LT6535_index<-file_index[file_index$Line=="6535",c(1,3,5,6)]
#Remove files relating to Myc
LT6535_index_no_myc<-subset(LT6535_index,Time!="Myc")
#Extract files names for LT6535 to read in and normalise
LT6535_file_names_no_myc<-as.vector(LT6535_index_no_myc[,1])
#Read files in for normalisation
data <- read.maimages(LT6535_file_names_no_myc, source="agilent",columns=list(G = "gMedianSignal", Gb = "gBGMeanSignal"))

##################### Remove unwanted probes ########################################

#after read data in to create EListRaw file, save as .csv file
#rearrange probes according to row and then col
#probes will now be ordered according to spot number to use to call individual spots
#only need to rearrange here is using spot correlation - not using here, so no need to do above

data_probes<-read.csv("data_probes.csv", header=TRUE)

#read in file with probes to be removed
control_probes <- read.table("ctrl_probes.txt")

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
#rename column names to sample names: load data frame to use for sorting
df_for_sorting_samples<-read.csv("D:/data_analysis/microarrays/FLAG_NmrA microarray/JP_raw_data_27-05-2015/JP_raw_data_27-05-2015/JP_data_files_txt_removed.csv")

names(data.q.no_controls.avg.frame)<-df_for_sorting_samples[match(names(data.q.no_controls.avg.frame),df_for_sorting_samples[,"Sample.Name"]), "Code"]
library(plyr)
names(data.q.no_controls.avg.frame) [5] <- "SystematicName"
rownames(data.q.no_controls.avg.frame) <- data.q.no_controls.avg.frame$SystematicName

data.q.no_controls.avg.frame <- data.q.no_controls.avg.frame[,-(1:5)] #remove unwanted columns

#calculate SE for each gene

TdTom_8hr_SE <- subset(data.q.no_controls.avg.frame, select = c("6535_TdTom_8h_rep1", "6535_TdTom_8h_rep2", "6535_TdTom_8h_rep3"))
TdTom_16hr_SE <- subset(data.q.no_controls.avg.frame, select = c("6535_TdTom_16h_rep1", "6535_TdTom_16h_rep2", "6535_TdTom_16h_rep3"))
TdTom_24hr_SE <- subset(data.q.no_controls.avg.frame, select = c("6535_TdTom_24h_rep1", "6535_TdTom_24h_rep2", "6535_TdTom_24h_rep3"))

FLAG_8hr_SE <- subset(data.q.no_controls.avg.frame, select = c("6535_FLAG2_8h_rep1", "6535_FLAG2_8h_rep2", "6535_FLAG2_8h_rep3"))
FLAG_16hr_SE <- subset(data.q.no_controls.avg.frame, select = c("6535_FLAG2_16h_rep1", "6535_FLAG2_16h_rep2", "6535_FLAG2_16h_rep3"))
FLAG_24hr_SE <- subset(data.q.no_controls.avg.frame, select = c("6535_FLAG2_24h_rep1", "6535_FLAG2_24h_rep2", "6535_FLAG2_24h_rep3"))


SE <- function(x){x/sqrt(3)}
TdTom_8hr_SE$means <- apply(TdTom_8hr_SE,1,mean)
TdTom_8hr_SE$SD <- apply(TdTom_8hr_SE,1,sd)
TdTom_8hr_SE$SE <- apply(TdTom_8hr_SE[5],2,SE)

TdTom_16hr_SE$means <- apply(TdTom_16hr_SE,1,mean)
TdTom_16hr_SE$SD <- apply(TdTom_16hr_SE,1,sd)
TdTom_16hr_SE$SE <- apply(TdTom_16hr_SE[5],2,SE)

TdTom_24hr_SE$means <- apply(TdTom_24hr_SE,1,mean)
TdTom_24hr_SE$SD <- apply(TdTom_24hr_SE,1,sd)
TdTom_24hr_SE$SE <- apply(TdTom_24hr_SE[5],2,SE)

FLAG_8hr_SE$means <- apply(FLAG_8hr_SE,1,mean)
FLAG_8hr_SE$SD <- apply(FLAG_8hr_SE,1,sd)
FLAG_8hr_SE$SE <- apply(FLAG_8hr_SE[5],2,SE)

FLAG_16hr_SE$means <- apply(FLAG_16hr_SE,1,mean)
FLAG_16hr_SE$SD <- apply(FLAG_16hr_SE,1,sd)
FLAG_16hr_SE$SE <- apply(FLAG_16hr_SE[5],2,SE)

FLAG_24hr_SE$means <- apply(FLAG_24hr_SE,1,mean)
FLAG_24hr_SE$SD <- apply(FLAG_24hr_SE,1,sd)
FLAG_24hr_SE$SE <- apply(FLAG_24hr_SE[5],2,SE)

SE_microarray <- cbind(TdTom_8hr_SE[,6],TdTom_16hr_SE[,6],TdTom_24hr_SE[,6],
                       FLAG_8hr_SE[,6],FLAG_16hr_SE[,6],FLAG_24hr_SE[,6])

rownames(SE_microarray) <- rownames(TdTom_8hr_SE)

colnames(SE_microarray) <- c("TdTom_8h","TdTom_16h","TdTom_24h","FLAG_8h","FLAG_16h","FLAG_24h")

write.csv(SE_microarray, file = "LT6535_SE_microarray.csv")

#calculate means over reps samples
#Use file names to subset normalised data 

data.q.no_controls.avg.frame$TdTom_8hr <- rowMeans(subset(data.q.no_controls.avg.frame, select = c("6535_TdTom_8h_rep1", "6535_TdTom_8h_rep2", "6535_TdTom_8h_rep3")), na.rm = TRUE)
data.q.no_controls.avg.frame$TdTom_16hr <- rowMeans(subset(data.q.no_controls.avg.frame, select = c("6535_TdTom_16h_rep1", "6535_TdTom_16h_rep2", "6535_TdTom_16h_rep3")), na.rm = TRUE)
data.q.no_controls.avg.frame$TdTom_24hr <- rowMeans(subset(data.q.no_controls.avg.frame, select = c("6535_TdTom_24h_rep1", "6535_TdTom_24h_rep2", "6535_TdTom_24h_rep3")), na.rm = TRUE)

data.q.no_controls.avg.frame$FLAG_8hr <- rowMeans(subset(data.q.no_controls.avg.frame, select = c("6535_FLAG2_8h_rep1", "6535_FLAG2_8h_rep2", "6535_FLAG2_8h_rep3")), na.rm = TRUE)
data.q.no_controls.avg.frame$FLAG_16h <- rowMeans(subset(data.q.no_controls.avg.frame, select = c("6535_FLAG2_16h_rep1", "6535_FLAG2_16h_rep2", "6535_FLAG2_16h_rep3")), na.rm = TRUE)
data.q.no_controls.avg.frame$FLAG_24hr <- rowMeans(subset(data.q.no_controls.avg.frame, select = c("6535_FLAG2_24h_rep1", "6535_FLAG2_24h_rep2", "6535_FLAG2_24h_rep3")), na.rm = TRUE)

# SUBSET MEANS TABLE, ORIGINAL TABLE
data.corrected <-subset(data.q.no_controls.avg.frame[,1:18])
data.corrected.means <- subset(data.q.no_controls.avg.frame[,19:24])

write.csv(data.corrected.means,file="LT6535 normalised data mean expression_2.csv")

######################### Plot graphs to see if normalisation is ok ######################################
par(mfrow=c(1,2))
boxplot(data.bg$E, las=2, col= c("pink", "pink","pink","pink","red","red","red","red",
                              "green","green","green", "green","magenta","magenta","magenta", "magenta"))
boxplot(data.q$E, las=2, col= c("pink", "pink","pink","pink","red","red","red","red",
                                "green","green","green", "green","magenta","magenta","magenta", "magenta"))







