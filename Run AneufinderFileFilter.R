###---------------Load required packages--------------###
###---------------------------------------------------###

library(AneuFinder)
#library(ggrepel)
library(colorspace)

#Set source with file that contains AneufinderFileFilter functions
source("C:/Users/t.ravesteyn/OneDrive - Hubrecht Institute/R workspace/Aneufinder/20210716_AneufinderFileFilter_v4,3.R")
#source("C:/Users/thoma/OneDrive - Hubrecht Institute/R workspace/AneufinderFileFilter/AneufinderFileFilter_Func.R")


###---------------Settings before start----------------###
###----------------------------------------------------###

#Set the working directory, this should be the folder that was used to store Aneufinder output
setwd("C:/Users/t.ravesteyn/OneDrive - Hubrecht Institute/R workspace/Aneufinder/Output/20210810_test")
#setwd("C:/Users/thoma/OneDrive - Hubrecht Institute/R workspace/Aneufinder/Output/20210721_TRscFACS_8")

#Set model from which output should be analyzed ("dnaCopy" / "hmm" / "edivisive") ! case-sensitive
selected.model <- "edivisive"

#Set project name (e.g. "My_filter_20210101")
project.name <- "filtertest4"


###------------ Filter settings --------------###
###-------------------------------------------###

#set filter on/off (TRUE/FALSE), if TRUE check values

#set minimum total read count per single cell
filter.total.read.count <- TRUE
min.total.read.count <- 25000

#set minimum amount of chromosome segments per single cell
filter.num.segments <- FALSE
min.num.segments <- 1

#set maximum level of spikiness
filter.spikiness <- TRUE
max.spikiness <- 0.25

#set minimum level of bhattacharyya distance
filter.bhattacharyya <- TRUE
min.bhattacharyya <- 0.7

#set max number of weighted copy number
filter.copy.number <- TRUE
max.filter.copy.number <- 3.0

#set filter to remove cells with perfect diploid state, set nXchr to 1 for male cells and n Xchr to 2 for female cells
filter.diploid <- FALSE
nXchr <- 1


###--------------- Obtain cell subsets ----------------###
###----------------------------------------------------###

#Copy model files of selected cells and removed
set.move.model.files <- T

#Copy model files of diploid cells
set.move.diploid.model.files <- FALSE


###---------------- Plotting settings -----------------###
###----------------------------------------------------###

#plot genomewide profiles
set.plot.genomewide <- T

#plot heterogeneity profiles
set.plot.heterogeneity <- T

#plot single profiles
set.plot.singles <- T

#save measures as .csv
set.measures <- T


###--------------Collect files to analyze--------------###
###----------------------------------------------------###

#puts the name of the files found in your working directory in a vector
sampleIDs <- list.files(pattern = "")
#sampleIDs <- sampleIDs[1]


#generate lists of the data found in your working directory
hmmFiles <- list()
dnaFiles <- list()
edivisiveFiles <- list()

#load Files from Aneufinder output
for(sample in sampleIDs) {
  hmmFiles[[sample]] <- list.files(paste0(sample, "/MODELS/method-HMM/"), full = TRUE)
  dnaFiles[[sample]] <- list.files(paste0(sample, "/MODELS/method-dnacopy/"), full = TRUE)
  edivisiveFiles[[sample]] <- list.files(paste0(sample, "/MODELS/method-edivisive/"), full = TRUE)
}

###--------------------Start script--------------------###
###----------------------------------------------------###

#Run to filter samples
AneufinderFileFilter(sampleIDs)

