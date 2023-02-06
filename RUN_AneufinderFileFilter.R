###--------------- AneufinderFileFilter ----------------###

#----Thomas van Ravesteyn
#----Kops Group
#----Hubrecht Institute


###--------------- Load required packages --------------###
###-----------------------------------------------------###

library(AneuFinder)
library(colorspace)

#Set source with file that contains AneufinderFileFilter functions
source("FUNC_AneufinderFileFilter.R")

###---------------Settings before start----------------###
###----------------------------------------------------###

#Set  input directory (end with "/"), folder should contain aneufinder output, 
#each sample in separate subdirectory
input_dir <- "Input/CC3A/"

#Set model from which output should be analyzed ("dnaCopy" / "hmm" / "edivisive") ! case-sensitive
selected.model <- "dnaCopy"

#Set project name, will be used as output folder (e.g. "My_project")
project.name <- "CC3A"


###------------ Filter QC control ------------###
###-------------------------------------------###

#set filter on/off (TRUE/FALSE), if TRUE check values

#set minimum total read count per single cell
filter.total.read.count <- T
min.total.read.count <- 25000

#set minimum amount of chromosome segments per single cell
filter.num.segments <- FALSE
min.num.segments <- 1

#set maximum level of spikiness
filter.spikiness <- T
max.spikiness <- 0.25

#set minimum level of bhattacharyya distance
filter.bhattacharyya <- F
min.bhattacharyya <- 0.7


###----------- Filter copy number ------------###
###-------------------------------------------###

#set min number of weighted average copy number
filter.min.copy.number <- FALSE
min.copy.number <- 2.5

#set max number of weighted average copy number
filter.max.copy.number <- T
max.copy.number <- 3

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
set.measures <- F


###--------------Collect files to analyze--------------###
###----------------------------------------------------###

#puts the name of the files found in your working directory in a vector
sampleIDs <- list.files(path = input_dir, pattern = "")
#sampleIDs <- sampleIDs[3]


#generate lists of the data found in your working directory
hmmFiles <- list()
dnaFiles <- list()
edivisiveFiles <- list()

#load Files from Aneufinder output
for(sample in sampleIDs) {
  hmmFiles[[sample]] <- list.files(paste0(input_dir, sample, "/MODELS/method-HMM/"), full = TRUE)
  dnaFiles[[sample]] <- list.files(paste0(input_dir, sample, "/MODELS/method-dnacopy/"), full = TRUE)
  edivisiveFiles[[sample]] <- list.files(paste0(input_dir, sample, "/MODELS/method-edivisive/"), full = TRUE)
}

###--------------------Start script--------------------###
###----------------------------------------------------###

#Run to filter samples
AneufinderFileFilter(sampleIDs)

