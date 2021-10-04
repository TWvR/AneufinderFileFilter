# AneufinderFileFilter

AneufinderFileFilter contains R scripts that allows for easy filtering of single-cell DNA sequencing output generated by the R package "Aneufinder". Filtering is based on QC data generated by Aneufinder and filtered model files can be reorganized to new directories. This script only reads the model files and doesn't make any changes to the orginal Aneufinder output and model files.

## Scripts

To use the code you need two included scripts:
1. __RUN_AneufinderFileFilter script__  
_Set input/output folders and choose correct settings_  

2. __FUNC_AneufinderFileFilter script__ 
_Contains the code that is used to filter the data_

## Required R packages

This script makes use of the following R packages:

1. __[Aneufinder](https://bioconductor.org/packages/release/bioc/html/AneuFinder.html)__
1. __[Colorspace](https://colorspace.r-forge.r-project.org/articles/colorspace.html)__

## Input / Output

__Input:__  
* Aneufinder Model files  

__Output:__

For each sample:
* .txt summary of used filter parameters
* directory with selected model files
* directory with excluded model files
* directory with perfect diploid model files
* .pdf genomewide karyotype plot based on selected/excluded files
* .pdf single karyotype plots based on selected/excluded files
* .pdf heterogeneity/aneuploidy plot based on selected/excluded files
* .csv with QC measurements for each file
* .csv with karyotype measurements for each chromosome
* .csv with karyotype measurements for whole genome

_many of the above are optional_


## Instructions for use

To get started I would advise users to make use of R studio and create a new project in R and name it 'AneufinderFileFilter'. Then download both scripts from GitHub and place these within the project folder of your new AneufinderFileFilter project.

In general there's no need to open and/or adjust the function script, this is only needed if you like to make adjustments to the code that performs the actual filtering  or the code by which the different plots are generated. You only need to make sure that the RUN script contains the correct source-path to the FUNC script.

The 'Run_AneufinderFileFilter'-script is subdivided in multiple sections to create a good overview of the different settings. Prior to each run you probably like to give your project a new name, assign the correct input folder and check the filtering and plotting settings.

After making all required adjustments, run the code line-by-line. The actual filtering is commenced at the end of the run script by running  __AneufinderFileFilter(sampleIDs)__. Soon thereafter you will be prompted to quickly check filter settings; if correct, please enter 'Y' to continue the script. 

## Options

__Available filtering options__

1. Filter Aneufinder model files generated via edivisive, dnaCopy or hmm.
1. Filter files based on total read count per cell, number of chromosome segements, spikiness and/or bhattacharyya distance.
1. Exclude model files with too high weighted average copy number
2. Exclude model files with a perfect diploid genome

__Obtain selected Aneufinder model files__

1. Copy selected model files to new folder
2. Copy model files from perfect diploid cells to new folder

__Plots__

1. PDF with summary statistics for included and excluded files
2. PDF with genomewide profile for selected files
3. PDF with single cell karyotype profiles for included or excluded files
4. PDF with heterogeneity profiles for selected model files
5. CSV file with measurement statistics for each model file

## Final comments

This script was built while being a novice in R programming. Hence, the actual code could have been propably much more efficient. Nevertheless I hope it can be used to your benefit.
If you have any questions or need help with running the script, please don't hesitate to send me a message.

Thomas van Ravesteyn

