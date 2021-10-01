# AneufinderFileFilter

AneufinderFileFilter contains R scripts that allows for easy filtering of single-cell DNA sequencing output generated by the R package "Aneufinder". Filtering is based on QC data generated by Aneufinder and filtered model files can be reorganized to new directories. 

## Scripts

To use the code you need two included scripts:
1. __Run_AneufinderFileFilter script:__ choose correct settings for filtering and possible output.
1. __The AneufinderFileFilter_Func script:__ contains the function that is used to filter the data

## Required R packages

This script makes use of the Aneufinder and colorspace R packages. So if you haven't installed them yet, please search for these in bioconductor and follow instructions to 

## Instructions for use

To get started I would advise users to make use of R studio and create a new project in R and name it 'AneufinderFileFilter'. Then download both scripts from GitHub and place these within the project folder of your new AneufinderFileFilter project.

In general there's no need to open and/or adjust the function script, this only needed if you like to make adjustments to the code that performs the actual filtering  or the code by which the different plots are generated.

The 'Run_AneufinderFileFilter'-script is subdivided in multiple sections to create a good overview of the different settings. First you like to give your project a name and set the model-type to match the previously used settings in Anuefinder (edivisive/dnaCopy/hmm).



