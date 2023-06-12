# ReefFishStability
# R code for the analysis of the paper by Benedetti-Cecchi et al (in review): Marine protected areas promote stability of reef fish communities under climate warming.

# USAGE: 
To run the analyses, first create a folder called ReefFishStability in your home directory, then download the scripts  and copy the folders named
DataManagement, MasterR and Run_analysis into the folder ReefFishStability.

## Description of files
 - DataManagement: includes an R script with the same name that does data celaning and generates file master.fish.dat for alpha stability analysis. This is required if you want to start the analysis from scratch. Otherwise, the necessary files are provided in folder Run_analaysis.
 - MasterR: includes the required scripts and aggregated data to run the full analysis, including plots. These scripts are called when running the code in Run_analysis.
 - Run_analysis: includes the code to run the various steps of the analysis. Please, note that some steps were originally performed on an HPC; the code provided here uses foreach loops, but in some cases, especieally hte jackknife simulation , the analysis is slow. This is noted when appropriate and the resulting files are also prpvided so the slow parts of the code can be skipped.
 
All analyses have been done in R version 4.1.3.
