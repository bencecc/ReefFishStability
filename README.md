# ReefFishStability
## R code for the analysis of the paper by Benedetti-Cecchi et al (in review): Marine protected areas promote stability of reef fish communities under climate warming.
### All analyses were done in R version 4.1.3 using librariers:
brms, broom, codyn, datawizard, DHARMa, doMC, FD, fishualize, foreach, ggeffects, ggstatsplot, knitr, lemon, lme4, lmerTest, mFD, mgcv, modelr, performance, piecewiseSEM, semEff, sjPlot, tidybayes, tidyverse, tidymv, tidyr.

## USAGE: 
To run the analyses, first create a folder called ReefFishStability in your home directory, then download the scripts  and copy the folders named
DataManagement, MasterR and Run_analysis into the folder ReefFishStability.

## Description of files
 - DataManagement: includes an R script with the same name that does data cleaning and generates the file called master.fish.dat for alpha stability analysis. This is required if you want to start the analysis from scratch. Otherwise, you can skip this part ond go directly to execute the code  provided in folder Run_analysis.
 - MasterR: includes the required scripts and aggregated data to run the full analysis, including plots. These scripts are called when running the code in Run_analysis.
 - Run_analysis: includes the code to run the various steps of the analysis. Please, note that some steps were originally performed on an HPC; the code provided here uses foreach loops, but in some cases, especially the jackknife simulation, the analysis is slow. This is noted when appropriate, but the slow parts of the code can be skipped since the resulting files are provided to run directly the downstream analysis.

