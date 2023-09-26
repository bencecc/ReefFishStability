# ReefFishStability
## R code for the analysis of the paper by Benedetti-Cecchi et al (in review): Marine protected areas promote stability of reef fish communities under climate warming.
### All analyses were done in R version 4.1.3 using libraries:
brms, broom, codyn, datawizard, DHARMa, doMC, FD, fishualize, foreach, ggeffects, ggstatsplot, igraph, knitr, lemon, lme4, lmerTest, mFD, mgcv, modelr, performance, piecewiseSEM, rnaturalearth, semEff, sf, sjPlot, SteinerNet, tidybayes, tidyverse, tidymv.

## USAGE: 
To run the analyses, first download the files on your computer, then:
1) create a folder in your home directory named Data_ReefFishStability and copy there the downloaded files from the same named folder;
2) create a folder in your workspace named ReefFishStability; the workspace folder should be in the home directory ~/workspace/ReefFishStability/;
3) copy all other folders (DataManagement, MasterR and Run_analysis) in ~/workspace/ReefFishStability/.
4) run code in the following order to matche the results as presented in the paper: alpha_stability.R, pSEM.R, sampling_coverage.R, thermal_affinity.R, gamma_stability.R, numm_models.R

## Description of files
 - Data_ReefFishStability: includes the raw data and derived data files for downstream analysis (this allows you to skip some slow computations that were originally performed on an HPC). 
 - DataManagement: includes R script DataManagement that aggregates the raw data by site and generates the file called master.fish.dat, which provides the input to perform  the alpha stability analysis, which, in turn, generates file site_fish_stab. These outputs are already provided in Data_ReefFishStability folder, so you can skip this part of the analysis and go directly to execute the code provided in folder Run_analysis. The DataManagement script includes the code to generate Fig. 2a (the map) and Supplementary Table 8, which reports information on data sources.
 - MasterR: includes the required scripts and aggregated data to run the full statistical and visualization analysis. These scripts are called when running the code in Run_analysis.
 - Run_analysis: includes the code to run the various steps of the analysis. Please, note that some steps were originally performed on an HPC; the code provided here uses foreach loops, but in some cases, especially the jackknife simulation, the analysis is slow. This is noted when appropriate. NOTE: the slow chunks of the code can be skipped since the files required to run the downstream analyses are provided in folder Data_ReefFishStability.

