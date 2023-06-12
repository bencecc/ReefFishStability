#### ---- DATA MANAGEMENT - CREATE FOLDERS FOR R SCRIPTS AND RAW DATA AND GENERATE -------- ####
#### ---- THE NECESSARY DATASETS TO RUN THE FULL ANALYSIS OF REEF FISH STABILITY ---------- ####

# Running the following lines of code will create the required directories to run the analysis
# in your home directory:

# create workspace directories for R scripts

# copy DataManagement.R here
if(!file.exists('~/workspace/ReefFishStability/DataManagement'))
	dir.create('~/workspace/ReefFishStability/DataManagement')

# copy alpha_stability.R, pSEM.R, thermal_affinity.R and and gamma_stability.R here 
if(!file.exists('~/workspace/ReefFishStability/RunAnalysis'))
	dir.create('~/workspace/ReefFishStability/RunAnalysis')

# copy all other scripts here
if(!file.exists('~/workspace/ReefFishStability/MasterR'))
	dir.create('~/workspace/ReefFishStability/MasterR')

# create directory for data
if(!file.exists('~/Data_ReefFishStability'))
	dir.create('~/Data_ReefFishStability')

# This folder should include the following data files:
# 1) master.fish.dat.RData - the raw fish abundance data
# 2) fish.traits.RData - fish traits
# 3) covars.RData - environmental covariates
# 4) master.ecoregion.fish.dat.RData
# 5) all distance matrices: .... dist.RData

# create directory for figures
if(!file.exists('~/Data_ReefFishStability/Figs'))
	dir.create('~/Data_ReefFishStability/Figs')

#### ---- GENERATE DATA FOR ALPHA STABILITY ANALYSIS ---- ####

# generate folder to save partial results. NOTE: if you change folder name here,
# you should also update folder name at line 199 of script "alpha_stab_foreach.r"
if(!file.exists('~/Data_ReefFishStability/alpha_stab_tmp'))
	dir.create('~/Data_ReefFishStability/alpha_stab_tmp')

# require libraries
require(dplyr)
require(tidyr)
require(tibble)
require(knitr)
require(codyn)
require(FD)
require(mFD)

# set number of cores
nc <- 15
require(foreach, quietly=T)
require(doMC, quietly=T)
registerDoMC(cores=nc)

# set working directory
setwd("~/workspace/ReefFishStability/")

rmarkdown::render("ReefFishStability.Rmd")

# load ,aster.fish.dat
#load("master.fish.dat.RData")

# load raw fish abundance data data
load("raw.fish.dat.RData")

master.fish.dat <- raw.fish.dat %>%
		group_by(ID,SITE_ID,LAT,LON,YEAR,MPA,MPA_STATUS,MPA_NAME,
				IUCN_CAT,NO_TAKE,STATUS_YR,STATUS,SAMP_AGE,N_YEARS,INI_YEAR,
				END_YEAR,N_REPS,N_DATES,SAMPLED_AREA,TRANSECT_SIZE,ECOREGION,
				ECO_CODE,PROVINCE,PROV_CODE,REALM,RLM_CODE,SPECIES) %>%
		summarise(abund=sum(abund), .groups="drop") 

#load("master.fish.dat.RData")
load("fish.traits.RData")
# load environmental covariates
load("covars.RData")

# load required functions
source("~/workspace/ReefFishStability/MasterR/gross_w_func.R")
source("~/workspace/ReefFishStability/MasterR/sync_det_funcs.R")
source("~/workspace/ReefFishStability/MasterR/func_div_run.R")
source("~/workspace/ReefFishStability/MasterR/alpha_stab_foreach.R")

#### ---- RUN ALPHA STABILITY ANALYSIS ---- ####
alpha_stab_foreach(df=master.fish.dat)

# collect results
setwd("~/Data_ReefFishStability/alpha_stab_tmp")
file_list <- as.list(list.files())

alpha_res <- foreach(i=1:length(file_list), .combine='rbind') %dopar% {
	mget(load(file_list[[i]]))[[1]]
}

alpha.stab.res <- alpha_res %>% distinct(LAT,LON, .keep_all=T)

#### ---- Combine with environmental data ---- ####
site_fish_stab <- alpha.stab.res %>%
		left_join(covars, by=c("LAT","LON")) %>%
		mutate(
				STAB=log(1/CV_TOT_ABUND),
				SP.STAB=log(1/CV_SP_ABUND),
				ASYNC=ASYNC_GROSS_W,
				SR=log(N_SPEC),
				FRIC=(mFRic),
				MEAN.ABUND=log(MEAN_TOT_ABUND),
				AREA=log(SAMPLED_AREA),
				MPA=recode(MPA, Yes="Protected", No="Unprotected"),
				MPA=factor(MPA, levels=c("Protected","Unprotected"))
		) %>%
		
		relocate(c(MPA,MPA_STATUS,MPA_NAME), .after=LON) %>%
		drop_na(FRIC) |> 
		distinct(LAT,LON, .keep_all=T)

site_fish_stab$ID <- as.factor(site_fish_stab$ID)
site_fish_stab$SITE_ID <- as.factor(site_fish_stab$SITE_ID)

# save alpha stability data
setwd("~/Data_ReefFishStability/")
save(site_fish_stab, file="site_fish_stab.RData")

# How many population timeseries, species, sites, MPAs, ecoregions?
n.pop.ts <- (master.fish.dat %>% group_by(ID,SITE_ID) %>% distinct(SPECIES) %>%
			ungroup() %>% summarise(n.ts=n()))$n.ts
n.sp <- (master.fish.dat %>% distinct(SPECIES) %>% summarise(n.sp=n()))$n.sp
n.mpa.sites <- (site_fish_stab %>% filter(MPA=="Protected") %>% distinct(LAT,LON) %>% summarise(n.mpa=n()))$n.mpa
n.open.sites <- (site_fish_stab %>% filter(MPA=="Unprotected") %>% distinct(LAT,LON) %>% summarise(n.open=n()))$n.open
n.mpa <- (site_fish_stab %>% distinct(MPA_NAME) %>% summarise(n.mpa=n()))$n.mpa
n.strict.mpa <- (site_fish_stab %>% filter(MPA=="Yes") %>% distinct(MPA_NAME) %>% summarise(n.str=n()))$n.str
n.ecoreg <- (site_fish_stab %>% distinct(ECOREGION) %>% summarise(n.ecor=n()))$n.ecor

n.pop.ts
n.sp
n.mpa.sites
n.open.sites
n.mpa
n.ecoreg

#### ---- Summary by Study ID (Table 6) ---- ####

tab6 <- master.fish.dat %>%
		#filter by sites used in alpha stability analysis
		select(-c(INI_YEAR,END_YEAR)) %>%
		mutate(ID=as.factor(ID), SITE_ID=as.factor(SITE_ID)) %>%
		right_join(site_fish_stab[,c("ID","SITE_ID","LAT","LON","INI_YEAR","END_YEAR")],
				by=c("ID","SITE_ID","LAT","LON")) %>%
		# rename IDs
		mutate(DATA_SOURCE=case_when(
						ID=="rls" ~ "Reef Life Survey",
						ID=="reefcheck" ~ "ReefCheck",
						ID=="biotime.271" ~ "BioTIME",
						ID=="biotime.359" ~ "BioTIME",
						ID=="biotime.365" ~ "BioTIME",
						ID=="biotime.436" ~ "BioTIME",
						ID=="biotime.438" ~ "BioTIME",
						ID=="heenan" ~ "Western Central Pacific",
						ID=="sbc_lter" ~ "BioTIME",
						ID=="channel_islands" ~ "Santa Barbara Channel",
						ID=="aims_gbr_rm" ~ "Great Barrier Reef Marine Park",
						ID=="aims_gbr_rap" ~ "Great Barrier Reef Marine Park",
						ID=="portugal" ~ "Portugal",
						ID=="au_atrc" ~ "Southern Australia",
						ID=="capcreus" ~ "Western Mediterranean",
						ID=="cpalos" ~ "Western Mediterranean",
						ID=="medes" ~ "Western Mediterranean",
						TRUE ~ as.character(ID)
						)) %>%
				# extract relevant information
				group_by(DATA_SOURCE) %>%
				mutate(N_MPAS=length(unique(MPA_NAME)),
						max.years=max(N_YEARS),
						ini.year=min(INI_YEAR),
						end.year=max(END_YEAR)) %>%
				group_by(DATA_SOURCE, MPA) %>%
				distinct(SITE_ID, .keep_all=T) %>%
				summarise(nsites=sum(n()), nmpas=unique(N_MPAS),
						sampling.interval=paste(unique(ini.year), unique(end.year), sep="-"),
						longest.ts=unique(max.years),
						.groups="drop") %>%
				pivot_wider(names_from=MPA, values_from=c(nsites,nmpas,
								sampling.interval, longest.ts), values_fill=NA) %>%
				mutate(nsites_Yes=ifelse(DATA_SOURCE=="BioTIME", 0, nsites_Yes),
						nmpas_No=ifelse(DATA_SOURCE=="BioTIME", 0, nmpas_No),
						nmpas_Yes=ifelse(DATA_SOURCE=="BioTIME", 0, nmpas_Yes)) %>%
				# reshape table
				select(-c(nmpas_No, sampling.interval_Yes, longest.ts_Yes)) %>%
				rename(N.MPAs=nmpas_Yes, N.MPA.SITES=nsites_Yes, N.OA.SITES=nsites_No, 
						SAMPLING.INTERVAL=sampling.interval_No,
						LONGEST.TIMESERIES=longest.ts_No) %>%
				relocate(c(N.MPAs, N.MPA.SITES), .before=N.OA.SITES)

kable(tab6)

# NOTE: Although the BioTIME database included only open areas (OA), these data were part
# of the sampling programs maintained in 6 marine protected areas (MPAs) and were therefore
# included in alpha and gamma stability analyses.


