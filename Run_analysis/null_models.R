#### ------------------------------ NULL MODELS OF GAMMA STABILITY ------------------------------------------ ####
#### ---- NOTE: this analysis is slow and was originally performed on an hpc cluster. The resulting file ---- ####
#### ---- cyclicshift.meta.dist.nullres.RData is provided in folder ~/Data_ReefFishStability and you may ---- ####
#### ---- want to go directly to the analysis at line 163 of this script. ----------------------------------- ####

# generate folder to save partial results. 
if(!file.exists("~/Data_ReefFishStability/mpa_null_dist_res"))
	dir.create("~/Data_ReefFishStability/mpa_null_dist_res")

# set working directory
setwd("~/Data_ReefFishStability")
#### ---- Load master.ecoregion.fish.dat, a subset of master.fish.dat that includes -------- ####
#### ---- ecoregions with suitable data for gamma stability analysis. Ecoregions        ---- ####
#### ---- where selected if they included at least 2 MPA and 2 OA sites sampled         ---- ####
#### ---- simultaneously for a minimum of 5 years. Selection was performed using        ---- ####
#### ---- an algorithm that maximised the number of sites for all possible combinations ---- ####
#### ---- of overlapping years ------------------------------------------------------------- ####
load("master.ecoregion.fish.dat.RData") 

# require libraries
require(dplyr)
require(codyn)
require(tidyr)
require(modelr)
require(ggplot2)
require(forcats)
require(ggstatsplot)
require(fishualize)

# set number of cores
nc <- 15
require(foreach, quietly=T)
require(doMC, quietly=T)
registerDoMC(cores=nc)

source("~/workspace/ReefFishStability/MasterR/stab_by_dist.R")
source("~/workspace/ReefFishStability/MasterR/meta_stab.R")
source("~/workspace/ReefFishStability/MasterR/gross_w_func.R")
source("~/workspace/ReefFishStability/MasterR/sync_det_funcs.R")
source("~/workspace/ReefFishStability/MasterR/stats_nullstab_dist.R")


# sync_det_funcs.R is needed for detrending with the
# three-term local variance, as in Leps et al. Ecography
source("~/workspace/ReefFishStability/MasterR/sync_det_funcs.R")

## function to perform cyclic-shift permutation as adapted from Hallett et al. (2014)
time_cyclic <- function(dat){
	dat <- data.frame(dat)
	rand.mat <- matrix(NA, nrow=nrow(dat), ncol=ncol(dat))
	rand.mat[,1] <- dat[,1]
	for (i in 2:ncol(dat)) rand.mat[,i] <- permute::shuffleSeries(as.vector(dat[,i]))
	rand.mat <- as.data.frame(rand.mat)
	colnames(rand.mat) <- colnames(dat)
	return(rand.mat)
}

# cyclic shift
cyclic_shift <- function(data, ...) {
	
	time.cyclic <- data %>%
			pivot_wider(names_from=c(SPECIES), values_from=abund, values_fill=0) %>%
			group_by(SITE_ID, MPA, MPA_NAME) %>%
			group_modify(~ as.data.frame(time_cyclic(.))) %>%
			pivot_longer(cols=which(colnames(.)==colnames(.)[5]):which(colnames(.)==colnames(.)[ncol(.)]),
					names_to=c("SPECIES"), values_to="abund")  %>%
			filter(abund>0) %>%
			arrange(SITE_ID, MPA, MPA_NAME, YEAR, SPECIES)
	time.cyclic$YEAR <- as.numeric(time.cyclic$YEAR)
	
	time.cyclic
}


id <- unique(master.ecoregion.fish.dat$ID)

for (i in id) {
	
	# select Ecoregion
	ecoreg.dat <- master.ecoregion.fish.dat %>% filter(ID%in%i) %>%
			group_by(SITE_ID, MPA, MPA_NAME, YEAR, SPECIES) %>%# remove possible duplicates
			summarise(abund=mean(abund), .groups="drop") %>% filter(abund>0)
	# load corresponding distance matrix (generated using script mpa_LeastCostDist.R)
	dist.dat.path <- paste("max", i, "dist.RData", sep=".")
	dist.dat <- mget(load(file.path('~/Data_ReefFishStability', dist.dat.path)))[[1]]
	
	# Define range of distances from 0 to the maximum distance between any two MPAs in ecoreg.dat
	mpa.sites <- (ecoreg.dat %>% 
				mutate(SITE_ID=as.character(SITE_ID)) %>%
				filter(MPA=="Yes") %>% distinct(SITE_ID) %>% dplyr::select(SITE_ID))$SITE_ID
	mpa.dist.df <- dist.dat[which(rownames(dist.dat)%in%mpa.sites),
					which(colnames(dist.dat)%in%mpa.sites)] %>%
			replace(., col(.) == row(.), NA) 
	target.range <- c(0, max(mpa.dist.df, na.rm=T))
	
	# run analysis
	rand.res <- foreach(k = 1:999, .combine="rbind",
					.packages=c("dplyr","codyn","tidyr","modelr","foreach"),
					.export=c("gross_w_func","cor_algo","phi_t3","var_t3","cov_t3","cor_t3",
							"time_cyclic","cyclic_shift","time_shuff_sync", "time_shuff_async")) %dopar% { 
				
				#cat('Doing iter ', k, ' of ', end, '\n', sep = '')
				# repeat randomization when any of the calculations produces only NAs
				MoveOn <- FALSE
				
				while(MoveOn==FALSE) {
					
					rand_dat <- do.call(cyclic_shift, list(ecoreg.dat)) %>% ungroup()
					
					# no need to reshape the distance matrix. The randomizing function maintains
					# the original matching between SITE_ID and MPA variables
					stab.rand.mpa <- (stab_by_dist(df=rand_dat, dist.dat=dist.dat, cut.range=target.range, mpa="Yes"))	
					stab.rand.oa <- try(stab_by_dist(df=rand_dat, dist.dat=dist.dat, cut.range=target.range, mpa="No"), silent=T)	
					
					if(!inherits(stab.rand.mpa, "try-error")&!inherits(stab.rand.oa, "try-error")) { 
						
						stab.rand=rbind(stab.rand.mpa, stab.rand.oa)
						stab.res <- data.frame(DATA.TYPE="Rand", stab.rand)
						MoveOn <- TRUE
					}
					
				}
				
				
				if(i==1) {
					
					stab.obs.mpa <- data.frame(DATA.TYPE="Obs", stab_by_dist(df=ecoreg.dat, dist.dat=dist.dat, cut.range=target.range, mpa="Yes"))				
					stab.obs.oa <- data.frame(DATA.TYPE="Obs", stab_by_dist(df=ecoreg.dat, dist.dat=dist.dat, cut.range=target.range, mpa="No"))				
					stab.obs=rbind(stab.obs.mpa, stab.obs.oa)
					stab.res <- rbind(stab.obs, stab.res)
					
				}	
				
				return(stab.res)
				
			}
	
	stats.res.tmp <- rand.res %>% group_by(MPA, DIST) %>%
			group_modify(~ stats_nullstab_dist(as.data.frame(.)))
	
	stats.res <- cbind(ID=i, stats.res.tmp)
	
	assign(paste(i, "_dist", floor(target.range[2]),  sep=""), value=stats.res, pos=1, inherits=T)
	outputName=paste(i, "_dist", floor(target.range[2]), ".RData",sep="")
	outputPath=file.path("~/Data_ReefFishStability/mpa_null_dist_res", outputName)
	save(list=paste(i, "_dist", floor(target.range[2]), sep=""), file=outputPath)
	
}


# collect and save results

meta_null_res <- foreach(i=1:length(file_list), .combine='rbind') %dopar% {
	mget(load(file_list[[i]]))[[1]]
}

setwd("~/Data_ReefFishStability/mpa_null_dist_res")

cyclicshift.meta.dist.nullres <- meta_null_res
save(cyclicshift.meta.dist.nullres,
		file="~/Data_ReefFishStability/cyclicshift.meta.dist.nullres.RData")

#### ---- PLOT NULL MODEL RESULTS ---- #### 

setwd("~/Data_ReefFishStability")
load("cyclicshift.meta.dist.nullres.RData")

# round distances and set maximum distance between MPA sites to 2 km for Portugal
cyclicshift.meta.dist.nullres <- cyclicshift.meta.dist.nullres %>%
		mutate(DIST=round(DIST,0)) %>%
		mutate(DIST=if_else(ID=="atleu", 2, DIST))
		
sel.metrics <- c(
		"loc_sp_comm_async_gross_w",
		"loc_comm_metacom_async_gross_w",
		"meta_pop_metacom_async_gross_w",
		"meta_sp_metacom_async_gross_w")

meta.dist.nullres <- cyclicshift.meta.dist.nullres  %>%
		filter(METRIC%in%sel.metrics) %>%
		mutate(
				ID=recode(ID,
						atleu="SEAS",bassian="BA",capehowe="CH",centsouthgbr="CSGBR",java="SSJS",
						lordhowe="LHNI",newcal="NEWC",northcalif="NC",northgbr="NGBR",nwmed="WM",
						southaust="SAG",southcalif="SCB",tweed="TM"),
				ID=factor(ID, 
						levels=rev(c("BA","CH","CSGBR","LHNI","NC","NEWC","NGBR",
										"SAG","SCB","SEAS","SSJS","TM","WM"))),
				METRIC=recode(METRIC,
						loc_sp_comm_async_gross_w="ASA",						
						loc_comm_metacom_async_gross_w="SCA",
						meta_pop_metacom_async_gross_w="SSA",
						meta_sp_metacom_async_gross_w="MPA"
				),
				METRIC=factor(METRIC,
						levels=c("SCA","SSA","ASA","MPA"))
		) %>%
		mutate(
				ID=factor(ID),
				ID=fct_reorder(ID, DIST, mean)
		)	

windows(width=8, height=10)
meta.dist.nullres %>%
		ggplot(aes(x = fct_reorder(ID, DIST, mean), y=OBSERVED,
						ymin=CI_LW95, ymax=CI_UP95)) +
		geom_pointrange(aes(col=MPA, group=MPA),
				shape = 18, size = 0.5,
				position=position_dodge(width = 0.8)) + 
		geom_hline(yintercept = 0, color = "grey60", linetype = 1, size = 1) +
		scale_color_fish_d(option = "Cirrhilabrus_solorensis", begin=0.95, end=0.05, direction=1) +
		xlab("") + 
		ylab("Asynchrony") + 
		theme_bw() +
		coord_flip() +
		facet_grid(DIST ~ METRIC, scales="free") +
		theme(
				legend.position="none",
				strip.background=element_rect(colour="grey90",
						fill="grey90"),
				strip.text=element_text(size=12),
				panel.spacing.x=unit(1.5,"line"),
				panel.border=element_rect(colour="grey60", size=0.2),
				panel.background = element_blank(),
				panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(), 
				#axis.line = element_line(colour = "black"),
				axis.text.y = element_text(size = 10, colour = "black"),
				axis.text.x.bottom = element_text(size = 8, colour = "black"),
				axis.title.x = element_text(size = 12, colour = "black"))

dev.off()
ggsave(file = "~/Data_ReefFishStability/Figs/cyclicshift.nullmeta.dist.pdf",
		dpi = 300, width = 8, height = 10, useDingbats=FALSE)


















