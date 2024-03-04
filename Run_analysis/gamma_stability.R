#### ------------------------------ GAMMA STABILITY ANALYSIS --------------------------------------- ####
#### ---- NOTE: this analysis is slow and was originally performed on an hpc cluster. The       ---- ####
#### ---- resulting file metastats.res.RData is provided in folder ~/Data_ReefFishStability     ---- ####
#### ---- and you may want to go directly to the analysis at line 123 of this script.           ---- ####

# generate folder to save partial results. 
if(!file.exists("~/Data_ReefFishStability/jnife_res"))
	dir.create("~/Data_ReefFishStability/jnife_res")

# set working directory
setwd("~/Data_ReefFishStability")
#### ---- Load master.ecoregion.fish.dat, a subset of master.fish.dat that includes -------- ####
#### ---- ecoregions with suitable data for gamma stability analysis. Ecoregions        ---- ####
#### ---- where selected if they included at least 2 MPA and 2 OA sites sampled         ---- ####
#### ---- simultaneously for a minimum of 5 years. Selection was performed using        ---- ####
#### ---- an algorithm that maximised the number of sites for all possible combinations ---- ####
#### ---- of overlapping years. ------------------------------------------------------------ ####
load("master.ecoregion.fish.dat.RData") 

# require libraries
require(tidyverse)
require(codyn)
require(broom)
require(brms)
require(tidybayes)
require(lemon)
require(fishualize)
require(ggpubr)
require(knitr)

# set number of cores
nc <- 15
require(foreach, quietly=T)
require(doMC, quietly=T)
registerDoMC(cores=nc)

source("~/workspace/ReefFishStability/MasterR/sites_jackknife_stab.R")
source("~/workspace/ReefFishStability/MasterR/stab_by_dist_sites.R")
source("~/workspace/ReefFishStability/MasterR/meta_stab.R")
source("~/workspace/ReefFishStability/MasterR/gross_w_func.R")
# sync_det_funcs.R is needed for detrending with the
# three-term local variance, as in Leps et al. Ecography
source("~/workspace/ReefFishStability/MasterR/sync_det_funcs.R")
source("~/workspace/ReefFishStability/MasterR/metacomm_sites_dist.R")

id <- unique(master.ecoregion.fish.dat$ID)

for (i in id) {
	
	# select Ecoregion
	ecoreg.dat <- master.ecoregion.fish.dat %>% filter(ID%in%i)
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
	mpa.jnife.est <- sites_jackknife_stab(df=ecoreg.dat, mpa="Yes",
			cut.range=target.range, dist.dat=dist.dat)
	oa.jnife.est <- sites_jackknife_stab(df=ecoreg.dat, mpa="No",
			cut.range=target.range, dist.dat=dist.dat)
	
	if(!is.null(mpa.jnife.est)&!is.null(oa.jnife.est)) {
		
		mpa.metastats.obs.res <- stab_by_dist_sites(df=ecoreg.dat, dist.dat=dist.dat, cut.range=target.range, mpa="Yes") %>%
				bind_rows()
		oa.metastats.obs.res <- stab_by_dist_sites(df=ecoreg.dat, dist.dat=dist.dat, cut.range=target.range, mpa="No") %>%
				bind_rows()
		metastats.obs.res <- rbind(mpa.metastats.obs.res, oa.metastats.obs.res)
		
		metastats.res.tmp <- metastats.obs.res
		metastats.res.tmp$ID <- rep(unique(ecoreg.dat$ID),2)
		metastats.res.tmp$PARM.TYPE <- "Obs"
		metastats.res.tmp <- metastats.res.tmp %>%
				relocate(ID, .before=MPA) %>% relocate(PARM.TYPE, .after=DIST)
		metastats.res.tmp1 <- metastats.res.tmp %>%
				gather('METRICS', 'response',
						which(colnames(.)==colnames(.)[10]):which(colnames(.)==colnames(.)[ncol(.)])) %>%
				pivot_wider(names_from=PARM.TYPE, values_from=response, values_fill=0)
		
		jnife.res <- rbind(mpa.jnife.est, oa.jnife.est) %>% rename(PARM.TYPE="Jnife.est") 
		jnife.resl <- jnife.res %>% gather('METRICS', 'response',
						which(colnames(.)==colnames(.)[6]):which(colnames(.)==colnames(.)[ncol(.)])) %>%
				pivot_wider(names_from=PARM.TYPE, values_from=response, values_fill=0)
		
		ecor.maxdist <- data.frame(ID=unique(ecoreg.dat$ID), MPANET_EXTENT=max(mpa.dist.df, na.rm=T))
		metastats.res <- metastats.res.tmp1 %>%
				left_join(jnife.resl, by=c("ID","MPA","DIST","METRICS")) %>%
				relocate(NSPEC, .before=Obs) %>% drop_na(Jnife.mean) %>%
				left_join(ecor.maxdist, by="ID") %>%
				relocate(MPANET_EXTENT, .after=SAMPLED_AREA)
		
		
		assign(paste(i, "_dist", round(target.range[2], 0),  sep=""), value=metastats.res, pos=1, inherits=T)
		outputName=paste(i, "_dist", round(target.range[2], 0), ".RData",sep="")
		outputPath=file.path("~/Data_ReefFishStability/jnife_res", outputName)
		save(list=paste(i, "_dist", floor(target.range[2]), sep=""), file=outputPath)
		
	}
		
}

# collect results
setwd("~/Data_ReefFishStability/jnife_res")
file_list <- as.list(list.files())

mpa_jnife_res <- foreach(i=1:length(file_list), .combine='rbind') %dopar% {
	mget(load(file_list[[i]]))[[1]]
}

metastats.res <- mpa_jnife_res %>% arrange(ID, DIST, MPA)
save(metastats.res,
		file="~/Data_ReefFishStability/metastats.res.RData")

#### --- RUN BAYESIAN META-ANALYSIS ON metastats.R ---- ####
setwd("~/Data_ReefFishStability")

# load data
load("~/Data_ReefFishStability/metastats.res.RData")

#### ---- The original analysis on HPC calculated Gamma stability and asynchorny measures at three distance ranges:       ---- ####
#### ---- of 0-20, 21-50, 51-100 kms and at the maximum distance between MPAs in an ecoregion. Statistics are included    ---- ####
#### ---- in the current metastats.res dataset. Below is the code to perform the maximum distance analysis (Fig. 5b.      ---- ####
#### ---- and Supplementary Fig. 12b.											  ---- ####
#### ---- line 480: Summary of characteristics of metacommunities - Supplementary Table 5                                 ---- ####
#### ---- line 513: Connectivity assessment - average metacommunity dissimilarity (also in Supplementary Table 5)         ---- ####
#### ---- CHECKS AND SENSITIVITY ANALYSIS:                                                                                ---- ####
#### ---- line 542: Connectivity assessment - global network connectivity - Supplementary Figs. 10 and 11                 ---- ####
#### ---- line 726: relationships between closenees centrality measured on from physically-derived graphs and degree      ---- ####
#### ---- centrality from biologically-derived graphs - Results for Supplementary Table 6                                 ---- ####
#### ---- line 757: relationhsips of stability and asynchony measures with MPA network characterstics (extent,            ---- #### 
#### ---- number of MPAs and number of sites) - Results for Supplementary Table 7                                         ---- ####
#### ---- line 801: Repeat metacommunity stability and asynchorny analysis at the 51-100 km distance range (Supplementary ---- ####
#### ---- Fig. 12)

#### ---- Fig. 12)

#### ---- MAXIMUM MPA DISTANCE ANALYSIS ---- ####
# filter DIST to select the maximum distance between any two MPAs in a metacommunity, which is equal to MPANET_EXTENT
metastats.res <- metastats.res %>% group_by(ID) %>%
		filter(DIST==MPANET_EXTENT) %>%
		mutate(DIST=round(DIST, 0))

# select metrics and prepare data for analysis
sel.metrics <- c("gamma_stab","alpha_stab","pop_stab","species_stab","loc_sp_comm_async_gross_w",
		"loc_comm_metacom_async_gross_w","meta_pop_metacom_async_gross_w",
		"meta_sp_metacom_async_gross_w")

dat.mpa <- metastats.res %>% filter(MPA=="Yes"&
						METRICS%in%sel.metrics) %>%
		rename(Obs.mpa=Obs, Nspec.mpa=NSPEC, Jmean.mpa=Jnife.mean,
				Jse.mpa=Jnife.se, Jvar.mpa=Jnife.var,
				SAMPLED_AREA_MPA=SAMPLED_AREA,
				MEAN_N_SITES_MPA=MEAN_N_SITES,
				MAX_N_SITES_MPA=MAX_N_SITES) %>%
		dplyr::select(-c("MPA"))

dat.oa <- metastats.res %>% filter(MPA=="No"&
						METRICS%in%sel.metrics) %>% 
		#select(-MPANET_EXTENT) %>%
		rename(Obs.oa=Obs, Nspec.oa=NSPEC, Jmean.oa=Jnife.mean,
				Jse.oa=Jnife.se, Jvar.oa=Jnife.var,
				SAMPLED_AREA_NOMPA=SAMPLED_AREA,
				MEAN_N_SITES_NOMPA=MEAN_N_SITES,
				MAX_N_SITES_NOMPA=MAX_N_SITES) %>%
		dplyr::select(-c("N_MPA"))

# calculate Hedge's g effect size; since the jackknife method uses the leave-one-out 
# approach on indivdual species, the total number of species in each ecoregion provides
# the sample size for calculating effect sizes
eff.size.df <- dat.mpa  %>% left_join(dat.oa,
				by=c("ID","DIST","MPANET_EXTENT","DIST_RANGE","METRICS")) %>%
			mutate(
				sampled.area=SAMPLED_AREA_MPA/SAMPLED_AREA_NOMPA,
				pooled.sd=sqrt(((Nspec.mpa-1)*Jvar.mpa+(Nspec.oa-1)*Jvar.oa)/(Nspec.mpa+Nspec.oa-2)),
				dObs=(Obs.mpa-Obs.oa)/pooled.sd,
				dJnife=(Jmean.mpa-Jmean.oa)/pooled.sd,
				se.dObs=sqrt((Nspec.mpa+Nspec.oa)/(Nspec.mpa*Nspec.oa)+(dObs^2/(2*(Nspec.mpa+Nspec.oa)))),
				se.dJnife=sqrt((Nspec.mpa+Nspec.oa)/(Nspec.mpa*Nspec.oa)+(dJnife^2/(2*(Nspec.mpa+Nspec.oa)))),
				gObs=(1-3/(4*(Nspec.mpa+Nspec.oa-2)-1))*dObs,
				gJnife=(1-3/(4*(Nspec.mpa+Nspec.oa-2)-1))*dJnife,
				se.gObs=sqrt((1-3/(4*(Nspec.mpa+Nspec.oa-2)-1))^2*se.dObs^2),
				se.gJnife=sqrt((1-3/(4*(Nspec.mpa+Nspec.oa-2)-1))^2*se.dJnife^2)
		) %>% 
		filter(!dObs%in%c(-Inf,Inf)) %>%
		drop_na() %>%
		mutate(ID=recode(ID,
						atleu="SEAS",bassian="BA",capehowe="CH",centsouthgbr="CSGBR",java="SSJS",
						lordhowe="LHNI",northcalif="NC",northgbr="NGBR",nwmed="WM",
						southaust="SAG",southcalif="SCB",tweed="TM")
		)

# meta-analysis - async stats
locsp <- brm(
		gObs | se(se.gObs) ~ 1 + (1|ID),
		data = eff.size.df %>%
				filter(METRICS=="loc_sp_comm_async_gross_w"), 
		prior = c(prior(normal(0, 1), class = Intercept),
				prior(cauchy(0, 1), class = sd)),
		iter = 4000, warmup = 1000, chains = 4,
		seed=237835)

spat <- brm(
		gObs | se(se.gObs) ~ 1 + (1|ID),
		data = eff.size.df %>%
				filter(METRICS=="loc_comm_metacom_async_gross_w"), 
		prior = c(prior(normal(0, 1), class = Intercept),
				prior(cauchy(0, 1), class = sd)),
		iter = 4000, warmup = 1000, chains = 4,
		seed=237835)

metapop <- brm(
		gObs | se(se.gObs) ~ 1 + (1|ID),
		data = eff.size.df %>%
				filter(METRICS=="meta_pop_metacom_async_gross_w"), 
		prior = c(prior(normal(0, 1), class = Intercept),
				prior(cauchy(0, 1), class = sd)),
		iter = 4000, warmup = 1000, chains = 4,
		seed=237835)

metasp <- brm(
		gObs | se(se.gObs) ~ 1 + (1|ID),
		data = eff.size.df %>%
				filter(METRICS=="meta_sp_metacom_async_gross_w"), 
		prior = c(prior(normal(0, 1), class = Intercept),
				prior(cauchy(0, 1), class = sd)),
		iter = 4000, warmup = 1000, chains = 4,
		seed=237835)

#### ---- Posterior distributions ---- ####

# this is average species asynchorny metric as in Lamy et al
post.locsp <- eff.size.df %>%
		filter(METRICS=="loc_sp_comm_async_gross_w") %>%
		# Get the posterior predictions
		add_epred_draws(object = locsp) %>%
		mutate(GM=fixef(locsp)[1,1], Q2.5=fixef(locsp)[1,3], Q97.5=fixef(locsp)[1,4])
post.locsp$METRICS <- "aver.sp"

post.spat <- eff.size.df %>%
		filter(METRICS=="loc_comm_metacom_async_gross_w") %>%
		# Get the posterior predictions
		add_epred_draws(object = spat) %>%
		mutate(GM=fixef(spat)[1,1], Q2.5=fixef(spat)[1,3], Q97.5=fixef(spat)[1,4])
post.spat$METRICS <- "spatial"

# this is spatial species asynchrony as in Lamy et al
post.metapop <- eff.size.df %>%
		filter(METRICS=="meta_pop_metacom_async_gross_w") %>%
		# Get the posterior predictions
		add_epred_draws(object = metapop) %>%
		mutate(GM=fixef(metapop)[1,1], Q2.5=fixef(metapop)[1,3], Q97.5=fixef(metapop)[1,4])
post.metapop$METRICS <- "spat.sp"

# this is metapopulation asynchorny as in Lamy et al
post.metasp <- eff.size.df %>%
		filter(METRICS=="meta_sp_metacom_async_gross_w") %>%
		# Get the posterior predictions
		add_epred_draws(object = metasp) %>%
		mutate(GM=fixef(metasp)[1,1], Q2.5=fixef(metasp)[1,3], Q97.5=fixef(metasp)[1,4])
post.metasp$METRICS <- "metapop"

posterior.async.ghedge <- rbind(post.locsp, post.spat, post.metapop, post.metasp)
save(posterior.async.ghedge, file="posterior.async.ghedge.RData")

#### ---- STABILITY MEASURES ---- ####

mg <- brm(
		gObs | se(se.gObs) ~ 1 +  (1|ID),
		data = eff.size.df %>%
				filter(METRICS=="gamma_stab"), 
		prior = c(prior(normal(0, 1), class = Intercept),
				prior(cauchy(0, 1), class = sd)),
		iter = 4000, warmup = 1000, chains = 4,
		seed=237835)

ma <- brm(
		gObs | se(se.gObs) ~ 1 + (1|ID),
		data = eff.size.df %>%
				filter(METRICS=="alpha_stab"), 
		prior = c(prior(normal(0, 1), class = Intercept),
				prior(cauchy(0, 1), class = sd)),
		iter = 4000, warmup = 1000, chains = 4,
		seed=237835)

ms <- brm(
		gObs | se(se.gObs) ~ 1 + (1|ID),
		data = eff.size.df %>%
				filter(METRICS=="species_stab"), 
		prior = c(prior(normal(0, 1), class = Intercept),
				prior(cauchy(0, 1), class = sd)),
		iter = 4000, warmup = 1000, chains = 4,
		seed=237835)

mp <- brm(
		gObs | se(se.gObs) ~ 1 + (1|ID),
		data = eff.size.df %>%
				filter(METRICS=="pop_stab"), 
		prior = c(prior(normal(0, 1), class = Intercept),
				prior(cauchy(0, 1), class = sd)),
		iter = 4000, warmup = 1000, chains = 4,
		seed=237835)

post.gamma <- eff.size.df %>%
		filter(METRICS=="gamma_stab") %>%
		# Get the posterior predictions
		add_epred_draws(object = mg) %>%
		mutate(GM=fixef(mg)[1,1], Q2.5=fixef(mg)[1,3], Q97.5=fixef(mg)[1,4])
post.gamma$METRICS <- "gamma"

post.alpha <- eff.size.df %>%
		filter(METRICS=="alpha_stab") %>%
		# Get the posterior predictions
		add_epred_draws(object = ma) %>%
		mutate(GM=fixef(ma)[1,1], Q2.5=fixef(ma)[1,3], Q97.5=fixef(ma)[1,4])
post.alpha$METRICS <- "alpha"

post.species <- eff.size.df %>%
		filter(METRICS=="species_stab") %>%
		# Get the posterior predictions
		add_epred_draws(object = ms) %>%
		mutate(GM=fixef(ms)[1,1], Q2.5=fixef(ma)[1,3], Q97.5=fixef(ms)[1,4])
post.species$METRICS <- "species"

post.population <- eff.size.df %>%
		filter(METRICS=="pop_stab") %>%
		# Get the posterior predictions
		add_epred_draws(object = mp) %>%
		mutate(GM=fixef(mp)[1,1], Q2.5=fixef(mp)[1,3], Q97.5=fixef(mp)[1,4])
post.population$METRICS <- "population"

posterior.stab.ghedge <- rbind(post.alpha, post.gamma, post.species, post.population)
save(posterior.stab.ghedge, file="posterior.stab.ghedge.RData")

#### ---- PLOT POSTERIOR DISTRUBUTIONS ---- ####

load("~/Data_ReefFishStability/posterior.async.ghedge.RData")
load("~/Data_ReefFishStability/posterior.stab.ghedge.RData")

# Recode metrics to match definitions in Lamy et al. 2019:
# GAS = Gamma stability
# AAS: average alpha stability
# ASS: average species stability
# MPS: Metapopulation stability

# SCA: Spatial community asynchrony
# SSA: Spatial species asynchrony
# ASA: Average species asynchrony
# MPA: Metapopulation asynchrony

# prepare data for plotting
stab.df <- posterior.stab.ghedge %>%
		mutate(
				METRICS=recode(METRICS, gamma="GAS", alpha="AAS",
						species="ASS", population="MPS"),
				METRICS=factor(METRICS,
						levels=c("GAS","AAS","ASS","MPS")),
				METRIC_TYPE="Stability"
		) %>%
		relocate(METRIC_TYPE, .after="METRICS")

async.df <- posterior.async.ghedge %>%
		mutate(
				METRICS=recode(METRICS, metapop="MPA", spat.sp="SSA",
						aver.sp="ASA", spatial="SCA"),
				METRICS=factor(METRICS,
						levels=c("SCA","SSA","ASA","MPA")),
				METRIC_TYPE="Asynchrony"
		) %>%
		relocate(METRIC_TYPE, .after="METRICS")

metastab.df <- rbind(stab.df, async.df) %>%
		arrange(DIST) %>%
		mutate(
				MPA.EXTENT=factor(round(DIST, 0)),
				MPA.EXTENT=fct_reorder(MPA.EXTENT, DIST, mean)
		) %>%
		ungroup()

point.est.df <- metastab.df %>%
		group_by(ID, DIST, METRIC_TYPE, METRICS,
				N_MPA, MEAN_N_SITES_MPA, MAX_N_SITES_MPA, 
				MEAN_N_SITES_NOMPA, MAX_N_SITES_NOMPA) %>%
		summarise(mepred=mean(.epred), .groups="drop") %>%
		mutate(MPA.EXTENT=factor(DIST),
				MPA.EXTENT=fct_reorder(MPA.EXTENT, DIST, mean)
		)

x.labs <- metastab.df %>%
		mutate(
				ID=factor(ID),
				ID=fct_reorder(ID, DIST, mean),
				#MPANET_EXTENT=floor(MPANET_EXTENT),
				SITE_RATIO=round(MEAN_N_SITES_MPA/MEAN_N_SITES_NOMPA,3),
				SRF=as.factor(SITE_RATIO)
		) %>%
		select(c(ID,DIST,MPA.EXTENT,SITE_RATIO, SRF)) %>%
		distinct(ID, .keep_all=T)

p.meta.stab <- metastab.df %>% filter(METRIC_TYPE=="Stability") %>%
		ggplot(aes(x = fct_reorder(MPA.EXTENT, DIST, mean), y =.epred,
						fill=MPA.EXTENT, col=MPA.EXTENT)) +
		stat_eye(side="right", alpha=0.8, size=0.1, stroke=0.1,
				point_size=0.5, point_color="grey40") +
		stat_pointinterval(.width = c(.66, .95),
				point_size=0.5, point_color="grey40", 
				interval_color="grey40") +
		scale_x_discrete(labels=x.labs$DIST) +
		geom_hline(yintercept = 0, col="red", linewidth=0.5, linetype=2) +
		scale_fill_fish_d(breaks=x.labs$DIST,
				option = "Cirrhilabrus_tonozukai", name="MPA extent (km)") +
		scale_color_fish_d(labels=x.labs$DIST,
				option = "Cirrhilabrus_tonozukai", name="MPA extent (km)") +
		labs(x="", y="") +
		theme_bw() +
		facet_rep_wrap(~ METRICS, scales="free",
				repeat.tick.labels = "none", ncol=2) +
		theme(
				panel.border=element_rect(colour="grey60", linewidth=0.2),
				strip.background=element_rect(colour="grey90",
						fill="grey90"),
				legend.position="none",
				panel.spacing = unit(0.4, "lines"),
				axis.ticks = element_line(color="grey60", linewidth=0.4),
				axis.text.y  = element_text(size=10, hjust = 0),
				axis.text.x = element_text(size=10, angle = 45, vjust=0.9, hjust=0.8),
				panel.grid = element_blank()
		
		)

#windows(width=8, height=6)
p.meta.stab
dev.off()
ggsave(file = "~/Data_ReefFishStability/Figs/Fig5b.pdf",
		dpi = 300, width = 8, height = 6, useDingbats=FALSE)

p.meta.async <- metastab.df %>% filter(METRIC_TYPE=="Asynchrony") %>%
		ggplot(aes(x = fct_reorder(MPA.EXTENT, DIST, mean), y =.epred,
						fill=MPA.EXTENT, col=MPA.EXTENT)) +
		stat_eye(side="right", alpha=0.8, size=0.1, stroke=0.1,
				point_size=0.5, point_color="grey40") +
		stat_pointinterval(.width = c(.66, .95),
				point_size=0.5, point_color="grey40", 
				interval_color="grey40") +
		scale_x_discrete(labels=x.labs$DIST) +
		geom_hline(yintercept = 0, col="red", linewidth=0.5, linetype=2) +
		scale_fill_fish_d(breaks=x.labs$DIST,
				option = "Cirrhilabrus_tonozukai", name="MPA extent (km)") +
		scale_color_fish_d(labels=x.labs$DIST,
				option = "Cirrhilabrus_tonozukai", name="MPA extent (km)") +
		labs(x="", y="") +
		theme_bw() +
		facet_rep_wrap(~ METRICS, scales="free",
				repeat.tick.labels = "none", ncol=2) +
		theme(
				panel.border=element_rect(colour="grey60", linewidth=0.2),
				strip.background=element_rect(colour="grey90",
						fill="grey90"),
				legend.position="none",
				panel.spacing = unit(0.4, "lines"),
				axis.ticks = element_line(color="grey60", linewidth=0.4),
				axis.text.y  = element_text(size=10, hjust = 0),
				axis.text.x = element_text(size=10, angle = 45, vjust=0.9, hjust=0.8),
				panel.grid = element_blank()
		
		)

#windows(width=8, height=6)
p.meta.async
dev.off()
ggsave(file = "~/Data_ReefFishStability/Figs/FigS12b.pdf",
		dpi = 300, width = 8, height = 6, useDingbats=FALSE)

# Summary of characteristics of metacommunities - Supplementary Table 5
n.sites <- metastab.df %>% group_by(ID) %>%
		filter(DIST==max(DIST)) %>%
		slice(METRIC_TYPE=1L) %>%
		select(ID, N_MPA, MAX_N_SITES_MPA, MAX_N_SITES_NOMPA) %>%
		mutate(N_SITES_MPA=(MAX_N_SITES_MPA),
				N_SITES_OA=(MAX_N_SITES_NOMPA)) %>%
		ungroup() %>%
		select(ID, N_MPA, N_SITES_MPA, N_SITES_OA)


supp.tab5 <- metastats.res %>%
		mutate(DIST=round(DIST, 0),
				DIST_RANGE=round(DIST_RANGE, 0),
				ID=recode(ID,
						atleu="SEAS",bassian="BA",capehowe="CH",centsouthgbr="CSGBR",java="SSJS",
						lordhowe="LHNI",northcalif="NC",northgbr="NGBR",nwmed="WM",
						southaust="SAG",southcalif="SCB",tweed="TM")
		) %>%
		group_by(ID, MPA, DIST, N_MPA, MEAN_N_SITES, MAX_N_SITES, SAMPLED_AREA,
				DIST_RANGE) %>%
		slice(1L) %>%
		ungroup() %>%
		select(ID, DIST, DIST_RANGE) %>%
		rename(MPA_SCALE=DIST, MAXIMUM_DISTANCE=DIST_RANGE) %>%
		distinct(ID, .keep_all=T) %>%
		right_join(n.sites, by="ID") %>%
		arrange(MPA_SCALE)

kable(supp.tab5) 
# NOTE: N_MPA for MPA = No indicates the number of MPAs within the distance range
# to OA sites; the actual number of MPAs is given by MPA=Yes. 

# add mean dissimilarity
diss.df <- master.ecoregion.fish.dat %>%
		filter(!ECOREGION%in%c("New Caledonia","Palawan/North Borneo")) %>%
		group_by(ECOREGION, SITE_ID, MPA, SPECIES) %>%
		summarise(abund=sum(abund), .groups="drop") %>%
		ungroup() %>%
		pivot_wider(names_from=SPECIES, values_from=abund, values_fill=0)

ecoreg <- unique(diss.df$ECOREGION)

metajacc.tmp <- foreach(i=1:length(ecoreg)) %dopar% {
	
	mpa.df <- diss.df %>% filter(MPA=="Yes"&ECOREGION%in%ecoreg[i]) %>%
			select(-c(ECOREGION,SITE_ID,MPA))
	mpa.diss <- vegan::vegdist(mpa.df[,which(colSums(mpa.df)>0)], method="jaccard", binary=T) 
	
	oa.df <- diss.df %>% filter(MPA=="No"&ECOREGION%in%ecoreg[i]) %>%
			select(-c(ECOREGION,SITE_ID,MPA))
	oa.diss <- vegan::vegdist(oa.df[,which(colSums(oa.df)>0)], method="jaccard", binary=T) 
	
	out <- rbind(data.frame(ECOREGION=ecoreg[i], MPA="No", MeanJacc=round(mean(mpa.diss),2)),
			data.frame(ECOREGION=ecoreg[i], MPA="Yes", MeanJacc=round(mean(oa.diss),2)))
	
}

metajacc.res <- do.call(rbind, metajacc.tmp)

kable(metajacc.res)

#### ---- NETWORK APPROACH TO CONNECTIVITY IN METACOMMUNITIES ---- ####

# load required libraries
require(igraph)

# load data 
setwd("~/Data_ReefFishStability")
load("master.ecoregion.fish.dat.RData")
load("~/Data_ReefFishStability/metastats.res.RData") # for information on network size
# source utility function to plot igraph networks with ggplot
source("~/workspace/ReefFishStability/MasterR/ggplot_igraph.R")

# run analysis; first generate a datatset with species as columns that includes the
# relevant information for plotting network graphs
diss.df.tmp <- master.ecoregion.fish.dat %>%
		filter(!ECOREGION%in%c("New Caledonia","Palawan/North Borneo")) %>%
		group_by(ECOREGION, SITE_ID, MPA, SPECIES) %>%
		summarise(abund=sum(abund), .groups="drop") %>%
		pivot_wider(names_from=SPECIES, values_from=abund, values_fill=0) %>%
		mutate(ECOREGION=case_when(
						ECOREGION=="Bassian" ~ "BA",
						ECOREGION=="Cape Howe" ~ "CH",
						ECOREGION=="Central and Southern Great Barrier Reef" ~ "CSGBR",
						ECOREGION=="Lord Howe and Norfolk Islands" ~ "LHNI",
						ECOREGION=="Northern California" ~ "NC",
						ECOREGION=="South Australian Gulfs" ~ "SAG",
						ECOREGION=="South European Atlantic Shelf" ~ "SEAS",
						ECOREGION=="Southern California Bight" ~ "SCB",
						ECOREGION=="Sunda Shelf/Java Sea" ~ "SSJS",
						ECOREGION=="Torres Strait Northern Great Barrier Reef" ~ "NGBR",
						ECOREGION=="Tweed-Moreton" ~ "TM",
						ECOREGION=="Western Mediterranean" ~ "WM",
						TRUE ~ as.character(ECOREGION)
				))

net.size.df <- metastats.res %>%
		group_by(ID) %>%
		filter(DIST==MPANET_EXTENT) %>%
		mutate(DIST=round(DIST, 0)) %>% ungroup() %>%
		select(ID, DIST, DIST_RANGE) %>%
		distinct(ID, .keep_all=T) %>%
		mutate(DIST_RANGE=round(DIST_RANGE, 0),
				ECOREGION=recode(ID,
						atleu="SEAS",bassian="BA",capehowe="CH",centsouthgbr="CSGBR",java="SSJS",
						lordhowe="LHNI",northcalif="NC",northgbr="NGBR",nwmed="WM",
						southaust="SAG",southcalif="SCB",tweed="TM")
		) %>%
		select(-ID)

diss.df <- diss.df.tmp %>%
		left_join(net.size.df, by="ECOREGION") %>%
		relocate(c("DIST","DIST_RANGE"), .after=MPA)

ecoreg <- unique(diss.df$ECOREGION)

# Network based on Jaccard dissimilarity
jacc.conn.net.tmp <- foreach(i=1:length(ecoreg)) %dopar% {
	
	ecoreg.dat <- diss.df %>% filter(ECOREGION%in%ecoreg[i])
	jacc.df <- ecoreg.dat %>% select(-c(ECOREGION,SITE_ID,MPA,DIST,DIST_RANGE))
	jacc.dist <- vegan::vegdist(jacc.df[,which(colSums(jacc.df)>0)], method="jaccard", binary=T) 
	jacc.mat <- as.matrix(jacc.dist)
	jacc.net.tmp <- graph.adjacency(jacc.mat, mode='undirected', weighted=TRUE)
	jacc.inet <- mst(jacc.net.tmp, weights = NULL, algorithm = NULL)
	jacc.inet.plot <- ggplot_igraph(plot.graph=jacc.inet, MPA=ecoreg.dat$MPA,
			mpa.extent=ecoreg.dat$DIST[1], net.extent=ecoreg.dat$DIST_RANGE[1], net.name=ecoreg[i])
}

jacc.conn.net <- bind_rows(jacc.conn.net.tmp) %>%
		mutate(IDlabel=as.factor(paste(ID, " (", net.extent, ")", sep="")),
		IDlabel=fct_reorder(IDlabel, mpa.extent, mean)
)

figs10 <- ggplot(data=jacc.conn.net , aes(x=V1, y=V2))+
		geom_segment(data=na.exclude(jacc.conn.net),
				aes(x=from.x,xend = to.x, y=from.y,yend = to.y),
				linewidth=1, colour="darkgrey") +
		geom_point(aes(color=MPA), alpha=1, size=2.5) +
		scale_color_fish_d(option = "Cirrhilabrus_solorensis", begin=0.05, end=0.95, direction=-1) +
		labs(x="", y="") +
		theme_bw() +
		facet_wrap(. ~ IDlabel, scales="free") +
		theme(
				legend.position="none",
				strip.background=element_rect(colour="grey90",
						fill="grey90"),
				plot.subtitle = element_text(size = 10),					
				panel.spacing.x=unit(0.1,"line"),
				panel.border=element_rect(colour="grey60", linewidth=0.2),
				panel.background = element_blank(),
				panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(),
				axis.text.x = element_blank(),
				axis.text.y = element_blank(),
				axis.ticks.x = element_blank(),
				axis.ticks.y = element_blank(),
				plot.margin = margin(0, 0.1, 0.1, 0.1, "cm")
		)

#windows(width=8, height=8)
figs10
dev.off()
ggsave(file = "~/Data_ReefFishStability/Figs/FigS10.pdf",
		dpi = 300, width = 160, height = 160, units="mm", device=cairo_pdf)

# networks based on physical distance
id <- (master.ecoregion.fish.dat %>%
			filter(!ECOREGION%in%c("New Caledonia","Palawan/North Borneo")) %>%
			distinct(ID))$ID

phy.conn.net.tmp <- foreach(i=1:length(id)) %dopar% {
	
	ecoreg.dat <- master.ecoregion.fish.dat %>% filter((ID%in%id[i])&
							(!ECOREGION%in%c("New Caledonia","Palawan/North Borneo"))) %>%
			group_by(ECOREGION, SITE_ID, MPA, SPECIES) %>%
			summarise(abund=sum(abund), .groups="drop") %>%
			pivot_wider(names_from=SPECIES, values_from=abund, values_fill=0)
	
	if(id[i]=="bassian") net.name = "BA"
	if(id[i]=="capehowe") net.name = "CH"
	if(id[i]=="centsouthgbr") net.name = "CSGBR"
	if(id[i]=="northcalif") net.name = "LHNI"
	if(id[i]=="lordhowe") net.name = "NC"
	if(id[i]=="southaust") net.name = "SAG"
	if(id[i]=="atleu") net.name = "SEAS"
	if(id[i]=="southcalif") net.name = "SCB"
	if(id[i]=="java") net.name = "SSJS"
	if(id[i]=="northgbr") net.name = "NGBR"
	if(id[i]=="tweed") net.name = "TM"
	if(id[i]=="nwmed") net.name =  "WM"	
	
	# filter net.size.df to retain maximum spatial distance among MPAs and
	# among any two sites in the whole network
	net.size <- net.size.df %>% filter(ECOREGION%in%net.name)
	
	# load Ecoregion distance matrix (generated using script mpa_LeastCostDist.R)
	dist.dat.path <- paste("max", id[i], "dist.RData", sep=".")
	dist.dat <- mget(load(file.path('~/Data_ReefFishStability', dist.dat.path)))[[1]]
	phy.mat <- as.matrix(dist.dat)
	phy.net.tmp <- graph.adjacency(phy.mat, mode='undirected', weighted=TRUE)
	phy.inet <- mst(phy.net.tmp, weights = NULL, algorithm = NULL)
	phy.inet.plot <- ggplot_igraph(plot.graph=phy.inet, MPA=ecoreg.dat$MPA,
			mpa.extent=net.size$DIST, net.extent=net.size$DIST_RANGE, net.name=net.name)	

}

phy.conn.net <- bind_rows(phy.conn.net.tmp) %>%
		mutate(IDlabel=as.factor(paste(ID, " (", net.extent, ")", sep="")),
				IDlabel=fct_reorder(IDlabel, mpa.extent, mean)
		)

figs11 <- ggplot(data=phy.conn.net , aes(x=V1, y=V2))+
		geom_segment(data=na.exclude(phy.conn.net),
				aes(x=from.x,xend = to.x, y=from.y,yend = to.y),
				linewidth=1, colour="darkgrey") +
		geom_point(aes(color=MPA), alpha=1, size=2.5) +
		scale_color_fish_d(option = "Cirrhilabrus_solorensis", begin=0.05, end=0.95, direction=-1) +
		labs(x="", y="") +
		theme_bw() +
		facet_wrap(. ~ IDlabel, scales="free") +
		theme(
				legend.position="none",
				strip.background=element_rect(colour="grey90",
						fill="grey90"),
				plot.subtitle = element_text(size = 10),					
				panel.spacing.x=unit(0.1,"line"),
				panel.border=element_rect(colour="grey60", linewidth=0.2),
				panel.background = element_blank(),
				panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(),
				axis.text.x = element_blank(),
				axis.text.y = element_blank(),
				axis.ticks.x = element_blank(),
				axis.ticks.y = element_blank(),
				plot.margin = margin(0, 0.1, 0.1, 0.1, "cm")
		)

#windows(width=8,height=8)
figs11
dev.off()
ggsave(file = "~/Data_ReefFishStability/Figs/FigS11.pdf",
		dpi = 300, width = 160, height = 160, units="mm", device=cairo_pdf)

# relationships between closenees centrality measured on from physically-derived graphs
# and degree centrality from biologically-derived graphs
jacc.cent.df.tmp <- jacc.conn.net %>%
		rename(jacc.closeness=cc, jacc.degree=dc) %>%
		select(c(ID, MPA, mpa.extent, net.extent, jacc.closeness, jacc.degree)) %>%
		arrange(ID)
phy.cent.df.tmp <- phy.conn.net %>%
		rename(phy.closeness=cc, phy.degree=dc) %>%
		select(c(ID, phy.closeness, phy.degree)) %>%
		arrange(ID)

cent.df <- as_tibble(cbind(jacc.cent.df.tmp,
						phy.closeness=phy.cent.df.tmp[,2],
						phy.degree=phy.cent.df.tmp[,3]))

# Results for Supplementary Table 6
lm.test.centrality <- cent.df %>%
		nest(data=-c("ID")) %>%
		mutate(fit = map(data, ~ lm(jacc.degree ~ phy.closeness, data = .x)),
				tidied = map(fit, tidy),
				augmented = map(fit, augment),
				glanced = map(fit, glance)
		) %>%
		unnest(tidied) %>% # replace tidied with glanced and rerun script omitting the last two lines to view R2.
		filter(term=="phy.closeness") %>%
		dplyr::rename(slope="estimate", x.var="term", t_stat="statistic") 

#### ---- RELATIONSIPS BETWEEN STABILITY AND ASYNCHONRY MEASURES AND   ---- ####
#### ---- MPA NETWORK CHARACTERISTICS (EXTENT, N. MPAs AND N. SITES    ---- ####
#### ---- ----------------------------------------------------------   ---- ####

fit.extent <- point.est.df %>%
		mutate(MPA.EXTENT=as.numeric(as.vector(MPA.EXTENT))) %>%
		nest(data=-c("METRIC_TYPE","METRICS")) %>%
		mutate(fit = map(data, ~ lm(mepred ~ MPA.EXTENT, data = .x)),
				tidied = map(fit, tidy),
				augmented = map(fit, augment),
				glanced = map(fit, glance)
		) %>%
		unnest(tidied) %>% # replace tidied with glanced and rerun script omitting the last two lines to view R2.
		filter(term=="MPA.EXTENT") %>%
		dplyr::rename(slope="estimate", x.var="term", t_stat="statistic") 

fit.nmpa <- point.est.df %>%
		nest(data=-c("METRIC_TYPE","METRICS")) %>%
		mutate(fit = map(data, ~ lm(mepred ~ N_MPA, data = .x)),
				tidied = map(fit, tidy),
				augmented = map(fit, augment),
				glanced = map(fit, glance)
		) %>%
		unnest(tidied) %>%
		filter(term=="N_MPA") %>%
		rename(slope="estimate", x.var="term", t_stat="statistic") 

fit.nsites <- point.est.df %>%
		mutate(NSITES=MEAN_N_SITES_MPA+MEAN_N_SITES_NOMPA) %>%
		nest(data=-c("METRIC_TYPE","METRICS")) %>%
		mutate(fit = map(data, ~ lm(mepred ~ NSITES, data = .x)),
				tidied = map(fit, tidy),
				augmented = map(fit, augment),
				glanced = map(fit, glance)
		) %>%
		unnest(tidied) %>%
		filter(term=="NSITES") %>%
		rename(slope="estimate", x.var="term", t_stat="statistic") 

# Results for Supplementary Table 7
fit.extent
fit.nmpa
fit.nsites

#### ---- REPEAT META-STABILITY ANALYSIS AT THE 51-100 KM DISTANCE RANGE ---- ####

# Clean and reload data
rm(list=ls())
setwd("~/Data_ReefFishStability")
load("~/Data_ReefFishStability/metastats.res.RData")

# fliter metacommunities with a distance range of at least 100 km, 
# drop MPANET_EXTENT (is equal to DIST)

metastats.res <- metastats.res %>% group_by(ID) %>%
		filter(DIST==100)

# select metrics and prepare data for analysis
sel.metrics <- c("gamma_stab","alpha_stab","pop_stab","species_stab","loc_sp_comm_async_gross_w",
		"loc_comm_metacom_async_gross_w","meta_pop_metacom_async_gross_w",
		"meta_sp_metacom_async_gross_w")

dat.mpa <- metastats.res %>% filter(MPA=="Yes"&
						METRICS%in%sel.metrics) %>%
		rename(Obs.mpa=Obs, Nspec.mpa=NSPEC, Jmean.mpa=Jnife.mean,
				Jse.mpa=Jnife.se, Jvar.mpa=Jnife.var,
				SAMPLED_AREA_MPA=SAMPLED_AREA,
				DIST_RANGE_MPA=MPANET_EXTENT,
				MEAN_N_SITES_MPA=MEAN_N_SITES,
				MAX_N_SITES_MPA=MAX_N_SITES) %>%
		dplyr::select(-c("MPA"))

dat.oa <- metastats.res %>% filter(MPA=="No"&
						METRICS%in%sel.metrics) %>% 
		#select(-MPANET_EXTENT) %>%
		rename(Obs.oa=Obs, Nspec.oa=NSPEC, Jmean.oa=Jnife.mean,
				Jse.oa=Jnife.se, Jvar.oa=Jnife.var,
				SAMPLED_AREA_NOMPA=SAMPLED_AREA,
				DIST_RANGE_NOMPA=DIST_RANGE,
				MEAN_N_SITES_NOMPA=MEAN_N_SITES,
				MAX_N_SITES_NOMPA=MAX_N_SITES) %>%
		dplyr::select(-c("N_MPA"))

# calculate Hedge's g effect size; since the jackknife procedure uses the leave-one-out
# of indivdual species, the total number of species in each ecoregion provides
# the sample size for calculating effect sizes
eff.size.df <- dat.mpa  %>% left_join(dat.oa,
				by=c("ID","DIST","METRICS")) %>%
		mutate(
				sampled.area=SAMPLED_AREA_MPA/SAMPLED_AREA_NOMPA,
				pooled.sd=sqrt(((Nspec.mpa-1)*Jvar.mpa+(Nspec.oa-1)*Jvar.oa)/(Nspec.mpa+Nspec.oa-2)),
				dObs=(Obs.mpa-Obs.oa)/pooled.sd,
				dJnife=(Jmean.mpa-Jmean.oa)/pooled.sd,
				se.dObs=sqrt((Nspec.mpa+Nspec.oa)/(Nspec.mpa*Nspec.oa)+(dObs^2/(2*(Nspec.mpa+Nspec.oa)))),
				se.dJnife=sqrt((Nspec.mpa+Nspec.oa)/(Nspec.mpa*Nspec.oa)+(dJnife^2/(2*(Nspec.mpa+Nspec.oa)))),
				gObs=(1-3/(4*(Nspec.mpa+Nspec.oa-2)-1))*dObs,
				gJnife=(1-3/(4*(Nspec.mpa+Nspec.oa-2)-1))*dJnife,
				se.gObs=sqrt((1-3/(4*(Nspec.mpa+Nspec.oa-2)-1))^2*se.dObs^2),
				se.gJnife=sqrt((1-3/(4*(Nspec.mpa+Nspec.oa-2)-1))^2*se.dJnife^2)
		) %>% 
		filter(!dObs%in%c(-Inf,Inf)) %>%
		drop_na() %>%
		mutate(ID=recode(ID,
						atleu="SEAS",bassian="BA",capehowe="CH",centsouthgbr="CSGBR",java="SSJS",
						lordhowe="LHNI",northcalif="NC",northgbr="NGBR",nwmed="WM",
						southaust="SAG",southcalif="SCB",tweed="TM")
		)

# meta-analysis - async stats
locsp <- brm(
		gObs | se(se.gObs) ~ 1 + (1|ID),
		data = eff.size.df %>%
				filter(METRICS=="loc_sp_comm_async_gross_w"), 
		prior = c(prior(normal(0, 1), class = Intercept),
				prior(cauchy(0, 1), class = sd)),
		iter = 4000, warmup = 1000, chains = 4,
		seed=237835)

spat <- brm(
		gObs | se(se.gObs) ~ 1 + (1|ID),
		data = eff.size.df %>%
				filter(METRICS=="loc_comm_metacom_async_gross_w"), 
		prior = c(prior(normal(0, 1), class = Intercept),
				prior(cauchy(0, 1), class = sd)),
		iter = 4000, warmup = 1000, chains = 4,
		seed=237835)

metapop <- brm(
		gObs | se(se.gObs) ~ 1 + (1|ID),
		data = eff.size.df %>%
				filter(METRICS=="meta_pop_metacom_async_gross_w"), 
		prior = c(prior(normal(0, 1), class = Intercept),
				prior(cauchy(0, 1), class = sd)),
		iter = 4000, warmup = 1000, chains = 4,
		seed=237835)

metasp <- brm(
		gObs | se(se.gObs) ~ 1 + (1|ID),
		data = eff.size.df %>%
				filter(METRICS=="meta_sp_metacom_async_gross_w"), 
		prior = c(prior(normal(0, 1), class = Intercept),
				prior(cauchy(0, 1), class = sd)),
		iter = 4000, warmup = 1000, chains = 4,
		seed=237835)

#### ---- Posterior distributions ---- ####

# this is average species asynchorny metric as in Lamy et al
post.locsp <- eff.size.df %>%
		filter(METRICS=="loc_sp_comm_async_gross_w") %>%
		# Get the posterior predictions
		add_epred_draws(object = locsp) %>%
		mutate(GM=fixef(locsp)[1,1], Q2.5=fixef(locsp)[1,3], Q97.5=fixef(locsp)[1,4])
post.locsp$METRICS <- "aver.sp"

post.spat <- eff.size.df %>%
		filter(METRICS=="loc_comm_metacom_async_gross_w") %>%
		# Get the posterior predictions
		add_epred_draws(object = spat) %>%
		mutate(GM=fixef(spat)[1,1], Q2.5=fixef(spat)[1,3], Q97.5=fixef(spat)[1,4])
post.spat$METRICS <- "spatial"

# this is spatial species asynchrony as in Lamy et al
post.metapop <- eff.size.df %>%
		filter(METRICS=="meta_pop_metacom_async_gross_w") %>%
		# Get the posterior predictions
		add_epred_draws(object = metapop) %>%
		mutate(GM=fixef(metapop)[1,1], Q2.5=fixef(metapop)[1,3], Q97.5=fixef(metapop)[1,4])
post.metapop$METRICS <- "spat.sp"

# this is metapopulation asynchorny as in Lamy et al
post.metasp <- eff.size.df %>%
		filter(METRICS=="meta_sp_metacom_async_gross_w") %>%
		# Get the posterior predictions
		add_epred_draws(object = metasp) %>%
		mutate(GM=fixef(metasp)[1,1], Q2.5=fixef(metasp)[1,3], Q97.5=fixef(metasp)[1,4])
post.metasp$METRICS <- "metapop"

posterior.async.ghedge.100 <- rbind(post.locsp, post.spat, post.metapop, post.metasp)
save(posterior.async.ghedge.100, file="posterior.async.ghedge.100.RData")

#### ---- STABILITY MEASURES ---- ####

mg <- brm(
		gObs | se(se.gObs) ~ 1 +  (1|ID),
		data = eff.size.df %>%
				filter(METRICS=="gamma_stab"), 
		prior = c(prior(normal(0, 1), class = Intercept),
				prior(cauchy(0, 1), class = sd)),
		iter = 4000, warmup = 1000, chains = 4,
		seed=237835)

ma <- brm(
		gObs | se(se.gObs) ~ 1 + (1|ID),
		data = eff.size.df %>%
				filter(METRICS=="alpha_stab"), 
		prior = c(prior(normal(0, 1), class = Intercept),
				prior(cauchy(0, 1), class = sd)),
		iter = 4000, warmup = 1000, chains = 4,
		seed=237835)

ms <- brm(
		gObs | se(se.gObs) ~ 1 + (1|ID),
		data = eff.size.df %>%
				filter(METRICS=="species_stab"), 
		prior = c(prior(normal(0, 1), class = Intercept),
				prior(cauchy(0, 1), class = sd)),
		iter = 4000, warmup = 1000, chains = 4,
		seed=237835)

mp <- brm(
		gObs | se(se.gObs) ~ 1 + (1|ID),
		data = eff.size.df %>%
				filter(METRICS=="pop_stab"), 
		prior = c(prior(normal(0, 1), class = Intercept),
				prior(cauchy(0, 1), class = sd)),
		iter = 4000, warmup = 1000, chains = 4,
		seed=237835)

post.gamma <- eff.size.df %>%
		filter(METRICS=="gamma_stab") %>%
		# Get the posterior predictions
		add_epred_draws(object = mg) %>%
		mutate(GM=fixef(mg)[1,1], Q2.5=fixef(mg)[1,3], Q97.5=fixef(mg)[1,4])
post.gamma$METRICS <- "gamma"

post.alpha <- eff.size.df %>%
		filter(METRICS=="alpha_stab") %>%
		# Get the posterior predictions
		add_epred_draws(object = ma) %>%
		mutate(GM=fixef(ma)[1,1], Q2.5=fixef(ma)[1,3], Q97.5=fixef(ma)[1,4])
post.alpha$METRICS <- "alpha"

post.species <- eff.size.df %>%
		filter(METRICS=="species_stab") %>%
		# Get the posterior predictions
		add_epred_draws(object = ms) %>%
		mutate(GM=fixef(ms)[1,1], Q2.5=fixef(ma)[1,3], Q97.5=fixef(ms)[1,4])
post.species$METRICS <- "species"

post.population <- eff.size.df %>%
		filter(METRICS=="pop_stab") %>%
		# Get the posterior predictions
		add_epred_draws(object = mp) %>%
		mutate(GM=fixef(mp)[1,1], Q2.5=fixef(mp)[1,3], Q97.5=fixef(mp)[1,4])
post.population$METRICS <- "population"

posterior.stab.ghedge.100 <- rbind(post.alpha, post.gamma, post.species, post.population)
save(posterior.stab.ghedge.100, file="posterior.stab.ghedge.100.RData")

#### ---- PLOT POSTERIOR DISTRIBUTIONS ---- ####
load("~/Data_ReefFishStability/posterior.async.ghedge.100.RData")
load("~/Data_ReefFishStability/posterior.stab.ghedge.100.RData")

# Recode metrics to match definitions in Lamy et al. 2019:
# GAS = Gamma stability
# AAS: average alpha stability
# ASS: average species stability
# MPS: Metapopulation stability

# SCA: Spatial community asynchrony
# SSA: Spatial species asynchrony
# ASA: Average species asynchrony
# MPA: Metapopulation asynchrony

# prepare data for plotting
stab.df <- posterior.stab.ghedge.100 %>%
		mutate(
				METRICS=recode(METRICS, gamma="GAS", alpha="AAS",
						species="ASS", population="MPS"),
				METRICS=factor(METRICS,
						levels=c("GAS","AAS","ASS","MPS")),
				METRIC_TYPE="Stability"
		) %>%
		relocate(METRIC_TYPE, .after="METRICS")

async.df <- posterior.async.ghedge.100 %>%
		mutate(
				METRICS=recode(METRICS, metapop="MPA", spat.sp="SSA",
						aver.sp="ASA", spatial="SCA"),
				METRICS=factor(METRICS,
						levels=c("SCA","SSA","ASA","MPA")),
				METRIC_TYPE="Asynchrony"
		) %>%
		relocate(METRIC_TYPE, .after="METRICS")

metastab.df <- rbind(stab.df, async.df) %>%
		arrange(DIST_RANGE_MPA) %>%
		mutate(
				MPA.EXTENT=factor(round(DIST_RANGE_MPA, 0)),
				MPA.EXTENT=fct_reorder(MPA.EXTENT, DIST_RANGE_MPA, mean)
		) %>%
		ungroup()

point.est.df <- metastab.df %>%
		group_by(ID, METRIC_TYPE, METRICS,
				N_MPA, MEAN_N_SITES_MPA, MAX_N_SITES_MPA, 
				MEAN_N_SITES_NOMPA, MAX_N_SITES_NOMPA) %>%
		summarise(mepred=mean(.epred), .groups="drop")

pal <- c("#FF9933", "#0000CC", "#FFCC00", "#6633CC", "#990066", "#3300CC")
ID.cols <- c("SCB"=pal[1], "CSGBR"=pal[2], "SAG"=pal[3], "BA"=pal[4], "SSJS"=pal[5], "CH"=pal[6])

p.meta.stab <- metastab.df %>% filter(METRIC_TYPE=="Stability") %>%
		ggplot(aes(x = fct_reorder(MPA.EXTENT, DIST_RANGE_MPA, mean), y =.epred,
						fill=ID, col=ID)) +
		stat_eye(side="right", alpha=0.8, size=0.1, stroke=0.1,
				point_size=0.5, point_color="grey40") +
		stat_pointinterval(.width = c(.66, .95),
				point_size=0.5, point_color="grey40", 
				interval_color="grey40") +
		scale_x_discrete(labels=c("SCB","CSGBR","SAG","BA","SSJS","CH")) +
		geom_hline(yintercept = 0, col="red", linewidth=0.5, linetype=2) +
		scale_fill_manual(values=ID.cols) +
		scale_color_manual(values=ID.cols) +
		labs(x="", y="") +
		theme_bw() +
		facet_rep_wrap(~ METRICS, scales="free",
				repeat.tick.labels = c("left"), ncol=2) +
		theme(
				panel.border=element_rect(colour="grey60", linewidth=0.2),
				strip.background=element_rect(colour="grey90",
						fill="grey90"),
				legend.position="none",
				panel.spacing = unit(0.4, "lines"),
				axis.ticks = element_line(color="grey60", linewidth=0.4),
				axis.text.y  = element_text(size=10, hjust = 0),
				axis.text.x = element_text(size=10, angle = 45, vjust=0.9, hjust=0.8),
				panel.grid = element_blank()
		
		)

#windows(width=8, height=6)
p.meta.stab
dev.off()
ggsave(file = "~/Data_ReefFishStability/Figs/FigS13_stability.pdf",
		dpi = 300, width = 8, height = 6, useDingbats=FALSE)

p.meta.async <- metastab.df %>% filter(METRIC_TYPE=="Asynchrony") %>%
		ggplot(aes(x = fct_reorder(MPA.EXTENT, DIST_RANGE_MPA, mean), y =.epred,
						fill=ID, col=ID)) +
		stat_eye(side="right", alpha=0.8, size=0.1, stroke=0.1,
				point_size=0.5, point_color="grey40") +
		stat_pointinterval(.width = c(.66, .95),
				point_size=0.5, point_color="grey40", 
				interval_color="grey40") +
		scale_x_discrete(labels=c("SCB","CSGBR","SAG","BA","SSJS","CH")) +
		geom_hline(yintercept = 0, col="red", linewidth=0.5, linetype=2) +
		scale_fill_manual(values=ID.cols) +
		scale_color_manual(values=ID.cols) +
		labs(x="", y="") +
		theme_bw() +
		facet_rep_wrap(~ METRICS, scales="free",
				repeat.tick.labels = c("left"), ncol=2) +
		theme(
				panel.border=element_rect(colour="grey60", linewidth=0.2),
				strip.background=element_rect(colour="grey90",
						fill="grey90"),
				legend.position="none",
				panel.spacing = unit(0.4, "lines"),
				axis.ticks = element_line(color="grey60", linewidth=0.4),
				axis.text.y  = element_text(size=10, hjust = 0),
				axis.text.x = element_text(size=10, angle = 45, vjust=0.9, hjust=0.8),
				panel.grid = element_blank()
		
		)

#windows(width=8, height=6)
p.meta.async
dev.off()
ggsave(file = "~/Data_ReefFishStability/Figs/FigS13_asynchrony.pdf",
		dpi = 300, width = 8, height = 6, useDingbats=FALSE)

#figS13 <- ggarrange(
#		p.meta.stab, p.meta.async,
#		ncol=1, nrow=2
#		)
#
#windows(width=7,height=10)
#figS13
#dev.off()
#ggsave(file = "~/Data_ReefFishStability/Figs/FigS13.pdf",
#		dpi = 300, width = 110, height = 180, units="mm", device=cairo_pdf)




























