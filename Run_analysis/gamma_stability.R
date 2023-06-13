#### ------------------------------ GAMMA STABILITY ANALYSIS --------------------------------------- ####
#### ---- NOTE: this analysis is slow and was originally performed on an hpc cluster. The       ---- ####
#### ---- resulting file metastats.res.RData is provided in folder ~/Data_ReefFishStability     ---- ####
#### ---- and you may want to go directly to the analysis at line 128 of this script.               ---- ####

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
require(dplyr)
require(codyn)
require(tidyr)
require(purrr)
require(broom)
require(brms)
require(tidybayes)
require(ggplot2)
require(lemon)
require(fishualize)
require(ggpubr)
require(forcats)
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

#### ---- The original analysis on HPC calculated gamma stability and asynchorny measures for three distance ranges in addition to  ---- ####
#### ---- the maximum distance among MPAs, which is the relevant measure here. The distance ranges were: 0-20, 21-50, 51-100 km     ---- ####
#### ---- and maximum distance between any two MPAs in an ecoregion. Statistics are included in the current metastats.res dataset,  ---- ####
#### ---- but are filtered here to select only those associated with the maximum distance analysis. ------------------------------------ ####

# round distances, make DIST equal to MPANET_EXTENT (they both indicate the
# maximum distance between two MPA sites; finlly, drop MPA_EXTENT.

metastats.res <- metastats.res %>% group_by(ID) %>%
		filter(DIST==max(DIST)) %>% mutate(DIST=MPANET_EXTENT) %>%
		mutate(DIST=round(DIST, 0)) %>%
		select(-MPANET_EXTENT)

# select metrics and prepare data for analysis
sel.metrics <- c("gamma_stab","alpha_stab","pop_stab","species_stab","loc_sp_comm_async_gross_w",
		"loc_comm_metacom_async_gross_w","meta_pop_metacom_async_gross_w",
		"meta_sp_metacom_async_gross_w")

dat.mpa <- metastats.res %>% filter(MPA=="Yes"&
						METRICS%in%sel.metrics) %>%
		rename(Obs.mpa=Obs, Nspec.mpa=NSPEC, Jmean.mpa=Jnife.mean,
				Jse.mpa=Jnife.se, Jvar.mpa=Jnife.var,
				SAMPLED_AREA_MPA=SAMPLED_AREA,
				DIST_RANGE_MPA=DIST_RANGE,
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
# method on indivdual species, the total number of species in each ecoregion provides
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
				MPA.EXTENT=factor(floor(DIST)),
				MPA.EXTENT=fct_reorder(MPA.EXTENT, DIST, mean)
		) %>%
		ungroup()

point.est.df <- metastab.df %>%
		group_by(ID, DIST, METRIC_TYPE, METRICS,
				N_MPA, MEAN_N_SITES_MPA, MAX_N_SITES_MPA, 
				MEAN_N_SITES_NOMPA, MAX_N_SITES_NOMPA,
				MPA.EXTENT) %>%
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
ggsave(file = "~/Data_ReefFishStability/Figs/FigS6b.pdf",
		dpi = 300, width = 8, height = 6, useDingbats=FALSE)

#### ---- RELATIONS BETWEEN STABILITY AND ASYNCHONRY MEASURES AND   ---- ####
#### ---- MPA NETWORK CHARACTERISTICS (EXTENT, N. MPAs AND N. SITES ---- ####
#### ---- --------------------------------------------------------  ---- ####
fit.extent <- point.est.df %>%
		nest(data=-c("METRIC_TYPE","METRICS")) %>%
		mutate(fit = map(data, ~ lm(mepred ~ MPA.EXTENT, data = .x)),
				tidied = map(fit, tidy),
				augmented = map(fit, augment),
				glanced = map(fit, glance)
		) %>%
		unnest(tidied) %>%
		filter(term=="MPA.EXTENT") %>%
		rename(slope="estimate", x.var="term", t_stat="statistic") 

#knitr::kable(fit.extent)

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

fit.extent
fit.nmpa
fit.nsites

#### ---- GENERATE SUPP TABLE 5 ---- ####
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



























