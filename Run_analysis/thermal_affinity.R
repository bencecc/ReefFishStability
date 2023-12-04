#### ---- Comparisons of fish thermal affinities for different trophic groups ---- ####
#### ---- by selecting species with a Species Temperature Index (STI) below   ---- #### 
#### ---- or equal-above maximum MHW intensity.                               ---- ####

require(tidyverse)
require(broom)
require(mgcv)
require(tidymv)
require(ggpubr)
require(fishualize)
require(sjPlot)
require(datawizard)
require(lme4)
require(lmerTest)
require(ggeffects)
require(performance)

setwd("~/Data_ReefFishStability/")

# load data
load("raw.fish.dat.RData")
load("site_fish_stab.RData")
load("fish.traits.dat.RData")
load("mhw.res.by.yr.RData")

master.fish.dat <- raw.fish.dat %>%
		group_by(ID,SITE_ID,LAT,LON,YEAR,MPA,MPA_STATUS,MPA_NAME,
				IUCN_CAT,NO_TAKE,STATUS_YR,STATUS,SAMP_AGE,N_YEARS,INI_YEAR,
				END_YEAR,N_REPS,N_DATES,SAMPLED_AREA,TRANSECT_SIZE,ECOREGION,
				ECO_CODE,PROVINCE,PROV_CODE,REALM,RLM_CODE,SPECIES) %>%
		summarise(abund=sum(abund), .groups="drop") 

# source functions
source("~/workspace/ReefFishStability/MasterR/sep_fit_plot.R") # separate model fits for MPAs and OAs
source("~/workspace/ReefFishStability/MasterR/site_plot_func.R") # fit interactions between covariates x MPA vs. OA
source("~/workspace/ReefFishStability/MasterR/annotate_stats.R") # annotate plots with statistics
source("~/workspace/ReefFishStability/MasterR/site_plot_gam.R") # plot results of GAMs
source("~/workspace/ReefFishStability/MasterR/sep_gam_plot.R") # plot results of GAMs fitted separately for MPAs and OAs

#### ---- SPECIES THERMAL AFFINITY VS MAX MHW INTENSITY - COMBINE DATA ---- ####
sub.mhw <- mhw.res.by.yr %>%
		filter(YEAR%in%master.fish.dat$YEAR) %>%
		dplyr::select(-c(ID,SITE_ID,MPA))

ta.df <- master.fish.dat %>%
		filter(SITE_ID%in%site_fish_stab$SITE_ID) %>%
		mutate(
				ID=as.factor(ID),
				SITE_ID=as.factor(SITE_ID),
				AREA=log(SAMPLED_AREA),
				MPA=recode(MPA, Yes="Protected", No="Unprotected"),
				MPA=factor(MPA, levels=c("Protected","Unprotected"))) %>%
		left_join(fish.traits.dat, by="SPECIES") %>%
		left_join(sub.mhw, by=c("LAT","LON","YEAR"))

#### ---- ANALYSIS AVERAGED OVER TIME ---- ####

# prepare data for temperature affinity (ta) analysis:
# define feeding categories, obtain total species abundances over time,
# distinguish species with ta above or below threshold,
# obtain total abundances by trophic and threshold categories and
# combine with site_fish_stab to restrict the analysis to those sites
# includedd in the alpha stability analysis

ta.mean <- ta.df %>%
		# consider removing site 680 (form Heenan dataset) that includes outlier fish abundances
		#filter(SITE_ID!=680) %>%
		mutate(
				FeedingType=case_when(FeedingType%in%c("browsing on substrate") ~ "Microphages",
						FeedingType%in%c("dplyr::selective plankton feeding",
								"filtering plankton","variable") ~ "Planktivores",
						FeedingType%in%c("grazing on aquatic plants") ~ "Grazers",
						FeedingType%in%c("hunting macrofauna (predator)","picking parasites off a host (cleaner)") ~ "Carnivores",
						FeedingType%in%c("other") ~ "NA",
						TRUE ~ as.character(FeedingType)
				),
				MPA=forcats::fct_relevel(MPA,
						c("Unprotected","Protected"))) %>%
		# Distinguish species based on their thermal affinity (above/below max MWH intensity)
		group_by(ID,SITE_ID,LAT,LON,MPA,SPECIES) %>%
		mutate(site.abund=sum(abund, na.rm=T)) %>%
		group_by(ID,SITE_ID,LAT,LON,MPA,FeedingType) %>%
		distinct(SPECIES, .keep_all=TRUE) %>%
		mutate(
				Thresh=case_when(
						sst_q95>=intensity_max_abs ~ "Above",
						sst_q95<intensity_max_abs ~ "Below"
				)
		) %>%
		drop_na(Thresh, FeedingType) %>%
		group_by(ID,SITE_ID,LAT,LON,MPA,Thresh,FeedingType,) %>%
		summarise(
				ABUND=sum(site.abund, na.rm=T),
				AREA=mean(AREA),
				.groups="drop") %>%
		dplyr::select(ID,SITE_ID,LAT,LON,MPA,FeedingType,Thresh,
				ABUND,AREA) %>%
		ungroup() %>%
		filter(FeedingType%in%c("Carnivores","Grazers","Microphages","Planktivores")) %>%
		right_join(site_fish_stab[,c("ID","SITE_ID","MHW")],
				by=c("ID","SITE_ID")) %>%
		mutate(
				ABUND=standardize(log(ABUND)),
				MHW=standardize(MHW),
				AREA=standardize(log(AREA)),
		) %>%
		drop_na(FeedingType)

# prepare data for analysis:
# obtain total species abundances for each year in each site
# obtain mean and standard deviations on log-transformed abundances

test.dat <- raw.fish.dat %>%
		# consider removing site 680 (form Heenan dataset) that includes outlier fish abundances
		#filter(SITE_ID!=680) %>%
		mutate(ID=as.factor(ID), SITE_ID=as.factor(SITE_ID)) %>%
		group_by(ID,SITE_ID,LAT,LON,YEAR,MPA,MPA_STATUS,MPA_NAME,
				IUCN_CAT,NO_TAKE,STATUS_YR,STATUS,SAMP_AGE,N_YEARS,INI_YEAR,
				END_YEAR,N_REPS,N_DATES,SAMPLED_AREA,TRANSECT_SIZE,ECOREGION,
				ECO_CODE,PROVINCE,PROV_CODE,REALM,RLM_CODE,SPECIES) %>%
		summarise(abund=sum(abund), .groups="drop") %>% # sum species abundances over transects
		group_by(ID,SITE_ID,LAT,LON,YEAR,MPA,SAMPLED_AREA) %>%
		summarise(abund=sum(abund), .groups="drop") %>% # total species abundances for each year in each site
		group_by(ID,SITE_ID,LAT,LON,MPA,SAMPLED_AREA) %>% # temporal mean and standard deviation
		summarise(MEAN.ABUND=mean(log(abund)), SD.ABUND=sd(log(abund)), AREA=mean(SAMPLED_AREA), .groups="drop") %>%
		select(-MPA) %>%
		right_join(site_fish_stab[,c("ID","SITE_ID","LAT","LON","MPA","MHW")],
				by=c("ID", "SITE_ID", "LAT", "LON")) %>%
		mutate(
				MEAN.ABUND=standardize(MEAN.ABUND),
				SD.ABUND=standardize(SD.ABUND),
				MHW=standardize(MHW),
				AREA=standardize(log(AREA)),
				MPA=forcats::fct_relevel(MPA,
						c("Unprotected","Protected"))
		) 		

#### ---- ANALYSIS USING lmer ---- ####

x.range <- (test.dat %>% summarise(range=range(MHW)))$range
fs9a <- sep_fit_plot(test.dat, resp="MEAN.ABUND", cov="MHW",
		x.range=x.range, y.lab=T, r2=T, plot=F)

fs9b <- sep_fit_plot(test.dat, resp="SD.ABUND", cov="MHW",
		x.range=x.range, y.lab=T, r2=T, plot=F)

figs9 <- ggarrange(
		fs9a[[1]], fs9b[[1]],
		ncol=2, nrow=1,
		align="hv",
		labels=c("a","b"),
		font.label = list(size = 10),
		label.x=0.02
)

#windows(width=6,height=4)
figs9
dev.off()
ggsave(file = "~/Data_ReefFishStability/Figs/FigS9.lmer.pdf",
		scale=0.9, dpi = 300, width = 95, height = 50, units="mm", device=cairo_pdf)

#### TROPHIC CATEGORIES ####

#### 2-way GAM ####
# fitted using tensor products 

ab.carn.gam <- gam(ABUND ~ MPA + #s(MHW,  bs="cs") +
				te(MHW, by=MPA, bs="cs") +
				s(ID, bs="re") + offset(AREA),
		data=ta.mean %>% filter(Thresh=="Above", FeedingType=="Carnivores"),
		na.action="na.omit")

#summary(ab.carn.gam)
#windows(height=5, width=6)
#gam.check(ab.carn.gam)	

bl.carn.gam <- gam(ABUND ~ MPA + #s(MHW,  bs="cs") +
				te(MHW, by=MPA, bs="cs") +
				s(ID, bs="re") + offset(AREA),
		data=ta.mean %>% filter(Thresh=="Below", FeedingType=="Carnivores"),
		na.action="na.omit")
#summary(bl.carn.gam)	

# Grazers
ab.gr.gam <- gam(ABUND ~ MPA + #s(MHW,  bs="cs") +
				te(MHW, by=MPA, bs="cs") +
				s(ID, bs="re") + offset(AREA),
		data=ta.mean %>% filter(Thresh=="Above", FeedingType=="Grazers"),
		na.action="na.omit")
#summary(ab.gr.gam)	

bl.gr.gam <- gam(ABUND ~ MPA + #s(MHW,  bs="cs") +
				te(MHW, by=MPA, bs="cs") +
				s(ID, bs="re") + offset(AREA),
		data=ta.mean %>% filter(Thresh=="Below", FeedingType=="Grazers"),
		na.action="na.omit")
#summary(bl.gr.gam)	

# Microphages
ab.mic.gam <- gam(ABUND ~ MPA + #s(MHW, bs="cs") +
				te(MHW, by=MPA, bs="cs") +
				s(ID, bs="re") + offset(AREA),
		data=ta.mean %>% filter(Thresh=="Above", FeedingType=="Microphages"),
		na.action="na.omit")	
#summary(ab.mic.gam)

bl.mic.gam <- gam(ABUND ~ MPA + #s(MHW, bs="cs") +
				te(MHW, by=MPA, bs="cs") +
				s(ID, bs="re") + offset(AREA),
		data=ta.mean %>% filter(Thresh=="Below", FeedingType=="Microphages"),
		na.action="na.omit")	
#summary(bl.mic.gam)

# Planktivores
ab.pl.gam <- gam(ABUND ~ MPA + #s(MHW,  bs="cs") +
				te(MHW, by=MPA, bs="cs") +
				s(ID, bs="re") + offset(AREA),
		data=ta.mean %>% filter(Thresh=="Above", FeedingType=="Planktivores"),
		na.action="na.omit")
#summary(ab.pl.gam)	

bl.pl.gam <- gam(ABUND ~ MPA + #s(MHW,  bs="cs") +
				te(MHW, by=MPA, bs="cs") +
				s(ID, bs="re") + offset(AREA),
		data=ta.mean %>% filter(Thresh=="Below", FeedingType=="Planktivores"),
		na.action="na.omit")
#summary(bl.pl.gam)	

#### ---- PLOT GAM RESULTS - COMPOSITE FIG 4---- ####

fig.4.3 <- site_plot_gam(ab.carn.gam, r2=T, y.lab=T, x.lab=F, plot=T)
fig.4.4 <- site_plot_gam(bl.carn.gam, r2=T, y.lab=F, x.lab=T, plot=T)

fig.4.5 <- site_plot_gam(ab.gr.gam, r2=T, y.lab=F, x.lab=F, plot=T)
fig.4.6 <- site_plot_gam(bl.gr.gam, r2=T, y.lab=F, x.lab=F, plot=T)

fig.4.7 <- site_plot_gam(ab.mic.gam, r2=T, y.lab=F, x.lab=F, plot=T)
fig.4.8 <- site_plot_gam(bl.mic.gam, r2=T, y.lab=F, x.lab=F, plot=T)

fig.4.9 <- site_plot_gam(ab.pl.gam, r2=T, y.lab=F, x.lab=F, plot=T)
fig.4.10 <- site_plot_gam(bl.pl.gam, r2=T, y.lab=F, x.lab=F, plot=T)

fig4 <- ggarrange(
		fig.4.3, fig.4.5, fig.4.7, fig.4.9,
		fig.4.4, fig.4.6, fig.4.8, fig.4.10,		
		ncol=4, nrow=2,
		align="hv",
		font.label = list(size = 10)
)


#windows(width=8,height=4)
fig4
dev.off()
ggsave(file = "~/Data_ReefFishStability/Figs/Fig4.pdf",
		scale=0.9, dpi = 300, width = 180, height = 100, units="mm")

#### ---- Extract fish silouettes ---- ####

gr <- fish.traits.dat %>% filter(FeedingType%in%"grazing on aquatic plants") %>% arrange(SPECIES)
pr <- fish.traits.dat %>% filter(FeedingType%in%"hunting macrofauna (predator)") %>% arrange(SPECIES)
pl <- fish.traits.dat %>% filter(FeedingType%in%"selective plankton feeding") %>% arrange(SPECIES)
mi <- fish.traits.dat %>% filter(FeedingType%in%"browsing on substrate") %>% arrange(SPECIES)

# see avilable silouettes
fishapes()

troph.sil <- ggplot(data.frame(x=1:10, y=1:10), aes(x=x, y=y)) +
		geom_point(col=NA) +
		# grazer
		add_fishape(family = "Siganidae",
				option = "Siganus_argenteus",
				xmin = 1, xmax = 3, ymin = 8, ymax = 10,
				#fill = fish(option = "Naso_lituratus", n = 4)[1],
				fill='chartreuse3',
				alpha = 1) +
		# predator
		add_fishape(family = "Carangidae",
				option = "Caranx_melampygus",
				xmin = 1, xmax = 3, ymin = 4, ymax = 6,
				#fill = fish(option = "Gramma_loreto", n = 54)[1],
				fill='mediumorchid2',
				alpha = 1) +
		# planktivore
		add_fishape(family = "Pomacentridae",
				option = "Chromis_iomelas",
				xmin = 4, xmax = 6, ymin = 8, ymax = 10,
				#fill = fish(option = "Gomphosus_varius", n = 4)[1],
				fill='cornflowerblue',
				alpha = 1) +
		# microphage
		add_fishape(family = "Acanthuridae",
				option = "Zebrasoma_scopas",
				xmin = 4, xmax = 6, ymin = 4, ymax = 6,
				#fill = fish(option = "Epinephelus_lanceolatus", n = 4)[1],
				fill='darkgoldenrod1',
				alpha = 1) +		
		
		theme_void()

#windows(width=20, height=20)
plot(troph.sil)		
dev.off()
ggsave(file = "~/Data_ReefFishStability/Figs/troph.sil_col.pdf",
		dpi = 300, width = 20, height = 20)











