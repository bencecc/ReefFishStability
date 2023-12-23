#### ---- Alpha stability and sampling completeness analysis ---- ####

# require libraries
require(tidyverse)
require(fishualize)
require(ggpubr)
require(lme4)
require(lmerTest)
require(ggeffects)
require(sjPlot)
require(iNEXT)
require(performance)
require(datawizard)

# load data
setwd("~/Data_ReefFishStability/")
load("site_fish_stab.RData")

# source functions
source("~/workspace/ReefFishStability/MasterR/sep_fit_plot.R") # separate model fits for MPAs and OAs
source("~/workspace/ReefFishStability/MasterR/site_plot_func.R") # fit interactions between covariates x MPA vs. OA
source("~/workspace/ReefFishStability/MasterR/onevar_fit_plot.R") # plot relationships with offset variable
source("~/workspace/ReefFishStability/MasterR/annotate_stats.R") # annotate plots with statistics
source("~/workspace/ReefFishStability/MasterR/plot.iNEXT.rar.R") # plotting function for iNEXT data to assess sampling completeness

# prepare data for analysis
test.dat <- site_fish_stab %>% 
		mutate(
				STAB=standardize(log(1/CV_TOT_ABUND)),
				SP.STAB=standardize(log(1/CV_SP_ABUND)),
				ASYNC=standardize(ASYNC_GROSS_W),
				FRIC=standardize(mFRic),
				MHW=standardize(MHW),
				REMOTENESS=standardize(REMOTENESS),
				AREA=standardize(log(SAMPLED_AREA)),
				MEAN.ABUND=standardize(log(MEAN_TOT_ABUND)),
				SD.ABUND=standardize(log(SD_TOT_ABUND)),
				MPA=fct_relevel(MPA,
						c("Unprotected","Protected"))
		) 		

#### ---- SEPARATE ANALYSIS MPAs/OAs ---- ####

# Use function sep_fit_plot with plot=F to tabulate the results of fitted models:
# the function returns a list (first element: plot; second element:stats
# for MPA; third element: stats for open areas (OA)); to run the functoin specify any
# covariate (cov argument); this only affects the plot; the full model is plotted
# as in Picewise Structural Equation Models - see below).

#### ALPHA STAB  ####

# ASYNC (analysis for plotting and to tabulate statistical results)
x.range <- (test.dat %>% summarise(range=range(ASYNC)))$range
alpha.async <- sep_fit_plot(df=test.dat, resp="STAB", cov="ASYNC", x.range=x.range, r2=T, plot=F)
alpha.async[[2]] # MPA
alpha.async[[3]] # OA
# Check model assumptions (repeat for other models)
#mpa.mod <- lmer(STAB ~ ASYNC + SP.STAB + FRIC + MHW + REMOTENESS + (1|ID) +
#				offset(AREA), data=test.dat %>% filter(MPA=="Protected"), REML=FALSE)
#oa.mod <- lmer(STAB ~ ASYNC + SP.STAB + FRIC + MHW + REMOTENESS + (1|ID) +
#				offset(AREA), data=test.dat %>% filter(MPA=="Unprotected"), REML=FALSE)
#check_model(mpa.mod)
#check_model(oa.mod)

# only for plotting (stats are the same as in alpha.async)
# SP.STAB
x.range <- (test.dat %>% summarise(range=range(SP.STAB)))$range
alpha.spstab <- sep_fit_plot(df=test.dat, resp="STAB", cov="SP.STAB", x.range=x.range, y.lab=F, plot=F)
# FRIC
x.range <- (test.dat %>% summarise(range=range(FRIC)))$range
alpha.fric <- sep_fit_plot(df=test.dat, resp="STAB", cov="FRIC", x.range=x.range, r2=T, plot=F)
# MHW
x.range <- (test.dat %>% summarise(range=range(MHW)))$range
alpha.mhw <- p.alpha.mhw <- sep_fit_plot(df=test.dat, resp="STAB", cov="MHW", x.range=x.range, y.lab=F, plot=F)
# REMOTENESS
x.range <- (test.dat %>% summarise(range=range(REMOTENESS)))$range
alpha.remot <- sep_fit_plot(df=test.dat, resp="STAB", cov="REMOTENESS", x.range=x.range, y.lab=F, plot=F)

#### SP.STAB ####
# FRIC (analysis for plotting and to tabulate statistical results)
x.range <- (test.dat %>% summarise(range=range(FRIC)))$range
spstab.fric <- sep_fit_plot(df=test.dat, resp="SP.STAB", cov="FRIC", x.range=x.range, r2=T, plot=F)
spstab.fric[[2]]
spstab.fric[[3]]

# only for plotting (stats are the same as in spstab.fric)
# MHW 
x.range <- (test.dat %>% summarise(range=range(MHW)))$range
spstab.mhw <- sep_fit_plot(df=test.dat, resp="SP.STAB", cov="MHW", x.range=x.range, r2=T, plot=F)
spstab.mhw[[2]]
spstab.mhw[[3]]
# REMOTENESS
x.range <- (test.dat %>% summarise(range=range(REMOTENESS)))$range
spstab.remot <- sep_fit_plot(df=test.dat, resp="SP.STAB", cov="REMOTENESS", x.range=x.range, y.lab=F, plot=F)

#### ASYNC ####
# FRIC (analysis for plotting and to tabulate statistical results)
x.range <- (test.dat %>% summarise(range=range(FRIC)))$range
async.fric <- sep_fit_plot(df=test.dat, resp="ASYNC", cov="FRIC", x.range=x.range, r2=T, plot=F)
async.fric[[2]]
async.fric[[3]]

# only for plotting (stats are the same as in spstab.fric)
# MHW
x.range <- (test.dat %>% summarise(range=range(MHW)))$range
async.mhw <- sep_fit_plot(df=test.dat, resp="ASYNC", cov="MHW", x.range=x.range, plot=F)
# REMOTENESS
x.range <- (test.dat %>% summarise(range=range(REMOTENESS)))$range
async.remot <- sep_fit_plot(df=test.dat, resp="ASYNC", cov="REMOTENESS", x.range=x.range, y.lab=F, plot=F)

#### FUNCIONAL RICHNESS ####
# FRIC (analysis for plotting and to tabulate statistical results)
x.range <- (test.dat %>% summarise(range=range(FRIC)))$range
fric.mhw <- sep_fit_plot(df=test.dat, resp="FRIC", cov="MHW", x.range=x.range, r2=T, plot=F)
fric.mhw[[2]]
fric.mhw[[3]]

# only for plotting (stats are the same as in spstab.fric)
fric.remot <- sep_fit_plot(df=test.dat, resp="FRIC", cov="REMOTENESS", x.range=x.range, r2=T, y.lab=F, plot=F)

# save Fig2
fig2bk <- ggarrange(
		alpha.async[[1]], alpha.spstab[[1]], alpha.mhw[[1]], alpha.remot[[1]],
		spstab.mhw[[1]], spstab.remot[[1]], async.mhw[[1]], async.remot[[1]],
		fric.mhw[[1]], fric.remot[[1]],
		ncol=4, nrow=3,
		align="hv",
		labels=c("b","c","d","e","f","g","h","i","j","k"),
		font.label = list(size = 10),
		label.x=0.1
)

#windows(width=12,height=8)
fig2bk
dev.off()
ggsave(file = "~/Data_ReefFishStability/Figs/Fig2b-k.pdf",
		dpi = 300, width = 180, height = 125, units="mm", device=cairo_pdf)

# save FigS1
figs1 <- ggarrange(
		alpha.fric[[1]], spstab.fric[[1]], async.fric[[1]],
		ncol=3, nrow=1,
		align="hv",
		labels=c("a","b","c"),
		font.label = list(size = 10),
		label.x=0.15
)

#windows(width=6,height=2)
figs1
dev.off()
ggsave(file = "~/Data_ReefFishStability/Figs/FigS1.pdf",
		dpi = 300, width = 150, height = 50, units="mm", device=cairo_pdf)

#### ----- 2-WAY MODELS ----- ####
# initialize predictions to apply site_plot_fun
x_range <- function(x) {
	xr <- range(x)
	seq(xr[1], xr[2], length.out=10)
}
# STAB
m.async.alpha.stab <- lmer(STAB ~ ASYNC*MPA + (1|ID) + offset(AREA), 
		data=test.dat, REML=FALSE)
tab_model(m.async.alpha.stab,  show.stat=T, show.ci=F, show.se=T)
# Check model assumptions (repeat for other models)
#check_model(m.async.alpha.stab)
fs2.1 <- site_plot_func(m.async.alpha.stab, r2=T)
m.sp.alpha.stab <- lmer(STAB ~ SP.STAB*MPA + (1|ID) + offset(AREA),
		data=test.dat, REML=FALSE)
tab_model(m.sp.alpha.stab,  show.stat=T, show.ci=F, show.se=T)	
fs2.2 <- site_plot_func(m.sp.alpha.stab, y.lab=F)
m.fric.alpha.stab <- lmer(STAB ~ FRIC*MPA + (1|ID) + offset(AREA),
		data=test.dat, REML=FALSE)
tab_model(m.fric.alpha.stab,  show.stat=T, show.ci=F, show.se=T)	
fs2.3 <- site_plot_func(m.fric.alpha.stab, y.lab=F)
m.mhw.alpha.stab <- lmer(STAB ~ MHW*MPA + (1|ID) + offset(AREA),
		data=test.dat, REML=FALSE)
tab_model(m.mhw.alpha.stab,  show.stat=T, show.ci=F, show.se=T)	
fs2.4 <- site_plot_func(m.mhw.alpha.stab)
m.remot.alpha.stab <- lmer(STAB ~ REMOTENESS*MPA + (1|ID) + offset(AREA),
		data=test.dat, REML=FALSE)
tab_model(m.remot.alpha.stab,  show.stat=T, show.ci=F, show.se=T)	
fs2.5 <- site_plot_func(m.remot.alpha.stab, y.lab=F)
# SP.STAB
m.fric.alpha.spstab <- lmer(SP.STAB ~ FRIC*MPA + (1|ID) + offset(AREA),
		data=test.dat, REML=FALSE)
tab_model(m.fric.alpha.spstab,  show.stat=T, show.ci=F, show.se=T)	
fs2.6 <- site_plot_func(m.fric.alpha.spstab, r2=T)
m.mhw.alpha.spstab <- lmer(SP.STAB ~ MHW*MPA + (1|ID) + offset(AREA),
		data=test.dat, REML=FALSE)
tab_model(m.mhw.alpha.spstab,  show.stat=T, show.ci=F, show.se=T)	
fs2.7 <- site_plot_func(m.mhw.alpha.spstab, y.lab=F)
m.remot.alpha.spstab <- lmer(STAB ~ REMOTENESS*MPA + (1|ID) + offset(AREA),
		data=test.dat, REML=FALSE)
tab_model(m.remot.alpha.spstab,  show.stat=T, show.ci=F, show.se=T)	
fs2.8 <- site_plot_func(m.remot.alpha.spstab, y.lab=F)
# ASYNC
m.fric.alpha.async <- lmer(ASYNC ~ FRIC*MPA + (1|ID) + offset(AREA),
		data=test.dat, REML=FALSE)
tab_model(m.fric.alpha.async,  show.stat=T, show.ci=F, show.se=T)	
fs2.9 <- site_plot_func(m.fric.alpha.async, r2=T)
m.mhw.alpha.async <- lmer(ASYNC ~ MHW*MPA + (1|ID) + offset(AREA),
		data=test.dat, REML=FALSE)
tab_model(m.mhw.alpha.async,  show.stat=T, show.ci=F, show.se=T)
fs2.10 <- site_plot_func(m.mhw.alpha.async, y.lab=F)
m.remot.alpha.async <- lmer(ASYNC ~ REMOTENESS*MPA + (1|ID) + offset(AREA),
		data=test.dat, REML=FALSE)
tab_model(m.remot.alpha.async,  show.stat=T, show.ci=F, show.se=T)	
fs2.11 <- site_plot_func(m.remot.alpha.async, y.lab=F)

# FUNCTIONAL RICHNESS
m.mhw.alpha.fric <- lmer(FRIC ~ MHW*MPA + (1|ID) + offset(AREA),
		data=test.dat, REML=FALSE)
tab_model(m.mhw.alpha.fric,  show.stat=T, show.ci=F, show.se=T)	
fs2.12 <- site_plot_func(m.mhw.alpha.fric, r2=T)
m.remot.alpha.fric <- lmer(FRIC ~ REMOTENESS*MPA + (1|ID) + offset(AREA),
		data=test.dat, REML=FALSE)
tab_model(m.remot.alpha.fric,  show.stat=T, show.ci=F, show.se=T)
fs2.13 <- site_plot_func(m.remot.alpha.fric, y.lab=F)

# FigS2
p.spacer <- ggplot() + theme(panel.border = element_blank())

figs2 <- ggarrange(
		fs2.1, fs2.2, fs2.3,
		fs2.4, fs2.5, p.spacer,
		fs2.6, fs2.7, fs2.8,
		fs2.9, fs2.10, fs2.11,
		fs2.12,fs2.13,
		ncol=3, nrow=5,
		align="hv",
		labels=c("a","b","c","d","e","","f","g","h","i","j","k","l","m"),
		font.label = list(size = 10),
		label.x=0.15
)

#windows(width=6,height=8)
figs2
dev.off()
ggsave(file = "~/Data_ReefFishStability/Figs/FigS2.pdf",
		dpi = 300, width = 150, height = 180, units="mm", device=cairo_pdf)

#### ---- Analysis with detrended ASYNC ---- ####
test.dat <- site_fish_stab %>% 
		mutate(
				STAB=standardize(log(1/CV_TOT_ABUND)),
				SP.STAB=standardize(log(1/CV_SP_ABUND)),
				ASYNC=standardize(ASYNC_GROSS_W_DETREG),
				FRIC=standardize(mFRic),
				MHW=standardize(MHW),
				REMOTENESS=standardize(REMOTENESS),
				AREA=standardize(log(SAMPLED_AREA)),
				MPA=forcats::fct_relevel(MPA,
						c("Unprotected","Protected"))
		) 		

#### ----- STAB ----- ####
x.range <- (test.dat %>% summarise(range=range(ASYNC)))$range
stab.async.det <- sep_fit_plot(df=test.dat, resp="STAB", cov="ASYNC", x.range=x.range, r2=T, plot=F)
stab.async.det[[2]]
stab.async.det[[3]]
#### ----- ASYNC ----- ####
# FRIC
x.range <- (test.dat %>% summarise(range=range(FRIC)))$range
async.fric.det <- sep_fit_plot(df=test.dat, resp="ASYNC", cov="FRIC", x.range=x.range, r2=T, plot=F)
async.fric.det[[2]]
async.fric.det[[3]]
# MHW
x.range <- (test.dat %>% summarise(range=range(MHW)))$range
async.mhw.det <- sep_fit_plot(df=test.dat, resp="ASYNC", cov="MHW", x.range=x.range, y.lab=F, plot=F)
async.mhw.det[[2]]
async.mhw.det[[3]]
# REMOTENESS
x.range <- (test.dat %>% summarise(range=range(REMOTENESS)))$range
async.remot.det <- sep_fit_plot(df=test.dat, resp="ASYNC", cov="REMOTENESS", x.range=x.range, y.lab=F, plot=F)
async.remot.det[[2]]
async.remot.det[[3]]

#### Analysis with  ASYNC  = LOREAU_DEMAZANCOURT ---- ####
test.dat <- site_fish_stab %>% 
		mutate(
				STAB=standardize(log(1/CV_TOT_ABUND)),
				SP.STAB=standardize(log(1/CV_SP_ABUND)),
				ASYNC=standardize(ASYNC_LOREAU_SQRT),
				FRIC=standardize(mFRic),
				MHW=standardize(MHW),
				REMOTENESS=standardize(REMOTENESS),
				AREA=standardize(log(SAMPLED_AREA)),
				MPA=forcats::fct_relevel(MPA,
						c("Unprotected","Protected"))
		) 		

#### ----- STAB ----- ####
x.range <- (test.dat %>% summarise(range=range(ASYNC)))$range
stab.async.lm <- sep_fit_plot(df=test.dat, resp="STAB", cov="ASYNC", x.range=x.range, r2=T, plot=F)
stab.async.lm[[2]]
stab.async.lm[[3]]

#### ----- ASYNC ----- ####
# FRIC
x.range <- (test.dat %>% summarise(range=range(FRIC)))$range
async.fric.lm <- sep_fit_plot(df=test.dat, resp="ASYNC", cov="FRIC", x.range=x.range, r2=T, plot=F)
async.fric.lm[[2]]
async.fric.lm[[3]]
# MHW
x.range <- (test.dat %>% summarise(range=range(MHW)))$range
async.mhw.lm <- sep_fit_plot(df=test.dat, resp="ASYNC", cov="MHW", x.range=x.range, y.lab=F, plot=F)
async.mhw.lm[[2]]
async.mhw.lm[[3]]
# REMOTENESS
x.range <- (test.dat %>% summarise(range=range(REMOTENESS)))$range
async.remot.lm <- sep_fit_plot(df=test.dat, resp="ASYNC", cov="REMOTENESS", x.range=x.range, y.lab=F, plot=F)
async.remot.lm[[2]]
async.remot.lm[[3]]

#### ---- PLOT FIG S3 ---- ####
figs3 <- ggarrange(
		stab.async.det[[1]], async.fric.det[[1]], async.mhw.det[[1]], async.remot.det[[1]],
		stab.async.lm[[1]], async.fric.lm[[1]], async.mhw.lm[[1]], async.remot.lm[[1]],
		ncol=4, nrow=2,
		align="hv",
		labels=c("a","b","c","d","e","f","g","h"),
		font.label = list(size = 10),
		label.x=0.15
)

#windows(width=8,height=4)
figs3
dev.off()
ggsave(file = "~/Data_ReefFishStability/Figs/FigS3.pdf",
		dpi = 300, width = 180, height = 90, units="mm", device=cairo_pdf)

#### Analysis with cumulative instead of mean MHW intensity  ---- ####
test.dat <- site_fish_stab %>% 
		mutate(
				STAB=standardize(log(1/CV_TOT_ABUND)),
				SP.STAB=standardize(log(1/CV_SP_ABUND)),
				ASYNC=standardize(ASYNC_LOREAU_SQRT),
				FRIC=standardize(mFRic),
				MHW=standardize(MHW.CUM),
				REMOTENESS=standardize(REMOTENESS),
				AREA=standardize(log(SAMPLED_AREA)),
				MPA=forcats::fct_relevel(MPA,
						c("Unprotected","Protected"))
		) 		

#### ALPHA STAB  ####
x.range <- (test.dat %>% summarise(range=range(ASYNC)))$range
alpha.async.mhwcum <- sep_fit_plot(df=test.dat, resp="STAB", cov="ASYNC", x.range=x.range, r2=T, plot=F)
# SP.STAB
x.range <- (test.dat %>% summarise(range=range(SP.STAB)))$range
alpha.spstab.mhwcum <- sep_fit_plot(df=test.dat, resp="STAB", cov="SP.STAB", x.range=x.range, y.lab=F, plot=F)
# FRIC
x.range <- (test.dat %>% summarise(range=range(FRIC)))$range
alpha.fric.mhwcum <- sep_fit_plot(df=test.dat, resp="STAB", cov="FRIC", x.range=x.range, y.lab=F, plot=F)
# MHW
x.range <- (test.dat %>% summarise(range=range(MHW)))$range
alpha.mhw.mhwcum <- p.alpha.mhw <- sep_fit_plot(df=test.dat, resp="STAB", cov="MHW", x.range=x.range, plot=F)
# REMOTENESS
x.range <- (test.dat %>% summarise(range=range(REMOTENESS)))$range
alpha.remot.mhwcum <- sep_fit_plot(df=test.dat, resp="STAB", cov="REMOTENESS", x.range=x.range, y.lab=F, plot=F)

#### SP.STAB ####
# FRIC (analysis for plotting and to tabulate statistical results)
x.range <- (test.dat %>% summarise(range=range(FRIC)))$range
spstab.fric.mhwcum <- sep_fit_plot(df=test.dat, resp="SP.STAB", cov="FRIC", x.range=x.range, r2=T, plot=F)
# MHW 
x.range <- (test.dat %>% summarise(range=range(MHW)))$range
spstab.mhwcum <- sep_fit_plot(df=test.dat, resp="SP.STAB", cov="MHW", x.range=x.range, r2=F, y.lab=F, plot=F)
# REMOTENESS
x.range <- (test.dat %>% summarise(range=range(REMOTENESS)))$range
spstab.remot.mhwcum <- sep_fit_plot(df=test.dat, resp="SP.STAB", cov="REMOTENESS", x.range=x.range, y.lab=F, plot=F)
#### ASYNC ####
# FRIC (analysis for plotting and to tabulate statistical results)
x.range <- (test.dat %>% summarise(range=range(FRIC)))$range
async.fric.mhwcum <- sep_fit_plot(df=test.dat, resp="ASYNC", cov="FRIC", x.range=x.range, r2=T, plot=F)
# MHW
x.range <- (test.dat %>% summarise(range=range(MHW)))$range
async.mhwcum <- sep_fit_plot(df=test.dat, resp="ASYNC", cov="MHW", x.range=x.range, y.lab=F, plot=F)
# REMOTENESS
x.range <- (test.dat %>% summarise(range=range(REMOTENESS)))$range
async.remot.mhwcum <- sep_fit_plot(df=test.dat, resp="ASYNC", cov="REMOTENESS", x.range=x.range, y.lab=F, plot=F)
#### FUNCIONAL RICHNESS ####
# FRIC 
x.range <- (test.dat %>% summarise(range=range(MHW)))$range
fric.mhwcum <- sep_fit_plot(df=test.dat, resp="FRIC", cov="MHW", x.range=x.range, r2=T, plot=F)
x.range <- (test.dat %>% summarise(range=range(REMOTENESS)))$range
fric.remot.mhwcum <- sep_fit_plot(df=test.dat, resp="FRIC", cov="REMOTENESS", x.range=x.range, r2=F, y.lab=F, plot=F)

# FigS4
p.spacer <- ggplot() + theme(panel.border = element_blank())

figs4 <- ggarrange(
		alpha.async.mhwcum[[1]], alpha.spstab.mhwcum[[1]], alpha.fric.mhwcum[[1]],
		alpha.mhw.mhwcum[[1]], alpha.remot.mhwcum[[1]], p.spacer,
		spstab.fric.mhwcum[[1]], spstab.mhwcum[[1]], spstab.remot.mhwcum[[1]],
		async.fric.mhwcum[[1]], async.mhwcum[[1]], async.remot.mhwcum[[1]],
		fric.mhwcum[[1]], fric.remot.mhwcum[[1]],
		ncol=3, nrow=5,
		align="hv",
		labels=c("a","b","c","d","e","","f","g","h","i","j","k","l","m"),
		font.label = list(size = 10),
		label.x=0.11
)

#windows(width=6,height=8)
figs4
dev.off()
ggsave(file = "~/Data_ReefFishStability/Figs/FigS4.pdf",
		dpi = 300, width = 150, height = 180, units="mm", device=cairo_pdf)

#### ---- Sensitivity analysis for the offset ---- ####
# Alternative analysis to control for sampling effort using log-response ratios
# obtained by dividing each response variable directly by the total area sampled
# at each site, instead of using an offset. The two approaches (offset vs. log-response
# ratios) are equivalent in gaussian linear models. This equivalence does not hold for
# generalized models (Crawley M.J. 2013. The R Book. Wiley, Second Edition, p. 415 and
# p. 566). The analysis below shows this equivalence (consistent results), which means
# that we could have examined log-response ratios directly instead of using the offset
# in the paper. We opted to present results based on the offset in the main text, since
# standardization improved data visualization. The use of an offset assumes that there
# is a significant linear relationhip between sampling effort and the response variable,
# with a fixed coefficient of 1. The analysis of log-response ratios does not assign
# any fixed coefficient to sampling effort, since it is part of the response variable.

# IMPORTANT: remove the offset from the models in function sep_fit_plot before running
# this analysis and remember to place the offset back once terminated.

#### ALPHA STAB  ####
alpha.dat <- site_fish_stab %>% 
		mutate(
				STAB=log(1/CV_TOT_ABUND/SAMPLED_AREA),
				SP.STAB=standardize(log(1/CV_SP_ABUND)),
				ASYNC=standardize(ASYNC_GROSS_W),
				FRIC=standardize(mFRic),
				SR=standardize(N_SPEC),
				MHW=standardize(MHW),
				REMOTENESS=standardize(REMOTENESS),
					MPA=fct_relevel(MPA,
						c("Unprotected","Protected"))
		) 		

x.range <- (alpha.dat %>% summarise(range=range(ASYNC)))$range
alpha1.async <- sep_fit_plot(df=alpha.dat, resp="STAB", cov="ASYNC", x.range=x.range, r2=T, plot=F)
x.range <- (alpha.dat %>% summarise(range=range(SP.STAB)))$range
alpha1.spstab <- sep_fit_plot(df=alpha.dat, resp="STAB", cov="SP.STAB", x.range=x.range, y.lab=F, plot=F)
x.range <- (alpha.dat %>% summarise(range=range(FRIC)))$range
alpha1.fric <- sep_fit_plot(df=alpha.dat, resp="STAB", cov="FRIC", x.range=x.range, y.lab=F, plot=F)
x.range <- (alpha.dat %>% summarise(range=range(MHW)))$range
alpha1.mhw <- sep_fit_plot(df=alpha.dat, resp="STAB", cov="MHW", x.range=x.range, plot=F)
x.range <- (alpha.dat %>% summarise(range=range(REMOTENESS)))$range
alpha1.remot <- sep_fit_plot(df=alpha.dat, resp="STAB", cov="REMOTENESS", x.range=x.range, y.lab=F, plot=F)

#### SP.STAB ####
sp.dat <- site_fish_stab %>% 
		mutate(
				SP.STAB=log(1/CV_SP_ABUND/SAMPLED_AREA),
				ASYNC=standardize(ASYNC_GROSS_W),
				FRIC=standardize(mFRic),
				SR=standardize(N_SPEC),
				MHW=standardize(MHW),
				REMOTENESS=standardize(REMOTENESS),
				MPA=fct_relevel(MPA,
						c("Unprotected","Protected"))
		) 		

x.range <- (sp.dat %>% summarise(range=range(FRIC)))$range
spstab1.fric <- sep_fit_plot(df=sp.dat, resp="SP.STAB", cov="FRIC", x.range=x.range, r2=T, plot=F)
x.range <- (sp.dat %>% summarise(range=range(MHW)))$range
spstab1.mhw <- sep_fit_plot(df=sp.dat, resp="SP.STAB", cov="MHW", x.range=x.range, y.lab=F, plot=F)
x.range <- (sp.dat %>% summarise(range=range(REMOTENESS)))$range
spstab1.remot <- sep_fit_plot(df=sp.dat, resp="SP.STAB", cov="REMOTENESS", x.range=x.range, y.lab=F, plot=F)

#### ASYNCHRONY ####
async.dat <- site_fish_stab %>% 
		mutate(
				ASYNC=log(ASYNC_LOREAU/SAMPLED_AREA), # use LM for asynchrony, which can be log-transfoermed (range 0-1)
				FRIC=standardize(mFRic),
				SR=standardize(N_SPEC),
				MHW=standardize(MHW),
				REMOTENESS=standardize(REMOTENESS),
				MPA=fct_relevel(MPA,
						c("Unprotected","Protected"))
		) 		

x.range <- (async.dat %>% summarise(range=range(FRIC)))$range
async1.fric <- sep_fit_plot(df=async.dat, resp="ASYNC", cov="FRIC", x.range=x.range, r2=T, plot=F)
x.range <- (async.dat %>% summarise(range=range(MHW)))$range
async1.mhw <- sep_fit_plot(df=async.dat, resp="ASYNC", cov="MHW", x.range=x.range, y.lab=F, plot=F)
x.range <- (async.dat %>% summarise(range=range(REMOTENESS)))$range
async1.remot <- sep_fit_plot(df=async.dat, resp="ASYNC", cov="REMOTENESS", x.range=x.range, y.lab=F, plot=F)

#### FUNCTIONAL RICHNESS ####
fric.dat <- site_fish_stab %>% 
		mutate(
				FRIC=log(mFRic/SAMPLED_AREA),
				SR=standardize(N_SPEC),
				MHW=standardize(MHW),
				REMOTENESS=standardize(REMOTENESS),
				MPA=fct_relevel(MPA,
						c("Unprotected","Protected"))
		) 		
x.range <- (fric.dat %>% summarise(range=range(MHW)))$range
fric1.mhw <- sep_fit_plot(df=fric.dat, resp="FRIC", cov="MHW", x.range=x.range, r2=T, plot=F)
x.range <- (fric.dat %>% summarise(range=range(REMOTENESS)))$range
fric1.remot <- sep_fit_plot(df=fric.dat, resp="FRIC", cov="REMOTENESS", x.range=x.range, r2=F, y.lab=F, plot=F)

# FigS5
p.spacer <- ggplot() + theme(panel.border = element_blank())

figs5 <- ggarrange(
		alpha1.async[[1]], alpha1.spstab[[1]], alpha1.fric[[1]],
		alpha1.mhw[[1]], alpha1.remot[[1]], p.spacer,
		spstab1.fric [[1]], spstab1.mhw[[1]], spstab1.remot[[1]],
		async1.fric[[1]], async1.mhw[[1]], async1.remot[[1]],
		fric1.mhw[[1]], fric1.remot[[1]],
		ncol=3, nrow=5,
		align="hv",
		labels=c("a","b","c","d","e","","f","g","h","i","j","k","l","m"),
		font.label = list(size = 10),
		label.x=0.11
)

windows(width=6,height=8)
figs5
dev.off()
ggsave(file = "~/Data_ReefFishStability/Figs/FigS5.pdf",
		dpi = 300, width = 150, height = 180, units="mm", device=cairo_pdf)

#### ---- Correlation among functional diversity indices ---- ####

cor.test(site_fish_stab$mFRic, site_fish_stab$mFDiv)
cor.test(site_fish_stab$mFRic, site_fish_stab$mFEv)
cor.test(site_fish_stab$mFRic, site_fish_stab$mFDis)
cor.test(site_fish_stab$mFRic, site_fish_stab$mFSpe)
cor.test(site_fish_stab$mFRic, site_fish_stab$mFOri)

cor.test(site_fish_stab$mFDiv, site_fish_stab$mFEv)
cor.test(site_fish_stab$mFDiv, site_fish_stab$mFDis)
cor.test(site_fish_stab$mFDiv, site_fish_stab$mFSpe)
cor.test(site_fish_stab$mFDiv, site_fish_stab$mFOri)

cor.test(site_fish_stab$mFEv, site_fish_stab$mFDis)
cor.test(site_fish_stab$mFEv, site_fish_stab$mFSpe)
cor.test(site_fish_stab$mFEv, site_fish_stab$mFOri)

cor.test(site_fish_stab$mFDis, site_fish_stab$mFSpe)
cor.test(site_fish_stab$mFDis, site_fish_stab$mFOri)

# Sensitivity analysis to assess the robustness of results to the varying levels of taxonomic
# scope used by different sampling programs (IDs). Most programs provide exhaustive species
# lists, but some use pre-determined target lists. Here, only sampling programs targeting
# more than 50 species are included and the relationships of alpha and species stability,
# asynchrony and functional richenss with MHWs is reassessed for MPAs and OAs.

# Calculate the total number of species targeted by each sampling program; use
# master.fish.dat, but filter for the sites used in the alpha stability analysis.
n.sp <- master.fish.dat %>%
		filter(SITE_ID%in%site_fish_stab$SITE_ID) %>%
		group_by(ID) %>%
		distinct(SPECIES, .keep_all=T) %>%
		mutate(n.sp=n()) %>%
		distinct(ID, .keep_all=T) %>% ungroup() %>%
		select(ID, n.sp)

# remove IDs with 50 species or less
sp.red <- n.sp %>% filter(n.sp>50) %>%
		mutate(ID=as.factor(ID))

# repeat analysis and extract main results (MHWs) 
red.dat <- site_fish_stab %>% 
		mutate(
				STAB=standardize(log(1/CV_TOT_ABUND)),
				SP.STAB=standardize(log(1/CV_SP_ABUND)),
				ASYNC=standardize(ASYNC_GROSS_W),
				FRIC=standardize(mFRic),
				MHW=standardize(MHW),
				REMOTENESS=standardize(REMOTENESS),
				AREA=standardize(log(SAMPLED_AREA)),
				MEAN.ABUND=standardize(log(MEAN_TOT_ABUND)),
				SD.ABUND=standardize(log(SD_TOT_ABUND)),
				MPA=fct_relevel(MPA,
						c("Unprotected","Protected"))
		) %>% filter(ID%in%sp.red$ID)	

x.range <- (red.dat %>% summarise(range=range(MHW)))$range
alpha.stab.mhw.red <- sep_fit_plot(df=red.dat, resp="STAB", cov="MHW", x.range=x.range, r2=T, y.lab=T, plot=F)
alpha.spstab.mhw.red <- sep_fit_plot(df=red.dat, resp="SP.STAB", cov="MHW", x.range=x.range, r2=T, y.lab=T, plot=F)
async.mhw.red <- sep_fit_plot(df=red.dat, resp="ASYNC", cov="MHW", x.range=x.range, r2=T, y.lab=T, plot=F)
fric.mhw.red <- sep_fit_plot(df=red.dat, resp="FRIC", cov="MHW", x.range=x.range, r2=T, y.lab=T, plot=F)

figs6 <- ggarrange(
		alpha.stab.mhw.red[[1]], alpha.spstab.mhw.red[[1]],
		async.mhw.red[[1]], fric.mhw.red[[1]],
		ncol=2, nrow=2,
		align="hv",
		labels=c("a","b","c","d"),
		font.label = list(size = 10),
		label.x=0.12
)

#windows(width=6,height=5)
figs6
dev.off()
ggsave(file = "~/Data_ReefFishStability/Figs/FigS6.pdf",
		dpi = 300, width = 110, height = 90, units="mm", device=cairo_pdf)

#### ---- Check sampling completeness by transect size ---- ####
# load data
setwd("~/Data_ReefFishStability/")
load("master.fish.dat.RData")

fish.sp.abund <- master.fish.dat %>% 
		filter(SITE_ID%in%site_fish_stab$SITE_ID) %>%
		#mutate(abund=exp(abund)-1) %>%
		mutate(TRANSECT_SIZE=case_when(
						TRANSECT_SIZE<=125 ~ 100,
						TRANSECT_SIZE>125&TRANSECT_SIZE<=200 ~ 180,
						TRANSECT_SIZE>200&TRANSECT_SIZE<500 ~ 250,
						TRANSECT_SIZE==500 ~ 500,
						TRANSECT_SIZE>500 ~ 2000, 
						TRUE ~ as.numeric(TRANSECT_SIZE))
		)
inext.res <- list()

for(i in 1:length(unique(fish.sp.abund$TRANSECT_SIZE))) {
	
	test.dat <- fish.sp.abund %>%
			filter(TRANSECT_SIZE==unique(fish.sp.abund$TRANSECT_SIZE)[i]) %>%
			group_by(SITE_ID,MPA,SPECIES) %>%
			summarise(abund=sum(abund), .groups="drop") %>%
			mutate(abund=1) %>%
			group_by(MPA) %>%
			pivot_wider(names_from="SITE_ID",
					values_from="abund", values_fill=0) %>%
			arrange(MPA,SPECIES)
	
	prot <- unique(test.dat$MPA)
	
	if(length(prot)==1) {
		
		if(prot=="No") {
			
			nompa.dat.tmp <- test.dat %>% ungroup() %>%
					dplyr::select(-c(MPA,SPECIES)) %>%
					mutate(total=rowSums(.)) %>% filter(total>0)
			nompa.dat <- c(ncol(nompa.dat.tmp)-1, nompa.dat.tmp$total)
			list.dat <- list(Unprotected=nompa.dat)
			
		}
		
		if(prot=="Yes") {
			
			mpa.dat.tmp <- test.dat %>% ungroup() %>%
					dplyr::select(-c(MPA,SPECIES)) %>%
					mutate(total=rowSums(.)) %>% filter(total>0)
			mpa.dat <- c(ncol(mpa.dat.tmp)-1, mpa.dat.tmp$total)
			list.dat <- list(Protected=mpa.dat)
			
		}
		
	}
	
	else {
		mpa.dat.tmp <- test.dat %>% ungroup() %>%
				filter(MPA=="Yes") %>% dplyr::select(-c(MPA,SPECIES)) %>%
				mutate(total=rowSums(.)) %>% filter(total>0)
		mpa.dat <- c(ncol(mpa.dat.tmp)-1, mpa.dat.tmp$total)
		
		nompa.dat.tmp <- test.dat %>% ungroup() %>%
				filter(MPA=="No") %>% dplyr::select(-c(MPA,SPECIES)) %>%
				mutate(total=rowSums(.)) %>% filter(total>0)
		nompa.dat <- c(ncol(nompa.dat.tmp)-1, nompa.dat.tmp$total)
		
		list.dat <- list(Protected=mpa.dat, Unprotected=nompa.dat)
		
	}
	
	inext.tmp <- try(iNEXT(list.dat, q=0, datatype="incidence_freq",
					nboot=100), silent=T)
	
	if (inherits(inext.tmp, "try-error")) {
		inext.plot <- NULL
	}
	
	else {
		
		inext.plot <- plot.iNEXT.rar(inext.tmp, plot.type="coverage")
		
	}
	
	inext.res[[i]] <- inext.plot
	
}

names(inext.res) <- unique(fish.sp.abund$TRANSECT_SIZE)
inext.div.mpa.plot <- inext.res

div.coverage <- ggarrange(
		inext.div.mpa.plot[[3]], inext.div.mpa.plot[[5]], inext.div.mpa.plot[[1]], 
		inext.div.mpa.plot[[4]], inext.div.mpa.plot[[2]],
		ncol=3, nrow=2,
		#align="hv",
		#labels=c("100","180","250","500","2000"),
		font.label = list(size = 10),
		#common.legend=TRUE,
		label.x=0.5,
		label.y=1.05
)


#windows(width=7, height=4)
div.coverage
dev.off()
ggsave(file = "~/Data_ReefFishStability/Figs/FigS7.pdf",
		dpi = 300, width = 7, height = 4, useDingbats=FALSE)

# evaluate proportion of transects in the 180 m^2 category.
size.180 <- fish.sp.abund %>%
		filter(TRANSECT_SIZE==180) %>%
		distinct(SITE_ID)
size.not180 <-  fish.sp.abund %>%
		filter(TRANSECT_SIZE!=180) %>%
		distinct(SITE_ID)

ratio.sites <- nrow(size.180)/nrow(size.not180) # 2.2%

# evaluate maximum coverage
max.cov <- NULL
for(i in 1:length(inext.div.mpa.plot)) {
	tmp.df <- inext.div.mpa.plot[[i]][[1]] %>%
			group_by(data.name) %>%
			summarise(max.cov=max(diversity))
	tmp.out <- cbind(
			transect.size=names(inext.div.mpa.plot)[i],
			tmp.df
	)
	max.cov <- rbind(max.cov, tmp.out)
}

#### -------------------------------------------------------------- ####

# how many MHWs before 2011 (end of climatology)?
load('mhw.res.by.yr.RData')

from.2012 <- site_fish_stab %>% filter(END_YEAR<2011)

ntot.mhws <- mhw.res.by.yr %>%
		summarise(n.mhws=sum(num_events, na.rm=T))
nmhw.pre.2011 <- mhw.res.by.yr %>%
		filter(SITE_ID%in%from.2012$SITE_ID) %>%
		summarise(n.mhws=sum(num_events, na.rm=T))
# 830/46976

# Repeat whole analysis excluding the 40 sites where sampling ended before 2011
# to generate FigS15.

test.dat <- site_fish_stab %>% 
		filter(END_YEAR>=2012) %>%  
		mutate(
				STAB=standardize(log(1/CV_TOT_ABUND)),
				SP.STAB=standardize(log(1/CV_SP_ABUND)),
				ASYNC=standardize(ASYNC_GROSS_W),
				FRIC=standardize(mFRic),
				SR=standardize(log(N_SPEC)),
				MHW=standardize(MHW),
				REMOTENESS=standardize(REMOTENESS),
				AREA=standardize(log(SAMPLED_AREA)),
				MEAN.ABUND=standardize(log(MEAN_TOT_ABUND)),
				SD.ABUND=standardize(log(SD_TOT_ABUND)),
				MPA=fct_relevel(MPA,
						c("Unprotected","Protected"))
		) 		

#### ALPHA STAB  ####

# ASYNC
x.range <- (test.dat %>% summarise(range=range(ASYNC)))$range
alpha2.async <- sep_fit_plot(df=test.dat, resp="STAB", cov="ASYNC", x.range=x.range, r2=T, plot=F)
# SP.STAB
x.range <- (test.dat %>% summarise(range=range(SP.STAB)))$range
alpha2.spstab <- sep_fit_plot(df=test.dat, resp="STAB", cov="SP.STAB", x.range=x.range, y.lab=F, plot=F)
# FRIC
x.range <- (test.dat %>% summarise(range=range(FRIC)))$range
alpha2.fric <- sep_fit_plot(df=test.dat, resp="STAB", cov="FRIC", x.range=x.range, y.lab=F, plot=F)
# MHW
x.range <- (test.dat %>% summarise(range=range(MHW)))$range
alpha2.mhw <- sep_fit_plot(df=test.dat, resp="STAB", cov="MHW", x.range=x.range, plot=F)
# REMOTENESS
x.range <- (test.dat %>% summarise(range=range(REMOTENESS)))$range
alpha2.remot <- sep_fit_plot(df=test.dat, resp="STAB", cov="REMOTENESS", x.range=x.range, y.lab=F, plot=F)

#### SP.STAB ####
# FRIC
x.range <- (test.dat %>% summarise(range=range(FRIC)))$range
spstab2.fric <- sep_fit_plot(df=test.dat, resp="SP.STAB", cov="FRIC", x.range=x.range, r2=T, plot=F)
# MHW 
x.range <- (test.dat %>% summarise(range=range(MHW)))$range
spstab2.mhw <- sep_fit_plot(df=test.dat, resp="SP.STAB", cov="MHW", x.range=x.range, y.lab=F, plot=F)

# REMOTENESS
x.range <- (test.dat %>% summarise(range=range(REMOTENESS)))$range
spstab2.remot <- sep_fit_plot(df=test.dat, resp="SP.STAB", cov="REMOTENESS", x.range=x.range, y.lab=F, plot=F)

#### ASYNC ####
# FRIC
x.range <- (test.dat %>% summarise(range=range(FRIC)))$range
async2.fric <- sep_fit_plot(df=test.dat, resp="ASYNC", cov="FRIC", x.range=x.range, r2=T, plot=F)
# MHW
x.range <- (test.dat %>% summarise(range=range(MHW)))$range
async2.mhw <- sep_fit_plot(df=test.dat, resp="ASYNC", cov="MHW", x.range=x.range, y.lab=F, plot=F)
# REMOTENESS
x.range <- (test.dat %>% summarise(range=range(REMOTENESS)))$range
async2.remot <- sep_fit_plot(df=test.dat, resp="ASYNC", cov="REMOTENESS", x.range=x.range, y.lab=F, plot=F)

#### FUNCIONAL RICHNESS ####
# FRIC
x.range <- (test.dat %>% summarise(range=range(FRIC)))$range
fric2.mhw <- sep_fit_plot(df=test.dat, resp="FRIC", cov="MHW", x.range=x.range, r2=T, plot=F)
fric2.remot <- sep_fit_plot(df=test.dat, resp="FRIC", cov="REMOTENESS", x.range=x.range, r2=F, y.lab=F, plot=F)

# FigS15
p.spacer <- ggplot() + theme(panel.border = element_blank())

figS15 <- ggarrange(
		alpha2.async[[1]], alpha2.spstab[[1]], alpha2.fric[[1]],
		alpha2.mhw[[1]], alpha2.remot[[1]], p.spacer,
		spstab2.fric [[1]], spstab2.mhw[[1]], spstab2.remot[[1]],
		async2.fric[[1]], async2.mhw[[1]], async2.remot[[1]],
		fric2.mhw[[1]], fric2.remot[[1]],
		ncol=3, nrow=5,
		align="hv",
		labels=c("a","b","c","d","e","","f","g","h","i","j","k","l","m"),
		font.label = list(size = 10),
		label.x=0.11)

#windows(width=6,height=8)
figS15
dev.off()
ggsave(file = "~/Data_ReefFishStability/Figs/FigS15.pdf",
		dpi = 300, width = 150, height = 180, units="mm", device=cairo_pdf)


#### ----------------------------------------------------------------------------------------------- ####





