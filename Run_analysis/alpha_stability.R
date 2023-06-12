#### ---- Alpha stability analysis ---- ####

# require libraries
require(dplyr)
require(tidyr)
require(ggplot2)
require(fishualize)
require(ggpubr)
require(lme4)
require(lmerTest)
require(ggeffects)
require(sjPlot)
require(performance)
require(datawizard)
require(forcats)

# load data
setwd("~/Data_ReefFishStability/")
load("site_fish_stab.RData")

# source functions
source("~/workspace/ReefFishStability/MasterR/sep_fit_plot.R") # separate model fits for MPAs and OAs
source("~/workspace/ReefFishStability/MasterR/site_plot_func.R") # fit interactions between covariates x MPA vs. OA
source("~/workspace/ReefFishStability/MasterR/annotate_stats.R") # annotate plots with statistics

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

p.spacer <- ggplot() +
		theme(panel.border = element_blank())

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

#### Analysis with detrended ASYNC ---- ####
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











