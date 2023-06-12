#### ---- Piecewise Sctrutural Equation Mmdels on alpha stability data ---- ####

# load required libraries
require(dplyr)
require(tidyr)
require(piecewiseSEM)
require(semEff)
require(lme4)
require(lmerTest)
require(forcats)
require(stringr)
require(datawizard)
require(ggplot2)
require(fishualize)

# load results of alpha stability analysis
setwd("~/Data_ReefFishStability/")
load("site_fish_stab.RData")

# prepare data for pSEM analysis
sem.dat <- site_fish_stab %>%
		mutate(
				STAB=standardize(log(1/CV_TOT_ABUND)),
				SP.STAB=standardize(log(1/CV_SP_ABUND)),
				ASYNC=standardize(ASYNC_GROSS_W),
				FRIC=standardize(mFRic),
				MHW=standardize(MHW),
				REMOTENESS=standardize(REMOTENESS),
				AREA=standardize(log(SAMPLED_AREA))
		) 		

sem.mpa <- sem.dat %>% filter(MPA=="Protected") 
sem.oa <- sem.dat %>% filter(MPA=="Unprotected") 

# analysis of MPA data
mod.mpa <- psem(
		STAB <- lmer(STAB ~ ASYNC + SP.STAB + FRIC + MHW + REMOTENESS + (1|ID),
				offset=sem.mpa[,"AREA"],
				data=sem.mpa,
				#control = contr.optim,
				REML=FALSE),
		ASYNC <- lmer(ASYNC ~ FRIC + MHW + REMOTENESS + (1|ID),
				offset=sem.mpa[,"AREA"],
				data=sem.mpa,
				#control = contr.optim,
				REML=FALSE),
		SP.STAB <- lmer(SP.STAB ~ FRIC + MHW + REMOTENESS + (1|ID),
				offset=sem.mpa[,"AREA"],
				data=sem.mpa,
				#control = contr.optim,
				REML=FALSE),
		FRIC <- lmer(FRIC ~ MHW + (1|ID), # REMOTENESS removed (p>0.6) to improve model fit 
				offset=sem.mpa[,"AREA"],
				data=sem.mpa,
				#control = contr.optim,
				REML=FALSE)
)

modmpa <- summary(mod.mpa)
modmpa
coefs(mod.mpa, standardize="scale")

mod.oa <-  psem(
		STAB <- lmer(STAB ~ ASYNC + SP.STAB + FRIC + MHW + REMOTENESS + (1|ID),
				offset=sem.oa[,"AREA"],
				data=sem.oa,
				#control = contr.optim,
				REML=FALSE),
		ASYNC <- lmer(ASYNC ~ FRIC + MHW + REMOTENESS + (1|ID),
				offset=sem.oa[,"AREA"],
				data=sem.oa,
				#control = contr.optim,
				REML=FALSE),
		SP.STAB <- lmer(SP.STAB ~ FRIC + MHW + REMOTENESS + (1|ID),
				offset=sem.oa[,"AREA"],
				data=sem.oa,
				#control = contr.optim,
				REML=FALSE),
		FRIC <- lmer(FRIC ~ MHW + REMOTENESS + (1|ID),
				offset=sem.oa[,"AREA"],
				data=sem.oa,
				#control = contr.optim,
				REML=FALSE)
)

modoa <- summary(mod.oa)
modoa
coefs(mod.oa, standardize="scale")

#### ---- USE semEff TO ESTIMATE DIRECT, INDIRECT AND TOTAL ---- ####
#### ---- EFFECTS AND COMPARE MPA AND OA MODELS ------------- ####

#### MPA
mod.mpa <-  list(
		STAB <- lmer(STAB ~ ASYNC + SP.STAB + FRIC + MHW + REMOTENESS + (1|ID) +
						offset(AREA),
				data=sem.mpa,
				#control = contr.optim,
				REML=FALSE),
		ASYNC <- lmer(ASYNC ~ FRIC + MHW + REMOTENESS + (1|ID) +
						offset(AREA),
				data=sem.mpa,
				#control = contr.optim,
				REML=FALSE),
		SP.STAB <- lmer(SP.STAB ~ FRIC + MHW + REMOTENESS + (1|ID) +
						offset(AREA),
				data=sem.mpa,
				#control = contr.optim,
				REML=FALSE),
		FRIC <- lmer(FRIC ~ MHW + (1|ID) + # REMOTENESS removed (p>0.6) to improve model fit 
						offset(AREA),
				data=sem.mpa,
				#control = contr.optim,
				REML=FALSE)
)

mpa.boot <- bootEff(mod.mpa, R = 10000, seed = 4444, ran.eff="ID", type="parametric", parallel = "snow")
save(mpa.boot, file="mpa.boot.RData")

#### OA
mod.oa <-  list(
		STAB <- lmer(STAB ~ ASYNC + SP.STAB + FRIC + MHW + REMOTENESS + (1|ID) +
						offset(AREA),
				data=sem.oa,
				#control = contr.optim,
				REML=FALSE),
		ASYNC <- lmer(ASYNC ~ FRIC + MHW + REMOTENESS + (1|ID) +
						offset(AREA),
				data=sem.oa,
				#control = contr.optim,
				REML=FALSE),
		SP.STAB <- lmer(SP.STAB ~ FRIC + MHW + REMOTENESS + (1|ID) +
						offset(AREA),
				data=sem.oa,
				#control = contr.optim,
				REML=FALSE),
		FRIC <- lmer(FRIC ~ MHW + REMOTENESS + (1|ID) +
						offset(AREA),
				data=sem.oa,
				#control = contr.optim,
				REML=FALSE)
)

oa.boot <- bootEff(mod.oa, R = 10000, seed = 4444, ran.eff="ID", type="parametric", parallel = "snow")
save(oa.boot,file="oa.boot.RData")

#### FOREST PLOTS

# function to reshape semEff output
shape_eff_sem <- function(eff.obj, var, plot=F, ...) {
	
	eff.obj <- eff.obj$Summary[[var]]
	
	dir.sel <- charmatch("DIRECT", eff.obj[,1])
	indir.sel <- charmatch("INDIRECT", eff.obj[,1])
	tot.sel <- charmatch("TOTAL", eff.obj[,1])
	
	var.dir <- str_trim(eff.obj[dir.sel:(indir.sel-1),2])
	est.dir <- as.numeric(eff.obj[dir.sel:(indir.sel-1),"Effect"])
	ci.low.dir <- as.numeric(eff.obj[dir.sel:(indir.sel-1),"Lower CI"])
	ci.high.dir <- as.numeric(eff.obj[dir.sel:(indir.sel-1),"Upper CI"])
	
	dir.eff <- data.frame(
			Effect=rep("Direct effects", length(var)),
			Var=var.dir,
			Coef=est.dir,
			CI.low=ci.low.dir,
			CI.upp=ci.high.dir
	)
	
	var.indir <- str_trim(eff.obj[indir.sel:(tot.sel-1),2])
	est.indir <- as.numeric(eff.obj[indir.sel:(tot.sel-1),"Effect"])
	ci.low.indir <- as.numeric(eff.obj[indir.sel:(tot.sel-1),"Lower CI"])
	ci.high.indir <- as.numeric(eff.obj[indir.sel:(tot.sel-1),"Upper CI"])
	
	indir.eff <- data.frame(
			Effect=rep("Indirect effects", length(var)),
			Var=var.indir,
			Coef=est.indir,
			CI.low=ci.low.indir,
			CI.upp=ci.high.indir
	)
	
	df <- rbind(dir.eff, indir.eff) %>% drop_na()
	df
	
}

mpa.eff <- semEff(mpa.boot)
stab.eff <- getTotEff(mpa.eff, "STAB")
stab.eff.boot <- getTotEff(mpa.eff, "STAB", type = "boot")

oa.eff <- semEff(oa.boot)
oa.stab.eff <- getTotEff(oa.eff, "STAB")
oa.stab.eff.boot <- getTotEff(oa.eff, "STAB", type = "boot")

mpa.forest.stab <- shape_eff_sem(mpa.eff, "STAB")
oa.forest.stab<- shape_eff_sem(oa.eff, "STAB")

stab.forest.dat <- data.frame(
				Metric="Alpha stability",
				MPA=rep(c("Protected","Unprotected"), each=nrow(mpa.forest.stab)),
				rbind(mpa.forest.stab, oa.forest.stab)) %>%
		mutate(Var=fct_relevel(Var,rev(c("ASYNC","SP.STAB","FRIC", #"SR",
										"MHW","REMOTENESS"))),
				MPA=fct_relevel(MPA,
						c("Unprotected","Protected")))

mpa.forest.async <- shape_eff_sem(mpa.eff, "ASYNC")
oa.forest.async <- shape_eff_sem(oa.eff, "ASYNC")

async.forest.dat <- data.frame(
				Metric="Asynchrony",
				MPA=c(rep("Protected", nrow(mpa.forest.async)),
						rep("Unprotected", each=nrow(oa.forest.async))),
				rbind(mpa.forest.async, oa.forest.async)) %>%
		mutate(Var=fct_relevel(Var,rev(c("FRIC","MHW", #"SR",
										"REMOTENESS"))),
				MPA=fct_relevel(MPA,
						c("Unprotected","Protected")))

mpa.forest.spstab <- shape_eff_sem(mpa.eff, "SP.STAB")
oa.forest.spstab <- shape_eff_sem(oa.eff, "SP.STAB")

spstab.forest.dat <- data.frame(
				Metric="Species stability",
				MPA=c(rep("Protected", nrow(mpa.forest.async)),
						rep("Unprotected", nrow(oa.forest.async))),
				rbind(mpa.forest.spstab, oa.forest.spstab)) %>%
		mutate(Var=fct_relevel(Var,rev(c("FRIC","MHW", #"SR",
										"REMOTENESS"))),
				MPA=fct_relevel(MPA,
						c("Unprotected","Protected")))

plot.semeff.dat <- rbind(stab.forest.dat, async.forest.dat,
				spstab.forest.dat) %>% 
		mutate(Var=recode(Var, REMOTENESS="REM"))

p.semeff <- ggplot(plot.semeff.dat,
				aes(x = Var, y = Coef, ymin=CI.low, ymax=CI.upp, group=MPA)) +
		geom_pointrange(aes(col=MPA, group=MPA), shape = 18, size=1, linewidth = 1,
				position=position_dodge(width = 0.8)) + 
		geom_hline(yintercept = 0, color = "red", linetype = 2, linewidth = 1) +
		scale_color_fish_d(option = "Cirrhilabrus_solorensis", begin=0.95, end=0.05, direction=1) +
		xlab(" ") + 
		ylab("Effect") + 
		theme_bw() +
		coord_flip() +
		facet_grid(Metric ~ Effect, scales="free_x") +
		theme(
				legend.position="none",
				strip.background=element_rect(colour="grey90",
						fill="grey90"),
				strip.text=element_text(size=12),
				panel.spacing.x=unit(1.5,"line"),
				panel.border=element_rect(colour="grey60", linewidth=0.2),
				panel.background = element_blank(),
				panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(), 
				axis.text.y = element_text(size = 10, colour = "black"),
				axis.text.x.bottom = element_text(size = 10, colour = "black"),
				axis.title.x = element_text(size = 12, colour = "black"))

windows(height=7, width=5)
p.semeff

dev.off()
ggsave(file = "~/Data_ReefFishStability/Figs/Fig3b.pdf",
		dpi = 300, width = 5, height = 7, useDingbats=FALSE)


#### ---- SEMS AND EFFECT SIZES WITH ALTERNATIVE ASYNCHRONY MEASURES ---- ####
# GROSS DETRENDED 
sem.dat <- site_fish_stab %>%
		mutate(
				STAB=standardize(log(1/CV_TOT_ABUND)),
				SP.STAB=standardize(log(1/CV_SP_ABUND)),
				ASYNC=standardize(ASYNC_GROSS_W_DETREG),
				FRIC=standardize(mFRic),
				MHW=standardize(MHW),
				REMOTENESS=standardize(REMOTENESS),
				AREA=standardize(log(SAMPLED_AREA))
		) 		

sem.mpa <- sem.dat %>% filter(MPA=="Protected") 
sem.oa <- sem.dat %>% filter(MPA=="Unprotected") 

#### MPA
mod.mpa <-  list(
		STAB <- lmer(STAB ~ ASYNC + SP.STAB + FRIC + MHW + REMOTENESS + (1|ID) +
						offset(AREA),
				data=sem.mpa,
				#control = contr.optim,
				REML=FALSE),
		ASYNC <- lmer(ASYNC ~ FRIC + MHW + REMOTENESS + (1|ID) +
						offset(AREA),
				data=sem.mpa,
				#control = contr.optim,
				REML=FALSE),
		SP.STAB <- lmer(SP.STAB ~ FRIC + MHW + REMOTENESS + (1|ID) +
						offset(AREA),
				data=sem.mpa,
				#control = contr.optim,
				REML=FALSE),
		FRIC <- lmer(FRIC ~ MHW + (1|ID) + # REMOTENESS removed (p>0.6) to improve model fit 
						offset(AREA),
				data=sem.mpa,
				#control = contr.optim,
				REML=FALSE)
)

mpa.boot.gross.det <- bootEff(mod.mpa, R = 10000, seed = 4444, ran.eff="ID", type="parametric", parallel = "snow")
save(mpa.boot.gross.det, file="mpa.boot.gross.det.RData")

#### OA
mod.oa <-  list(
		STAB <- lmer(STAB ~ ASYNC + SP.STAB + FRIC + MHW + REMOTENESS + (1|ID) +
						offset(AREA),
				data=sem.oa,
				#control = contr.optim,
				REML=FALSE),
		ASYNC <- lmer(ASYNC ~ FRIC + MHW + REMOTENESS + (1|ID) +
						offset(AREA),
				data=sem.oa,
				#control = contr.optim,
				REML=FALSE),
		SP.STAB <- lmer(SP.STAB ~ FRIC + MHW + REMOTENESS + (1|ID) +
						offset(AREA),
				data=sem.oa,
				#control = contr.optim,
				REML=FALSE),
		FRIC <- lmer(FRIC ~ MHW + REMOTENESS + (1|ID) +
						offset(AREA),
				data=sem.oa,
				#control = contr.optim,
				REML=FALSE)
)

oa.boot.gross.det <- bootEff(mod.oa, R = 10000, seed = 4444, ran.eff="ID", type="parametric", parallel = "snow")
save(oa.boot.gross.det, file="oa.boot.gross.det.RData")

mpa.eff <- semEff(mpa.boot.gross.det)
oa.eff <- semEff(oa.boot.gross.det)

mpa.forest.stab <- shape_eff_sem(mpa.eff, "STAB")
oa.forest.stab<- shape_eff_sem(oa.eff, "STAB")

stab.forest.dat <- data.frame(
				Metric="Alpha stability",
				MPA=rep(c("Protected","Unprotected"), each=nrow(mpa.forest.stab)),
				rbind(mpa.forest.stab, oa.forest.stab)) %>%
		mutate(Var=fct_relevel(Var,rev(c("ASYNC","SP.STAB","FRIC",
										"MHW","REMOTENESS"))),
				MPA=fct_relevel(MPA,
						c("Unprotected","Protected")))

mpa.forest.async <- shape_eff_sem(mpa.eff, "ASYNC")
oa.forest.async <- shape_eff_sem(oa.eff, "ASYNC")

async.forest.dat <- data.frame(
				Metric="Asynchrony",
				MPA=c(rep("Protected", nrow(mpa.forest.async)),
						rep("Unprotected", each=nrow(oa.forest.async))),
				rbind(mpa.forest.async, oa.forest.async)) %>%
		mutate(Var=fct_relevel(Var,rev(c("FRIC","MHW",
										"REMOTENESS"))),
				MPA=fct_relevel(MPA,
						c("Unprotected","Protected")))

mpa.forest.spstab <- shape_eff_sem(mpa.eff, "SP.STAB")
oa.forest.spstab <- shape_eff_sem(oa.eff, "SP.STAB")

spstab.forest.dat <- data.frame(
				Metric="Species stability",
				MPA=c(rep("Protected", nrow(mpa.forest.async)),
						rep("Unprotected", each=nrow(oa.forest.async))),
				rbind(mpa.forest.spstab, oa.forest.spstab)) %>%
		mutate(Var=fct_relevel(Var,rev(c("FRIC","MHW",
										"REMOTENESS"))),
				MPA=fct_relevel(MPA,
						c("Unprotected","Protected")))

plot.semeff.dat <- rbind(stab.forest.dat, async.forest.dat,
				spstab.forest.dat) %>% 
		mutate(Var=recode(Var, REMOTENESS="REM"))

p.semeff.gross.det <- ggplot(plot.semeff.dat,
				aes(x = Var, y = Coef, ymin=CI.low, ymax=CI.upp, group=MPA)) +
		geom_pointrange(aes(col=MPA, group=MPA), shape = 18, size = 1,
				position=position_dodge(width = 0.8)) + 
		geom_hline(yintercept = 0, color = "red", linetype = 2, linewidth = 1) +
		scale_color_fish_d(option = "Cirrhilabrus_solorensis", begin=0.95, end=0.05, direction=1) +
		xlab(" ") + 
		ylab("Effect") + 
		theme_bw() +
		coord_flip() +
		facet_grid(Metric ~ Effect, scales="free_x") +
		theme(
				legend.position="none",
				strip.background=element_rect(colour="grey90",
						fill="grey90"),
				strip.text=element_text(size=12),
				panel.spacing.x=unit(1.5,"line"),
				panel.border=element_rect(colour="grey60", linewidth=0.2),
				panel.background = element_blank(),
				panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(), 
				#axis.line = element_line(colour = "black"),
				axis.text.y = element_text(size = 10, colour = "black"),
				axis.text.x.bottom = element_text(size= 10, colour = "black"),
				axis.title.x = element_text(size = 12, colour = "black"))

windows(height=8, width=6)
p.semeff.gross.det

dev.off()
ggsave(file = "~/Data_ReefFishStability/Figs/FigS4a.pdf",
		dpi = 300, width = 5, height = 7, useDingbats=FALSE)

# SQRT(LOREAU_DEMAZANCOURT)
sem.dat <- site_fish_stab %>%
		mutate(
				STAB=standardize(log(1/CV_TOT_ABUND)),
				SP.STAB=standardize(log(1/CV_SP_ABUND)),
				ASYNC=standardize(ASYNC_LOREAU_SQRT),
				FRIC=standardize(mFRic),
				MHW=standardize(MHW),
				REMOTENESS=standardize(REMOTENESS),
				AREA=standardize(log(SAMPLED_AREA))
		) 		

sem.mpa <- sem.dat %>% filter(MPA=="Protected") 
sem.oa <- sem.dat %>% filter(MPA=="Unprotected") 

#### MPA
mod.mpa <-  list(
		STAB <- lmer(STAB ~ ASYNC + SP.STAB + FRIC + MHW + REMOTENESS + (1|ID) +
						offset(AREA),
				data=sem.mpa,
				#control = contr.optim,
				REML=FALSE),
		ASYNC <- lmer(ASYNC ~ FRIC + MHW + REMOTENESS + (1|ID) +
						offset(AREA),
				data=sem.mpa,
				#control = contr.optim,
				REML=FALSE),
		SP.STAB <- lmer(SP.STAB ~ FRIC + MHW + REMOTENESS + (1|ID) +
						offset(AREA),
				data=sem.mpa,
				#control = contr.optim,
				REML=FALSE),
		FRIC <- lmer(FRIC ~ MHW + (1|ID) + # REMOTENESS removed (p>0.6) to improve model fit 
						offset(AREA),
				data=sem.mpa,
				#control = contr.optim,
				REML=FALSE)
)

mpa.boot.loreau.sqrt <- bootEff(mod.mpa, R = 10000, seed = 4444, ran.eff="ID", type="parametric", parallel = "snow")
save(mpa.boot.loreau.sqrt, file="mpa.boot.loreau.sqrt.RData")

#### OA
mod.oa <-  list(
		STAB <- lmer(STAB ~ ASYNC + SP.STAB + FRIC + MHW + REMOTENESS + (1|ID) +
						offset(AREA),
				data=sem.oa,
				#control = contr.optim,
				REML=FALSE),
		ASYNC <- lmer(ASYNC ~ FRIC + MHW + REMOTENESS + (1|ID) +
						offset(AREA),
				data=sem.oa,
				#control = contr.optim,
				REML=FALSE),
		SP.STAB <- lmer(SP.STAB ~ FRIC + MHW + REMOTENESS + (1|ID) +
						offset(AREA),
				data=sem.oa,
				#control = contr.optim,
				REML=FALSE),
		FRIC <- lmer(FRIC ~ MHW + REMOTENESS + (1|ID) +
						offset(AREA),
				data=sem.oa,
				#control = contr.optim,
				REML=FALSE)
)

oa.boot.loreau.sqrt <- bootEff(mod.oa, R = 10000, seed = 4444, ran.eff="ID", type="parametric", parallel = "snow")
save(oa.boot.loreau.sqrt, file="oa.boot.loreau.sqrt.RData")

mpa.eff <- semEff(mpa.boot.loreau.sqrt)
oa.eff <- semEff(oa.boot.loreau.sqrt)

mpa.forest.stab <- shape_eff_sem(mpa.eff, "STAB")
oa.forest.stab<- shape_eff_sem(oa.eff, "STAB")

stab.forest.dat <- data.frame(
				Metric="Alpha stability",
				MPA=rep(c("Protected","Unprotected"), each=nrow(mpa.forest.stab)),
				rbind(mpa.forest.stab, oa.forest.stab)) %>%
		mutate(Var=fct_relevel(Var,rev(c("ASYNC","SP.STAB","FRIC",
										"MHW","REMOTENESS"))),
				MPA=fct_relevel(MPA,
						c("Unprotected","Protected")))

mpa.forest.async <- shape_eff_sem(mpa.eff, "ASYNC")
oa.forest.async <- shape_eff_sem(oa.eff, "ASYNC")

async.forest.dat <- data.frame(
				Metric="Asynchrony",
				MPA=c(rep("Protected", nrow(mpa.forest.async)),
						rep("Unprotected", each=nrow(oa.forest.async))),
				rbind(mpa.forest.async, oa.forest.async)) %>%
		mutate(Var=fct_relevel(Var,rev(c("FRIC","MHW",
										"REMOTENESS"))),
				MPA=fct_relevel(MPA,
						c("Unprotected","Protected")))

mpa.forest.spstab <- shape_eff_sem(mpa.eff, "SP.STAB")
oa.forest.spstab <- shape_eff_sem(oa.eff, "SP.STAB")

spstab.forest.dat <- data.frame(
				Metric="Species stability",
				MPA=c(rep("Protected", nrow(mpa.forest.async)),
						rep("Unprotected", each=nrow(oa.forest.async))),
				rbind(mpa.forest.spstab, oa.forest.spstab)) %>%
		mutate(Var=fct_relevel(Var,rev(c("FRIC","MHW",
										"REMOTENESS"))),
				MPA=fct_relevel(MPA,
						c("Unprotected","Protected")))

plot.semeff.dat <- rbind(stab.forest.dat, async.forest.dat,
				spstab.forest.dat) %>% 
		mutate(Var=recode(Var, REMOTENESS="REM"))

p.semeff.loreau.sqrt <- ggplot(plot.semeff.dat,
				aes(x = Var, y = Coef, ymin=CI.low, ymax=CI.upp, group=MPA)) +
		geom_pointrange(aes(col=MPA, group=MPA), shape = 18, size = 1,
				position=position_dodge(width = 0.8)) + 
		geom_hline(yintercept = 0, color = "red", linetype = 2, linewidth = 1) +
		scale_color_fish_d(option = "Cirrhilabrus_solorensis", begin=0.95, end=0.05, direction=1) +
		xlab(" ") + 
		ylab("Effect") + 
		theme_bw() +
		coord_flip() +
		facet_grid(Metric ~ Effect, scales="free_x") +
		theme(
				legend.position="none",
				strip.background=element_rect(colour="grey90",
						fill="grey90"),
				strip.text=element_text(size=12),
				panel.spacing.x=unit(1.5,"line"),
				panel.border=element_rect(colour="grey60", linewidth=0.2),
				panel.background = element_blank(),
				panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(), 
				#axis.line = element_line(colour = "black"),
				axis.text.y = element_text(size = 10, colour = "black"),
				axis.text.x.bottom = element_text(size= 10, colour = "black"),
				axis.title.x = element_text(size = 12, colour = "black"))

windows(height=8, width=6)
p.semeff.loreau.sqrt

dev.off()
ggsave(file = "~/Data_ReefFishStability/Figs/FigS4b.pdf",
		dpi = 300, width = 5, height = 7, useDingbats=FALSE)






