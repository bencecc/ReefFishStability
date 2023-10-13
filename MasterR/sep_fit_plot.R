# Utility function to plot results from lmer models
# examining bivariate relations of stability separately
# between levels of protection (MPAs vs OAs). Use in
# combinaiton with function annotate_stats. 
##############################################################

sep_fit_plot <- function(df,
		resp=c("STAB","ASYNC","SP.STAB","SR","FRIC"),
		cov=NULL, x.range, x.lab=T, y.lab=T, y.lim=NULL, r2=FALSE, plot=F) {
	
	df <- df %>% 
			mutate(
					MPA=forcats::fct_relevel(MPA,
							c("Unprotected","Protected"))
			) 		
	
	mpa.dat <- df %>% dplyr::filter(MPA=="Protected")
	oa.dat <- df %>% dplyr::filter(MPA=="Unprotected")

	if(resp=="STAB") {
		mpa.mod <- lmer(STAB ~ ASYNC + SP.STAB + FRIC + MHW + REMOTENESS + (1|ID) +
						offset(AREA), data=mpa.dat, REML=FALSE)
		oa.mod <- lmer(STAB ~ ASYNC + SP.STAB + FRIC + MHW + REMOTENESS + (1|ID) +
						offset(AREA), data=oa.dat, REML=FALSE)
	}
	
	if(resp=="ASYNC"){
		mpa.mod <- lmer(ASYNC ~ FRIC + MHW + REMOTENESS + (1|ID) +
						offset(AREA), data=mpa.dat, REML=FALSE)
		oa.mod <- lmer(ASYNC ~ FRIC + MHW + REMOTENESS + (1|ID) +
						offset(AREA), data=oa.dat, REML=FALSE)		
	}
	
	if(resp=="SP.STAB") {
		mpa.mod <- lmer(SP.STAB ~ FRIC + MHW + REMOTENESS + (1|ID) +
						offset(AREA), data=mpa.dat, REML=FALSE)
		oa.mod <- lmer(SP.STAB ~ FRIC + MHW + REMOTENESS + (1|ID) +
						offset(AREA), data=oa.dat, REML=FALSE)
	}
	
#	if(resp=="SR") {
#		mpa.mod <- lmer(SR ~ MHW + REMOTENESS + (1|ID) +
#						offset(AREA), data=mpa.dat, REML=FALSE)
#		oa.mod <- lmer(SR ~ MHW + REMOTENESS + (1|ID) +
#						offset(AREA), data=oa.dat, REML=FALSE)
#	}
	
	if(resp=="FRIC") {
		mpa.mod <- lmer(FRIC ~ MHW + REMOTENESS + (1|ID) +
						offset(AREA), data=mpa.dat, REML=FALSE)
		oa.mod <- lmer(FRIC ~ MHW + REMOTENESS + (1|ID) +
						offset(AREA), data=oa.dat, REML=FALSE)
	}
	
	if(resp=="MEAN.ABUND") {
		mpa.mod <- lmer(MEAN.ABUND ~ MHW + (1|ID) + offset(AREA),
				data=mpa.dat, REML=FALSE)
		oa.mod <- lmer(MEAN.ABUND ~ MHW + (1|ID) + offset(AREA),
				data=oa.dat, REML=FALSE)
	}
	
	if(resp=="SD.ABUND") {
		mpa.mod <- lmer(SD.ABUND ~ MHW + (1|ID) + offset(AREA),
				data=mpa.dat, REML=FALSE)
		oa.mod <- lmer(SD.ABUND ~ MHW + (1|ID) + offset(AREA),
				data=oa.dat, REML=FALSE)
		
	}
	
	model.eff.mpa <- ggpredict(mpa.mod,
					terms=paste(cov, "[x.range]", sep=" "),
					ci.lvl=0.95) %>%
			rename(Cov=x, Pred=predicted)
	
	model.eff.oa <- ggpredict(oa.mod,
					terms=paste(cov, "[x.range]", sep=" "),
					ci.lvl=0.95) %>%
			rename(Cov=x, Pred=predicted)
	
	mycol <- fish(5, option="Cirrhilabrus_solorensis",
			begin=0.95, end=0.05, direction=1)
	
	if(y.lab) {
		
		y.lab <- switch(
				resp,
				"STAB"=(expression(paste(alpha, " stability", sep=" "))),
				"ASYNC"="Species asynchrony",
				"SP.STAB"="Species stability",
				"SR"="Species richness",
				"FRIC"="Functional richness",
				"abund"="Abundance",
				"MEAN.ABUND"="Abundance",
				"SD.ABUND"="Standard deviation",
				"STAB.FT"="Stability",
				"cti"="Upper thermal niche (°C)"
		)
		
	} else {y.lab <- NULL}
	
	if(x.lab) {
		
		x.lab <- switch(
				cov,
				"MHW"=(expression(paste("MHW intensity (", degree, "C)"))),
				"REMOTENESS"="Remoteness",
				"STAB"="Stability",
				"ASYNC"="Asynchrony",
				"SP.STAB"="Species stability",
				"SR"="Species richness",
				"FRIC"="Functional richness",
				"abund"="Abundance",
				"MEAN.ABUND"="Abundance",
				"SD.ABUND"="Standard deviation",
				"STAB.FT"="Stability",
				"cti"="Upper thermal niche (°C)"
		)
	} else {x.lab <- NULL}
	
	mpa.stats <- annotate_stats(mod=mpa.mod, cov=cov, r2=r2)
	oa.stats <- annotate_stats(mod=oa.mod, cov=cov, r2=r2)
	
	xrange <- range(df[,cov])
	yrange <- range(df[,resp])
	y.shift <- yrange[2] + 0.2*diff(yrange)
	
	x1 <- xrange[1] + 0.02*diff(xrange)
	x2 <- xrange[1] + 0.55*diff(xrange)
	y1 <- y.shift - 0.01 * diff(yrange)
	y2 <-  y.shift - 0.11*diff(yrange)
	
	p <- ggplot()+
			geom_point(data=df,
					aes(x=.data[[cov]], y=.data[[resp]], shape=.data[["MPA"]], col=.data[["MPA"]]),
					alpha=0.2, size=1.2)+
			geom_line(data=model.eff.mpa, aes(x=Cov, y=Pred), color=mycol[5], linewidth=1) +
			geom_line(data=model.eff.oa, aes(x=Cov, y=Pred), color=mycol[1], linewidth=1) +
			annotate("text", x = x1, y = y1, hjust=0, vjust=1, label = mpa.stats[[1]], parse=T, size=2, col=mycol[5]) +
			annotate("text", x = x1, y = y2, hjust=0, vjust=0.7, label = mpa.stats[[2]], parse=T, size=2, col=mycol[5]) +
			annotate("text", x = x2, y = y1, hjust=0, vjust=1, label = oa.stats[[1]], parse=T, size=2, col=mycol[1]) +
			annotate("text", x = x2, y = y2, hjust=0, vjust=0.7, label = oa.stats[[2]], parse=T, size=2, col=mycol[1]) +
			scale_size_manual(values=c(1,1))+	
#			geom_ribbon(data=model.eff.oa, aes(x=Cov, y=Pred, ymin=conf.low,ymax=conf.high),
#					fill=mycol[1], alpha=0.15)+
#			geom_ribbon(data=model.eff.mpa, aes(x=Cov, y=Pred, ymin=conf.low,ymax=conf.high),
#					fill=mycol[5], alpha=0.15)+
			scale_fill_fish_d(option = "Cirrhilabrus_solorensis", begin=0.95, end=0.05, direction=1) +
			scale_color_fish_d(option = "Cirrhilabrus_solorensis", begin=0.95, end=0.05, direction=1) +
			set_theme(base=theme_bw(),
					legend.pos="none")+
			theme(
					#panel.border = element_blank(),
					panel.background = element_blank(),
					panel.grid.major = element_blank(), 
					panel.grid.minor = element_blank(),
					axis.text.y = element_text(size = 8, colour = "black"),
					axis.text.x = element_text(size = 8, colour = "black"),
					axis.title.x = element_text(size = 10, colour = "black"),
					axis.title.y = element_text(size = 10, colour = "black", vjust=-0.5)
					#plot.margin = margin(0.15,0.1,0.1,0.1, "cm")
					#axis.line = element_line(colour = "black")
			) +
			labs(y=y.lab, x=x.lab)
	
	if(!is.null(y.lim)) {
		
		p <- p + coord_cartesian(ylim=y.lim)
	}
	
	if(plot) {
		windows(height=3, width=3)
		plot(p)
	}
	
	else {
		
		mpa.tab <- tab_model(mpa.mod,  show.stat=T, show.ci=F, show.se=T)	
		oa.tab <- tab_model(oa.mod,  show.stat=T, show.ci=F, show.se=T)
		return(list(p, mpa.tab, oa.tab))	
	}
	
}


