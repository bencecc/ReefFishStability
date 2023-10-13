# Utility function to plot results from one-way lmer models
##############################################################

onevar_fit_plot <- function(model, x.lab=T, y.lab=T,
		symb.shape=21, r2=FALSE, plot=F, ...) {
	
	mf <- model@frame
	fm <- formula(model)
	str <- strsplit(as.character(fm), split=c('* ~'))
	resp <- str[[2]][1]
	cov <- attr(terms(model), 'term.labels')[1]
	
	
	#### ---- PREDICITONS FOR FITTING ---- ####
	# To fit lines only within data range
	
	x_range <- function(x) {
		xr <- range(x)
		seq(xr[1], xr[2], length.out=10)
	}
	
# To limit x range for predictions to observed values seprarately for protected and unprotected zones; 
# currently set for full range
#	x.range.mpa <- (mf %>% filter(MPA=="Protected") %>% summarise(range=range(.data[[cov]])))$range
#	x.range.oa <- (mf %>% filter(MPA=="Unprotected") %>% summarise(range=range(.data[[cov]])))$range
#	x.mpa <- seq(x.range.mpa[1], x.range.mpa[2], length.out=10)
#	x.oa <- seq(x.range.oa[1], x.range.oa[2], length.out=10)
	
	model.eff <- ggpredict(model,
					terms=paste(cov, "[x_range]", sep=" "),
					ci.lvl=0.95) %>%
			rename(Cov=x, Pred=predicted)
	
	minmax.data <- mf %>% group_by(ID) %>%
			summarise(
					cov.min=min(.data[[cov]], na.rm=T),
					cov.max=max(.data[[cov]], na.rm=T), .groups = "drop")
	
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
				"AREA"="Area",
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
	
	mod.stats <- annotate_stats(mod=model, cov=cov, r2=r2)
		
	xrange <- range(mf[,cov])
	yrange <- range(mf[,resp])
	y.shift <- yrange[2] + 0.1*diff(yrange)
	
	x1 <- xrange[1] + 0.02*diff(xrange)
	x2 <- xrange[1] + 0.55*diff(xrange)
	y1 <- y.shift - 0.01 * diff(yrange)
	y2 <-  y.shift - 0.12*diff(yrange)
	
	p <- ggplot(data=mf, aes(x=.data[[cov]], y=.data[[resp]]))+
			geom_point(shape=symb.shape, col="#46A9BE", alpha=0.1, size=0.8)+
			geom_line(data=model.eff, aes(x=Cov, y=Pred), color="#46A9BE", linewidth=1) +
			annotate("text", x = x1, y = y1, hjust=0, vjust=1, label = mod.stats[[1]], parse=T, size=2, col="#46A9BE") +
			annotate("text", x = x2, y = y1, hjust=0, vjust=0.7, label = mod.stats[[2]], parse=T, size=2, col="#46A9BE") +
			scale_size_manual(values=c(1,1))+	
			geom_ribbon(data=model.eff, aes(x=Cov, y=Pred, ymin=conf.low,ymax=conf.high),
					fill="#46A9BE", alpha=0.5)+
			set_theme(base=theme_bw(),
					legend.pos="none")+
			theme(
					panel.background = element_blank(),
					panel.grid.major = element_blank(), 
					panel.grid.minor = element_blank(),
					axis.text.y = element_text(size = 8, colour = "black"),
					axis.text.x = element_text(size = 8, colour = "black"),
					axis.title.x = element_text(size = 10, colour = "black"),
					axis.title.y = element_text(size = 10, colour = "black", vjust=-0.5)
			) +
			labs(y=y.lab, x=x.lab)
	
	
	if(plot) {
		windows(height=3, width=3)
		plot(p)
	}
	
	else {
		
		mod.tab <- tab_model(model,  show.stat=T, show.ci=F, show.se=T)	
		return(list(p, mod.tab))	
	}
	
}


