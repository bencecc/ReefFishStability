# Utility function to plot results from two-way gam models
# examining interactions between stability and level of
# protection (MPAs vs OAs). Use in combinaiton with function
# annotate_stats. The funciton returns prediciotns if plot=F.
###############################################################################

# require(mgcv) # for fitting gams
# require(tidymv) # for predictions

site_plot_gam <- function(gam.fit, x.lab=T, y.lab=T, r2=F, plot=F, ...) {
	
	is.gamm <- FALSE 
	
	if(any(class(gam.fit)%in%"gamm")) {
		gam.fit <- gam.fit$gam
		is.gamm <- T
	}
	
	mf <- gam.fit$model
	cov <- attr(terms(gam.fit), 'term.labels')[2]
	fm <- formula(gam.fit)
	str <- strsplit(as.character(fm), split=c('* MPA'))
	resp <- str[[2]][1]
	
	sm <- summary(gam.fit)
	coef.table <- sm$s.table
	mpa.level.order <- paste("MPA", levels(mf$MPA), sep="")
	model.sources <- unlist(strsplit(as.character(rownames(coef.table)), split=c(':')))
	cov.int <- mpa.level.order[which(mpa.level.order%in%model.sources)]
	
	x_range <- function(x) {
		xr <- range(x)
		seq(xr[1], xr[2], length.out=10)
	}
	
	# To limit x range for predictions to observed values seprarately for protected and unprotected zones; 
	# currently set for full range
	x.range.mpa <- (mf %>% filter(MPA=="Protected") %>% summarise(range=range(.data[[cov]])))$range
	x.range.oa <- (mf %>% filter(MPA=="Unprotected") %>% summarise(range=range(.data[[cov]])))$range
	x.mpa <- seq(x.range.mpa[1], x.range.mpa[2], length.out=10)
	x.oa <- seq(x.range.oa[1], x.range.oa[2], length.out=10)
	
	# Prediction from model
	if(is.gamm) {
	
		mpa.model.eff <- tidygam::predict_gam(gam.fit, values=list(MHW=x.mpa)) %>%
				rename(Cov=!!sym(cov), Pred=!!sym(resp), conf.low=lower_ci, conf.high=upper_ci) %>%
				filter(MPA=="Protected")
		
		oa.model.eff <- tidygam::predict_gam(gam.fit, values=list(MHW=x.oa)) %>%
				rename(Cov=!!sym(cov), Pred=!!sym(resp), conf.low=lower_ci, conf.high=upper_ci) %>%
				filter(MPA=="Unprotected")
	
	} else {
		
		mpa.model.eff <- tidygam::predict_gam(gam.fit, exclude_terms=c("s(ID)"), values=list(MHW=x.mpa)) %>%
				rename(Cov=!!sym(cov), Pred=!!sym(resp), conf.low=lower_ci, conf.high=upper_ci) %>%
				filter(MPA=="Protected")
		
		oa.model.eff <- tidygam::predict_gam(gam.fit, exclude_terms=c("s(ID)"), values=list(MHW=x.oa)) %>%
				rename(Cov=!!sym(cov), Pred=!!sym(resp), conf.low=lower_ci, conf.high=upper_ci) %>%
				filter(MPA=="Unprotected")
		
	}

	model.eff <- rbind(
			cbind(mpa.model.eff),
			cbind(oa.model.eff)
	)
	
	# Extract edf, p-values and R^2 to add to plots
	row.mpa <- which(grepl("MPAProtected", rownames(coef.table))==T)
	row.oa <- which(grepl("MPAUnprotected", rownames(coef.table))==T)
	
	star <- function(pval) {
		if(pval<0.001) {
			star <- "***"
		} else if(pval<0.01){
			star <- "**"
		} else if(pval<0.05){
			star <- "*"
		} else{star <- " "}
		
	}
	
	if(r2) {
		r2 <- round(sm$r.sq, 2)
		stat.r2 <- sprintf("italic(R^2) == %.3f", r2)
	} else { stat.r2 <- ""}
	
	mpa.stats <- list(paste("edf = ", round(coef.table[row.mpa,1],2),
					star(coef.table[row.mpa,4]), sep=""), stat.r2)
	oa.stats <- list(paste("edf = ", round(coef.table[row.oa,1], 2),
					star(coef.table[row.oa,4]), sep=""))
	
	# set y maximum limit and position for text statistics
	y.max <- max(
			max(model.eff[,"Pred"]+model.eff[,"conf.high"]),
			max(mf[,resp]))
	ylim.high <- y.max+ 0.4*y.max
	
	xrange <- range(mf[,cov])
	yrange <- range(mf[,resp])
	y.shift <- yrange[2] + 0.2*diff(yrange)
	
	x1 <- xrange[1] + 0.02*diff(xrange)
	x2 <- xrange[1] + 0.55*diff(xrange)
	y1 <- Inf #y.shift - 0.01 * diff(yrange)
	y2 <-  Inf #y.shift - 0.12*diff(yrange)
	
	# Graphical parameters
	mycol <- fish(5, option="Cirrhilabrus_solorensis",
			begin=0.95, end=0.05, direction=1)
	
	if(y.lab) {
		
		y.lab <- switch(
				resp,
				"STAB"=(expression(paste(alpha, " stability", sep=" "))),
				"ASYNC"="Species asynchrony",
				"SP.STAB"="Species stability",
				"SR"="Species richness",
				"FR"="Functional richness",
				"abund"="Abundance",
				"ABUND"="Abundance",
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
	
	p <- ggplot()+
			geom_point(data=mf,
					aes_string(x=cov, y=resp, shape="MPA", col="MPA"),
					alpha=0.1, size=0.8)+
			geom_path(data=model.eff, aes(x=Cov, y=Pred, color=MPA), linewidth=0.3) +
			annotate("text", x = x1, y = y1, hjust=0, vjust=1.8, label = mpa.stats[[1]], parse=F, size=2.5, col=mycol[5]) +
			annotate("text", x = x2, y = y1, hjust=0, vjust=1.8, label = oa.stats[[1]], parse=F, size=2.5, col=mycol[1]) +
			annotate("text", x = x1, y = y2, hjust=0, vjust=2.5, label = mpa.stats[[2]], parse=T, size=2.5, col="magenta") +
			
			scale_size_manual(values=c(1,1))+	
			geom_ribbon(data=model.eff, aes(x=Cov, y=Pred, ymin=conf.low,ymax=conf.high, fill=MPA, col=MPA),
					alpha=0.25, linetype=2, linewidth=0.2)+
			scale_fill_fish_d(option = "Cirrhilabrus_solorensis", begin=0.95, end=0.05, direction=1) +
			scale_color_fish_d(option = "Cirrhilabrus_solorensis", begin=0.95, end=0.05, direction=1) +
			#theme_sjplot(base_family="arial") +
			set_theme(base=theme_bw(),
					legend.pos="none") +
			ylim(NA,ylim.high) +
			theme(
					legend.position = "none",
					panel.background = element_blank(),
					panel.grid.major = element_blank(), 
					panel.grid.minor = element_blank(),
					axis.text.y = element_text(size = 8, colour = "black"),
					axis.text.x = element_text(size = 8, colour = "black"),
					axis.title.x = element_text(size = 10, colour = "black"),
					axis.title.y = element_text(size = 10, colour = "black", vjust=-0.5),
					plot.margin = margin(0.1,0.1,0,0, "cm")
			) +
			labs(y=y.lab, x=x.lab)
	
	if(plot) {
		#windows(height=3, width=3)
		plot(p)
	}
	
	else {
		
		out <- model.eff %>%
				mutate(edf = case_when(
								
								MPA=="Protected" ~ mpa.stats[[1]],
								MPA=="Unprotected" ~ oa.stats[[1]]),
						
						r2 = case_when(
								
								MPA=="Protected" ~ mpa.stats[[2]],
								MPA=="Unprotected" ~ ""
						
						),
						
						x1=x1,
						x2=x2,
						y1=y1,
						y2=y2,
						ylim.high=ylim.high
				)
		return(out)
		
	}	
	
}
