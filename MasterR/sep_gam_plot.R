# Utility function to plot results from two-way gam models
# examining interactions between stability and level of
# protection (MPAs vs OAs). Use in combinaiton with function
# annotate_stats. The funciton returns prediciotns if plot=F.
###############################################################################

# require(mgcv) # for fitting gams
# require(tidymv) # for predictions

sep_gam_plot <- function(df,
		resp=c("ABUND","MEAN.ABUND","SD.ABUND"),
		cov="MHW", x.range, x.lab=T, y.lab=T, r2=T, plot=F, ...) {
	
	df <- df %>% 
			mutate(
					MPA=forcats::fct_relevel(MPA,
							c("Unprotected","Protected"))
			) 		
	
	# sepratae MPA and OA data
	mpa.dat <- df %>% dplyr::filter(MPA=="Protected")
	oa.dat <- df %>% dplyr::filter(MPA=="Unprotected")
	
	# fit models
	fm <- formula(paste(resp, "~ s(", cov, ", bs='cs') + s(ID, bs='re') + offset(AREA)"))
	mpa.mod <- gam(fm, data=mpa.dat, na.action="na.omit")
	oa.mod <- gam(fm, data=oa.dat, na.action="na.omit")
	
	# extract attributes and model summaries
	mf <- rbind(cbind(MPA="Protected", mpa.mod$model),
			cbind(MPA="Unprotected", oa.mod$model))	
	cov <- attr(terms(mpa.mod), 'term.labels')[1]
	str <- strsplit(as.character(fm), split=c('~'))
	resp <- str[[2]][1]
	mpa.sm <- summary(mpa.mod)
	mpa.coef.table <- mpa.sm$s.table
	oa.sm <- summary(oa.mod)
	oa.coef.table <- oa.sm$s.table
	
	# range of data for plotting
	x_range <- function(x) {
		xr <- range(x)
		seq(xr[1], xr[2], length.out=10)
	}
	
	# To limit x range for predictions to observed values seprarately for protected and unprotected zones; 
	# currently set for full range
	x.range.mpa <- (mpa.mod$model %>% summarise(range=range(.data[[cov]])))$range
	x.range.oa <- (oa.mod$model %>% summarise(range=range(.data[[cov]])))$range
	x.mpa <- seq(x.range.mpa[1], x.range.mpa[2], length.out=10)
	x.oa <- seq(x.range.oa[1], x.range.oa[2], length.out=10)
	
	# Prediction from model
	mpa.model.eff <- tidygam::predict_gam(mpa.mod, exclude_terms=c("s(ID)"), values=list(cov=x.mpa)) %>%
			rename(Cov=!!sym(cov), Pred=!!sym(resp), conf.low=lower_ci, conf.high=upper_ci)
	
	oa.model.eff <- tidygam::predict_gam(oa.mod, exclude_terms=c("s(ID)"), values=list(cov=x.oa)) %>%
			rename(Cov=!!sym(cov), Pred=!!sym(resp), conf.low=lower_ci, conf.high=upper_ci)
	
	model.eff <- rbind(
			cbind(MPA="Protected", mpa.model.eff),
			cbind(MPA="Unprotected",oa.model.eff)
	)
	
	# Extract edf, p-values and R^2 to add to plots
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
		mpa.r2 <- round(mpa.sm$r.sq, 2)
		mpa.stat.r2 <- sprintf("italic(R^2) == %.3f", mpa.r2)
		oa.r2 <- round(oa.sm$r.sq, 2)
		oa.stat.r2 <- sprintf("italic(R^2) == %.3f", oa.r2)
	
	} else { 
		mpa.stat.r2 <- ""
		oa.stat.r2 <- ""
	}
	
	mpa.stats <- list(paste("edf = ", round(mpa.coef.table[1,1],2),
					star(mpa.coef.table[1,4]), sep=""), mpa.stat.r2)
	oa.stats <- list(paste("edf = ", round(oa.coef.table[1,1], 2),
					star(oa.coef.table[1,4]), sep=""), oa.stat.r2)
	
	# set y maximum limit and position for text statistics
	y.max <- max(
			max(model.eff[,"Pred"]+model.eff[,"conf.high"]),
			max(c(mpa.mod$model[,resp], oa.mod$model[,resp])))
	ylim.high <- y.max+ 0.4*y.max
	
	xrange <- range(c(mpa.mod$model[,cov], oa.mod$model[,cov]))
	yrange <- range(c(mpa.mod$model[,resp], oa.mod$model[,resp]))
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
				"ABUND"="Abundance",
				"MEAN.ABUND"="Abundance",
				"SD.ABUND"="Standard deviation",
		)
		
	} else {y.lab <- NULL}
	
	if(x.lab) {
		
		x.lab <- switch(
				cov,
				"MHW"=(expression(paste("MHW intensity (", degree, "C)"))),
				"REMOTENESS"="Remoteness",
				"MEAN.ABUND"="Abundance",
				"ABUND"="Abundance",
				"SD.ABUND"="Standard deviation",
		)
	} else {x.lab <- NULL}
	
	p <- ggplot()+
			geom_point(data=mf,
					aes_string(x=cov, y=resp, shape="MPA", col="MPA"),
					alpha=0.1, size=0.8)+
			geom_path(data=model.eff, aes(x=Cov, y=Pred, color=MPA), linewidth=0.3) +
			annotate("text", x = x1, y = y1, hjust=0, vjust=1.8, label = mpa.stats[[1]], parse=F, size=2.5, col=mycol[5]) +
			annotate("text", x = x2, y = y1, hjust=0, vjust=1.8, label = oa.stats[[1]], parse=F, size=2.5, col=mycol[1]) +
			annotate("text", x = x1, y = y2, hjust=0, vjust=2.5, label = mpa.stats[[2]], parse=T, size=2.5, col=mycol[5]) +
			annotate("text", x = x2, y = y2, hjust=0, vjust=2.5, label = oa.stats[[2]], parse=T, size=2.5, col=mycol[1]) +
			scale_size_manual(values=c(1,1))+	
			geom_ribbon(data=model.eff, aes(x=Cov, y=Pred, ymin=conf.low,ymax=conf.high, fill=MPA, col=MPA),
					alpha=0.25, linetype=2, linewidth=0.2)+
			scale_fill_fish_d(option = "Cirrhilabrus_solorensis", begin=0.95, end=0.05, direction=-1) +
			scale_color_fish_d(option = "Cirrhilabrus_solorensis", begin=0.95, end=0.05, direction=-1) +
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
	} else {
		
		return(list(p, mpa.sm, oa.sm))	
	}	
	
}
