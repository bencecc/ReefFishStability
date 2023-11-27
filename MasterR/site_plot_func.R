# Utility function to plot results from two-way lmer models
# examining interactions between stability and level of
# protection (MPAs vs OAs). Use in combinaiton with function
# annotate_stats. pred.type argument allows to specify plots
# of mean prototypical trajectoris (fixed effect: "fe") or
# individual trajectoris by site (random effect: "re") 
##############################################################

site_plot_func <- function(model, x.var=NULL, pred.type="fe",
		off.val=NULL, x.lab=T, y.lab=T, r2=F, y.lim=NULL, plot=F, ...) {
	
	# derive response variable and covariate from the model
	mf <- model@frame
	fm <- formula(model)
	str <- strsplit(as.character(fm), split=c('* MPA'))
	resp <- str[[2]][1]
	
	#reff <- ranef(model)
	
	if(is.null(x.var)) {
		cov <- attr(terms(model), 'term.labels')[1]
	}
	
	else { 
		cov <- x.var
	}
	
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
	
# Using ggpredict - ok both with/without offset; fitted lines beyond data range
	if(!is.null(off.val)) {
		
		model.eff <- ggpredict(model,
						terms=c(paste(cov, "[x_range]", sep=" "), "MPA"),
						condition=c(AREA=off.val),
						ci.lvl=0.95) %>%
				rename(Cov=x, MPA=group, Pred=predicted)
	}
	
	else {
		
		model.eff <- ggpredict(model,
						terms=c(paste(cov, "[x_range]", sep=" "), "MPA"),
						ci.lvl=0.95) %>%
				rename(Cov=x, MPA=group, Pred=predicted)
	}
	
	model.eff.mpa <- model.eff %>%
			filter(MPA=="Protected")
	
	model.eff.oa <- model.eff %>%
			filter(MPA=="Unprotected")
	
	if(pred.type=="re") {
		
		ran.eff <- ranef(model)$ID
		
		# only random intercept
		if(ncol(ran.eff)==1) {
			# Refit ID level-1 regression using model slope
			fix.eff.tmp <- coef(model)$ID %>%
					tibble::rownames_to_column("ID")				
			cov.terms <- grep(cov, colnames(fix.eff.tmp))
			int.term <- colnames(fix.eff.tmp)[cov.terms[2]]
			
			fix.eff <- fix.eff.tmp %>%
					select(ID, all_of(cov), !!sym(int.term)) %>%
					rename(Slope=!!sym(cov),
							MPAeff=!!sym(int.term))
			pred.dat <- mf %>% group_by(MPA,ID) %>%
					left_join(fix.eff, by="ID") %>%
					mutate(Slope = case_when(MPA=="Protected" ~ Slope + MPAeff,
									TRUE ~ Slope)) %>%			
					mutate(Int=mean(!!sym(resp)-!!sym(cov)*Slope)) %>%
					summarise(Int=mean(Int),
							Slope=mean(Slope),
							cov.min=min(.data[[cov]]),
							cov.max=max(.data[[cov]]), .groups = "drop",
							pred1=Int+Slope*cov.min,
							pred2=Int+Slope*cov.max)
			
		}
		
		# random intercept and slope
		else {
			
			if(any(colnames(mf)=="offset(AREA)")) {
				
				mf <- mf %>% 
						rename(AREA="offset(AREA)")
			}
			
			minmax.data <- mf %>% group_by(MPA,ID) %>%
					summarise(AREA=mean(AREA),
							cov.min=min(.data[[cov]], na.rm=T),
							cov.max=max(.data[[cov]], na.rm=T), .groups = "drop") %>%
					select(-"MPA")
			
			pred.dat <- mf %>% nest(data = -c(ID, MPA)) %>% 
					mutate(
							fit = map(data, ~ lm(!!sym(resp) ~ !!sym(cov) + offset(AREA), data = .x)),
							tidied = map(fit, tidy)
					#glanced = map(fit, glance),
					#augmented = map(fit, augment)
					
					) %>%
					#mutate(new.dat=data.frame(SR=c(.$data$cov.min,.$data$cov.max)))
					unnest(tidied) %>%
					select(c("MPA","ID","term","estimate")) %>%
					pivot_wider(names_from=term, values_from=estimate) %>%
					rename(Int="(Intercept)", Slope=!!sym(cov)) %>%
					#mutate(term=as.factor(term)) %>%
					#group_by(MPA,ID) %>%
					left_join(minmax.data, by=c("ID")) %>%				
					mutate(pred1=Int+AREA+Slope*cov.min,
							pred2=Int+AREA+Slope*cov.max)
			
		}
		
		model.eff.mpa.re <- pred.dat %>%
				filter(MPA=="Protected")
		model.eff.oa.re <- pred.dat %>%
				filter(MPA=="Unprotected")
		
	}
	
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
	
	# identify the term indicating the interaction, whether it is MPAProtected or MPAUnprotected to assign
	# model statiscs correctly
	mpa.level.order <- paste("MPA", levels(mf$MPA), sep="")
	model.sources <- unlist(strsplit(as.character(rownames(summary(model)$coefficients)), split=c(':')))
	cov.int <- mpa.level.order[which(mpa.level.order%in%model.sources)]
	
	if(cov.int=="MPAProtected") {
		
		mpa.stats <- annotate_stats(mod=model, cov=paste(cov, ":", cov.int, sep=""), r2=r2)
		oa.stats <- annotate_stats(mod=model, cov=cov, r2=F)
		
	} else if(cov.int=="MPAUnprotected") {
		
		mpa.stats <- annotate_stats(mod=model, cov=cov, r2=r2)
		oa.stats <- annotate_stats(mod=model, cov=paste(cov, ":", cov.int, sep=""), r2=F)
		
	}
	
	xrange <- range(mf[,cov])
	yrange <- range(mf[,resp])
	y.shift <- yrange[2] + 0.2*diff(yrange)
	
	x1 <- xrange[1] + 0.02*diff(xrange)
	x2 <- xrange[1] + 0.55*diff(xrange)
	y1 <- y.shift - 0.01 * diff(yrange)
	y2 <-  y.shift - 0.12*diff(yrange)
	
	p <- ggplot()+
			geom_point(data=mf,
					aes(x=.data[[cov]], y=.data[[resp]], shape=.data[["MPA"]], col=.data[["MPA"]]),
					alpha=0.2, size=1.2)+
			geom_line(data=model.eff.oa, aes(x=Cov, y=Pred), color=mycol[1], linewidth=1) +
			geom_line(data=model.eff.mpa, aes(x=Cov, y=Pred), color=mycol[5], linewidth=1) +
			annotate("text", x = x1, y = y1, hjust=0, vjust=1, label = mpa.stats[[1]], parse=T, size=2, col=mycol[5]) +
			annotate("text", x = x2, y = y1, hjust=0, vjust=1, label = oa.stats[[1]], parse=T, size=2, col=mycol[1]) +
			annotate("text", x = x1, y = y2, hjust=0, vjust=0.7, label = mpa.stats[[2]], parse=T, size=2, col="magenta") +	
			scale_size_manual(values=c(1,1))+	
#			scale_linetype_manual(values=c("Protected"="solid",
#							"Unprotected"="dashed"))+
#			geom_ribbon(data=model.eff.oa, aes(x=Cov, y=Pred, ymin=conf.low,ymax=conf.high),
#					fill=mycol[1], alpha=0.15)+
#			geom_ribbon(data=model.eff.mpa, aes(x=Cov, y=Pred, ymin=conf.low,ymax=conf.high),
#					fill=mycol[5], alpha=0.15)+
			scale_fill_fish_d(option = "Cirrhilabrus_solorensis", begin=0.95, end=0.05, direction=1) +
			scale_color_fish_d(option = "Cirrhilabrus_solorensis", begin=0.95, end=0.05, direction=1) +
			theme_sjplot(base_family="arial") +
			set_theme(base=theme_bw(),
					legend.pos="none") +
			theme(
					panel.background = element_blank(),
					panel.grid.major = element_blank(), 
					panel.grid.minor = element_blank(),
					axis.text.y = element_text(size = 8, colour = "black"),
					axis.text.x = element_text(size = 8, colour = "black"),
					axis.title.x = element_text(size = 10, colour = "black"),
					axis.title.y = element_text(size = 10, colour = "black", vjust=-0.5),
					plot.margin = margin(0.15,0.1,0.2,0.1, "cm")
					#axis.line = element_line(colour = "black")
			) +
			labs(y=y.lab, x=x.lab)
	
	if(!is.null(y.lim)) {
		
		p <- p + coord_cartesian(ylim=y.lim)
	}
	
	if(pred.type=="re") {
		
		p <- p + 
				geom_segment(data=model.eff.oa.re, aes_string(x="cov.min", xend="cov.max", y="pred1", yend="pred2", group="ID"),
						color=mycol[1], linetype=1, size=0.2) +
				geom_segment(data=model.eff.mpa.re, aes_string(x="cov.min", xend="cov.max", y="pred1", yend="pred2", group="ID"),
						color=mycol[5], linetype=1, size=0.2)
		
	}
	
	if(plot) {
		#windows(height=3, width=3)
		plot(p)
	}
	
	else return(p)
	
}


