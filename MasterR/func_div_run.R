# Compute functional diversity indices using individual site data 
###############################################################################
#require(ape)
#require(FD)

func_div_run <- function(data, fish.traits.data, ...) {
	
	dat <- data %>% ungroup() %>% dplyr::select(YEAR, SPECIES, abund) %>%
			filter(abund>0) %>% arrange(YEAR, SPECIES)
	
	if(length(unique(dat$SPECIES))>2) {
		
		traits <- fish.traits.data %>% filter(SPECIES%in%dat$SPECIES) %>%
				dplyr::select(-c("SPECIES"))
		gw.diss <- gowdis(traits)
		
		# remove species that have more than three missing traits
		if(any(is.na(gw.diss))) {
			
			na.find <- which(apply(traits, 1, function(x)
										length(which(is.na(x))))
							>3)
			
			bad.spp <- fish.traits.data %>% filter(SPECIES%in%dat$SPECIES) %>%
					slice(na.find) %>% dplyr::select(SPECIES)
			dat <- dat %>% filter(!SPECIES%in%bad.spp$SPECIES)
			traits <- fish.traits.data %>% filter(SPECIES%in%dat$SPECIES) %>%
					dplyr::select(-SPECIES)
			gw.diss <- gowdis(traits)
		}
		
		pco.res <- try(pcoa(gw.diss, correction="cailliez"), silent=TRUE)
		if(inherits(pco.res, "try-error")) {
			pco.res <- pcoa(gw.diss, correction="none")
		}
			
		# Use first three cordinates to speed-up computation
		if(ncol(pco.res$vectors)>3) {
			coords <- as.matrix(pco.res$vectors[,1:3])
			perc.axes <- sum(pco.res$values$Rel_corr_eig[1:3])
			
			if(is.null(perc.axes)) {
				perc.axes <- sum(pco.res$values$Relative_eig[1:3])
			}
		} else {
			coords <- as.matrix(pco.res$vectors)
			perc.axes <- sum(pco.res$values$Rel_corr_eig)
			
			if(is.null(perc.axes)) {
				perc.axes <- sum(pco.res$values$Relative_eig)
			}		
		}
		
		rownames(coords) <- unique(dat$SPECIES) 
		
		# For analysis by year
		weights <- as.matrix(dat %>% group_by(YEAR, SPECIES) %>%
						summarise(abund=mean(abund, na.rm=T), .groups="drop") %>%
						pivot_wider(names_from=SPECIES, values_from=abund, values_fill=0) %>%
						dplyr::select(-YEAR))
		rownames(weights) <- unique(dat$YEAR)
		
		# For cumulative analysis over years
		#	weights <- as.matrix(t(dat$abund))
		#	colnames(weights) <- dat$SPECIES
		#	rownames(weights) <- "abund"
		
		# Computation of FRic requires that the number of species in each year
		# is greater than the number of functional axes + 1. This is done by
		# adjusting the number of functional axes to be no greater than the 
		# number of species - 2.
		
		FDs <- NULL
		
		for(k in 1:nrow(weights)) {
			
			wg <- t(data.frame(wg=weights[k,]) %>% filter(wg>0))
			
			if(ncol(wg)>2) {
				coord.nm <- ifelse(ncol(wg)>ncol(coords)+1, ncol(coords), ncol(wg)-1)
				axs <- cbind(coords[, 1:coord.nm])
				if(ncol(axs)==1) {colnames(axs) <- "Axis.1"}
				FDs <- rbind(FDs, 
						alpha.fd.multidim(
										sp_faxes_coord = axs,
										asb_sp_w = as.matrix(wg),
										ind_vect = c("fric","fdiv","feve","fdis","fspe","fori"),
										details_returned = FALSE,
										verb=FALSE)[[1]] %>%
								dplyr::select(c("fric","fdiv","feve","fdis","fspe","fori")) %>%
								rename(FRic=fric,FDIv=fdiv, FEve=feve, FDis=fdis, FSpe=fspe, FOri=fori)
				
				)
				
			}
		}
		
#		# Computation of functional diversity measures without correcting
#		# the number of functional axes with respect to the number of species
#		# in each year 
#		FDs <- alpha.fd.multidim(
#						sp_faxes_coord = coords,
#						asb_sp_w = weights,
#						ind_vect = c("fric","fdiv","feve","fdis","fspe","fori"),
#						details_returned = FALSE,
#						verb=FALSE)[[1]] %>%
#				dplyr::select(c("fric","fdiv","feve","fdis","fspe","fori")) %>%
#				rename(FRic=fric,FDIv=fdiv, FEve=feve, FDis=fdis, FSpe=fspe, FOri=fori)
		
# 		to use old Villager multidimFD function; NOTE: this requires selecting
# 		the appropriate columns for the apply function below
#		FDs <- multidimFD(coord=coords, weight=weights, verb=FALSE
		## 			Uncomment for plotting
		##			folder_plot='~/Lavori/MPA_timeseries/plot.tmp',
		##			plot_pool=TRUE, nm_asb_plot=row.names(weights),
		##			Faxes_plot=c("Axis.1","Axis.2")
#		)
#		col.sel <- which(colnames(FDs)%in%"FRic")
#		apply(FDs[,col.sel:(col.sel+5)]					
		
# 		to use functions in FD package
#		trait.mat <- as.matrix(gw.diss) #as.matrix(traits)
#		rownames(trait.mat) <- unique(dat$SPECIES)
#		fd.tmp <- dbFD(trait.mat, weights, corr="cailliez")
		
		if(!is.null(FDs)) {
			
			fd.means <- apply(FDs, 2, function(x) mean(x, na.rm=T))
			
		}	else {
			
			perc.axes <- NA
			fd.means <- rep(NA, 6)
		}
		
	} else {
		
		perc.axes <- NA
		fd.means <- rep(NA, 6)
	}
	
	names(fd.means) <- c("mFRic","mFDiv","mFEve","mFDis","mFSpe","mFOri")
	out <- cbind(perc.pcoa=perc.axes,as.data.frame(t(fd.means)))
	return(out)
}





