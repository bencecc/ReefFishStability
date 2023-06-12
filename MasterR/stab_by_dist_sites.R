# Function to compute stability measures as a function of distance from a set of
# sites; this can be applied to any geographic region (e.g., ECOREGION) with multiple
# MPA and/or OA sites. Distances among all sites are obtained with function
# mpa_LeastCostDist and the resulting matrix provides input to argument dist.dat.
# The dist.dat is subsetted into two distance data frames seaprately for MPAs and OAs
# to provide distance ranges between each site and all other sites for each level of
# protection. The function cycles over the columns of these distance data frames
# (thus, separately for MPA and OA levels) to identify all sites within the current
# distance range (set by cut.range) from the focal site (the column name). Stability
# and asynchrony metrics are then computed among the sites within the distance range
# (through function meta_stab). This is repeated across all sites (columns of the
# distance data frames) and values are then averaged for each metric. This generates
# a unique mean value of each metric for MPA or OA protection levels within the given
# distance range. A distance range is examined if it includes at least one new MPA
# site compared to the MPA sites that were present in the previous range, the maximum
# number of protected and unprotected sites contributing to the stability metrics and
# the maximum distance among sites.
######################################################################################

stab_by_dist_sites <- function(df, dist.dat, cut.range, mpa=c("Yes","No"), ...) {
	
# Split into one or two distance matrices depending on mpa argument:
# A distance matrix between each MPA and all other sites. This ù
# matrix is used to identify sites for comparison within given distance
# ranges; it is also used to determine if the current distance range
# includes at least one new MPA site compared to the previous distance range.
	
	mpa.sites <- (df %>% 
				mutate(SITE_ID=as.character(SITE_ID)) %>%
				filter(MPA=="Yes") %>% distinct(SITE_ID) %>% dplyr::select(SITE_ID))$SITE_ID
	mpa.dist <- dist.dat[which(rownames(dist.dat)%in%mpa.sites),
			which(colnames(dist.dat)%in%mpa.sites)] # only MPA sites in cols and rows
	# Rownames as a variable; this is used to verify that expanding
	# range includes at least one new MPA site
	mpa.site.sel <- mpa.dist %>% 
			tibble::rownames_to_column(var="SITE_ID") 
	
	sites=mpa.sites
	site.dist=dist.dat[,
			which(colnames(dist.dat)%in%mpa.sites)] # colnames MPA sites, rows include all sites (MPAs and OAs)
	site.sel <- site.dist %>% 
			tibble::rownames_to_column(var="SITE_ID") 
	
	# Cycle over sites; for each focal site, select all sites within a range of distances
	# between 0 and the upper level of cut.range. Proceed with computations only if there
	# are at least two sites within a given range
	stab.res <- foreach(i = 1:ncol(site.dist),
					.packages=c("dplyr","codyn","tidyr","modelr"),
					.export=c("meta_stab", "gross_w_func","cor_algo",
							"phi_t3","var_t3","cov_t3","cor_t3")) %dopar% { 													
					
					# Check sites within current distance range from focal site
					site.i <- site.sel %>% dplyr::select(c(SITE_ID,sites[i])) %>%
							filter(if_any(where(is.numeric), ~ .x>=0&.x<=cut.range[2])|
											SITE_ID==sites[i]) %>%
							dplyr::select("SITE_ID")
					
					# sites in the previous spatial range determined as < cut.range[1]
					site.previous.range <- site.sel %>% dplyr::select(c(SITE_ID,sites[i])) %>%
							filter(if_any(where(is.numeric), ~ .x>=0&.x<cut.range[1])|
											SITE_ID==mpa.sites[i]) %>%
							dplyr::select("SITE_ID")
					
					# assess whetehr the larger range includes at least one new MPA site
					added.mpa.sites <- length(which(sites%in%
											site.i[which(!site.i$SITE_ID%in%site.previous.range$SITE_ID),"SITE_ID"]))
					
					# proceed only if there is at least one new MPA site compared to the previous range;
					if(added.mpa.sites>0&length(site.i$SITE_ID)>1) {
						
					# number of sites
					n.sites <- (df %>% ungroup() %>% filter(MPA==mpa&SITE_ID%in%site.i$SITE_ID) %>%
							dplyr::select(SITE_ID) %>% distinct() %>% summarise(n.sites=n(), .groups="drop"))$n.sites
					
					# number of MPAs
					n.mpa <- (df %>% ungroup() %>% filter(MPA==mpa&SITE_ID%in%site.i$SITE_ID) %>%
								dplyr::select(MPA_NAME) %>% distinct() %>% summarise(n.mpa=n(), .groups="drop"))$n.mpa
										
					# yearly mean sampled area
					sa <- (df %>% ungroup() %>% filter(MPA==mpa&SITE_ID%in%site.i$SITE_ID) %>%
								summarise(sa=mean(SAMPLED_AREA), .groups="drop"))$sa

					dr <- max(site.dist)
					
					# Move on if there are at least two sites within the cut.range[2] distance.
					if(n.sites>1) {
						
						stab.sites <- meta_stab(df %>% filter(MPA==mpa&SITE_ID%in%site.i$SITE_ID))
											
						out <- list(
								stab.sites,
								n.mpa,
								n.sites,
								sa,
								dr
								)
						out
					}
					
				}
				
			}
	
	check.empty.list <- plyr::compact(stab.res)
	
	# If list is not empty, average across sites to obtain a mean value for protected
	# and unprotected conditions and extract max number of MPAs and protected and
	# unprotected sites involved
	
	if(length(check.empty.list)>0) {
		
		n.mpa <- max(unlist(lapply(stab.res, function(x) x[[2]])))
		n.sites <- max(unlist(lapply(stab.res, function(x) x[[3]])))
		sa <- mean(unlist(lapply(stab.res, function(x) x[[4]])))
		dr <- max(unlist(lapply(stab.res, function(x) x[[5]])))
		
		stab.out <- lapply(stab.res, function(x) x[[1]]) %>% bind_rows() %>%
				summarise_if(is.numeric, mean)
		
		obs_stab <- data.frame(MPA=mpa, DIST=cut.range[2], N_MPA=n.mpa, MAX_N_SITES=n.sites,
				SAMPLED_AREA=sa, DIST_RANGE=dr, stab.out) %>%
				rename(MEAN_N_SITES=N_SITES) %>% relocate(MEAN_N_SITES, .before=MAX_N_SITES)
		
		# OUTPUT OF ORIGINAL STAB METRICS 
		return(obs_stab)
		
	}
	
}




