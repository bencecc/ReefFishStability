# Function to compute robust estimates and standard errors of stability metrics
# using a jackknife procedure (leave-one-out) of species in a a metacommunity
# framework. The function repeats calculations for different spatial ranges for
# MPA and OA sites separately.
#################################################################################

sites_jackknife_stab <- function(df, mpa=c("Yes","No"), cut.range=target.range, dist.dat, ...) {
	
# Split into one or two distance matrices depending on mpa argument:
# A distance matrix between each MPA and all other sites. This
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
	
	# Subset by mpa condition(MPA or OA sites) 
	sub.dat <- df %>% filter(MPA==mpa) %>% mutate(SITE_ID=as.character(SITE_ID))
	
	# Proceed only if there are sites 
	#if(!is.null(sub.dat)&added.mpa.sites>0) {
	if(!is.null(sub.dat)) {
		
		sp <- unique(sub.dat$SPECIES)
		n.spec <- length(sp)
		
		jack_res <- foreach(i = 1:n.spec, .combine="rbind",
						.packages=c("dplyr","codyn","tidyr","modelr","foreach"),
						.export=c("meta_stab","gross_w_func","cor_algo","phi_t3","var_t3","cov_t3","cor_t3")) %dopar%
				{ 
					
					# remove species i
					jack.dat.tmp <- sub.dat %>% filter(SPECIES!=sp[i])
					
					# Computation of stability metrics for sites (either within or outside MPAs
					# as per mpa argument) within target range of distance.
					stab.res <- NULL
					
					for(j in 1:ncol(site.dist)) { 
						
						# Select sites within the specified distance range {0, cut.range[2]} from the focal site  
						site.j <- site.sel %>% dplyr::select(c(SITE_ID,sites[j])) %>%
								filter(if_any(where(is.numeric), ~ .x>=0&.x<=cut.range[2])|
												SITE_ID==sites[j]) %>%
								dplyr::select("SITE_ID")
						
						site.previous.range <- site.sel %>% dplyr::select(c(SITE_ID,sites[j])) %>%
								filter(if_any(where(is.numeric), ~ .x>=0&.x<cut.range[1])|
												SITE_ID==mpa.sites[j]) %>%
								dplyr::select("SITE_ID")
						
						# assess whether the larger range includes at least one new MPA site
						added.mpa.sites <- length(which(sites%in%
												site.j[which(!site.j$SITE_ID%in%site.previous.range$SITE_ID),"SITE_ID"]))
						
						# Move on only if there is at least one new MPA site compared to the previous range and
						# if there are at least two sites within the distance range
						if(added.mpa.sites>0&length(site.j$SITE_ID)>1) {
							
							# This is the dataset for computation with sites (either MPA or OA)
							# within selected distance from focal site sites[j]
							jack.dat <- jack.dat.tmp %>% ungroup() %>%
									filter(SITE_ID%in%site.j$SITE_ID) %>% 
									dplyr::select(SITE_ID, MPA, YEAR, SPECIES, abund)
							
							# check there are at least two sites remaining after removing species i
							# to proceed with the computation (although unlikely, removing species
							# i may have determined the disappearance of a site)
							nsites <- jack.dat %>% group_by(SITE_ID) %>% summarise(sum.abund=sum(abund)) %>%
									filter(sum.abund>0) %>% summarise(n.sites=length(unique(SITE_ID)),
											.groups="drop")
							
							# proceed with the computaion if there are at least two sites
							if(nsites$n.sites>1) {
								
								stab.res <- rbind(stab.res, meta_stab(jack.dat))
								
							}
							
						}
						
						
					}
					
					if(!is.null(stab.res)) {
						
						out <- data.frame(
								t(apply(stab.res[,3:ncol(stab.res)], 2, function(x) mean(x, na.rm=T)))
						)
						
					} else {
						
						out <- NULL
						
					}
					
					out
				}
		
		if(!is.null(jack_res)) {
			
			jnife.mean <- data.frame(Jnife.est="Jnife.mean",
					t(apply(jack_res, 2, function(x) mean(x, na.rm=T))))
			jnife.se <- data.frame(Jnife.est="Jnife.se",
					t(apply(jack_res, 2, function(x) sqrt(var(x, na.rm=T)*(n.spec-1)*(1-1/n.spec)))))
			jnife.var <- data.frame(Jnife.est="Jnife.var",
					t(apply(jack_res , 2, function(x) var(x, na.rm=T)*(n.spec-1)*(1-1/n.spec)))) 
			
			# (NSPEC-1)*(1-1/NSPEC)) = (NSPEC-1)^2/NSPEC is a correction term in the variance formula
			# (see Jackknife_Millar_etal function and corresponding paper)
			
			jnife.est <- data.frame(
					ID=rep(unique(df$ID), 3),
					MPA=rep(mpa, 3),
					DIST=rep(cut.range[2], 3),
					NSPEC=rep(n.spec, 3),
					rbind(jnife.mean, jnife.se, jnife.var)
			)
			
			return(jnife.est)			
			
		}
		
	}
	
}

