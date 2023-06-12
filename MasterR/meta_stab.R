#### ---- Function to perform spatiotemporal stability analysis (e.g., Wilcox et al. ---- ####
#### ---- Ecol Lett 2017 and Lamy et al. Ecology 2019)                               ---- ####
#### ---- Requires columns SITE_ID, MPA, YEAR, SPECIES, abund                         ---- ####
#### ---- ########################################################################## ---- ####

meta_stab <- function(dat, ...) {
		
	dat$SITE_ID <- as.factor(as.vector(dat$SITE_ID))
	
	dat <- dat %>% arrange(SITE_ID, YEAR, SPECIES) %>%
			group_by(SITE_ID, MPA, YEAR, SPECIES) %>%
			summarise(abund=mean(abund), .groups="drop") %>%
			pivot_wider(names_from=SPECIES, values_from=abund, values_fill=0) %>%
			gather('SPECIES', 'abund',
					which(colnames(.)==colnames(.)[4]):which(colnames(.)==colnames(.)[ncol(.)])) %>%
			ungroup()
	
	gamma_dat <- dat %>% group_by(YEAR) %>%
			summarise(sum_abund=sum(abund), nsites=length(unique(SITE_ID)),
					mean_abund=sum(abund)/length(unique(SITE_ID)),
					.groups="drop")
	
	#### ---------------------- STABILITY METRICS ---------------------- ####
# gamma stability
	mu_gamma <- mean(gamma_dat$sum_abund)
	sd_gamma <- sd(gamma_dat$sum_abund)
	gamma_stab <- mu_gamma/sd_gamma
	
# alpha stability
	site_dat_abund <- dat %>%
			group_by(SITE_ID, YEAR) %>% summarise(site_abund=sum(abund),
					.groups="drop") %>%
			group_by(SITE_ID) %>%
			summarise(mean_site_abund=mean(site_abund),
					sd_site_abund=sd(site_abund), .groups="drop")
	alpha_stab <- mu_gamma/sum(site_dat_abund$sd_site_abund)
	
# species stability
	species_dat_abund <- dat %>%
			group_by(SITE_ID, SPECIES) %>% summarise(mean_sp_abund=mean(abund),
					sd_sp_abund=sd(abund), tot_abund=sum(abund),
					.groups="drop") %>% filter(mean_sp_abund>0)
	sum_sd_sp_site <- species_dat_abund %>% group_by(SITE_ID) %>%
			summarise(sum_sp_sd=sum(sd_sp_abund), nspec=n(), .groups="drop")
	cv_sp_site <- sum_sd_sp_site$sum_sp_sd/site_dat_abund$mean_site_abund
	cv_sp <- sum(site_dat_abund$mean_site_abund/mu_gamma*cv_sp_site)
	scaled_species_stab <- 1/cv_sp
	
# population stability
	mean_sp_abund <- dat %>%
			group_by(SPECIES, YEAR) %>% summarise(tot_sp_abund=sum(abund),
					.groups="drop") %>%
			pivot_wider(names_from=SPECIES, values_from=tot_sp_abund, values_fill=0) %>%
			gather('SPECIES', 'abund',
					which(colnames(.)==colnames(.)[2]):which(colnames(.)==colnames(.)[ncol(.)])) %>%
			group_by(SPECIES) %>% summarise(mean_sp_abund=mean(abund),
					sd_sp_abund=sd(abund), .groups="drop")
	
	cv_sp <- mean_sp_abund$sd_sp_abund/mean_sp_abund$mean_sp_abund
	cv_pop <- sum(mean_sp_abund$mean_sp_abund/mu_gamma*cv_sp)
	scaled_population_stab <- 1/cv_pop
	
	#### ---- Data aggregations to compute synchrony metrics ---- ####
	#### ---- based on variances/covariances -------------------- ####
	
# Year x Site total species abundance
#	vcov_dat <-  dat %>% arrange(SITE_ID, YEAR) %>%
#			group_by(SITE_ID, YEAR) %>% summarise(tot_abund=sum(abund),
#					.groups="drop") %>%
#			pivot_wider(names_from=SITE_ID, values_from=tot_abund, values_fill=0) %>%
#			dplyr::select(-YEAR)
	
# SD of summed species abundances within sites
#	sd_sum_sp <- dat %>% 
#			group_by(SPECIES, YEAR) %>% summarise(tot_abund=sum(abund),
#					.groups="drop") %>%
#			group_by(SPECIES) %>% summarise(sd_sum_sp=sd(tot_abund),
#					.groups="drop")
	
# temporal SD of summed species abundances for each site
#	sd_sum_site <- dat %>%
#			group_by(YEAR, SITE_ID) %>% summarise(tot_abund=sum(abund)) %>%
#			group_by(SITE_ID) %>% summarise(sd_sum_site=sd(tot_abund),
#					.groups="drop")
	
# temporal SD of individual species abundances within each site
#	sd_sp_site <- dat %>% 
#			group_by(SPECIES, SITE_ID) %>% summarise(sd_sp_site=sd(abund),
#					.groups="drop") %>%
#			pivot_wider(names_from=SPECIES, values_from=sd_sp_site, values_fill=0) %>%
#			ungroup() %>% dplyr::select(-SITE_ID)
	
	#### --------------- SYNCHRONY METRICS ---------------------- ####
	
# Within-Site species synchrony: from species to individual sites
# [phi_SC_L (phiP->C,L)  in Lamy et al)].
		
	# Unweighted SYNC_GROSS
	# NOTE: gross.sync_func automatically removes columns (species) with all zeros 
#	ssync_gross <- dat %>% pivot_wider(names_from=SPECIES, values_from=abund, values_fill=0) %>%
#			dplyr::select(-c(MPA)) %>% group_by(SITE_ID) %>%
#			group_modify(~ gross_w_func(., weight=FALSE))
	
	# Weighted SYNC_GROSS	
	ssync_gw <- dat %>% pivot_wider(names_from=SPECIES, values_from=abund, values_fill=0) %>%
			dplyr::select(-c(MPA)) %>% group_by(SITE_ID) %>%
			group_modify(~ gross_w_func(., weight=TRUE))
	
	# Detrended SYNC_GROSS - three-term local variance in Leps et al. Ecography
	# NOTE: cor_algo function has been modified to remove empty columns
	# (zero species abundances) and SITE_ID column
#	ssync_gw_det <- dat %>%
#			pivot_wider(names_from=SPECIES, values_from=abund, values_fill=0) %>%
#			dplyr::select(-c(YEAR, MPA)) %>% group_by(SITE_ID) %>%
#			#dplyr::select(where(~ is.numeric(.x) && colSums(data.frame(.x)) > 0)) %>%
#			group_modify(~ data.frame(sync=cor_algo(as.data.frame(.x),
#									method=cor_t3, weighted=TRUE, rm.col=TRUE)))
	
	# Detrended SYNC_GROSS - standard method of linear detrending
	sync_gross.sync_detreg <- dat %>% pivot_wider(names_from=SPECIES, 
					values_from=abund, values_fill=0) %>%
			dplyr::select(-c(MPA)) %>% group_by(SITE_ID) %>%
			group_modify(~ gross_w_func(., weight=TRUE, det=TRUE)) %>%
			ungroup()
	
	# Computation of LOREAU sync from function synchrony
	ssync_loreau <- synchrony(df=dat,
			time.var='YEAR', species.var='SPECIES',
			abundance.var='abund', replicate.var='SITE_ID', metric='Loreau')
	
	# Detrended SYNC LOREAU - three-term local variance in Leps et al. Ecography
#	ssync_loreau_det <- dat %>% pivot_wider(names_from=SPECIES,
#					values_from=abund, values_fill=0) %>%
#			dplyr::select(-c(MPA, YEAR)) %>% group_by(SITE_ID) %>%
#			group_modify(~ data.frame(sync=phi_t3(as.data.frame(.))))
	
	# detrended SYNC LOREAU - standard method of linear detrending
	ssync_loreau_detreg <- dat %>% group_by(SITE_ID) %>% 
			arrange(SITE_ID, YEAR, SPECIES) %>%
			dplyr::select(SITE_ID, YEAR, SPECIES, abund) %>%
			pivot_wider(names_from=SPECIES, values_from=abund, values_fill=0) %>%
			mutate(across(which(colnames(.)==colnames(.)[2]):
									which(colnames(.)==colnames(.)[ncol(.) - 1]),
							~ resid(lm(.x ~ YEAR)))) %>%
			gather('SPECIES', 'abund',
					which(colnames(.)==colnames(.)[3]):which(colnames(.)==colnames(.)[ncol(.)])) %>%
			arrange(SITE_ID, YEAR, SPECIES) %>% ungroup() %>%
			synchrony(time.var="YEAR", abundance.var="abund",
					species.var="SPECIES", replicate.var='SITE_ID', metric='Loreau')
	
#	scaled_ssync_gross <- (site_dat_abund$mean_site_abund/mu_gamma)*ssync_gross$gross.sync
	scaled_ssync_gw <- (site_dat_abund$mean_site_abund/mu_gamma)*ssync_gw$gross.sync
#	scaled_ssync_gw_det <- (site_dat_abund$mean_site_abund/mu_gamma)*ssync_gw_det$sync
	scaled_ssync_gw_detreg <- (site_dat_abund$mean_site_abund/mu_gamma)*sync_gross.sync_detreg$gross.sync
	scaled_ssync_loreau <- (site_dat_abund$mean_site_abund/mu_gamma)*ssync_loreau$synchrony
#	scaled_ssync_loreau_det <- (site_dat_abund$mean_site_abund/mu_gamma)*ssync_loreau_det$sync
	scaled_ssync_loreau_detreg <- (site_dat_abund$mean_site_abund/mu_gamma)*ssync_loreau_detreg$synchrony
	
#	avg_sp_sync_gross <- sum(scaled_ssync_gross)
	avg_sp_sync_gw <- sum(scaled_ssync_gw)
#	avg_sp_sync_gw_det <- sum(scaled_ssync_gw_det, na.rm=T)
	avg_sp_sync_gw_detreg <- sum(scaled_ssync_gw_detreg, na.rm=T)
	avg_sp_sync_loreau <- sum(scaled_ssync_loreau)
#	avg_sp_sync_loreau_det <- sum(scaled_ssync_loreau_det)
	avg_sp_sync_loreau_detreg <- sum(scaled_ssync_loreau_detreg)
	
#	sp_sync_var <- sd_sum_site$sd_sum_site/rowSums(sd_sp_site)
#	weight_sp <- rowSums(sd_sp_site)/sum(sd_sp_site)
#	avg_sp_sync_var <- sum(weight_sp*sp_sync_var)
	
# Metacommunity assemblage-level spatial synchrony: from sites to the metacommunity
# [phi_C_LR (phyC,L->R) in Lamy et al.]
	site_abund <- dat %>%
			group_by(SITE_ID, YEAR) %>%
			summarise(site_abund=sum(abund), .groups="drop") %>%
			filter(site_abund>0)
	
	# Unweighted SYNC_GROSS
#	spatial_sync_gross <- site_abund %>% 
#		pivot_wider(names_from=SITE_ID, values_from=site_abund, values_fill=0) %>%
#		group_modify(~ gross_w_func(., weight=FALSE))	
	
#	spatial_synch_gross <- site_abund %>%
#			synchrony(time.var='YEAR', species.var='SITE_ID',
#					abundance.var='site_abund', metric='Gross')
	
	# Weighted SYNC_GROSS
	spatial_sync_gw <- site_abund %>% 
			pivot_wider(names_from=SITE_ID, values_from=site_abund, values_fill=0) %>%
			group_modify(~ gross_w_func(., weight=TRUE))
	
	# Detrended SYNC_GROSS - three-term local variance in Leps et al. Ecography
#	spatial_sync_gw_det <- site_abund %>% 
#			pivot_wider(names_from=SITE_ID, values_from=site_abund, values_fill=0) %>%
#			dplyr::select(-YEAR) %>%
#			group_modify(~ data.frame(sync=cor_algo(as.data.frame(.x), method=cor_t3, weighted=TRUE)))
	
	# Detrended SYNC_GROSS - standard method of linear detrending
	spatial_sync_gw_detreg <- site_abund %>% 
			pivot_wider(names_from=SITE_ID, values_from=site_abund, values_fill=0) %>%
			group_modify(~ gross_w_func(., weight=TRUE, det=TRUE))
	
	# computation of LOREAU sync from function synchrony
	spatial_sync_loreau <- site_abund %>%
			synchrony(time.var='YEAR', species.var='SITE_ID',
					abundance.var='site_abund', metric='Loreau')
	
	# Detrended SYNC LOREAU - three-term local variance in Leps et al. Ecography
#	spatial_sync_loreau_det <- site_abund %>%
#			pivot_wider(names_from=SITE_ID, values_from=site_abund, values_fill=0) %>%
#			dplyr::select(-YEAR) %>% 
#			group_modify(~ data.frame(sync=phi_t3(as.data.frame(.))))
	
	# detrended SYNC LOREAU - standard method of linear detrending
	spatial_sync_loreau_detreg <- site_abund %>%
			pivot_wider(names_from=SITE_ID, values_from=site_abund, values_fill=0) %>%
			mutate(across(which(colnames(.)==colnames(.)[2]):
									which(colnames(.)==colnames(.)[ncol(.)]),
							~ resid(lm(.x ~ YEAR)))) %>%
			gather('SITE_ID', 'site_abund',
					which(colnames(.)==colnames(.)[2]):which(colnames(.)==colnames(.)[ncol(.)])) %>%
			ungroup() %>%
			synchrony(time.var="YEAR", abundance.var="site_abund",
					species.var="SITE_ID", metric='Loreau')
	
	## Year x Site variance/covariance in total species abundance 
#	vcov <- cov(vcov_dat)
	## Metacommunity assemblage-level spatial synchrony based on sqrt of Loreau & de Mazancourt
	## as in Lamy et al. Ecology.
#	spatial_sync_var <- sqrt(sum(vcov))/((sum(sd_sum_site$sd_sum_site)))
	
# population synchrony: from population to metacommunity
# [phi_S_LR (phiP,L->R in Table 1) in Lamy et al.]
	
	# Unweighted SYNC_GROSS	
#	pop_sync_gross_tmp <- dat %>% pivot_wider(names_from=SITE_ID, values_from=abund, values_fill=0) %>%
#		dplyr::select(-c(MPA)) %>%  group_by(SPECIES) %>%  #filter(SPECIES%in%"Aetapcus maculatus") %>% 
#		group_modify(~ gross_w_func(., weight=FALSE))
	
#	pop_synch_gross_tmp <- synchrony(df=dat,
#			time.var='YEAR', species.var='SITE_ID',
#			abundance.var='abund', replicate.var='SPECIES', metric='Gross')

	# Weighted SYNC_GROSS
	pop_sync_gw_tmp <- dat %>% pivot_wider(names_from=SITE_ID, values_from=abund, values_fill=0) %>%
			dplyr::select(-c(MPA)) %>% group_by(SPECIES) %>%
			group_modify(~ gross_w_func(., weight=TRUE)) %>%
			ungroup()
	
	# Detrended SYNC_GROSS - three-term local variance in Leps et al. Ecography
	#### ---- To run cor_algo for detrended and weighted gross synchrony, first select species
	#### ---- that have at least one occurrence in a year; the function stops when a species
	#### ---- does not have at least one occurrence in one year in a given site. Selection of
	#### ---- species is done using mean_sp_abund as per the code below.
#	dat_pop <- dat %>% ungroup() %>%
#				pivot_wider(names_from=SITE_ID, values_from=abund, values_fill=0) %>%
#				dplyr::select(-c(MPA, YEAR)) #%>% #group_by(SPECIES) %>%
#	
#	pop_sync_gw_det_tmp <- NULL
#		
#		for(i in unique(dat_pop$SPECIES)) {
#			dat.pop.tmp <- dat_pop[which(dat_pop$SPECIES==i),] %>% dplyr::select(-SPECIES)
#			dat.pop.i <- dat.pop.tmp[,which(colSums(dat.pop.tmp)>0)] %>%
#					group_modify(~ data.frame(sync=cor_algo(as.data.frame(.x), method=cor_t3,
#											weighted=TRUE)))
#			pop_sync_gw_det_tmp <- rbind(pop_sync_gw_det_tmp, dat.pop.i)		
#					
#		}
		
	# Detrended SYNC_GROSS - standard method of linear detrending
	pop_sync_gw_detreg_tmp <- dat %>% pivot_wider(names_from=SITE_ID,
					values_from=abund, values_fill=0) %>%
			dplyr::select(-c(MPA)) %>% group_by(SPECIES) %>%
			group_modify(~ gross_w_func(., weight=TRUE, det=TRUE)) %>%
			ungroup()
	
	# Computation of LOREAU sync from function synchrony
	pop_sync_loreau_tmp <- synchrony(df=dat %>%
					pivot_wider(names_from=SPECIES, values_from=abund,
							values_fill=0) %>%
					gather('SPECIES', 'abund',
							which(colnames(.)==colnames(.)[4]):which(colnames(.)==colnames(.)[ncol(.)])),
			time.var='YEAR', species.var='SITE_ID',
			abundance.var='abund', replicate.var='SPECIES', metric='Loreau')
	
	# Detrended SYNC LOREAU - three-term local variance in Leps et al. Ecography
#	pop_sync_loreau_tmp_det <- dat %>% pivot_wider(names_from=SITE_ID, values_from=abund,
#					values_fill=0) %>%
#			dplyr::select(-c(MPA, YEAR)) %>% group_by(SPECIES) %>%
#			group_modify(~ data.frame(sync=phi_t3(as.data.frame(.))))
	
	# detrended SYNC LOREAU - standard method of linear detrending
	pop_sync_loreau_tmp_detreg <- dat %>% group_by(SPECIES) %>%
			dplyr::select(SITE_ID, YEAR, SPECIES, abund) %>%
			pivot_wider(names_from=SITE_ID, values_from=abund, values_fill=0) %>%
			arrange(SPECIES, YEAR) %>%
			mutate(across(which(colnames(.)==colnames(.)[2]):
									which(colnames(.)==colnames(.)[ncol(.) - 1]),
							~ resid(lm(.x ~ YEAR)))) %>%
			gather('SITE_ID', 'abund',
					which(colnames(.)==colnames(.)[3]):which(colnames(.)==colnames(.)[ncol(.)])) %>%
			arrange(SITE_ID, YEAR, SPECIES) %>% ungroup() %>%
			synchrony(time.var="YEAR", abundance.var="abund",
					species.var="SITE_ID", replicate.var='SPECIES', metric='Loreau')
			
#	pop_sync_gross <- pop_sync_gross_tmp %>% filter(SPECIES%in%unique(mean_sp_abund$SPECIES))		
	pop_sync_gw <- pop_sync_gw_tmp %>% filter(SPECIES%in%unique(mean_sp_abund$SPECIES))
	pop_sync_loreau <- pop_sync_loreau_tmp %>% filter(SPECIES%in%unique(mean_sp_abund$SPECIES))	
#	pop_sync_loreau_det <- pop_sync_loreau_tmp_det %>% filter(SPECIES%in%unique(mean_sp_abund$SPECIES))	
	pop_sync_loreau_detreg <- pop_sync_loreau_tmp_detreg %>% filter(SPECIES%in%unique(mean_sp_abund$SPECIES))	
	
#	scaled_pop_sync_gross <- (mean_sp_abund$mean_sp_abund/mu_gamma)*pop_sync_gross$gross.sync
	scaled_pop_sync_gw <- (mean_sp_abund$mean_sp_abund/mu_gamma)*pop_sync_gw$gross.sync
#	scaled_pop_sync_gw_det <- (mean_sp_abund$mean_sp_abund/mu_gamma)*pop_sync_gw_det_tmp$sync
	scaled_pop_sync_gw_detreg <- (mean_sp_abund$mean_sp_abund/mu_gamma)*as.numeric(pop_sync_gw_detreg_tmp$gross.sync)
	scaled_pop_sync_loreau <- (mean_sp_abund$mean_sp_abund/mu_gamma)*pop_sync_loreau$synchrony
#	scaled_pop_sync_loreau_det <- (mean_sp_abund$mean_sp_abund/mu_gamma)*pop_sync_loreau_det$sync
	scaled_pop_sync_loreau_detreg <- (mean_sp_abund$mean_sp_abund/mu_gamma)*pop_sync_loreau_detreg$synchrony
	
#	avg_pop_sync_gross <- sum(scaled_pop_sync_gross, na.rm=T)
	avg_pop_sync_gw <- sum(scaled_pop_sync_gw, na.rm=T)
#	avg_pop_sync_gw_det <- sum(scaled_pop_sync_gw_det, na.rm=T)
	avg_pop_sync_gw_detreg <- sum(scaled_pop_sync_gw_detreg, na.rm=T)
		
	avg_pop_sync_loreau <- sum(scaled_pop_sync_loreau)
#	avg_pop_sync_loreau_det <- sum(scaled_pop_sync_loreau_det)
	avg_pop_sync_loreau_detreg <- sum(scaled_pop_sync_loreau_detreg)
	
#	pop_sync_var <- sd_sum_sp$sd_sum_sp/colSums(sd_sp_site)
#	weight_pop <- colSums(sd_sp_site)/sum(sd_sp_site) 
#	avg_pop_sync_var <- sum(weight_pop*pop_sync_var)
	
# From the metapopulation to the metacommunity
# [phi_SC_R in (phiP->C,R in Table 1) Lamy et al.]
	sp_sum_site <- dat %>% 
			group_by(SPECIES, YEAR) %>%
			summarise(tot_sp_abund=sum(abund), .groups="drop") 
	
	# Unweighted SYNC_GROSS
#	mm_sync_gross <- sp_sum_site %>%
#			pivot_wider(names_from=SPECIES, values_from=tot_sp_abund, values_fill=0) %>%
#			group_modify(~ gross_w_func(., weight=FALSE))
#	mm_sync_gross <- synchrony(sp_sum_site, time.var='YEAR', species.var='SPECIES',
#			abundance.var='tot_sp_abund', metric='Gross')
	
	# Weighted SYNC_GROSS	
	mm_sync_gw <- sp_sum_site %>%
			pivot_wider(names_from=SPECIES, values_from=tot_sp_abund, values_fill=0) %>%
			group_modify(~ gross_w_func(., weight=TRUE))
	
	# Ddetrended SYNC_GROSS - three-term local variance in Leps et al. Ecography
#	mm_sync_gw_det <- sp_sum_site %>%
#			pivot_wider(names_from=SPECIES, values_from=tot_sp_abund, values_fill=0) %>%
#			dplyr::select(-YEAR) %>%
#			group_modify(~ data.frame(sync=cor_algo(as.data.frame(.x),
#									method=cor_t3, weighted=TRUE)))
	
	# Detrended SYNC_GROSS - standard method of linear detrending
	mm_sync_gw_detreg <- sp_sum_site %>%
			pivot_wider(names_from=SPECIES, values_from=tot_sp_abund, values_fill=0) %>%
			group_modify(~ gross_w_func(., weight=TRUE, det=TRUE))
	
	# Computation of LOREAU sync from function synchrony
	mm_sync_loreau <- synchrony(sp_sum_site, time.var='YEAR', species.var='SPECIES',
			abundance.var='tot_sp_abund', metric='Loreau')
	
	# Detrended SYNC LOREAU - three-term local variance in Leps et al. Ecography
#	mm_sync_loreau_det <- sp_sum_site %>%
#			pivot_wider(names_from=SPECIES, values_from=tot_sp_abund, values_fill=0) %>%
#			dplyr::select(-YEAR) %>%
#			group_modify(~ data.frame(sync=phi_t3(as.data.frame(.))))
	
	# detrended SYNC LOREAU - standard method of linear detrending
	mm_sync_loreau_detreg <- sp_sum_site %>%
			pivot_wider(names_from=SPECIES, values_from=tot_sp_abund, values_fill=0) %>%
			mutate(across(which(colnames(.)==colnames(.)[2]):
									which(colnames(.)==colnames(.)[ncol(.)]),
							~ resid(lm(.x ~ YEAR)))) %>%
			gather('SPECIES', 'tot_sp_abund',
					which(colnames(.)==colnames(.)[2]):which(colnames(.)==colnames(.)[ncol(.)])) %>%
			synchrony(time.var="YEAR", abundance.var="tot_sp_abund",
					species.var="SPECIES", metric='Loreau')
	
	# using sqrt of Loreau and de Mazancourt
#	mm_sync_var <- sqrt(sum(vcov))/sum(sd_sum_sp$sd_sum_sp)
	
	#### ----------------------- METACOMMUNITY DIVERSITY INDICES ---------------------- ####
	#### ---- n_spec is equal to effective number of species or Hill numbers (ens) ---- ####
	#### ---- with q=0; shannon_ens is ens with q=1; simp_dom is ens with q=2;     ---- ####
	#### ---- simp_div is the gini-simpson index and it is equal to Hurlbert's PIE ---- ####
	tot_sp_abund <- sum(dat$abund, na.rm=T)
	sp_rel_abund <- dat %>% group_by(SPECIES) %>%
			summarise(sp_rel_abund=sum(abund, na.rm=T)/tot_sp_abund)
	
	n_spec <- sp_rel_abund %>% drop_na(sp_rel_abund) %>%
			filter(sp_rel_abund>0) %>% summarise(n_spec=n())
	shannon <- -sum(sp_rel_abund[,2]*log(sp_rel_abund[,2]))
	shannon_ens <- exp(shannon)
	simp_div <- 1-sum(sp_rel_abund[,2]^2)
	simp_dom <- 1/sum(sp_rel_abund[,2]^2)
	simp_evenn <- simp_dom/nrow(sp_rel_abund)
	
	# Hurlbert rarefacion
#	raref.div <- dat %>% group_by(SPECIES) %>%
#		summarise(tot.abund=sum(floor(exp(abund)))) %>%
#		pivot_wider(names_from=SPECIES, values_from=tot.abund, values_fill=0) %>%
#		group_modify(~ data.frame(simp=vegan::rarefy(as.data.frame(.), floor(tot_sp_abund/n_spec$n_spec)))) %>%
#		ungroup() 
	
	# beta diversity - multiplicative partitioning: gamma/mean(alpha) using simp_div
	alpha.beta.div <- dat %>% group_by(SITE_ID, SPECIES) %>%
			summarise(tot.abund=sum(abund), .groups="drop") %>%
			pivot_wider(names_from=SPECIES, values_from=tot.abund, values_fill=0) %>%
			group_by(SITE_ID) %>%
			group_modify(~ data.frame(simp=vegan::diversity(as.data.frame(.), "inv"))) %>%
			ungroup() %>%
			summarise(beta=simp_dom/mean(simp), mean.alpha=mean(simp))
			
# loreau metrics are returned as 1- metric - gross metrics as synchrony*-1
# thus, all metrics are returned as asynchrony measures
	out_df <- data.frame(
			MPA=dat[1,"MPA"],
			N_SITES=length(unique(dat$SITE_ID)),
			gamma_stab=gamma_stab,
			alpha_stab=alpha_stab,
			species_stab=scaled_species_stab,
			pop_stab=scaled_population_stab,
			
#			loc_sp_comm_async_gross=avg_sp_sync_gross*-1,
			loc_sp_comm_async_gross_w=avg_sp_sync_gw*-1,
#			loc_sp_comm_async_gross_w_det=avg_sp_sync_gw_det*-1,
			loc_sp_comm_async_gross_w_detreg=avg_sp_sync_gw_detreg*-1,
#			loc_sp_comm_async_sqrtlor=1-avg_sp_sync_var,
			loc_sp_comm_async_loreau=1-avg_sp_sync_loreau,
#			loc_sp_comm_async_loreau_det=1-avg_sp_sync_loreau_det,
			loc_sp_comm_async_loreau_detreg=1-avg_sp_sync_loreau_detreg,
			
#			loc_comm_metacom_async_gross=spatial_sync_gross$gross.sync*-1,
			loc_comm_metacom_async_gross_w=spatial_sync_gw$gross.sync*-1,
#			loc_comm_metacom_async_gross_w_det=spatial_sync_gw_det$sync*-1,
			loc_comm_metacom_async_gross_w_detreg=spatial_sync_gw_detreg$gross.sync*-1,
#			loc_comm_metacomm_async_sqrtlor=1-spatial_sync_var,
			loc_comm_metacom_async_loreau=1-spatial_sync_loreau,
#			loc_comm_metacom_async_loreau_det=1-spatial_sync_loreau_det$sync,
			loc_comm_metacom_async_loreau_detreg=1-spatial_sync_loreau_detreg,
			
#			meta_pop_metacom_async_gross=avg_pop_sync_gross*-1,
			meta_pop_metacom_async_gross_w=avg_pop_sync_gw*-1,
#			meta_pop_metacom_async_gross_w_det=avg_pop_sync_gw_det*-1,
			meta_pop_metacom_async_gross_w_detreg=avg_pop_sync_gw_detreg*-1,
#			meta_pop_metacom_async_sqrtlor=1-avg_pop_sync_var,
			meta_pop_metacom_async_loreau=1-avg_pop_sync_loreau,
#			meta_pop_metacom_async_loreau_det=1-avg_pop_sync_loreau_det,
			meta_pop_metacom_async_loreau_detreg=1-avg_pop_sync_loreau_detreg,
			
#			meta_sp_metacom_async_gross=mm_sync_gross$gross.sync*-1,
			meta_sp_metacom_async_gross_w=mm_sync_gw$gross.sync*-1,
#			meta_sp_metacom_async_gross_w_det=mm_sync_gw_det$sync*-1,
			meta_sp_metacom_async_gross_w_detreg=mm_sync_gw_detreg$gross.sync*-1,
#			meta_sp_metacom_async_sqrtlor=1-mm_sync_var,
			meta_sp_metacom_async_loreau=1-mm_sync_loreau,
#			meta_sp_metacom_async_loreau_det=1-mm_sync_loreau_det$sync,
			meta_sp_metacom_async_loreau_detreg=1-mm_sync_loreau_detreg,
					
			N_SPEC=n_spec$n_spec,
			SHANNON=shannon,
			SHANNON_ENS=shannon_ens,
			SIMP_DIV=simp_div,
			SIMP_DOM=simp_dom,
			SIMP_EVENN=simp_evenn,
			MEAN.ALPHA=alpha.beta.div$mean.alpha,
			BETA.DIV=alpha.beta.div$beta
	)
	
	out_df
	
}






