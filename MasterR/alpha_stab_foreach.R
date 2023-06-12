#### ---- Function to perform alpha (within site) stabilty analysis on workstation ---- ####
#### ---- using the combined datatset of fish abundances across studies. See       ---- ####
#### ---- function alpha_stab.R in folder cluster for same function to run on HPC  ---- ####

alpha_stab_foreach <- function(df, ...) {
		
	ids <- unique(df$ID)
	
	for(j in 1:length(ids)) {
		
		study_id <- ids[j]
		
		dat <- df %>% filter(ID%in%study_id)
		
		site.id <- unique(dat$SITE_ID)
		
		
		out <- foreach(i=1:length(site.id),
						.packages=c("dplyr","codyn","tidyr","FD","mFD"),
						.combine="rbind") %dopar% { 
					
					cat('Doing SITE_ID ', i, ' of ',
							length(site.id),
							' in ID ', study_id, '\n', sep = '')
					
					sp.sel <- dat %>% filter(SITE_ID%in%site.id[i]) %>% group_by(SPECIES) %>%
							summarise(sum.abund=sum(abund)) %>% filter(sum.abund>0)
					
					dat_sub <- dat %>% filter(SITE_ID%in%site.id[i]) %>% arrange(YEAR)
					
					sel <- dat_sub %>% group_by(YEAR) %>% distinct(SPECIES) %>%
							summarise(nspec=n())
					
					if(!any(sel$nspec<2)&nrow(sel)>4) {
						
						# add zeros in long format
						dat_sub_sub <- dat_sub %>% group_by(SPECIES) %>%
								pivot_wider(names_from=SPECIES, values_from=abund, values_fill=0) %>%
								gather('SPECIES', 'abund',
										which(colnames(.)==colnames(.)[27]):
												which(colnames(.)==colnames(.)[ncol(.)])) %>%
								arrange(YEAR, SPECIES)
						
						# compute mean and standard deviation of total and species abundances
						mean_sd_tot_abund <- dat_sub_sub %>% group_by(YEAR) %>%
								summarise(sabund=sum(abund)) %>%
								summarise(mabund=mean(sabund), sdabund=sd(sabund))
						mean_sd_sp_abund <- dat_sub_sub %>% group_by(SPECIES) %>%
								summarise(mabund=mean(abund), sdabund=sd(abund)) %>%
								summarise(m_sp_abund=mean(mabund), sd_sp_abund=mean(sdabund))
						
						# compute diversity indices
						tot_sp_abund <- sum(dat_sub_sub$abund, na.rm=T)
						sp_rel_abund <- dat_sub_sub %>% group_by(SPECIES) %>%
								summarise(sp_rel_abund=sum(abund, na.rm=T)/tot_sp_abund)
						
						shannon <- -sum(sp_rel_abund[,2]*log(sp_rel_abund[,2]))
						simp_div <- 1-sum(sp_rel_abund[,2]^2)
						simp_dom <- 1/sum(sp_rel_abund[,2]^2)
						simp_evenn <- simp_dom/nrow(sp_rel_abund)
						
						# weighted SYNC_GROSS
						sync_gross_w <- dat_sub_sub %>% ungroup() %>% arrange(YEAR, SPECIES) %>%
								dplyr::select(YEAR, SPECIES, abund) %>%
								pivot_wider(names_from=SPECIES, values_from=abund, values_fill=0) %>%
								group_modify(~ gross_w_func(., weight=TRUE), det=FALSE)
						
						# detrended SYNC_GROSS - three-term local variance in Leps et al. Ecography
#						sync_gross_w_det <- dat_sub_sub %>% ungroup() %>% arrange(YEAR, SPECIES) %>%
#						dplyr::select(YEAR, SPECIES, abund) %>%
#						pivot_wider(names_from=SPECIES, values_from=abund, values_fill=0) %>%
#						dplyr::select(-YEAR) %>%
#						group_modify(~ cor_algo(as.data.frame(.), method=cor_t3, weighted=TRUE))
						
						# detrended SYNC_GROSS - standard method of linear detrending
						sync_gross_w_detreg <- dat_sub_sub %>% ungroup() %>% arrange(YEAR, SPECIES) %>%
								dplyr::select(YEAR, SPECIES, abund) %>%
								pivot_wider(names_from=SPECIES, values_from=abund, values_fill=0) %>%
								group_modify(~ gross_w_func(., weight=TRUE, det=TRUE))
						
						# unweighted SYNC_GROSS
#						sync_gross <- dat_sub_sub %>% ungroup() %>% arrange(YEAR, SPECIES) %>%
#						dplyr::select(YEAR, SPECIES, abund) %>%
#						pivot_wider(names_from=SPECIES, values_from=abund, values_fill=0) %>%
#						group_modify(~ gross_w_func(., weight=FALSE, det=FALSE))
						##				sync_gross <- dat_sub_sub %>% arrange(YEAR, SPECIES) %>% 
						##						synchrony(time.var="YEAR", abundance.var="abund",
						##						species.var="SPECIES", metric='Gross')
						
						# computation of LOREAU sync from function synchrony
						sync_loreau <- dat_sub_sub %>% arrange(YEAR, SPECIES) %>% 
								synchrony(time.var="YEAR", abundance.var="abund",
										species.var="SPECIES", metric='Loreau')
						
						# detrended SYNC LOREAU - three-term local variance in Leps et al. Ecography
#						sync_loreau_det <- dat_sub_sub %>% ungroup() %>% arrange(YEAR, SPECIES) %>%
#						dplyr::select(YEAR, SPECIES, abund) %>%
#						pivot_wider(names_from=SPECIES, values_from=abund, values_fill=0) %>%
#						dplyr::select(-YEAR) %>%
#						group_modify(~ phi_t3(as.data.frame(.)))
						
						# detrended SYNC LOREAU - standard method of linear detrending
						sync_loreau_detreg <- dat_sub_sub %>% ungroup() %>% arrange(YEAR, SPECIES) %>%
								dplyr::select(YEAR, SPECIES, abund) %>%
								pivot_wider(names_from=SPECIES, values_from=abund, values_fill=0) %>%
								mutate(across(which(colnames(.)==colnames(.)[2]):
														which(colnames(.)==colnames(.)[ncol(.)]),
												~ resid(lm(.x ~ YEAR)))) %>%
								gather('SPECIES', 'abund',
										which(colnames(.)==colnames(.)[2]):which(colnames(.)==colnames(.)[ncol(.)])) %>%
								arrange(YEAR, SPECIES) %>% 
								synchrony(time.var="YEAR", abundance.var="abund",
										species.var="SPECIES", metric='Loreau')
						
						# computation of synchrony as sqrt of Loreau & de Mazancourt; code below is
						# equal to sqrt(sync_loreau)
#						sd_sum_sp <- dat_sub_sub %>% group_by(YEAR) %>% summarise(sum_abund=sum(abund)) %>%
#						ungroup() %>% summarise(sd_sum_sp=sd(sum_abund))
#						sd_sp <- dat_sub_sub %>% group_by(SPECIES) %>% summarise(sd_sp=sd(abund))
#						sync_var <- sd_sum_sp/sum(sd_sp$sd_sp)
						
						# assemblage-level CV
						cv_tot_abund <- dat_sub_sub %>% group_by(YEAR) %>%
								summarise(tot_abund=sum(abund, na.rm=T)) %>%
								summarise(cv_tot_abund=sd(tot_abund, na.rm=T)/mean(tot_abund, na.rm=T))
						
						# species-level CV weighted by relative species abundance
						cv_sp_abund <- dat_sub_sub %>% group_by(SPECIES) %>%
								summarise(cv_tmp=sd(abund, na.rm=T)/mean(abund, na.rm=T)) %>% ungroup() %>%
								summarise(cv_sp_abund=stats::weighted.mean(cv_tmp, sp_rel_abund$sp_rel_abund, na.rm=T))	
						
						#functional diversity
						func_div <- func_div_run(data=dat_sub_sub, fish.traits.data=fish.traits)
						
						sync_dat_out <- data.frame(
								ID=dat_sub_sub[1,"ID"],
								SITE_ID=dat_sub[1, "SITE_ID"],
								LAT=as.numeric(unique(dat_sub_sub$LAT)),
								LON=as.numeric(unique(dat_sub_sub$LON)),
								MPA=dat_sub_sub[1, "MPA"],
								MPA_STATUS=dat_sub_sub[1, "MPA_STATUS"],
								MPA_NAME=dat_sub_sub[1, "MPA_NAME"],
								IUCN_CAT=dat_sub_sub[1, "IUCN_CAT"],
								NO_TAKE=dat_sub_sub[1, "NO_TAKE"],
								STATUS_YR=dat_sub_sub[1, "STATUS_YR"],
								STATUS=dat_sub_sub[1, "STATUS"],
								SAMP_AGE=dat_sub_sub[1, "SAMP_AGE"],
								N_YEARS=nrow(unique(dat_sub_sub[, "YEAR"])),
								INI_YEAR=unique(min(dat_sub_sub[, "YEAR"])),
								END_YEAR=unique(max(dat_sub_sub[, "YEAR"])),
								N_REPS=dat_sub_sub[1, "N_REPS"],
								N_DATES=dat_sub_sub[1, "N_DATES"],
								SAMPLED_AREA=dat_sub_sub[1, "SAMPLED_AREA"],
								TRANSECT_SIZE=dat_sub_sub[1, "TRANSECT_SIZE"],
								ECOREGION=dat_sub_sub[1, "ECOREGION"],
								ECO_CODE=dat_sub_sub[1, "ECO_CODE"],
								PROVINCE=dat_sub_sub[1, "PROVINCE"],
								PROV_CODE=dat_sub_sub[1, "PROV_CODE"],
								REALM=dat_sub_sub[1, "REALM"],
								RLM_CODE=dat_sub_sub[1, "RLM_CODE"],
								AVG_DELTA_YEAR=mean(diff(unlist(unique(dat_sub_sub[, "YEAR"])))),
								ASYNC_GROSS_W=as.numeric(sync_gross_w)*-1,
								#ASYNC_GROSS_W_DET=sync_gross_w_det*-1,
								ASYNC_GROSS_W_DETREG=as.numeric(sync_gross_w_detreg)*-1,
								#ASYNC_GROSS=as.numeric(sync_gross)*-1,
								ASYNC_LOREAU=1-as.numeric(sync_loreau),
								#ASYNC_LOREAU_DET=1-sync_loreau_det,
								ASYNC_LOREAU_DETREG=1-sync_loreau_detreg,
								ASYNC_LOREAU_SQRT=1-sqrt(as.numeric(sync_loreau)),
								CV_TOT_ABUND=as.numeric(cv_tot_abund),
								CV_SP_ABUND=as.numeric(cv_sp_abund),
								N_SPEC=nrow(sp_rel_abund),
								SHANNON_DIV = shannon,
								SIMP_DIV = simp_div,
								SIMP_DOM = simp_dom,
								SIMP_EVENN = simp_evenn,
								MEAN_TOT_ABUND = as.vector(unlist(mean_sd_tot_abund$mabund)),
								MEAN_SP_ABUND = as.vector(unlist(mean_sd_sp_abund$m_sp_abund)),
								SD_TOT_ABUND = as.vector(unlist(mean_sd_tot_abund$sdabund)),
								SD_SP_ABUND = as.vector(unlist(mean_sd_sp_abund$sd_sp_abund)),
								func_div
						)
						
						
					}
					
					else {
						
						sync_dat_out <-  NULL
						
					}
					
					return(sync_dat_out)
					
				}
		
		assign(paste("id_", study_id, sep=""), value=out, pos=1, inherits=T)
		outputName=paste("id_", study_id, ".RData", sep="")
		outputPath=file.path('~/Data_ReefFishStability/alpha_stab_tmp', outputName)
		save(list=paste("id_", study_id, sep=""), file=outputPath)
		
	}
	
}


