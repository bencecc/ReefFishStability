# Function to summarise results from null stability analysis by distance
###############################################################################

stats_nullstab_dist <- function(df, ...) {
	
	obs <- df[1,]
	perm <- df[2:nrow(df),]
	n_perm <- nrow(perm)
	
	mean_perm <- NULL
	pval_greater <- NULL
	pval_lower <- NULL
	ci_lw95 <- NULL
	ci_up95 <- NULL
	se <- NULL
	eff_size <- NULL
	perc_change <- NULL
	
	for (i in 6:ncol(perm)) {
		
		pval_greater <- rbind(pval_greater, sum(perm[,i]>=obs[,i], na.rm=T)/(n_perm+1))
		pval_lower <- rbind(pval_lower, sum(perm[,i]<=obs[,i], na.rm=T)/(n_perm+1))
		mean_perm <- rbind(mean_perm, mean(perm[,i], na.rm=T))
		ci_lw95 <- rbind(ci_lw95, quantile(perm[,i], 0.025, na.rm=T))
		ci_up95 <- rbind(ci_up95, quantile(perm[,i], 0.975, na.rm=T))
		se <- rbind(se, sd(perm[,i], na.rm=T))
		eff_size <- rbind(eff_size, (obs[,i]-mean(perm[,i], na.rm=T))/sd(perm[,i], na.rm=T))
		perc_change <- rbind(perc_change, (obs[,i]-mean(perm[,i], na.rm=T))/obs[,i]*100)
		
	}
	
	stab_res <- data.frame(
			df[1:nrow(pval_greater),2:5],
			METRIC=rownames(t(obs[1,6:ncol(obs)])),
			OBSERVED=as.vector(t(obs[1,6:ncol(obs)])),
			PERMUTED=mean_perm,
			CI_LW95=as.vector(ci_lw95), CI_UP95=as.vector(ci_up95),
			SE=as.vector(se),
			PVAL_GREATER=pval_greater, PVAL_LOWER=pval_lower,
			EFF_SIZE=eff_size,
			PERC_CHANGE=perc_change)
	
	return(stab_res)
}



