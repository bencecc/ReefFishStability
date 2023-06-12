# Annotate coeffients, significance and R2 on ggplot
# Use in combinaiton with functions sep_fit_plot and
# site_plot_func to fit lmer models to stability data
# seprately for MPAs and OAs or two-way models that include
# interactions with MPAs vs OAs, respectively
# Author: Lisandro Bendetti-Cecchi
###############################################################################

annotate_stats <- function(mod, cov, r2=r2, ...) {
	
	coef <- (summary(mod))$coefficients
	pval <- coef[which(rownames(coef)==cov), ncol(coef)]
	slope <- round(coef[which(rownames(coef)==cov), "Estimate"],3)
	
	if(pval<0.001) {
		star <- "***"
	} else if(pval<0.01){
		star <- "**"
	} else if(pval<0.05){
		star <- "*"
	} else{star <- " "}
	
	stat.coef <- sprintf("italic(beta) == %.3f~'%s'", slope, star)
	
	if(r2) {
		
		r2 <- round(r2_nakagawa(mod)[1][[1]], 2)
		stat.r2 <- sprintf("italic(R^2) == %.3f", r2)
		
	} else(stat.r2 <- "")
	
	list(stat.coef, stat.r2)
}


