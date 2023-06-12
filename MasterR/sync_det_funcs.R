# Functions to compute synchony measures (Gross and Loreau and De Mazancourt)
# using the three-term local variance # detrending approach
# (Leps et al Ecography 42: 1728–1741, 2019)
###############################################################################

# Function to calculate different versions of the eta synchrony measure
# (Gross et al. 2014)

cor_algo <- function(x, method, weighted = FALSE, rm.col=FALSE) {
	
	#### ---- Modifications to the original function introduced ---- ####
	#### ---- for the fish stability analysis ---------------------- #### 
	if(any(colnames(x)%in%c("SITE_ID","SPECIES"))) {
		col.rm <- which(colnames(x)%in%c("SITE_ID","SPECIES"))
		x <- x[,-col.rm]
	}
	
	if(rm.col) {
		if(any(colSums(x)==0)) x <- x[,-(which(colSums(x)==0))]
	}
	
	#### ------------------- STOP MODIFICATIONS -------------------- ####
	
	n <- ncol(x)
	cor_sp <- vector("numeric", length = n)
	for (i in 1:n) {
		xi <- as.numeric(x[, i])
		xnoti <- rowSums(x[,-i, drop = FALSE])
		cor_sp[i] <- method(xi, xnoti)
	}
	if (weighted) {
		rel_ab <- colSums(x) / sum(x)
		sync <- sum(cor_sp * rel_ab)
	} else {
		sync <- mean(cor_sp)
	}
	return(sync)
}

# Loreau and De Mazancoutr phi and detrended phi
phi <- function (x) {
	species.sd = apply(x, MARGIN = 2, FUN = sd)
	community.var = var(rowSums(x))
	return(community.var / sum(species.sd, na.rm = TRUE) ^ 2)
}

phi_t3 <- function (x) {
	species.var = apply(x, MARGIN = 2, FUN = var_t3)
	species.sd = sqrt(species.var)
	community.var = var_t3(rowSums(x))
	return(community.var / sum(species.sd, na.rm = TRUE) ^ 2)
}

# Variance, covariance, and correlations based on three term
# local quadrat variance (t3)
var_t3 <- function(x) {
	ws = 1
	n <- length(x) # sample size
	neg <- rep(c(1, -1), ws + 1)[1:(ws + 2)]
	multip  <- c(1, seq(ws + 1, 1)) * neg
	xsq <- vector("numeric", n - ws - 2)
	for (i in 1:(n - ws - 1)) {
		# calculate the "within brackets" component of the 3 term local variance for
		# xi, then the actual 3 term local quadratic variance of xi
		xisq <- x[i:(i + ws + 1)] * multip
		xsq[i] <- sum(xisq)^2 / 6
	}
	v <- sum(xsq) / (n - 2) # three term local quadratic variance of x
	return(v)
}

cov_t3 <- function(x, y) {
	cov <- (var_t3(x + y) - var_t3(x) - var_t3(y)) / 2
	return(cov)
}

cor_t3 <- function(x, y) {
	cor <- cov_t3(x, y) / sqrt(var_t3(x) * var_t3(y))
	return(cor)
}
