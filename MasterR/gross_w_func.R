# Function to compute the weighthed version of Gross synchrony across
# multiple columns (species, sites). The input file must be a data frame
# or matrix with columns as species (or any other aggregation) and rows
# as Years (or time more generally)
###############################################################################

gross_w_func <- function(data, weight=TRUE, det=FALSE, ...) {
	
	if(colnames(data)[1]!="YEAR") stop("The first column must be YEAR")
	
	#data <- as.data.frame(data)
	#sp.data <- as.data.frame(data[,2:ncol(data)])
	sp.data <- data %>% dplyr::select(-YEAR) %>%
			dplyr:: select(where(~ any(. != 0)))
	
	if(ncol(sp.data)>1) {
		
		year <- as.vector(unlist(data[,"YEAR"]))
		
		#if(any(colSums(sp.data)==0)) sp.data <- sp.data[,-(which(colSums(sp.data)==0))]
		
		rel_abund <- as.vector(unlist(apply(sp.data, 2,
								function(x) sum(x)/sum(sp.data))))
		
		# the first column of data must be "YEAR" for detrending
		if(det) {
			
			sp.data <- apply(sp.data, 2, function(x) resid(lm(x ~ year)))
			
		}
		
		tmp_cor_abund <- NULL
		
		for(i in 1:ncol(sp.data)) {
			
			focal_sp <- data.frame(sp=sp.data[,i])
			
			# exclude species with no temporal variation in abundance
			if(nrow(unique(focal_sp))>1) { 
				
				other_sp <- rowSums(as.data.frame(sp.data[,-i]))
				cor_dat <- cbind(focal_sp, other_sp)
				
				if(weight) {
					tmp_cor_abund <- c(tmp_cor_abund, rel_abund[i]*cor(focal_sp, other_sp))
				}
				
				else {
					tmp_cor_abund <- c(tmp_cor_abund, cor(focal_sp, other_sp))
				}
				
			}	
			
		}
		
		if(!all(is.na(tmp_cor_abund))) {
			
			if(weight) {
				data.frame(gross.sync=sum(tmp_cor_abund, na.rm=T))
			}
			
			else {
				data.frame(gross.sync=mean(tmp_cor_abund, na.rm=T))
			}
		}
		
		else {
			data.frame(gross.sync=NA)
		}
		
	}
	
	else {
		data.frame(gross.sync=NA)
	}

} 

