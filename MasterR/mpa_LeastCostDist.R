# Function to aggregate sites as a function of least-cost distance from MPA sites
# for subsequent analysis of meta-stability; the value is a distance matrix with
# columns names as SITE_ID filtered by MPA sites from the input datatset and
# rows are all sites at different distances from the column (MPA) site. Distances
# are computed using package gdistance after applying a Mollweide projection
# as in function site_by_dist. Arguments are: the input dataset (df), a
# projected map of the world (rast object from terra package), and its 
# resolution (res)
################################################################################

mpa_LeastCostDist <- function(df, proj.world.rast, res=0.3, onlyMPA=FALSE, ...)

{
	
	dat <- df %>%
			group_by(ECOREGION, PROVINCE, REALM, SITE_ID, LAT, LON, MPA, YEAR, SPECIES) %>%
			summarise(abund=mean(abund), .groups="drop") %>% filter(abund>0) %>%
			distinct(LAT,LON,.keep_all=TRUE) %>%
			ungroup()
	
	sf.dat <- st_as_sf(dat, coords = c("LON", "LAT"), 
			crs = "+proj=longlat +datum=WGS84")
	
	# projection
	my.proj <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84"
	# rasterize and project data
	r.dat.templ <- rast(vect(sf.dat), res=res, crs="+proj=longlat +datum=WGS84")
	r.dat <- rasterize(vect(sf.dat), r.dat.templ, "MPA")
	r.dat <- terra::project(r.dat, my.proj)
	
	# crop world to 2x data extent
	r.land <- terra::crop(proj.world.rast, 2*ext(r.dat))
	# change raster values to 0 (land) and 1 (sea)
	r.land <- is.na(r.land)
	#change from terra to raster for use with gdistance
	r.land <- raster(r.land)
	
	# transition layer
	tr <- transition(r.land, transitionFunction=function(x) mean(x), directions = 8)
	tr <- geoCorrection(tr, type="c")
	
	# project data points before computing distances
	sf.dat <- st_transform(sf.dat, my.proj)
	
	dist.mat <- costDistance(tr, st_coordinates(sf.dat),
			st_coordinates(sf.dat))/1000
	colnames(dist.mat) <- rownames(dist.mat) <- dat$SITE_ID
	
	# select columns that are MPA sites
	if(onlyMPA) {
		mpa.id <- dat %>% filter(MPA=="Yes") %>% distinct(SITE_ID)
		out <- dist.mat[,colnames(dist.mat)%in%as.vector(mpa.id$SITE_ID)]
	} else {
		out <- dist.mat
		
	}
	
	return(as.data.frame(out))
	
}


