# Assemble fish traits data generated on hpc and match names with fish_names
# 
###############################################################################

require(dplyr)
require(tidyr)
require(rfishbase)
require(stringr)
require(worrms)

#require(cati)
# load functions for trait analysis
source('~/workspace/BiodGlob/MPA_timeseries/Analysis/Traits/fish_traits_extract.R')
source('~/workspace/BiodGlob/MPA_timeseries/Analysis/Traits/fish_traits_utils.R')

# load data
setwd('~/Lavori/MPA_timeseries/')
load("fish_names.RData")
load("fish.traits.names.RData")
load("fish_traits_res.RData")
load("ccope_names.RData")
load("cpalos_names.RData")
# NOTE: cpalos_names already includes ccope_names

#out <- NULL
#for(i in 1:nrow(cpalos_names)) {
#traits <- fish_traits(cpalos_names$SPECIES[i])
#out <- rbind(out, traits)
#}
#
#cpalos_ccope_traits <- out
#
#fish_names <- rbind(fish_names, cpalos_names) %>% distinct(SPECIES)
#fish.traits.names <- rbind(fish.traits.names, cpalos_ccope_traits$SPECIES) %>% distinct(SPECIES)
##save(fish_names, file="fish_names.RData")
##save(fish.traits.names, file="fish.traits.names.RData")


# rename some taxa to family (to be named back to original for
# subsequent analyses
too.bad <- unlist(lapply(fish_names$SPECIES, function(x) {
					str.tmp <- str_locate(x, "id")
					if(any(!is.na(str.tmp))) {out <- TRUE}
					else {
						out=FALSE
					}
				}))

fam.tmp <- fish_names[too.bad,1]		


fish.traits.names <- fish_names %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="baitfish, unidentified", "baitfish")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Cottid sp", "Cottidae")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Bothid sp", "Botthidae")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Atherinid spp.", "Atherinidae")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Gobiid spp.", "Gobidae")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Blenniid spp.", "Blennidae")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Tripterygiid spp.", "Tripterygiidae")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Platycephalid spp.", "Platycephalidae")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Monacanthid spp.", "Monacanthidae")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Clinid spp.", "Clinidae")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Ophichthid spp.", "Ophichthidae")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Clupeid spp.", "Clupeidae")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Scorpaenid spp.", "Scorpaenidae")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Acanthurid spp.", "Acanthuridae")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Pomacentrid spp.", "Pomacentridae")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Mugilid spp.", "Mugilidae")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Creediid spp.", "Creediidae")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Callionymid spp.", "Callionymidae")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Carangid spp.", "Carangidae")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Belonid spp.", "Belonidae")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Morid spp.", "Moridae")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Pomacanthid spp.", "Pomacanthidae")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Muraenid spp.", "Muraenidae")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Lutjanid spp.", "Lutjanidae")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Syngnathid spp.", "Syngnathidae")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Sphyraenid spp.", "Sphyraenidae")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Pseudochromid spp.", "Pseudochromidae")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Sphyraenid spp.", "Sphyraenidae")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="baitfish, unidentified", "baitfish"))

#save(fish.traits.names, file="fish.traits.names.RData")

# revert species names in fish_traits_res to original names
traits.tmp <- fish_traits_res %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="baitfish", "baitfish, unidentified")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Cottidae", "Cottid sp")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Botthidae", "Bothid sp")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Atherinidae", "Atherinid spp.")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Gobidae", "Gobiid spp.")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Blennidae", "Blenniid spp.")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Tripterygiidae", "Tripterygiid spp.")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Platycephalidae", "Platycephalid spp.")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Monacanthidae", "Monacanthid spp.")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Clinidae", "Clinid spp.")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Ophichthidae", "Ophichthid spp.")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Clupeidae", "Clupeid spp.")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Scorpaenidae", "Scorpaenid spp.")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES== "Acanthuridae", "Acanthurid spp.")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Pomacentridae", "Pomacentrid spp.")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Mugilidae", "Mugilid spp.")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Creediidae", "Creediid spp.")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Callionymidae", "Callionymid spp.")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Carangidae", "Carangid spp.")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Belonidae", "Belonid spp.")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Moridae", "Morid spp.")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Pomacanthidae", "Pomacanthid spp.")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Muraenidae", "Muraenid spp.")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Lutjanidae", "Lutjanid spp.")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Syngnathidae", "Syngnathid spp.")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Sphyraenidae", "Sphyraenid spp.")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Pseudochromidae", "Pseudochromid spp.")) %>%
		mutate(SPECIES=replace(SPECIES, SPECIES=="Sphyraenidae", "Sphyraenid spp.")) %>%
		distinct(SPECIES, .keep_all=TRUE)

#miss.sp <- fish_names %>% filter(!SPECIES%in%traits.tmp$SPECIES)
#which(fish_names$SPECIES=="Gobidae")

fish.traits.dat.tmp <- rbind(fish_traits_res, traits.tmp, cpalos_ccope_traits) %>% distinct(SPECIES, .keep_all=TRUE)
#miss.sp <- fish.traits.dat %>% filter(!SPECIES%in%fish_names$SPECIES)

load("fish_sti_res.RData")
load("medes_sti_res.RData")
load("cpalos_sti_res.RData") #also includes ccope sti res

sti.tmp <- rbind(fish_sti_res, medes_sti_res, cpalos_sti_res) %>% distinct(SPECIES, .keep_all=TRUE)

fish.traits.dat <- fish.traits.dat.tmp %>% left_join(sti.tmp, by="SPECIES")
miss.sp <- fish_names %>% filter(!SPECIES%in%fish.traits.dat$SPECIES)

# some FoodTroph values are wrong
ft.wr <- which(fish.traits.dat$FoodTroph>5)
# so far these are for Acanthurus nigrofuscus and Plectropomus leopardus;
# values are replaced manually with data provided from fihsbase
# (mean and se): 2.0 (0.0) and 4.42 (0.71), respectively 
fish.traits.dat[13, 3] <- 2.0
fish.traits.dat[13, 4] <- 0.0
fish.traits.dat[1036, 3] <- 4.42
fish.traits.dat[1036, 4] <- 0.71

setwd('~/Lavori/MPA_timeseries/')
#save(fish.traits.dat, file="fish.traits.dat.RData")

# fish.traits.dat is the master database originating from fishbase
# other fish.traits is modified with selected traits for the analysis

# MODIFICATION 1: FeedingType
# combine fish.traits.dat with cluster feeding mode from Parraviciniet al. PlOS ONE
# and associate FeedingType categories with clusters on the basis of the most frequent
# matching
load("~/Lavori/MPA_timeseries/trophic_guilds_medoid.RData")

ft.len <- fish.traits.dat %>% group_by(FeedingType) %>% summarise(n=n())
mod.pred <- mod_predators %>% rename(SPECIES=species)
feed.cat <- fish.traits.dat %>% left_join(mod.pred, by="SPECIES") %>%
		dplyr::select(SPECIES, FeedingType, cluster) %>% drop_na() %>%
		group_by(FeedingType) %>%
		summarise(clust.sel=as.numeric(names(table(cluster)[which.max(table(cluster))])))
# the follwing feeding categories are defined from this analysis:
# browsing on substrate -> microphages
# selective plankton feeding + filtering plankton + variable -> planktivores
# grazing on aquatic plants -> grazers
# hunting macrofauna (predator) + picking parasites off a host (cleaner) -> carnivores
# others -> NA

# Select five traits:
# Maximum Length - continuous (body size)
# Trophic position - continuous (trophic)
# Thermal affinity - continuous (physiological)
# Vulnerability - continuous (life-history); integrates many life-history aspects and may be redundant with other traits
# Gregariousness - ordered factor (behaviour)
# Water position - ordered factor (behaviour)
# NOT USED, but potentially useful
# Diet; carnivores are not distinguished among corallivore, pishivore etc
# Hab.comp - continuous (habitat use) - Many NAs

load("fish.traits.dat.RData")

fish.traits <- fish.traits.dat %>%
		mutate(FeedingType=replace(FeedingType, FeedingType=="browsing on substrate", "microphages")) %>%
		mutate(FeedingType=replace(FeedingType, FeedingType%in%c("selective plankton feeding",
								"filtering plankton", "variable"), "planktivores")) %>%
		mutate(FeedingType=replace(FeedingType, FeedingType=="grazing on aquatic plants", "grazers")) %>%
		mutate(FeedingType=replace(FeedingType, FeedingType%in%c("hunting macrofauna (predator)",
								"picking parasites off a host (cleaner)"), "carnivores")) %>%
		mutate(FeedingType=replace(FeedingType, FeedingType=="other", NA)) %>%
		mutate(Greg=cut(greg.score, breaks=c(1, 1.6, 2.5, 3), labels=c("solitary","pairing","schooling"),
						ordered_result=TRUE)) %>%
		mutate(DemersPelag=replace(DemersPelag, DemersPelag=="reef-associated", "benthic")) %>%
		mutate(DemersPelag=replace(DemersPelag, DemersPelag=="bathydemersal", "demersal")) %>%
		mutate(DemersPelag=replace(DemersPelag, DemersPelag=="benthopelagic", "palagic.reef")) %>%
		mutate(DemersPelag=replace(DemersPelag, DemersPelag%in%c("pelagic-neritic","pelagic-oceanic","bethypelagic"),
						"palagic")) %>%
		mutate(Water.pos=factor(DemersPelag, ordered=TRUE, levels=c("benthic","demersal","pelagic.reef","pelagic"))) %>%
		mutate(Thermal.affinity=as.numeric(sst_mean_TN)) %>%
		mutate(Diet=factor(FeedingType)) %>%
		rename(Troph.pos=FoodTroph) %>%
		rename(Hab.comp=habitat.compl) %>%
		dplyr::select(SPECIES, FeedingType, Length, Troph.pos, Thermal.affinity, Vulnerability, Greg, Water.pos)
		
#save(fish.traits, file="fish.traits.RData")		
		













