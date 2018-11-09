# EMODnet-thermal-traits
Deriving, summarising and visualising thermal affinities for European marine species
# Overview
The aim is to derive thermal affinities for all European marine species, by matching occurrence records from [OBIS](http://www.iobis.org) to gridded temperature products. These species-level thermal affinities are then used to produce assemblage-level averages on a 0.5º grid covering European seas, separately for major functional groups (e.g. benthos, zooplankton, fish). Finally these gridded assemblage-level averages are compared to current and projected future sea temperatures to identify areas of high climate vulnerability for each functional group.
# R session info, required packages
```R
# Load required packages
library(tidyverse)
library(robis)
library(worrms)
```
# European species list
The European species list was derived on the basis of occurrence records in [OBIS](iobis.org). Specifically, we used the `checklist` function in the [`robis`](https://github.com/iobis/robis) package to produce a summary of all species with occurrence records since 2000 within the European region (defined by the rectangle 45ºW–70ºE, 26ºN–90ºN, see http://www.eurobis.org/geo_scope), and limiting the results to taxa identified at the species level, as follows:
```R
eur_spp <- checklist(
  geometry = "POLYGON((-45 26,-45 90,70 90,70 26,-45 26))", year = 2000:2018,
  fields = c("species", "records", "worms_id", "obisid", "rank_name"))
```
This is a large query, so the results are provided here as a csv file that can be read in directly:
```R
eur_spp <- read_csv("eur_spp.csv")
# get rid of excess column and restrict to species
eur_spp <- eur_spp %>% dplyr::select(-X1) %>%
	filter(rank_name == "Species")
# needs some grouping
eur_spp <- eur_spp %>% group_by(worms_id) %>%
	summarise(species = first(species), records = sum(records))
length(unique(eur_spp$worms_id))
#' 11184 species
```
A second species list was derived directly from the OBIS database by Pieter Provoost, this is for the same region but with no time limit (i.e. records from all years are considered). It was generated using the SQL query
```SQL
select p.tname as scientificname, p.worms_id, rt.rank_name as rank, p.species, st.worms_id as worms_id_species, count(*)
from explore.points p
left join explore.taxon st on st.id = p.species_id
left join explore.taxon rt on rt.id = p.valid_id
where ST_Within(geom, 'POLYGON((-45 26,-45 90,70 90,70 26,-45 26))'::geography::geometry)
group by p.tname, p.worms_id, rt.rank_name, p.species, st.worms_id
order by count(*) desc
```
This is also available as a csv file:
```R
eur_spp_full <- read_csv2("checklist_europe.csv")
# restrict to species or below only
eur_spp_full <- eur_spp_full %>%
  filter(!is.na(worms_id_species))
# group to species level
eur_spp_full <- eur_spp_full %>% group_by(worms_id_species) %>%
	summarise(species = first(species), records = sum(count))
length(unique(eur_spp_full$worms_id_species))
#' 26569
```
The two species lists are tidied up a bit and combined:
```R
eur_spp <- eur_spp %>% rename(AphiaID = worms_id)
eur_spp_full <- eur_spp_full %>% rename(AphiaID = worms_id_species)

# check (lack of) overlap
sum(!(eur_spp$AphiaID %in% eur_spp_full$AphiaID))
#' 1 species

#' add `since2000` as variable to `eur_spp_full` to identify those species observed since 2000
eur_spp_full <- eur_spp_full %>% mutate(since2000 = AphiaID %in% eur_spp$AphiaID)
sum(eur_spp_full$since2000)
#' 11183
```
# Adding functional groups
We then obtained functional groups for these species using attributes from [WoRMS](http://www.marinespecies.org) using the [`worrms`](https://github.com/ropensci/worrms) package and a function written specifically for this purpose available from https://github.com/tomjwebb/WoRMS-functional-groups. NB - this will take up to a couple of hours to run for all species:
```R
species_attr <- eur_spp_full %>%
  group_by(AphiaID) %>%
  do(get_worms_fgrp(AphiaID = .$AphiaID))
```
Count the number of species within each function group for each life stage:
```R
count(species_attr, adult)
# adult     n
#<fctr> <int>
# 1          algae   386
# 2        benthos 16054
# 3          birds   163
# 4     epibenthos     8
# 5           fish  2539
# 6   hyperbenthos     5
# 7          macro     5
# 8     macroalgae   985
# 9   macrobenthos    38
# 10        mammals    53
# 11           meso     1
# 12         nekton   307
# 13        neuston     4
# 14 not applicable    18
# 15  phytoplankton  1489
# 16       plankton     3
# 17       reptiles    13
# 18    zooplankton  1414
# 19           <NA>  3084
```
Create an unknown Fg for adult species with FG
```R
species_attr$adult <- factor(species_attr$adult, levels = levels(addNA(species_attr$adult)),
  labels = c(levels(species_attr$adult), "unknown"), exclude = NULL)
```
Find the number of OBIS records for each functional gorup, especially unknown 
```R
# join FG data to dat with number of obis records
species_attr <- left_join(species_attr, eur_spp_full, by = "AphiaID")
# calculate summary of obis records
record_group_sum <- species_attr %>% group_by(adult) %>% summarise(min.recs = min(records),
                                                             max.recs = max(records),
                                                             mean.recs = mean(records),
                                                             total.recs = sum(records),
                                                             num.species = n())
record_sum <- species_attr %>% summarise(min.recs = min(records),
                                         max.recs = max(records),
                                         mean.recs = mean(records),
                                         total.recs = sum(records),
                                         num.species = n())

species_attr$adult[species_attr$adult == "not applicable"]<- "unknown"
```
Tidy up the dataset a bit
```R
new_checklist_species_fgrps <- species_attr %>% dplyr::select(AphiaID, adult) %>%
  rename(functional_group = adult)
```
Subset the checklist to just include species in OBIS records - this should be all species but a few do not return OBIS records for some reason, these need to be excluded here because the temperature matching functions at present assume that all species are present in OBIS.
```
new_checklist_obis_test <- new_checklist_species_fgrps %>% 
  dplyr::select(AphiaID) %>% 
  arrange(AphiaID) %>% group_by(AphiaID) %>%
  do(checklist(aphiaid = .$AphiaID)) %>% summarise(records = sum(records))
# 26338 species
```
Add functional group data to this set of species with records available in OBIS:
```R
new_checklist_species_fgrps <- left_join(new_checklist_obis_test, new_checklist_species_fgrps, by = "AphiaID")
```
Subset by functional groups to enable each to be processed separately below:
```R
levels(new_checklist_species_fgrps$functional_group)
# "algae"          "benthos"        "birds"          "epibenthos"     "fish"           "hyperbenthos"   "macro"
# "macroalgae"     "macrobenthos"  "mammals"        "meso"           "nekton"         "neuston"     "phytoplankton"
# "plankton"       "reptiles"       "unknown"   "zooplankton" 

# No. species listed before and after 
algae_sp <- new_checklist_species_fgrps %>% filter(functional_group == "algae") # 386 species / 384
benthos_sp <- new_checklist_species_fgrps %>% filter(functional_group == "benthos") # 16054 species/ 15965
birds_sp <- new_checklist_species_fgrps %>% filter(functional_group == "birds") # 163 / 163
epibenthos_sp <- new_checklist_species_fgrps %>% filter(functional_group == "epibenthos") #8 / 8
fish_sp <- new_checklist_species_fgrps %>% filter(functional_group == "fish") #2539 / 2535
hyperbenthos_sp <- new_checklist_species_fgrps %>% filter(functional_group == "hyperbenthos") #5/5
macro_sp <- new_checklist_species_fgrps %>% filter(functional_group == "macro") #5 / 5
macroalgae_sp <- new_checklist_species_fgrps %>% filter(functional_group == "macroalgae") #985 / 979
macrobenthos_sp <- new_checklist_species_fgrps %>% filter(functional_group == "macrobenthos") #38 / 38
mammals_sp <- new_checklist_species_fgrps %>% filter(functional_group == "mammals") #53 / 53
meso_sp <- new_checklist_species_fgrps %>% filter(functional_group == "meso") #1 / 1
nekton_sp <- new_checklist_species_fgrps %>% filter(functional_group == "nekton") #307  / 303
neuston_sp <- new_checklist_species_fgrps %>% filter(functional_group == "neuston") #4 / 4
phytoplankton_sp <- new_checklist_species_fgrps %>% filter(functional_group == "phytoplankton") #1489 / 1486
plankton_sp <- new_checklist_species_fgrps %>% filter(functional_group == "plankton") #3 / 3
reptiles_sp <- new_checklist_species_fgrps %>% filter(functional_group == "reptiles") #13 / 13
unknown_sp <- new_checklist_species_fgrps %>% filter(functional_group == "unknown") #3084 / 2994
zooplankton_sp <- new_checklist_species_fgrps %>% filter(functional_group == "zooplankton") #1414 / 1399
```

# Temperature Matching
Thermal affinities were derived from all European marine species using a suite of functions described and justified in full here: https://github.com/tomjwebb/marine_thermal_traits. Briefly, for each species, all global occurrences were extracted from OBIS, these were then matched in space (latitude, longitude) to a global temperature climatology (SST and SBT from [Bio-ORACLE](http://www.bio-oracle.org) SST and SBT, obtained using the [`sdmpredictors`](https://github.com/lifewatch/sdmpredictors) package), and in space, time (month, year), and sample depth to a global gridded temperature product (Institute of Atmospheric Physics, [IAP gridded temperature](http://159.226.119.60/cheng/)). Temperatures were then summarised to provide a range of thermal affinity measures for each species (including mean, min, max, 5th and 95% quantiles, for each temperature measure). The full set of functions is provided and documented [elsewhere](https://github.com/tomjwebb/marine_thermal_traits) but specifically for this implementation, we processed species from each functional group in turn, and we developed a parallelised workflow to allow for more efficient application of the general functions over large species lists. This was run on a local server using the approach documented here:
```R
# load packages for parallel processing
library(parallel)
library(multidplyr)
```
# Cluster management
```R
# create cluster
cluster20 <- create_cluster(20)

# register functions and data with cluster20
cluster_assign_value(cluster20, 'get_temp_summ_by_sp', get_temp_summ_by_sp)
cluster_assign_value(cluster20, 'get_bio_oracle_t', get_bio_oracle_t)
cluster_assign_value(cluster20,'get_ersem_gridded_t', get_ersem_gridded_t)
cluster_assign_value(cluster20, 'get_iap_gridded_t', get_iap_gridded_t)
cluster_assign_value(cluster20, 'get_obis_recs', get_obis_recs)
cluster_assign_value(cluster20, 'save_full_recs', save_full_recs)
cluster_assign_value(cluster20, 't_summary', t_summary)
cluster_assign_value(cluster20, 'use_defaults', use_defaults)
cluster_assign_value(cluster20, 'use_ersem', use_ersem)

# attach the required packages to each node
# create a vector of packages to load into each cluster
packages <- c("tidyverse", "lubridate", "robis", "raster", "sdmpredictors", "naniar", "ncdf4", "multidplyr", "parallel")
cluster_library(cluster20, packages)
```
Then run for each functional group in turn. Algae:
```R
# algae_sp t_matching
# add algae_sp data to the cluster
cluster_assign_value(cluster20, 'algae_sp', algae_sp)
# partition the algae_sp data
algae_sp <- algae_sp %>% partition(AphiaID, cluster = cluster20)
algae_sp
# run function on test data # 
algae_sp_t_matching <- algae_sp %>% do(get_temp_summ_by_sp(sp_id = .$AphiaID)) %>% collect()

write_csv(algae_sp_t_matching, "~/EMODnet_Thermal_Traits/t_matching/algae_sp_t_matching.csv")
```
Benthos:
```R
# benthos_sp t_matching
# add benthos_sp data to the cluster
cluster_assign_value(cluster20, 'benthos_sp', benthos_sp)
# partition the benthos_sp data
benthos_sp <- benthos_sp %>% partition(AphiaID, cluster = cluster20)
benthos_sp
# run function on test data # 
benthos_sp_t_matching <- benthos_sp %>% do(get_temp_summ_by_sp(sp_id = .$AphiaID)) %>% collect()

write_csv(benthos_sp_t_matching, "~/EMODnet_Thermal_Traits/t_matching/benthos_sp_t_matching.csv")
```
Birds:
```R
# birds_sp t_matching
# add birds_sp data to the cluster
cluster_assign_value(cluster20, 'birds_sp', birds_sp)
# partition the birds_sp data
birds_sp <- birds_sp %>% partition(AphiaID, cluster = cluster20)
birds_sp
# run function on test data # 
birds_sp_t_matching <- birds_sp %>% do(get_temp_summ_by_sp(sp_id = .$AphiaID)) %>% collect()

write_csv(birds_sp_t_matching, "~/EMODnet_Thermal_Traits/t_matching/birds_sp_t_matching.csv")
```
Epibenthos:
```R
# epibenthos_sp t_matching
# add epibenthos_sp data to the cluster
cluster_assign_value(cluster20, 'epibenthos_sp', epibenthos_sp)
# partition the epibenthos_sp data
epibenthos_sp <- epibenthos_sp %>% partition(AphiaID, cluster = cluster20)
epibenthos_sp
# run function on test data # 
epibenthos_sp_t_matching <- epibenthos_sp %>% do(get_temp_summ_by_sp(sp_id = .$AphiaID)) %>% collect()

write_csv(epibenthos_sp_t_matching, "~/EMODnet_Thermal_Traits/t_matching/epibenthos_sp_t_matching.csv")
```
Fish:
```R
# fish_sp t_matching
# add fish_sp data to the cluster
cluster_assign_value(cluster20, 'fish_sp', fish_sp)
# partition the fish_sp data
fish_sp <- fish_sp %>% partition(AphiaID, cluster = cluster20)
fish_sp
# run function on test data # 
fish_sp_t_matching <- fish_sp %>% do(get_temp_summ_by_sp(sp_id = .$AphiaID)) %>% collect()

write_csv(fish_sp_t_matching, "~/EMODnet_Thermal_Traits/t_matching/fish_sp_t_matching.csv")
```
Hyperbenthos:
```R
# hyperbenthos_sp t_matching
# add hyperbenthos_sp data to the cluster
cluster_assign_value(cluster20, 'hyperbenthos_sp', hyperbenthos_sp)
# partition the hyperbenthos_sp data
hyperbenthos_sp <- hyperbenthos_sp %>% partition(AphiaID, cluster = cluster20)
hyperbenthos_sp
# run function on test data # 
hyperbenthos_sp_t_matching <- hyperbenthos_sp %>% do(get_temp_summ_by_sp(sp_id = .$AphiaID)) %>% collect()

write_csv(hyperbenthos_sp_t_matching, "~/EMODnet_Thermal_Traits/t_matching/hyperbenthos_sp_t_matching.csv")
```
Macro:
```
# macro_sp t_matching
# add macro_sp data to the cluster
cluster_assign_value(cluster20, 'macro_sp', macro_sp)
# partition the macro_sp data
macro_sp <- macro_sp %>% partition(AphiaID, cluster = cluster20)
macro_sp
# run function on test data # 
macro_sp_t_matching <- macro_sp %>% do(get_temp_summ_by_sp(sp_id = .$AphiaID)) %>% collect()

write_csv(macro_sp_t_matching, "~/EMODnet_Thermal_Traits/t_matching/macro_sp_t_matching.csv")
```
Macroalgae:
```
# macroalgae_sp t_matching
# add macroalgae_sp data to the cluster
cluster_assign_value(cluster20, 'macroalgae_sp', macroalgae_sp)
# partition the macroalgae_sp data
macroalgae_sp <- macroalgae_sp %>% partition(AphiaID, cluster = cluster20)
macroalgae_sp
# run function on test data # 
macroalgae_sp_t_matching <- macroalgae_sp %>% do(get_temp_summ_by_sp(sp_id = .$AphiaID)) %>% collect()

write_csv(macroalgae_sp_t_matching, "~/EMODnet_Thermal_Traits/t_matching/macroalgae_sp_t_matching.csv")
```
Macrobenthos
```R
# macrobenthos_sp t_matching
# add macrobenthos_sp data to the cluster
cluster_assign_value(cluster20, 'macrobenthos_sp', macrobenthos_sp)
# partition the macrobenthos_sp data
macrobenthos_sp <- macrobenthos_sp %>% partition(AphiaID, cluster = cluster20)
macrobenthos_sp
# run function on test data # 
macrobenthos_sp_t_matching <- macrobenthos_sp %>% do(get_temp_summ_by_sp(sp_id = .$AphiaID)) %>% collect()

write_csv(macrobenthos_sp_t_matching, "~/EMODnet_Thermal_Traits/t_matching/macrobenthos_sp_t_matching.csv")
```
Mammals:
```R
# mammals_sp t_matching
# add mammals_sp data to the cluster
cluster_assign_value(cluster20, 'mammals_sp', mammals_sp)
# partition the mammals_sp data
mammals_sp <- mammals_sp %>% partition(AphiaID, cluster = cluster20)
mammals_sp
# run function on test data # 
mammals_sp_t_matching <- mammals_sp %>% do(get_temp_summ_by_sp(sp_id = .$AphiaID)) %>% collect()

write_csv(mammals_sp_t_matching, "~/EMODnet_Thermal_Traits/t_matching/mammals_sp_t_matching.csv")
```
Meso:
```R
# meso_sp t_matching
# add meso_sp data to the cluster
cluster_assign_value(cluster20, 'meso_sp', meso_sp)
# partition the meso_sp data
meso_sp <- meso_sp %>% partition(AphiaID, cluster = cluster20)
meso_sp
# run function on test data # 
meso_sp_t_matching <- meso_sp %>% do(get_temp_summ_by_sp(sp_id = .$AphiaID)) %>% collect()

write_csv(meso_sp_t_matching, "~/EMODnet_Thermal_Traits/t_matching/meso_sp_t_matching.csv")
```
Nekton:
```R
# nekton_sp t_matching
# add nekton_sp data to the cluster
cluster_assign_value(cluster20, 'nekton_sp', nekton_sp)
# partition the nekton_sp data
nekton_sp <- nekton_sp %>% partition(AphiaID, cluster = cluster20)
nekton_sp
# run function on test data # 
nekton_sp_t_matching <- nekton_sp %>% do(get_temp_summ_by_sp(sp_id = .$AphiaID)) %>% collect()

write_csv(nekton_sp_t_matching, "~/EMODnet_Thermal_Traits/t_matching/nekton_sp_t_matching.csv")
```
Neuston:
```R
# neuston_sp t_matching
# add neuston_sp data to the cluster
cluster_assign_value(cluster20, 'neuston_sp', neuston_sp)
# partition the neuston_sp data
neuston_sp <- neuston_sp %>% partition(AphiaID, cluster = cluster20)
neuston_sp
# run function on test data # 
neuston_sp_t_matching <- neuston_sp %>% do(get_temp_summ_by_sp(sp_id = .$AphiaID)) %>% collect()

write_csv(neuston_sp_t_matching, "~/EMODnet_Thermal_Traits/t_matching/neuston_sp_t_matching.csv")
```
Phytoplankton
```R
# phytoplankton_sp t_matching
# add phytoplankton_sp data to the cluster
cluster_assign_value(cluster20, 'phytoplankton_sp', phytoplankton_sp)
# partition the phytoplankton_sp data
phytoplankton_sp <- phytoplankton_sp %>% partition(AphiaID, cluster = cluster20)
phytoplankton_sp
# run function on test data # 
phytoplankton_sp_t_matching <- phytoplankton_sp %>% do(get_temp_summ_by_sp(sp_id = .$AphiaID)) %>% collect()

write_csv(phytoplankton_sp_t_matching, "~/EMODnet_Thermal_Traits/t_matching/phytoplankton_sp_t_matching.csv")
```
Plankton:
```R
# plankton_sp t_matching
# add plankton_sp data to the cluster
cluster_assign_value(cluster20, 'plankton_sp', plankton_sp)
# partition the plankton_sp data
plankton_sp <- plankton_sp %>% partition(AphiaID, cluster = cluster20)
plankton_sp
# run function on test data # 
plankton_sp_t_matching <- plankton_sp %>% do(get_temp_summ_by_sp(sp_id = .$AphiaID)) %>% collect()

write_csv(plankton_sp_t_matching, "~/EMODnet_Thermal_Traits/t_matching/plankton_sp_t_matching.csv")
```
Reptiles:
```R
# reptiles_sp t_matching
# add reptiles_sp data to the cluster
cluster_assign_value(cluster20, 'reptiles_sp', reptiles_sp)
# partition the reptiles_sp data
reptiles_sp <- reptiles_sp %>% partition(AphiaID, cluster = cluster20)
reptiles_sp
# run function on test data # 
reptiles_sp_t_matching <- reptiles_sp %>% do(get_temp_summ_by_sp(sp_id = .$AphiaID)) %>% collect()

write_csv(reptiles_sp_t_matching, "~/EMODnet_Thermal_Traits/t_matching/reptiles_sp_t_matching.csv")
```
Unknown:
```R
# unknown_sp t_matching
# add unknown_sp data to the cluster
cluster_assign_value(cluster20, 'unknown_sp', unknown_sp)
# partition the unknown_sp data
unknown_sp <- unknown_sp %>% partition(AphiaID, cluster = cluster20)
unknown_sp
# run function on test data # 
unknown_sp_t_matching <- unknown_sp %>% do(get_temp_summ_by_sp(sp_id = .$AphiaID)) %>% collect()

write_csv(unknown_sp_t_matching, "~/EMODnet_Thermal_Traits/t_matching/unknown_sp_t_matching.csv")
```
Zooplankton:
```R
# zooplankton_sp t_matching
# add zooplankton_sp data to the cluster
cluster_assign_value(cluster20, 'zooplankton_sp', zooplankton_sp)
# partition the zooplankton_sp data
zooplankton_sp <- zooplankton_sp %>% partition(AphiaID, cluster = cluster20)
zooplankton_sp
# run function on test data # 
zooplankton_sp_t_matching <- zooplankton_sp %>% do(get_temp_summ_by_sp(sp_id = .$AphiaID)) %>% collect()

write_csv(zooplankton_sp_t_matching, "~/EMODnet_Thermal_Traits/t_matching/zooplankton_sp_t_matching.csv")
```
Bind all groups together:
```R
all_species_t_matching <- bind_rows(algae_sp_t_matching, benthos_sp_t_matching, birds_sp_t_matching,
                                    epibenthos_sp_t_matching, fish_sp_t_matching, hyperbenthos_sp_t_matching,
                                    macro_sp_t_matching, macroalgae_sp_t_matching, macrobenthos_sp_t_matching,
                                    mammals_sp_t_matching, meso_sp_t_matching, nekton_sp_t_matching,
                                    neuston_sp_t_matching, phytoplankton_sp_t_matching, plankton_sp_t_matching,
                                    reptiles_sp_t_matching, unknown_sp_t_matching, zooplankton_sp_t_matching)
```
Join to Checked_species_functinal_groups:
```R
all_species_fgrps_t_matched <- left_join(new_checklist_species_fgrps, all_species_t_matching, by = "AphiaID")
```
Export full dataset if required:
```R
write_csv(all_species_fgrps_t_matched, "~/all_species_fgrps_t_matched.csv")
```
# Assemblage-level thermal summaries and maps
The next step is to derive gridded products for European seas. To do this, we set up a 0.5º grid, and for each cell extracted a list of species with recorded occurrences since 2000. These species were then matched to our thermal affinity database, allowing summaries of 'community thermal affinity' for each functional group in each grid cell. These summaries include, for example, mean thermal affinity for each temperature measure, across all species within a group occurring within a grid cell, weighted by the number of occurrences for each species within that cell.
```R
# load packages
library(raster)
library(mregions)
library(tidyverse)
library(sf)
library(robis)
```
Get a sensible base map of European seas, using FAO fishing regions
```R
fao <- mr_shp(key = "MarineRegions:fao")
# subset to those covering Europe
fao_eur <- subset(fao,
  name %in% c("Atlantic, Northeast", "Mediterranean and Black Sea"))
fao_eur <- aggregate(fao_eur, dissolve = TRUE)

# plot if required
plot(fao_eur)
```
Set up a raster based on EMODNet 'europe' extent (45W-70E, 26N-90N). For now set resolution = 0.5; could refine to 0.1 if needed
```R
eur_r <- raster(extent(c(-45, 70, 26, 90)),
  resolution = 0.5, crs = crs("+proj=longlat +datum=WGS84"))
# need to fill with values, but these are arbitrary
values(eur_r) <- 0
```
Mask to fao areas and convert to spatial polygons
```R
eur_r_mask <- mask(eur_r, fao_eur)
eur_sp_poly <- rasterToPolygons(eur_r_mask)
nrow(eur_sp_poly)
# 14777 grid squares
```
Function to get lon and lat of SW corner of each cell
```R
get_extent <- function(x){
  cbind(extent(x)[c(1, 3)])
}
i <- 1:nrow(eur_sp_poly)
lat_lon <- as_tibble(t(sapply(i, function(i){get_extent(eur_sp_poly[i,])}))) %>%
  rename(lon = V1, lat = V2)
```
Convert to sf to allow more efficient working
```R
eur_sf <- st_as_sf(eur_sp_poly)
# add lon and lat
eur_sf <- bind_cols(eur_sf, lat_lon)
# add rownames as a cell ID variable
eur_sf <- eur_sf %>% rownames_to_column %>%
  rename(cell_id = rowname)
# tidy up by removing 'layer'
eur_sf <- eur_sf %>% dplyr::select(-layer)
```
Get obis checklist for a cell - records since 2000
```R
cell_checklist <- checklist(
  geometry = st_as_text(eur_sf$geometry[12000]),
  year = 2000:2018,
  fields = c("species", "records", "worms_id", "rank_name")) %>%
  filter(rank_name == "Species")
```
Now add back in the species-level thermal affinities - either using the object created above (renamed here for convenience):
```R
spp_t_matched <- all_species_fgrps_t_matched
```
Or you can read in this file directly from the csv provided here:
```R
spp_t_matched <- read_csv("all_species_fgrps_t_matched.csv")
```
Then check functional groups, combine some (e.g. all sub-categories of benthos -> benthos), and filter out species with no functional group information, or poorly specified functional group:
```R
spp_t_matched %>% group_by(functional_group) %>% count()

# combine benthic functional groups, and filter to remove unknown or poorly-specified groups:
spp_t_matched <- spp_t_matched %>%
  mutate(fg = case_when(
    functional_group %in% c("benthos", "epibenthos", "macrobenthos", "hyperbenthos") ~ "benthos",
    functional_group %in% c("macro", "meso", "neuston", "unknown", "plankton") ~ "NA",
    TRUE ~ functional_group
    )) %>%
  filter(fg != "NA")
```
Now you can get the average temperature affinity of any given grid square, separately for the different functional groups, using the `get_cell_t_affin` function in the `R` folder. For cells with no species in this returns an empty tibble, e.g.:
```R
# empty example:
get_cell_t_affin(cell_geom = slice(eur_sf, 1)$geometry,
  c_id = 1, t_dat = spp_t_matched)
```
For cells with species from one or more functional groups, you get a tibble with a row for each functional group, which contains some summary information (`cell_id` to identify the grid square, functional group as `fg`, the number of species of that group within the cell `n_sp`, the total number of OBIS records in the cell for that group `n_rec`, and then a range of assemblage-level temperature affinity summaries. For a working example:
```R
# working example:
get_cell_t_affin(cell_geom = slice(eur_sf, 12000)$geometry,
  c_id = 12000, t_dat = spp_t_matched)
```
This can be run over multiple grid cells like this:
```R
# test over a few cells
test_cell_affins <- sample_n(eur_sf, 10) %>% group_by(cell_id) %>%
  do(get_cell_t_affin(cell_geom = .$geometry, c_id = .$cell_id))
```
To run over all the grid cells in our European data takes approx 30-45 mins to run:
```R
t_affin_grid <- eur_sf %>% group_by(cell_id) %>%
  do(get_cell_t_affin(cell_geom = .$geometry, c_id = .$cell_id))
```
Alternatively you can read in the ouput directly as a csv provided here:
```R
t_affin_grid <- read_csv("t_affin_by_grid_fg.csv")
```
Tidy up and summarise a bit, to give number and proportion of grid cells each functional group occurs in, as well as median and maximum number of species and occurrence records per cell:
```R
t_affin_grid %>% ungroup %>% group_by(fg) %>% summarise(
  n_cell = n(), p_cell = round(n()/14777, 2),
  med_n_sp = median(n_sp), max_n_sp = max(n_sp),
  med_n_rec = median(n_rec), max_n_rec = max(n_rec)
)

# >   fg            n_cell p_cell med_n_sp max_n_sp med_n_rec max_n_rec
# >   <chr>          <int>  <dbl>    <dbl>    <dbl>     <dbl>     <dbl>
# > 1 benthos         2182   0.15        3      984       9       35077
# > 2 birds            948   0.06        3       66      16       19977
# > 3 fish            2486   0.17        3      212      25       25174
# > 4 macroalgae      2161   0.15        3      182       5        4952
# > 5 mammals         1377   0.09        1       16       3         939
# > 6 nekton           799   0.05        2       40       9         591
# > 7 phytoplankton   2224   0.15        3      182       4        9735
# > 8 reptiles         410   0.03        1        3       3         189
# > 9 zooplankton     3950   0.27        4      241       8.5      5267
```
# Match gridded temperature affinities to current and future temperature
The next step is to compare the gridded, functional group-level temperature affinities to current and projected future temperature. To do this we use [`sdmpredictors`](https://github.com/lifewatch/sdmpredictors) to access current and future temperature data from [Bio-ORACLE](http://www.bio-oracle.org). First, load `sdmpredictors` and download required temperature layers to a sensible directory (here, a folder called `biooracle` is created within your current working directory):
```R
library(sdmpredictors)
# set path to store layers
bo_path <- paste0(file.path(getwd()), "/biooracle")
# create directory if needed
if(!dir.exists(file.path(bo_path))){
  dir.create(path = file.path(bo_path),
    recursive = FALSE, showWarnings = FALSE)}
```
The temperature layers downloaded here are mean and max SST and SBT, both current and future scenarios - here the only future scenario downloaded is RCP8.5, for 2050, but other RCP scenarios can be loaded in the same way, as can projections to 2100:
```R
bo_t_dat <- load_layers(layercodes = c(
  "BO_sstmean", "BO_sstmax",
  "BO2_RCP85_2050_tempmean_ss", "BO2_RCP85_2050_tempmax_ss",
  "BO2_tempmean_bdmean", "BO2_tempmax_bdmean",
  "BO2_RCP85_2050_tempmean_bdmean", "BO2_RCP85_2050_tempmax_bdmean"
  ),
  equalarea = FALSE, datadir = "biooracle")
```
These layers are then restricted to Europe and resampled to 0.5 degrees:
```R
eur_bo_t <- raster::resample(bo_t_dat, eur_r, method = "bilinear")
```
These can be plotted if desired:
```R
plot(eur_bo_t, col = RColorBrewer::brewer.pal(9, "Reds"))
```
The next step is to populate `eur_sf` with this temperature data. This requires getting the lon and lat of each gridsquare, which is done by adding a small amount to the SW corner coordinates to ensure the point is within the correct square, then extracting the relevant value from the temperature rasters for each grid cell, and finally adding these back into `eur_sf`:
```R
points <- eur_sf %>% st_set_geometry(NULL) %>% dplyr::select(lon, lat) %>%
  mutate(lon = lon + 0.01, lat = lat + 0.01)
points <- SpatialPoints(points, lonlatproj)

bo_sst <- raster::extract(eur_bo_t[[1]], points)
bo_sst_max <- raster::extract(eur_bo_t[[2]], points)
bo_sst_rcp85_2050 <- raster::extract(eur_bo_t[[3]], points)
bo_sst_max_rcp85_2050 <- raster::extract(eur_bo_t[[4]], points)
bo_sbt <- raster::extract(eur_bo_t[[5]], points)
bo_sbt_max <- raster::extract(eur_bo_t[[6]], points)
bo_sbt_rcp85_2050 <- raster::extract(eur_bo_t[[7]], points)
bo_sbt_max_rcp85_2050 <- raster::extract(eur_bo_t[[8]], points)

eur_sf <- eur_sf %>% bind_cols(
  bo_sst = bo_sst,
  bo_sst_max = bo_sst_max,
  bo_sst_rcp85_2050 = bo_sst_rcp85_2050,
  bo_sst_max_rcp85_2050 = bo_sst_max_rcp85_2050,
  bo_sbt = bo_sbt,
  bo_sbt_max = bo_sbt_max,
  bo_sbt_rcp85_2050 = bo_sbt_rcp85_2050,
  bo_sbt_max_rcp85_2050 = bo_sbt_max_rcp85_2050)
 ```
# Creating maps
Creating the maps requires loading a couple more packages:
```R
library(viridis)
library(maps)
```
You can also create the European base map here (the plotting function will do this if not, but it may be quicker just to create it once rather than with each plot):
```R
eur_dat <- map_data("world") %>% filter(
  long >= -45 & long <= 70 & lat >= 26 & lat <= 90)
```
You can then create maps using the `do_gridded_t_affin_plot` function provided in the `R` folder here. This allows you to plot a single variable for a specified functional group - for instance to map the temperature affinity of benthos, derived using mean sea bottom temperature:
```R
do_gridded_t_affin_plot(fgrp = "benthos", tvar = "bo_sbt_mean", mapdat = eur_dat)
```
This map is found in `bo_sbt_mean benthos.png`. Any combination of functional group and temperature variable available in `t_affin_grid` can be entered into the function. In addition, you can map the difference between functional group-level temperature affinity and one of the environmental temperature layers obtained from Bio-ORACLE. For instance to plot the difference between mean zooplankton thermal affinity based on mean SST, and expected max SST in 2050 under RCP 8.5:
```R
do_gridded_t_affin_plot(fgrp = "zooplankton", tvar = c("bo_sst_mean", "bo_sst_max_rcp85_2050"),
  mapdat = eur_dat)
```
This produces the map in `bo_sst_mean-bo_sst_max_rcp85_2050 zooplankton.png`. You can also do some basic QC by plotting, for instance, number of species per grid square for a given group, to see if extreme temperature affinity values can be explained by a lack of data. For instance for phytoplankton:
```R
do_gridded_t_affin_plot(fgrp = "phytoplankton", tvar = "n_sp", mapdat = eur_dat)
```
which shows the low species richness throughout most of the region (see `n_sp phytoplankton.png`).
