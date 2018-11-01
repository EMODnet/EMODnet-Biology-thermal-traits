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
This is a large query, so the results are provided in the `data` folder here as a csv file that can be read in directly:
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
This is also available as a csv file in the `data` folder:
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
We then obtained functional groups for these species using attributes from [WoRMS](http://www.marinespecies.org) using the [`worrms`](https://github.com/ropensci/worrms) package and a function written specifically for this purpose available from https://github.com/tomjwebb/WoRMS-functional-groups:
```R
species_attr <- eur_spp_full %>%
  group_by(AphiaID) %>%
  do(get_worms_fgrp(AphiaID = .$AphiaID))
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
For an example using benthos data
```R
# add benthos_spp data to the cluster
cluster_assign_value(cluster20, 'benthos_spp', benthos_spp)
# partition the benthos_spp data
benthos_spp <- benthos_spp %>% partition(AphiaID, cluster = cluster20)
benthos_spp

# run t matching functions #
benthos_t_matching <- benthos_spp %>% do(get_temp_summ_by_sp(sp_id = .$AphiaID)) %>% collect()

# write output
write_csv(benthos_t_matching, "/benthos_ad_t_matching.csv")
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
Set up a raster based on EMODNet 'europe' extent. For now set resolution = 0.5; could refine to 0.1 if needed
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
