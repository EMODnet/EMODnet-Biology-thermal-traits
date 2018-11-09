get_cell_t_affin <- function(cell_geom, c_id, t_dat = spp_t_matched){

  # Function to get obis checklist for a cell - records since 2000
  # and then match to t-matched data and get summaries

  # get checklist
  cell_cl <- checklist(geometry = st_as_text(cell_geom),
    year = 2000:2018,
    fields = c("species", "records", "worms_id", "rank_name")) %>%
    filter(rank_name == "Species")
	
  if(nrow(cell_cl) > 0){
    cell_cl <- left_join(cell_cl,
      dplyr::select(t_dat, AphiaID, fg, iap_t_mean:bo_sbt_mean,
        iap_t_min:bo_sbt_min, iap_t_max:bo_sbt_max, iap_t_q5:bo_sbt_q5,
        iap_t_q95:bo_sbt_q95),
      by = c("worms_id" = "AphiaID")) %>%
      filter(!is.na(fg)) %>%
      group_by(fg) %>%
      summarise(n_sp = n(), n_rec = sum(records),
        iap_t_mean = weighted.mean(iap_t_mean, w = records),
        iap_t_q5 = weighted.mean(iap_t_q5, w = records),
        iap_t_q95 = weighted.mean(iap_t_q95, w = records),
        bo_sst_mean = weighted.mean(bo_sst_mean, w = records),
        bo_sst_q5 = weighted.mean(bo_sst_q5, w = records),
        bo_sst_q95 = weighted.mean(bo_sst_q95, w = records),
        bo_sbt_mean = weighted.mean(bo_sbt_mean, w = records),
        bo_sbt_q5 = weighted.mean(bo_sbt_q5, w = records),
        bo_sbt_q95 = weighted.mean(bo_sbt_q95, w = records)
      ) %>%
      mutate(cell_id = c_id) %>%
      dplyr::select(cell_id, everything())
      # some extraneous variables seem to come in here - species, records, etc.
    }

  # add an option to save cell-level checklists to file?

  # return summary data	
  cell_cl
  
}
