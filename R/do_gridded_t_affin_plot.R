do_gridded_t_affin_plot <- function(fgrp,
  sp_dat = eur_sf, t_dat = t_affin_grid, tvar,
  mapdat = NULL, print_plot = TRUE, save_plot = FALSE){

  # Function to plot gridded temperature affinities of a given functional group
  # Will plot either raw affinity (if tvar is a single variable)
  # Or difference between affinity and gridded temperature (if tvar is two vars)
  
  # create europe map data if needed
  if(is.null(mapdat)){
  	mapdat <- map_data("world") %>% filter(
      long >= -45 & long <= 70 & lat >= 26 & lat <= 90)
  }
  
  theme_black_map = function(base_size = 12, base_family = "") {
    # a theme for a black background map
    theme_grey(base_size = base_size, base_family = base_family) %+replace%
    
    theme(
      # Specify axis options
      axis.line = element_blank(),  
      axis.text.x = element_blank(),  
      axis.text.y = element_blank(),  
      axis.ticks = element_blank(),  
      axis.title.x = element_blank(),  
      axis.title.y = element_blank(),  
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"),  
      legend.key = element_rect(color = "black",  fill = "black"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*0.8, color = "white"),  
      legend.title = element_text(
        size = base_size*0.8, face = "bold", hjust = 0, color = "white"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.background = element_rect(fill = "black", color  =  NA),  
      panel.border = element_rect(fill = NA, color = "grey50"),  
      panel.grid.major = element_line(color = "grey35"),  
      panel.grid.minor = element_line(color = "grey20"),  
      panel.spacing = unit(0.5, "lines"),   
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),  
      strip.text.x = element_text(size = base_size*0.8, color = "white"),  
      strip.text.y = element_text(
        size = base_size*0.8, color = "white",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),  
      plot.title = element_text(size = base_size*1.2, color = "white"),  
      plot.margin = unit(rep(1, 4), "lines")
      )
  }
  
  fg_grid_dat <- left_join(
    sp_dat, filter(t_dat, fg == fgrp), by = "cell_id") %>%
    filter(!is.na(fg))
  
  if(length(tvar) == 1){
    # add single temperature variable to plot if length(tvar) is 1
    t_temp <- fg_grid_dat %>% st_set_geometry(NULL) %>%
      dplyr::select(tvar) %>%
      as_tibble() %>%
      rename(t_plot = tvar)
    fg_grid_dat <- fg_grid_dat %>% bind_cols(t_temp)

    # generate a plot title	
    plot_tit <- paste0(tvar, ", ", fgrp)
    # create the plot
    t_plot <- ggplot(fg_grid_dat) +
      geom_raster(aes(x = lon, y = lat, fill = t_plot)) +
      scale_fill_viridis() +
      geom_polygon(data = mapdat,
        aes(x = long, y = lat, group = group), fill = "grey25") +
      theme_black_map() +
      ggtitle(plot_tit)
      
  } else if(length(tvar) == 2){
    # add temperature difference to plot if length(tvar) is 2
    t_temp <- fg_grid_dat %>% st_set_geometry(NULL) %>%
      dplyr::select(tvar)
    t_temp <- mutate(t_temp, t_plot = t_temp[,1] - t_temp[,2]) %>%
      as_tibble()
    fg_grid_dat <- fg_grid_dat %>% bind_cols(dplyr::select(t_temp, t_plot))

    # generate plot title
    plot_tit <- paste0(tvar[1], " - ", tvar[2], ", ", fgrp)
    # create the plot
    t_plot <- ggplot(fg_grid_dat) +
      geom_raster(aes(x = lon, y = lat, fill = t_plot)) +
      scale_fill_gradient2() +
      geom_polygon(data = mapdat,
        aes(x = long, y = lat, group = group), fill = "grey25") +
      theme_black_map() +
      ggtitle(plot_tit)
      }
    
  # print and/or save plot if requested
  if(print_plot == TRUE){print(t_plot)}
  if(save_plot == TRUE){
    ggsave(filename = paste0(plot_tit, ".png"), plot = t_plot)}	
}
