#########################################################################
#####             Functions to align multiple ggplots              ######
#########################################################################


######## align legends (written by Bobby) #######
#Extract the legend from a ggplot object with this function:
g_legend <- function(a.gplot){
  tmp    <- ggplot_gtable(ggplot_build(a.gplot))
  leg    <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if(length(leg) == 0){
    return(NULL)
  }
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#Then use the following function on a list of legends to return a single grob that will contain all the legends aligned:
align.legends.fun  <- function(legend.list){
  aligned.leg <- legend.list[[1]]
  for(i in 2:length(legend.list)){
    leg1        <- legend.list[[i]]$grobs[[1]]
    leg_a       <- gtable_add_rows(aligned.leg, pos = nrow(aligned.leg) - 1, heights = sum(leg1$heights))
    leg_final   <- gtable_add_grob(leg_a, leg1, t = nrow(leg_a) - 1, l = 3)
    aligned.leg <- leg_final
  }
  return(aligned.leg)
}



######## align plots and add tiles ##########
addTiles_ggplot <- function(ggplot_list, tile_data, order_samples,  heights, widths,
                            text_size = 20,
                            categorical   = c('Set1', 'Set2', 'Set3', 'Accent'),
                            numerical     = c('RdPu', 'Purples', 'Blues', 'Greens', 'YlOrBr'),
                            legend_col_n  = 1){
  
  #Parameter information:
  # *ggplot_list = List with ggplots to which the tiles should be added
  # *tile_data = dataframe with samples in first column and all tiles that should be included as columns (column names is used as label)
  # *order_samples = samples in order they are supposed to be plotted
  # *heights and weights = heights and weights for grid plot
  # *text_size = size for text in plots
  # *categorical and numerical = color palletes for tile plots
  # *legend_col_n = value from 1 to 3 to classify in how many columns the legend should be split
  
  #libraries
  require(grid)
  require(gridExtra)
  
  #extract classes of variables and assign color palettes
  color_palettes <- lapply(2:ncol(tile_data), function(x){
    name  <- colnames(tile_data)[x]
    class <- class(tile_data[,x])
    n     <- ifelse(class == 'numeric', NA, length(unique(as.character(tile_data[,x]))))
    data.frame(name, class, n)
  })
  color_palettes <- Reduce(rbind, color_palettes)
  color_palettes$palette <- NA
  color_palettes$palette[color_palettes$class != 'numeric'] <- rep(categorical, ceiling(sum(color_palettes$class != 'numeric') / length(categorical)))[1:sum(color_palettes$class != 'numeric')]
  color_palettes$palette[color_palettes$class == 'numeric'] <- rep(numerical, ceiling(sum(color_palettes$class == 'numeric') / length(numerical)))[1:sum(color_palettes$class == 'numeric')]
  
  
  #create tile plot for each column of tile_data
  tile_list <- lapply(2:ncol(tile_data), function(x){
    name          <- colnames(tile_data)[x]
    plot_data     <- data.frame(sample = tile_data[,1], value = tile_data[,x])
    plot_data[,1] <- factor(plot_data[,1], levels = order_samples)
    
    
    
    if(color_palettes$class[color_palettes$name == name] == 'numeric'){
      
      colours <- brewer.pal(n = 8, color_palettes$palette[color_palettes$name == name])
      
      ggplot(plot_data, aes(x = sample, y = 1, fill = value)) + 
        geom_tile() +
        scale_y_continuous(expand = c(0,0)) +
        scale_x_discrete(expand = c(0,0)) +
        scale_fill_gradient(name = name, low = colours[2], high = colours[6], na.value = '#f0f0f0') +
        ylab(name) +
        theme_bw() +
        theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),
              axis.title.x = element_blank(), text = element_text(size = text_size))
    } else {
      if(name == 'response'){
        ggplot(plot_data, aes(x = sample, y = 1, fill = value)) + 
          geom_tile() +
          scale_y_continuous(expand = c(0,0)) +
          scale_x_discrete(expand = c(0,0)) +
          # scale_fill_manual(name = name, na.value = '#f0f0f0') +
          ylab(name) +
          theme_bw() +
          theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),
                axis.title.x = element_blank(), text = element_text(size = text_size))
      }else {
        colours <- rep(brewer.pal(n=10, color_palettes$palette[color_palettes$name == name]), ceiling(color_palettes$n[color_palettes$name == name] / 8))[1:color_palettes$n[color_palettes$name == name]]
        
        ggplot(plot_data, aes(x = sample, y = 1, fill = value)) + 
          geom_tile() +
          scale_y_continuous(expand = c(0,0)) +
          scale_x_discrete(expand = c(0,0)) +
          scale_fill_manual(name = name, values = colours, na.value = '#f0f0f0') +
          ylab(name) +
          theme_bw() +
          theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),
                axis.title.x = element_blank(), text = element_text(size = text_size))
      }
    }
  })
  
  ### combine plots ###
  plot_list <- c(ggplot_list, tile_list)
  
  #extract legends
  legend_list     <- lapply(plot_list, function(x){ g_legend(x) }) 
  legend_list     <- plyr::compact(legend_list)
  breaks_legend    <- split(1:length(legend_list), sort((1:length(legend_list))%%legend_col_n))
  aligned.legends <- lapply(1:legend_col_n, function(x){
    align.legends.fun(legend.list = legend_list[breaks_legend[[x]]])
  })
  
  #delete legends
  plot_list_withoutLegends <- lapply(1:length(plot_list), function(x){
    p <- plot_list[[x]]
    p <- p + theme(legend.position = 'none')
    return(p)
  })
  
  #convert into ggplotGrob object
  plot_list_ggplotGrop <- lapply(1:length(plot_list_withoutLegends), function(x){
    p <- plot_list_withoutLegends[[x]]
    if(x == 1){
      p <- ggplotGrob(p + theme(plot.margin = unit(c(1, 0.5, 0, 0.5), "cm")))
    } else if(x == length(plot_list_withoutLegends)) {
      p <- ggplotGrob(p + theme(plot.margin = unit(c(0, 0.5, 1, 0.5), "cm")))
    } else{
      p <- ggplotGrob(p + theme(plot.margin = unit(c(0, 0.5, 0, 0.5), "cm")))
    }
    return(p)
  })
  
  #align widths 
  all_widths <- lapply(plot_list_ggplotGrop, function(x) {x$widths})
  plot_list_alignedWiths <- lapply(plot_list_ggplotGrop, function(x){
    x$widths <- do.call(unit.pmax, all_widths)
    return(x)
  })
  
  
  #set heights
  g <- do.call(gtable_rbind, plot_list_alignedWiths)
  id_panels_h <- unique(g$layout[g$layout$name=="panel", "t"])
  g$heights[id_panels_h] <- grid::unit(heights, "null")
  
  
  #plot
  if(legend_col_n == 1){
    grid.arrange(g, aligned.legends[[1]], ncol = 2, widths = widths)
  } 
  
  if(legend_col_n == 2){
    grid.arrange(g, aligned.legends[[1]], aligned.legends[[2]], ncol = 3, widths = widths)
  } 
  
  if(legend_col_n == 3){
    grid.arrange(g, aligned.legends[[1]], aligned.legends[[2]], aligned.legends[[3]], ncol = 3, widths = widths)
  } 
}
