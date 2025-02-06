
# This function creates a plot to show changing accuracy (mae or rsq) with
# proportions on x

prpn_plot <- function(dat, prpn, samp, samp_name, fills, cls, y, x, lgnd = F) {
  # Create the plot
  p <- ggplot(data = dat, 
              aes(x = get(prpn), y = mean, colour = get(samp), fill = get(samp))) +
    geom_ribbon(aes(ymin = mean - ci, ymax = mean + ci), alpha = 0.5, colour = NA) +
    geom_point() + 
    geom_line() +
    scale_fill_manual(values = c(fills[1], fills[2]), name = samp_name) +
    scale_colour_manual(values = c(cls[1], cls[2]), name = samp_name) +
    theme(panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid = element_blank(),
          plot.margin = unit(c(0.25, 0.25, 0.25, .5), 'cm'),
          axis.title.y = element_text(size = 15, colour = 'black', vjust = 3),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 15, colour = 'black'),
          legend.text = element_text(size = 15, colour = 'black'),
          legend.title = element_text(size = 15, colour = 'black'),
          legend.key.spacing.y = unit(0.2, 'cm')) +
    labs(y = y, x = '')
  # Remove legend if specified
  if(isFALSE(lgnd)) {
    p <- p +
      theme(legend.position = 'none')
  }
  # Return the plot
  return(p)
}

