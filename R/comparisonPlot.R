x3pComparisonPlot <- function(reference,
                              target,
                              plotNames = c("x3p1","x3p2","Element-wise Average","x3p1 diff.","x3p2 diff."),
                              unit = "Norm.",
                              cutoffThresh = 1){

  # for the life of me, I can't replicate what facet_wrap does to put everything
  # on the exact same color scale (limits and all) when I try to plot each
  # surface matrix individually - no matter what I set for limits, values, etc.
  # in the call to scale_fill_gradientn. So the code below is a wonky workaround
  # where I plot all of the surface matrices using facet_wrap, then I "build"
  # the plot to get at the data that ggplot uses to when it creates the plot.
  # this is a data frame containing x,y, hex code fill, and alpha values - so it
  # no longer corresponds directly to the surface values. However, using
  # scale_*_identity() tells ggplot to directly interpret the column elements
  # *as* the appropriate value at that pixel (e.g., an element of the alpha
  # column being 0 means make that pixel clear.) This may come back to shoot me
  # in the foot if I ever try to further develop upon these functions, but oh
  # well.

  # to keep a consistent color scheme across all plots, combine the values from
  # the reference and target and give these to the values argument in each
  # ggplot scale_colour_gradientn call
  refTargCombined <- purrr::map2_dfr(.x = list(reference,target),
                                     .y = plotNames[1:2],
                                     function(x3p,name){

                                       x3p$header.info$incrementX <- 1
                                       x3p$header.info$incrementY <- 1

                                       x3p %>%
                                         x3ptools::x3p_to_df() %>%
                                         #perform some transformations on the
                                         #x,y values so that the plot is
                                         #representative of the actual surface
                                         #matrix (i.e., element [1,1] of the
                                         #surface matrix is in the top-left
                                         #corner)
                                         dplyr::mutate(xnew = max(y) - y,
                                                       ynew = max(x) - x,
                                                       # ,value = .data$value - median(.data$value,na.rm = TRUE)
                                         ) %>%
                                         dplyr::select(-c(x,y)) %>%
                                         dplyr::rename(x=xnew,
                                                       y=ynew) %>%
                                         mutate(x3pName = name)
                                     }) %>%
    mutate(alpha = 1)

  # pixelwise average and absolute difference between the two scans
  refTargAverage <- bind_rows(reference %>%
                                x3pToDF() %>%
                                mutate(value = value),
                              target %>%
                                x3pToDF() %>%
                                mutate(value = value)) %>%
    group_by(x,y) %>%
    summarise(valueDiff = abs(diff(value)),
              value = mean(value),
              .groups = "drop")

  # add to the df created above data on the pixelwise average between ref and
  # targ scans. alpha-blend pixels where the pixelwise difference is large.
  surfaceMat_df <- refTargAverage %>%
    mutate(alpha = ifelse(valueDiff <= cutoffThresh,1,0),
           x3pName = plotNames[3]) %>%
    select(-valueDiff)

  refTargCombined <- bind_rows(refTargCombined,
                               surfaceMat_df)

  # add to the df created above the pixelwise difference between the reference
  # and average scans. alpha-blend pixels
  surfaceMat_df <- reference %>%
    x3pToDF() %>%
    left_join(refTargAverage %>%
                rename(aveValue = value),
              by = c("x","y")) %>%
    mutate(value = ifelse(is.na(aveValue),NA,value)) %>%
    mutate(alpha = ifelse(valueDiff > cutoffThresh,1,0),
           x3pName = plotNames[4]) %>%
    select(-c(valueDiff,aveValue))

  refTargCombined <- bind_rows(refTargCombined,
                               surfaceMat_df)

  surfaceMat_df <- target %>%
    x3pToDF() %>%
    left_join(
      refTargAverage %>%
        rename(aveValue = value),
      by = c("x","y")) %>%
    mutate(value = ifelse(is.na(aveValue),NA,value)) %>%
    mutate(alpha = ifelse(valueDiff > cutoffThresh,1,0),
           x3pName = plotNames[5]) %>%
    select(-c(valueDiff,aveValue))

  refTargCombined <- bind_rows(refTargCombined,
                               surfaceMat_df)


  x3pPlt <- refTargCombined %>%
    mutate(x3pName = factor(x3pName,levels = plotNames)) %>%
    ggplot(aes(x=x,y=y)) +
    # geom_raster(fill = "gray80") +
    geom_raster(aes(fill=value,alpha = alpha)) +
    ggplot2::scale_fill_gradientn(colours = c("#2d004b","#542788","#8073ac","#b2abd2","#d8daeb","#f7f7f7","#fee0b6","#fdb863","#e08214","#b35806","#7f3b08"),
                                  values = scales::rescale(quantile(refTargCombined$value,c(0,.01,.025,.1,.25,.5,.75,0.9,.975,.99,1),na.rm = TRUE)),
                                  breaks = function(lims){
                                    dat <- quantile(refTargCombined$value,c(0,.01,.25,.5,.75,.99,1),na.rm = TRUE)

                                    dat <- dat %>%
                                      setNames(paste0(names(dat),"\n[",round(dat,1),"]"))

                                    return(dat)
                                  },
                                  oob = scales::oob_keep,
                                  limits = range(refTargCombined$value),
                                  na.value = "gray65") +
    # labs(fill = expression("Rel. Height [Norm.]")) +
    labs(fill = paste0("Rel. Height\n[",unit,"]")) +
    ggplot2::guides(fill = ggplot2::guide_colourbar(barwidth = grid::unit(2,"in"),
                                                    label.theme = ggplot2::element_text(size = 6),
                                                    title.theme = ggplot2::element_text(size = 8),
                                                    frame.colour = "black",
                                                    ticks.colour = "black"),
                    colour = 'none',
                    alpha = "none") +
    facet_wrap(~ x3pName) +
    ggplot2::scale_alpha_identity() +
    ggplot2::theme(legend.direction = "horizontal")

  plt <- ggplot2::ggplot_build(x3pPlt)

  pltLegend <- cowplot::get_legend(plt$plot)

  ret <- plt$data[[1]] %>%
    group_by(PANEL) %>%
    group_split() %>%
    map2(plotNames,
         function(dat,title){

           # dat_notGray <- dat %>%
           #   filter(!str_detect(fill,"gray"))

           # return(length(unique(dat_notGray$fill)))

           ggplot(data = dat,
                  aes(x=x,y=y)) +
             geom_raster(aes(fill=fill,alpha=alpha)) +
             theme_minimal() +
             scale_fill_identity() +
             scale_alpha_identity(limits = c(0,1)) +
             ggplot2::coord_fixed(expand = FALSE) +
             ggplot2::theme_minimal() +
             ggplot2::theme(
               legend.position = "none",
               axis.title.x = ggplot2::element_blank(),
               axis.text.x = ggplot2::element_blank(),
               axis.ticks.x = ggplot2::element_blank(),
               axis.title.y = ggplot2::element_blank(),
               axis.text.y = ggplot2::element_blank(),
               axis.ticks.y = ggplot2::element_blank(),
               panel.grid.major = ggplot2::element_blank(),
               panel.grid.minor = ggplot2::element_blank(),
               panel.background = ggplot2::element_rect(fill = "gray80")
               ,plot.title = element_text(size = 7,hjust = .5)
               # ,legend.position = "bottom"
             ) +
             labs(title = title)

         })

  ret <- c(ret,list(pltLegend))  %>%
    set_names(c(plotNames,"legend"))

  return(ret)

}

