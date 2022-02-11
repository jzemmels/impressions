x3pComparisonPlot <- function(reference,target,x3pNames = c("x3p1","x3p2")){

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
                                     .y = x3pNames,
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
    mutate(
      # value = value - median(value,na.rm = TRUE),
      # alpha = ifelse(valueDiff <= quantile(valueDiff,.5,na.rm = TRUE),1,0),
      alpha = ifelse(valueDiff <= 1,1,0),
      # alpha = scales::rescale(sqrt(valueDiff),to = c(1,0)),
      x3pName = "Pixelwise Average") %>%
    select(-valueDiff)

  refTargCombined <- bind_rows(refTargCombined,
                               surfaceMat_df)

  # add to the df created above the pixelwise difference between the reference
  # and average scans. alpha-blend pixels
  surfaceMat_df <- reference %>%
    x3pToDF() %>%
    mutate(
      value = value,
      # value = value - median(value,na.rm = TRUE)
    ) %>%
    left_join(refTargAverage %>%
                rename(aveValue = value),
              by = c("x","y")) %>%
    mutate(
      # value = value - median(value,na.rm = TRUE),
      value = ifelse(is.na(aveValue),NA,value),
      # valueDiff = abs(value - aveValue)
    ) %>%
    mutate(
      # alpha = ifelse(valueDiff > quantile(valueDiff,.5,na.rm = TRUE),1,0),
      alpha = ifelse(valueDiff > 1,1,0),
      # alpha = scales::rescale(sqrt(valueDiff),to = c(0,1)),
      x3pName = paste0(x3pNames[1]," diff")) %>%
    select(-c(valueDiff,aveValue))

  refTargCombined <- bind_rows(refTargCombined,
                               surfaceMat_df)

  surfaceMat_df <- target %>%
    x3pToDF() %>%
    mutate(
      # value = value - median(value,na.rm = TRUE)
    ) %>%
    left_join(
      refTargAverage %>%
        rename(aveValue = value),
      by = c("x","y")) %>%
    mutate(
      # value = value - median(value,na.rm = TRUE),
      value = ifelse(is.na(aveValue),NA,value),
      # valueDiff = abs(value - aveValue)
    ) %>%
    mutate(
      # alpha = ifelse(valueDiff > quantile(valueDiff,.5,na.rm = TRUE),1,0),
      alpha = ifelse(valueDiff > 1,1,0),
      #alpha = scales::rescale(sqrt(valueDiff),to = c(0,1)),
      x3pName = paste0(x3pNames[2]," diff")) %>%
    select(-c(valueDiff,aveValue))

  refTargCombined <- bind_rows(refTargCombined,
                               surfaceMat_df)


  x3pPlt <- refTargCombined %>%
    mutate(x3pName = factor(x3pName,levels = c(x3pNames,"Pixelwise Average",paste0(x3pNames," diff")))) %>%
    ggplot(aes(x=x,y=y)) +
    # geom_raster(fill = "gray80") +
    geom_raster(aes(fill=value,alpha = alpha)) +
    ggplot2::scale_fill_gradientn(colours = c("#2d004b","#542788","#8073ac","#b2abd2","#d8daeb","#f7f7f7","#fee0b6","#fdb863","#e08214","#b35806","#7f3b08"),
                                  values = scales::rescale(quantile(refTargCombined$value,c(0,.01,.025,.1,.25,.5,.75,0.9,.975,.99,1),na.rm = TRUE)),
                                  breaks = function(lims){
                                    dat <- quantile(refTargCombined$value,c(0,.01,.25,.5,.75,.99,1),na.rm = TRUE)

                                    dat <- dat %>%
                                      setNames(paste0(names(dat)," [",round(dat,1),"]"))

                                    return(dat)
                                  },
                                  oob = scales::oob_keep,
                                  limits = range(refTargCombined$value),
                                  na.value = "gray65") +
    labs(fill = expression("Rel. Height [Norm.]")) +
    ggplot2::guides(fill = ggplot2::guide_colourbar(barheight = grid::unit(3,"in"),
                                                    label.theme = ggplot2::element_text(size = 8),
                                                    title.theme = ggplot2::element_text(size = 10),
                                                    frame.colour = "black",
                                                    ticks.colour = "black"),
                    colour = 'none',
                    alpha = "none") +
    facet_wrap(~ x3pName)

  plt <- ggplot2::ggplot_build(x3pPlt)

  pltLegend <- cowplot::get_legend(plt$plot)

  ret <- plt$data[[1]] %>%
    group_by(PANEL) %>%
    group_split() %>%
    map2(c(x3pNames,"Pixelwise Average",paste0(x3pNames," diff")),
         function(dat,title){

           # dat_notGray <- dat %>%
           #   filter(!str_detect(fill,"gray"))

           # return(length(unique(dat_notGray$fill)))

           ggplot(data = dat,
                  aes(x=x,y=y)) +
             geom_raster(aes(fill=fill,alpha=alpha)) +
             theme_minimal() +
             theme(legend.position = "none") +
             scale_fill_identity() +
             scale_alpha_identity() +
             ggplot2::coord_fixed(expand = FALSE) +
             ggplot2::theme_minimal() +
             ggplot2::theme(
               axis.title.x = ggplot2::element_blank(),
               axis.text.x = ggplot2::element_blank(),
               axis.ticks.x = ggplot2::element_blank(),
               axis.title.y = ggplot2::element_blank(),
               axis.text.y = ggplot2::element_blank(),
               axis.ticks.y = ggplot2::element_blank(),
               panel.grid.major = ggplot2::element_blank(),
               panel.grid.minor = ggplot2::element_blank(),
               panel.background = ggplot2::element_rect(fill = "gray80")
               ,plot.title = element_text(size = 9,hjust = .5)
               # ,legend.position = "bottom"
             ) +
             labs(title = title)

         })

  ret <- c(ret,list(pltLegend))  %>%
    set_names(c(x3pNames,"Pixelwise Average",paste0(x3pNames," diff"),"legend"))

  return(ret)

}

