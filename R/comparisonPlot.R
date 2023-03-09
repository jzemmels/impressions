# this function allows one to slip a ggplot layer underneath the elements of an
# existing ggplot
`-.gg` <- function(plot, layer) {
  if (missing(layer)) {
    stop("Cannot use `-.gg()` with a single argument. Did you accidentally put - on a new line?")
  }
  if (!ggplot2::is.ggplot(plot)) {
    stop('Need a plot on the left side')
  }
  plot$layers = c(layer, plot$layers)
  plot
}

#' Plot filtered element-wise average and differences for two cartridge case
#' scans
#' @name x3p_comparisonPlot
#'
#' @param x3p1 an x3p object
#' @param x3p2 another x3p object
#' @param threshold the default filtering threshold. Defaults to a scalar (1
#'   micron = 1e-6 meters), but can also be set to a scalar-valued function that
#'   takes x3p1 and x3p2 as arguments. For example, threshold =
#'   impressions::x3p_sd will use the joint standard deviation of x3p1 and x3p2
#'   as the threshold.
#' @param plotLabels a character vector of five elements that will display as
#'   labels on each plot
#' @param labelSize font size for the plot labels
#' @param label_x horizontal location of the plot labels
#' @param label_y vertical location of the plot labels
#' @param type return the five plots in a single, faceted plot or as a list
#' @param legendLength length of the plot legend. Passed to the barwidth
#'   argument of the ggplot2::guide_colorbar function
#' @param legendUnit unit of measurement for the surface values (typically
#'   either "micron" or "Norm.")
#' @param legendHoriz horizontal location of the legend. Passed to the
#'   patchwork::inset_element function
#' @param legendQuantiles quantile values to be shown as legend labels. Passed
#'   within the breaks argument of the ggplot2::scale_fill_gradientn function
#'
#' @seealso dplyr::guide_colorbar
#' @importFrom dplyr select
#'
#' @examples
#' data("K013sA1","K013sA2")
#'
#' compData <- cmcR::comparison_allTogether(reference = K013sA1,
#'                                          target = K013sA2,
#'                                          theta = 3,numCells = c(1,1),
#'                                          maxMissingProp = .99,
#'                                          sideLengthMultiplier = 1.1,
#'                                          returnX3Ps = TRUE)
#'
#' x3p_comparisonPlot(x3p1 = compData$cellHeightValues[[1]],x3p2 = compData$alignedTargetCell[[1]])
#'
#' @export
x3p_comparisonPlot <- function(x3p1,
                               x3p2,
                               threshold = 1,
                               plotLabels = c("x3p1","x3p2",
                                             "Element-wise Average",
                                             "x3p1 diff.","x3p2 diff."),
                               labelSize = 4,
                               label_x = ncol(x3p1$surface.matrix)/2,
                               label_y = nrow(x3p1$surface.matrix)/2,
                               type = "faceted",
                               legendLength = grid::unit(3,"in"),
                               legendUnit = "Norm.",
                               legendHoriz = -1.2,
                               legendQuantiles = c(0,.01,.25,.5,.75,.99,1)){

  stopifnot(type %in% c("faceted","list"))

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

  if(is.function(threshold)){
    cutoffThresh <- threshold(x3p1,x3p2)
  }
  else{
    cutoffThresh <- threshold
  }

  x3pAveraged <- x3p_filter(x3p = x3p_elemAverage(x3p1,x3p2),
                            cond = function(x,y,thresh) abs(y) <= thresh,
                            y = c({x3p1$surface.matrix - x3p2$surface.matrix}),
                            thresh = cutoffThresh)

  x3p1Differences <- x3p_filter(x3p = x3p1,
                                cond = function(x,y,thresh) abs(x - y) > thresh,
                                y = c(x3p2$surface.matrix),
                                thresh = cutoffThresh)

  x3p2Differences <- x3p_filter(x3p = x3p2,
                                cond = function(x,y,thresh) abs(x - y) > thresh,
                                y = c(x3p1$surface.matrix),
                                thresh = cutoffThresh)

  # to keep a consistent color scheme across all plots, combine the values from
  # x3p1 and x3p2 and give these to the values argument in each ggplot
  # scale_colour_gradientn call
  x3pCombined <- purrr::map2_dfr(.x = list(x3p1,x3p2,
                                           x3pAveraged,
                                           x3p1Differences,x3p2Differences),
                                 .y = plotLabels,
                                 function(x3p,name){

                                   x3p %>%
                                     x3p_to_dataFrame(preserveResolution = FALSE) %>%
                                     dplyr::mutate(x3pName = name,
                                                   alpha = ifelse(!is.na(value),1,0))
                                 })%>%
    dplyr::mutate(x3pName = factor(x3pName,levels = plotLabels))

  x3pPlt <- x3pCombined %>%
    ggplot2::ggplot(ggplot2::aes(x=x,y=y)) +
    ggplot2::geom_raster(ggplot2::aes(fill=value,alpha = alpha)) +
    ggplot2::scale_fill_gradientn(colours = c("#2d004b","#542788","#8073ac","#b2abd2","#d8daeb","#f7f7f7","#fee0b6","#fdb863","#e08214","#b35806","#7f3b08"),
                                  values = scales::rescale(quantile(x3pCombined$value,c(0,.01,.025,.1,.25,.5,.75,0.9,.975,.99,1),na.rm = TRUE)),
                                  breaks = function(lims){
                                    dat <- quantile(x3pCombined$value,legendQuantiles,na.rm = TRUE)

                                    dat <- dat %>%
                                      setNames(paste0(names(dat),"\n[",round(dat,1),"]"))

                                    return(dat)
                                  },
                                  oob = scales::oob_keep,
                                  limits = range(x3pCombined$value),
                                  na.value = "gray65") +
    ggplot2::labs(fill = paste0("Rel. Height\n[",legendUnit,"]")) +
    ggplot2::guides(fill = ggplot2::guide_colourbar(barwidth = legendLength,
                                                    label.theme = ggplot2::element_text(size = 6),
                                                    title.theme = ggplot2::element_text(size = 8),
                                                    frame.colour = "black",
                                                    ticks.colour = "black"),
                    colour = 'none',
                    alpha = "none") +
    ggplot2::facet_wrap(~ x3pName) +
    ggplot2::scale_alpha_identity() +
    ggplot2::theme(legend.direction = "horizontal")

  # next, we build the ggplot and return a list of individual plots that are all
  # on the same color scale

  plt <- ggplot2::ggplot_build(x3pPlt)

  pltLegend <- plt$plot %>%
    cowplot::get_legend() %>%
    cowplot::plot_grid()

  ret <- plt$data[[1]] %>%
    dplyr::group_by(PANEL) %>%
    dplyr::group_split() %>%
    purrr::map(function(dat){

      ggplot2::ggplot(data = dat,
                      ggplot2::aes(x=x,y=y)) +
        ggplot2::geom_raster(ggplot2::aes(fill=fill,alpha=alpha)) +
        ggplot2::theme_minimal() +
        ggplot2::scale_fill_identity() +
        ggplot2::scale_alpha_identity(limits = c(0,1)) +
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
          panel.background = ggplot2::element_rect(fill = "gray65"),
          plot.title = ggplot2::element_text(size = 7,hjust = .5)
        )

    })

  patchComparisonPlts <- c(ret,list(pltLegend))

  # next, we use a labeling algorithm to identify the borders of the filtered
  # regions and plot these borders on top of the filtered plots created in the
  # ret object

  averageBinarized <- x3pAveraged %>%
    x3p_to_dataFrame(preserveResolution = FALSE) %>%
    dplyr::mutate(value = (abs(c({x3p1$surface.matrix - x3p2$surface.matrix})) > cutoffThresh))

  outline <- filterBoundaries(averageBinarized)

  combinedValues <-  x3p1 %>%
    impressions::x3p_to_dataFrame() %>%
    dplyr::rename(refValue = value) %>%
    dplyr::left_join(x3p2 %>%
                       impressions::x3p_to_dataFrame() %>%
                       dplyr::rename(targValue = value),
                     by = c("x","y"))

  topLeft <- patchComparisonPlts[[1]] +
    # cowplot::theme_nothing() +
    ggplot2::theme(plot.margin = ggplot2::margin(0,0,5,0)) +
    ggplot2::geom_raster(data = combinedValues %>%
                           dplyr::filter(is.na(refValue) & !is.na(targValue)),
                         fill = "gray40") +
    ggplot2::annotate(x = label_x,
                      y = label_y,
                      geom = "text",
                      size = labelSize,
                      label = plotLabels[1])


  bottomLeft <-patchComparisonPlts[[2]] +
    # cowplot::theme_nothing() +
    ggplot2::theme(plot.margin = ggplot2::margin(-20,-100,30,-100)) +
    ggplot2::geom_raster(data = combinedValues %>%
                           dplyr::filter(!is.na(refValue) & is.na(targValue)),
                         fill = "gray40") +
    ggplot2::annotate(x = label_x,
                      y = label_y,
                      geom = "text",
                      size = labelSize,
                      label = plotLabels[2])

  middle <- patchComparisonPlts[[3]] -
    ggplot2::geom_raster(data = averageBinarized %>%
                           dplyr::filter(!is.na(value)),
                         ggplot2::aes(x=x,y=y),fill="gray80",
                         inherit.aes = FALSE) +
    ggplot2::annotate(x = label_x,
                      y = label_y,
                      geom = "text",
                      size = labelSize,
                      label = plotLabels[3]) +
    ggplot2::theme(plot.margin = ggplot2::margin(0,25,0,25)) +
    ggplot2::geom_path(data = outline,  color = "grey40",
                       ggplot2::aes(x=long,y=lat,group=group),
                       colour = "gray40",
                       inherit.aes = FALSE,
                       size = .2)

  topRight <- patchComparisonPlts[[4]] -
    ggplot2::geom_raster(data = averageBinarized %>%
                           dplyr::filter(!is.na(value)),
                         ggplot2::aes(x=x,y=y),fill="gray80",
                         inherit.aes = FALSE) +
    ggplot2::annotate(x = label_x,
                      y = label_y,
                      geom = "text",
                      size = labelSize,
                      label = plotLabels[4]) +
    ggplot2::theme(plot.margin = ggplot2::margin(0,0,5,0))+
    ggplot2::geom_path(data = outline,  color = "grey40",
                       ggplot2::aes(x=long,y=lat,group=group),
                       colour = "gray40",
                       inherit.aes = FALSE,
                       size = .1)

  bottomRight <- patchComparisonPlts[[5]] -
    ggplot2::geom_raster(data = averageBinarized %>%
                           dplyr::filter(!is.na(value)),
                         ggplot2::aes(x=x,y=y),fill="gray80",
                         inherit.aes = FALSE) +
    ggplot2::annotate(x = label_x,
                      y = label_y,
                      geom = "text",
                      size = labelSize,
                      label = plotLabels[5]) +
    ggplot2::theme(plot.margin = ggplot2::margin(-20,-100,30,-100))+
    ggplot2::geom_path(data = outline, color = "grey40",
                       ggplot2::aes(x=long,y=lat,group=group),
                       colour = "gray40",
                       inherit.aes = FALSE,
                       size = .1)

  if(type == "faceted"){

    design <- "ACCD\nBCCE"

    return(patchwork::wrap_plots(topLeft,bottomLeft,middle,topRight,bottomRight,design = design) +
             patchwork::inset_element(pltLegend,left = legendHoriz,bottom = 0,right = legendHoriz,top = 0,
                                      on_top = FALSE,align_to = 'full'))

  }
  if(type == "list")
  {
    return(list(topLeft,bottomLeft,middle,topRight,bottomRight,pltLegend) %>%
             purrr::set_names(c(plotLabels,"legend")))
  }

}


