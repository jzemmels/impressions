#'Plot a list of x3ps
#'@name x3pPlot
#'
#'@description Plots the surface matrices in a list of x3p objects. Either
#'  creates one plot faceted by surface matrix or creates individual plots per
#'  surface matrix and returns them in a list.
#'
#'@param ... one or more x3p objects
#'@param x3pNames character vector containing names of each x3p object
#'@param type dictates whether one plot faceted by surface matrix or a list of
#'  plots per surface matrix is returned. The faceted plot will have a
#'  consistent height scale across all surface matrices.
#'@param legend.quantiles vector of quantiles to be shown as tick marks on
#'  legend plot
#'@param height.quantiles vector of quantiles associated with each color defined
#'  in the height.colors argument
#'@param height.colors vector of colors to be passed to scale_fill_gradientn
#'  that dictates the height value colorscale
#'@param na.value color to be used for NA values (passed to
#'  scale_fill_gradientn)
#'@param legendLength length of the plot legend. Passed to the barwidth
#'   argument of the ggplot2::guide_colorbar function
#'@return A ggplot object or list of ggplot objects showing the surface matrix
#'  height values.
#' @examples
#'data("K013sA1","K013sA2")
#'
#' x3pPlot(K013sA1,K013sA2,x3pNames = c("Scan A","Scan B"))
#'@export
#'
#'@importFrom stats setNames median quantile
#'@importFrom rlang .data

x3pPlot <- function(...,
                    x3pNames = NULL,
                    type = "faceted",
                    legend.quantiles = c(0,.01,.25,.5,.75,.99,1),
                    height.quantiles = c(0,.01,.025,.1,.25,.5,.75,0.9,.975,.99,1),
                    height.colors = rev(c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b')),
                    na.value = "gray65",
                    legendLength = grid::unit(3,"in")){

  x3pList <- list(...)

  if(is.null(x3pNames)){
    x3pList <- setNames(x3pList,paste0("x3p",1:length(x3pList)))
  }
  else{
    x3pList <- x3pList %>% purrr::set_names(x3pNames)
  }

  if(type == "faceted"){
    surfaceMat_df <- purrr::pmap_dfr(.l = list(x3pList,
                                               names(x3pList)),
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
                                         dplyr::mutate(xnew = max(.data$y) - .data$y,
                                                       ynew = max(.data$x) - .data$x,
                                                       value = .data$value - median(.data$value,na.rm = TRUE)) %>%
                                         dplyr::select(-c(.data$x,.data$y)) %>%
                                         dplyr::rename(x=.data$xnew,
                                                       y=.data$ynew) %>%
                                         dplyr::mutate(x3p = rep(name,times = nrow(.)))
                                     }) %>%
      dplyr::mutate(x3p = factor(.data$x3p,levels = names(x3pList)))

    plts <- surfaceMat_df %>%
      ggplot2::ggplot(ggplot2::aes(x = .data$x,y = .data$y)) +
      ggplot2::geom_raster(ggplot2::aes(fill = .data$value))  +
      ggplot2::scale_fill_gradientn(colours = height.colors,
                                    values = scales::rescale(quantile(surfaceMat_df$value,height.quantiles,na.rm = TRUE)),
                                    breaks = function(lims){
                                      dat <- quantile(surfaceMat_df$value,legend.quantiles,na.rm = TRUE)

                                      dat <- dat %>%
                                        setNames(paste0(names(dat)," [",round(dat,3),"]"))

                                      return(dat)
                                    },
                                    na.value = na.value) +
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
        panel.background = ggplot2::element_blank()) +
      ggplot2::guides(fill = ggplot2::guide_colourbar(barheight = legendLength,
                                                      label.theme = ggplot2::element_text(size = 8),
                                                      title.theme = ggplot2::element_text(size = 10),
                                                      frame.colour = "black",
                                                      ticks.colour = "black"),
                      colour = 'none') +
      ggplot2::labs(fill = expression("Rel. Height ["*mu*"m]")) +
      ggplot2::facet_wrap(~ x3p)
  }
  else if(type == "list"){
    plts <- purrr::pmap(.l = list(x3pList,
                                  names(x3pList)),
                        function(x3p,name){

                          surfaceMat_df <- x3p %>%
                            x3ptools::x3p_to_df() %>%
                            #perform some transformations on the
                            #x,y values so that the plot is
                            #representative of the actual surface
                            #matrix (i.e., element [1,1] of the
                            #surface matrix is in the top-left
                            #corner)
                            dplyr::mutate(xnew = max(.data$y) - .data$y,
                                          ynew = max(.data$x) - .data$x,
                                          value = .data$value - median(.data$value,na.rm = TRUE)) %>%
                            dplyr::select(-c(.data$x,.data$y)) %>%
                            dplyr::rename(x=.data$xnew,
                                          y=.data$ynew) %>%
                            dplyr::mutate(value = .data$value - median(.data$value,na.rm = TRUE)) %>%
                            dplyr::mutate(x3p = rep(name,times = nrow(.)))

                          surfaceMat_df %>%
                            ggplot2::ggplot(ggplot2::aes(x = .data$x,y = .data$y)) +
                            ggplot2::geom_raster(ggplot2::aes(fill = .data$value))  +
                            ggplot2::scale_fill_gradientn(colours = height.colors,
                                                          values = scales::rescale(quantile(surfaceMat_df$value,height.quantiles,na.rm = TRUE)),
                                                          breaks = function(lims){
                                                            dat <- quantile(surfaceMat_df$value,legend.quantiles,na.rm = TRUE)

                                                            dat <- dat %>%
                                                              setNames(paste0(names(dat)," [",round(dat,3),"]"))

                                                            return(dat)
                                                          },
                                                          na.value = na.value) +
                            ggplot2::theme_minimal() +
                            ggplot2::coord_fixed(expand = FALSE) +
                            ggplot2::theme(
                              axis.title.x = ggplot2::element_blank(),
                              axis.text.x = ggplot2::element_blank(),
                              axis.ticks.x = ggplot2::element_blank(),
                              axis.title.y = ggplot2::element_blank(),
                              axis.text.y = ggplot2::element_blank(),
                              axis.ticks.y = ggplot2::element_blank(),
                              panel.grid.major = ggplot2::element_blank(),
                              panel.grid.minor = ggplot2::element_blank(),
                              panel.background = ggplot2::element_blank(),
                              plot.title = ggplot2::element_text(hjust = .5,
                                                                 size = 11)) +
                            ggplot2::guides(fill = ggplot2::guide_colourbar(barheight = legendLength,
                                                                            label.theme = ggplot2::element_text(size = 8),
                                                                            title.theme = ggplot2::element_text(size = 10),
                                                                            frame.colour = "black",
                                                                            ticks.colour = "black"),
                                            colour =  'none') +
                            ggplot2::labs(fill = expression("Rel. Height ["*mu*"m]")) +
                            ggplot2::labs(title = name)
                        })
  }
  return(plts)
}
