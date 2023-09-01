#'Plot a set of x3ps
#'@name x3pPlot
#'
#'@description Plots the surface matrices in a list of x3p objects. Either
#'  creates one plot faceted by surface matrix or creates individual plots per
#'  surface matrix and returns them in a list.
#'
#'@param ... one or more x3p objects
#'@param x3p.names character vector containing names of each x3p object
#'@param output dictates whether one plot faceted by surface matrix or a list of
#'  plots per surface matrix is returned. The faceted plot will have a
#'  consistent height scale across all surface matrices.
#'@param height.colors vector of colors to be passed to scale_fill_gradientn
#'  that dictates the height value colorscale
#'@param height.quantiles vector of quantiles associated with each color defined
#'  in the height.colors argument
#'@param na.value color to be used for NA values (passed to
#'  scale_fill_gradientn)
#'@param legend.quantiles vector of quantiles to be shown as tick marks on
#'  legend plot
#'@param legend.length length of the plot legend. Passed to the barheight
#'   argument of the ggplot2::guide_colorbar function
#'@return A ggplot object or list of ggplot objects showing the surface matrix
#'  height values.
#' @examples
#'data("K013sA1","K013sA2")
#'
#' x3pPlot(K013sA1,K013sA2,x3p.names = c("Scan A","Scan B"))
#'@export
#'
#'@importFrom stats setNames median quantile
#'@importFrom rlang .data

x3pPlot <- function(...,
                    x3p.names = NULL,
                    output = "faceted",
                    height.colors = c('#2d004b','#542788','#8073ac','#b2abd2','#d8daeb','#f7f7f7','#fee0b6','#fdb863','#e08214','#b35806','#7f3b08'),
                    height.quantiles = c(0,.01,.025,.1,.25,.5,.75,0.9,.975,.99,1),
                    legend.quantiles = c(0,.01,.25,.5,.75,.99,1),
                    na.value = "gray65",
                    legend.length = grid::unit(3,"in")){

  # input can either be a sequence of x3p objects  of arbitrary length OR a list
  # of x3p objects. We want the x3p objects bundled in a list either way
  x3pList <- list(...)

  # check if we have a list containing a list containing x3ps, or just a list of
  # x3ps. If the former, then flatten the list by one level
  if(!all(purrr::map_chr(x3pList,class) == "x3p")){

    x3pList <- purrr::flatten(x3pList)

  }

  stopifnot("Input should be either a sequence of x3p objects, 'x3p1,x3p2,...', or a list of x3p objects, 'list(x3p1,x3p2,...)'" = all(purrr::map_chr(x3pList,class) == "x3p"))

  if(is.null(names(x3pList))){
    if(is.null(x3p.names)){
      x3pList <- setNames(x3pList,paste0("x3p",1:length(x3pList)))
    }
    else{
      x3pList <- x3pList %>% purrr::set_names(x3p.names)
    }

  }

  if(output == "faceted"){
    surfaceMat_df <- purrr::pmap_dfr(.l = list(x3pList,
                                               names(x3pList)),
                                     function(x3p,name){

                                       x3p$header.info$incrementX <- 1
                                       x3p$header.info$incrementY <- 1
                                       x3p$mask <- NULL

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
                                                       value = value - median(value,na.rm = TRUE)) %>%
                                         dplyr::select(-c(x,y)) %>%
                                         dplyr::rename(x=xnew,
                                                       y=ynew) %>%
                                         dplyr::mutate(x3p = rep(name,times = nrow(.)))
                                     }) %>%
      dplyr::mutate(x3p = factor(x3p,levels = names(x3pList)))

    plts <- surfaceMat_df %>%
      ggplot2::ggplot(ggplot2::aes(x = x,y = y)) +
      ggplot2::geom_raster(ggplot2::aes(fill = value))  +
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
      ggplot2::guides(fill = ggplot2::guide_colourbar(barheight = legend.length,
                                                      label.theme = ggplot2::element_text(size = 8),
                                                      title.theme = ggplot2::element_text(size = 10),
                                                      frame.colour = "black",
                                                      ticks.colour = "black"),
                      colour = 'none') +
      ggplot2::labs(fill = expression("Rel. Height ["*mu*"m]")) +
      ggplot2::facet_wrap(~ x3p)
  }
  else if(output == "list"){
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
                            dplyr::mutate(xnew = max(y) - y,
                                          ynew = max(x) - x,
                                          value = value - median(value,na.rm = TRUE)) %>%
                            dplyr::select(-c(x,y)) %>%
                            dplyr::rename(x=xnew,
                                          y=ynew) %>%
                            dplyr::mutate(x3p = rep(name,times = nrow(.)))

                          surfaceMat_df %>%
                            ggplot2::ggplot(ggplot2::aes(x = x,y = y)) +
                            ggplot2::geom_raster(ggplot2::aes(fill = value))  +
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
                              panel.background = ggplot2::element_blank()) +
                            ggplot2::guides(fill = ggplot2::guide_colourbar(barheight = legend.length,
                                                                            label.theme = ggplot2::element_text(size = 8),
                                                                            title.theme = ggplot2::element_text(size = 10),
                                                                            frame.colour = "black",
                                                                            ticks.colour = "black"),
                                            colour =  'none') +
                            ggplot2::labs(fill = expression("Rel. Height ["*mu*"m]")) +
                            ggplot2::facet_wrap(~ x3p)
                        })
  }
  return(plts)
}

#' @rdname x3pPlot
#' @export

x3p_plot <- x3pPlot
