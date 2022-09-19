#'Transform x3p surface matrix values
#'
#'Transform the values of a x3p surface matrix
#'
#'@param ... one or more x3p objects
#'@param x3p an x3p object
#'@param x3p1 an x3p object
#'@param x3p2 another x3p object
#'@param cond a Boolean function whose first argument is `x`.
#'@param replacement value to replace each element for which cond returns FALSE
#'@param na.rm logical. Should missing values be removed?
#'@param croppingThresh minimum number of non-NA pixels that need to be in a
#'  row/column for it to not be cropped out of the surface matrix
#'@param mask_vals a hexidecimal color value that corresponds to indices in a
#'  mask whose indices are to be replaced in the associated surface matrix
#'@param preserveResolution a Boolean dictating whether the scan resolution is
#'  preserved in the returned data frame. If FALSE, then the x,y data frame
#'  columns will be integer-valued. Otherwise, the difference between
#'  consecutive x,y values will be equal to the scan resolution.
#'
#'@description
#'`x3p_elemAverage()` calculate the element-wise average between two surface
#'matrices (of the same dimensions)
#'
#'`x3p_sd()` Calculate the standard deviation of x3p surface values
#'
#'`x3p_filter()` replace values of a surface matrix based on an element-wise
#'conditional function
#'
#'`x3p_cropWS()` Crop rows/columns of missing values around an x3p
#'
#'`x3p_to_dataFrame()` Convert an x3p object to a data frame
#'
#' @examples
#' data("K013sA1","K013sA2")
#'
#' # calculates the sd for a single x3p's surface values
#' x3p_sd(K013sA1)
#'
#' # calculates the sd for the joint surface values for two x3ps
#' x3p_sd(K013sA1,K013sA2)
#'
#' # calculate optimal alignment between the two x3ps
#' K013sA2_aligned <- cmcR::comparison_allTogether(K013sA1,K013sA2,theta = -3,
#'                                                 returnX3Ps = TRUE,numCells = c(1,1),
#'                                                 maxMissingProp = .99)$alignedTargetCell[[1]]
#'
#' averaged <- x3p_elemAverage(K013sA1,K013sA2_aligned)
#'
#' # this will replace values that are larger (in magnitude) than one standard
#' # deviation of the input x3p's surface values with NA:
#' filtered1 <- x3p_filter(K013sA1,
#'                         cond = function(x,thresh) x < thresh,
#'                         thresh = x3p_sd(K013sA1))
#'
#' # this will replace all surface matrix values between -1 and 1 with 0
#' filtered2 <- x3p_filter(K013sA1,cond = function(x) abs(x) > 1,replacement = 0)
#'
#' # exaggerated cropping for the sake of an example
#' cropped <- x3p_cropWS(K013sA1,croppingThresh = 100)
#'
#' x3pPlot(K013sA1,K013sA2_aligned,averaged,filtered1,filtered2,cropped,
#'         x3pNames = c("K01sA1","K013sA2 Aligned","Averaged","Filtered1","Filtered2","Cropped"))
#'
#'@rdname x3pOperations
#'@export
x3p_elemAverage <- function(x3p1,x3p2){

  stopifnot(all(dim(x3p1$surface.matrix) == dim(x3p2$surface.matrix)))

  if(abs(x3p1$header.info$incrementX - x3p2$header.info$incrementX) > 1e-8 |
     abs(x3p1$header.info$incrementY - x3p2$header.info$incrementY) > 1e-8){

    warning("x3p resolutions are not equal - may lead to unexpected output")

  }

  ret <- x3p1

  ret$surface.matrix <- (x3p1$surface.matrix + x3p2$surface.matrix)/2

  return(ret)

}

#' @rdname x3pOperations
#' @export
x3p_sd <- function(...,na.rm = TRUE){

  list(...) %>%
    purrr::map(~ c(.$surface.matrix)) %>%
    unlist() %>%
    sd(na.rm = na.rm)

}

#'@rdname x3pOperations
#'@export
x3p_filter <- function(x3p,cond,replacement = NA,...){

  filteredDF <- x3p %>%
    x3p_to_dataFrame(preserveResolution = FALSE) %>%
    dplyr::mutate(value = ifelse(do.call(cond,args = alist(x = value,...)),
                                 value,
                                 replacement))

  x3p$surface.matrix <- matrix(filteredDF$value,
                               nrow = nrow(x3p$surface.matrix),
                               ncol = ncol(x3p$surface.matrix))

  return(x3p)

}

#'@rdname x3pOperations
#'@export
x3p_cropWS <- function (x3p,
                        croppingThresh = 1){

  surfaceMat <- x3p$surface.matrix

  # count the number of NAs in each row/column
  colSum <-
    surfaceMat[(nrow(surfaceMat)/2 - .5 * nrow(surfaceMat)):(
      nrow(surfaceMat)/2 + .5 * nrow(surfaceMat)), ] %>%
    is.na() %>%
    magrittr::not() %>%
    colSums()

  rowSum <-
    surfaceMat[, (ncol(surfaceMat)/2 - .5 * ncol(surfaceMat)):(
      ncol(surfaceMat)/2 + .5 * ncol(surfaceMat))] %>%
    is.na() %>%
    magrittr::not() %>%
    rowSums()

  # crop to rows/columns that have at least croppingThresh non-missing values
  surfaceMatCropped <-
    surfaceMat[min(which(rowSum >= croppingThresh)):max(
      which(rowSum >= croppingThresh)),
      min(which(colSum >= croppingThresh)):max(
        which(colSum >= croppingThresh))]

  x3p$surface.matrix <- surfaceMatCropped
  x3p$header.info$sizeX <- nrow(surfaceMatCropped)
  x3p$header.info$sizeY <- ncol(surfaceMatCropped)

  return(x3p)
}

#'@rdname x3pOperations
#'@export
x3p_to_dataFrame <- function(x3p,preserveResolution = FALSE){

  if(!preserveResolution){

    x3p$header.info$incrementX <- 1
    x3p$header.info$incrementY <- 1

  }

  ret <- x3p %>%
    x3ptools::x3p_to_df() %>%
    #perform some transformations on the x,y values so that the plot is
    #representative of the actual surface matrix (i.e., element [1,1] of the
    #surface matrix is in the top-left corner)
    dplyr::mutate(xnew = max(y) - y,
                  ynew = max(x) - x) %>%
    dplyr::select(-c(x,y)) %>%
    dplyr::rename(x=xnew,
                  y=ynew)

  return(ret)
}
