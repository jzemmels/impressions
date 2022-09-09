#' Calculate the element-wise average between two x3ps
#' @name x3p_elemAverage
#'
#' @description Calculates the element-wise average between the surface values
#'   of two x3ps (of the same dimension)
#'
#' @param x3p1 an x3p object
#' @param x3p2 another x3p object
#'
#' @return An x3p object containing the element-wise average between the surface
#'   values of x3p1 and x3p1
#'
#' @export
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

#' Calculate the standard deviation of x3p surface values
#' @name x3p_sd
#'
#' @description Calculates the standard deviation of all surface values of one
#'   or more x3p objects
#'
#' @param ... one or more x3p objects
#' @param na.rm logical. Should missing values be removed?
#'
#' @return The standard deviation of the input surface matrices
#'
#' @examples
#' # calculates the sd for a single x3p's surface values
#' x3p_sd(x3p)
#'
#' # calculates the sd for the joint surface values for two x3ps
#' x3p_sd(x3p1,x3p2)
#'
#' @export
x3p_sd <- function(...,na.rm = TRUE){

  list(...) %>%
    purrr::map(~ c(.$surface.matrix)) %>%
    unlist() %>%
    sd(na.rm = na.rm)

}

#' Filter a surface matrix by an element-wise conditional function
#' @name x3p_filter
#'
#' @description Filters the elements of an x3p surface matrix by replacing
#'   values based on an element-wise conditional (Boolean) function
#'
#' @param x3p an x3p object
#' @param cond a Boolean function whose first argument is `x`.
#' @param replacement value to replace each element for which cond returns FALSE
#' @param ... additional arguments for the cond function
#'
#' @note The value returned by the cond function should be the same length as
#'   the total number of elements in the surface matrix (e.g., use vectorized
#'   operations)
#'
#' @return An x3p object containing a filtered version of the input x3p's
#'   surface matrix
#'
#' @examples
#'
#' # this will replace values that are larger than one standard deviation of
#' # the input x3p's surface values with NA:
#' x3p_filter(x3p,cond = function(x,thresh) x < thresh,thresh = x3p_sd(x3p))
#'
#' # this will replace all surface matrix values between -1 and 1 with 0
#' x3p_filter(x3p,cond = function(x) abs(x) > 1,replacement = 0)
#'
#' @export
x3p_filter <- function(x3p,cond,replacement = NA,...){

  filteredDF <- x3p %>%
    x3pToDF(preserveResolution = FALSE) %>%
    dplyr::mutate(value = ifelse(do.call(cond,args = alist(x = value,...)),
                                 value,
                                 replacement))

  x3p$surface.matrix <- matrix(filteredDF$value,
                               nrow = nrow(x3p$surface.matrix),
                               ncol = ncol(x3p$surface.matrix))

  return(x3p)
  # %>%
  #   x3ptools::df_to_x3p() %>%
  #   x3ptools::x3p_flip_y() %>%
  #   x3ptools::x3p_rotate(angle = 270)

}

#' Crop rows/columns of missing values around an x3p
#'
#' @name preProcess_cropWS
#'
#' @param x3p an x3p object containing a surface matrix
#' @param croppingThresh minimum number of non-NA pixels that need to be in a
#'   row/column for it to not be cropped out of the surface matrix
#' @param crop
#'
#' @return a surface matrix with outer rows/columns removed depending on
#'   croppingThresh
#'
#' @examples
#' \dontrun{
#' raw_x3p <- x3ptools::read_x3p("path/to/file.x3p") %>%
#'   x3ptools::sample_x3p(m = 2)
#'
#' raw_x3p$surface.matrix <- raw_x3p$surface.matrix %>%
#'   cmcR::preProcess_ransacLevel() %>%
#'   cmcR::preProcess_levelBF() %>%
#'   cmcR::preProcess_cropWS(croppingThresh = 2)
#' }
#'
#' @export

x3p_cropWS <- function (x3p,
                        croppingThresh = 1,
                        # croppingProp = 0.5,
                        # robust = FALSE,
                        ...){

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

#' Replace elements of a surface matrix based on a mask
#'
#' @name x3p_delete
#'
#' @param x3p an x3p object
#' @param mask_vals a hexidecimal color value that corresponds to indices in a
#'   mask whose indices are to be replaced in the associated surface
#'   matrix
#' @param replacement value used replace indices in the surface matrix
#'
#' @export

x3p_delete <- function(x3p, mask_vals, replacement = NA) {

  idx <- which(t(as.matrix(x3p$mask)) %in% mask_vals)
  x3p$surface.matrix[idx] <- replacement
  x3p$mask <- NULL

  x3p <- x3p_cropWS(x3p)

  return(x3p)
}
