# Helper functions in this chunk

# this function is mainly intended to transform an x3p object to a data.frame
# for easy ggplot-ing
x3pToDF <- function(x3p,preserveResolution = FALSE){

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

# combine a named list of x3ps into similarities, differences, or missing values

# For type = "missing" returns an x3p object containing a surface matrix with
# labeled pixels. The labels correspond to which x3p object has a missing value
# in that pixel index.

combineX3Ps <- function(x3p1,x3p2,
                        x3pNames = c("x3p1","x3p2"),
                        type = "difference",
                        quantileThresh = .5){

  stopifnot("Surface matrix dimensions should be equal" =
              all(dim(x3p1$surface.matrix) == dim(x3p2$surface.matrix)))

  #TODO: make this work with a list of > 2 x3ps
  # if(is.null(x3pList)){
  #   x3pList <- set_names(x3pList,paste0("x3p",1:length(x3pList)))
  # }

  if (type == "difference"){

    difference <- x3p1

    difference$surface.matrix <- x3p1$surface.matrix - x3p2$surface.matrix

    return(difference)

  }
  if(type == "missing"){

    missingVals <- x3p1

    missingVals$surface.matrix <- matrix(NA,nrow = nrow(missingVals$surface.matrix),
                                         ncol = ncol(missingVals$surface.matrix))

    missingVals$surface.matrix[!is.na(x3p1$surface.matrix) & !is.na(x3p2$surface.matrix)] <- "Neither"

    missingVals$surface.matrix[!is.na(x3p1$surface.matrix) & is.na(x3p2$surface.matrix)] <- x3pNames[1]

    missingVals$surface.matrix[is.na(x3p1$surface.matrix) & !is.na(x3p2$surface.matrix)] <- x3pNames[2]

    missingVals$surface.matrix[is.na(x3p1$surface.matrix) & is.na(x3p2$surface.matrix)] <- "Both"

    return(missingVals)

  }

}
