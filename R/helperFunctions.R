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

targetCellCorners <- function(alignedTargetCell,cellIndex,theta,cmcClassif,target){

  targetScanRows <- alignedTargetCell$cmcR.info$regionIndices[c(3)] + alignedTargetCell$cmcR.info$regionRows - 1
  targetScanCols <- alignedTargetCell$cmcR.info$regionIndices[c(1)] + alignedTargetCell$cmcR.info$regionCols - 1

  rotatedMask <- cmcR:::rotateSurfaceMatrix(target$surface.matrix,theta)

  rowPad <- 0
  colPad <- 0

  if(targetScanRows[1] <= 0){

    rowPad <- abs(targetScanRows[1]) + 1

    rotatedMask <- rbind(matrix(NA,nrow = rowPad,ncol = ncol(rotatedMask)),
                         rotatedMask)

    targetScanRows <- targetScanRows + rowPad
  }

  if(targetScanCols[1] <= 0){

    colPad <- abs(targetScanCols[1]) + 1

    rotatedMask <- cbind(matrix(NA,nrow = nrow(rotatedMask),ncol = colPad),
                         rotatedMask)

    targetScanCols <- targetScanCols + colPad
  }

  if(targetScanRows[2] > nrow(rotatedMask)){

    rowPad <- targetScanRows[2] - nrow(rotatedMask)

    rotatedMask <- rbind(rotatedMask,
                         matrix(NA,nrow = rowPad,ncol = ncol(rotatedMask)))

  }

  if(targetScanCols[2] > ncol(rotatedMask)){

    colPad <- targetScanCols[2] - ncol(rotatedMask)

    rotatedMask <- cbind(rotatedMask,
                         matrix(NA,nrow = nrow(rotatedMask),ncol = colPad))

  }

  rotatedMask[targetScanRows[1]:targetScanRows[2],targetScanCols[1]:targetScanCols[2]] <- 100

  rotatedMask <- cmcR:::rotateSurfaceMatrix_noCrop(rotatedMask,theta = -1*theta)
  #make a copy that isn't going to have the target cell indices added so that we
  #know exactly how many rows/cols we need to translate everything to get back to
  #the original scan indices
  rotatedMaskCopy <- rotatedMask#cmcR:::rotateSurfaceMatrix_noCrop(rotatedMaskCopy,theta = 0)#-1*(-30))

  rotatedMaskCopy[rotatedMaskCopy == 100] <- NA

  newColPad <- 0
  newRowPad <- 0
  if(theta != 0){

    # the cells that are rotated may have been "shifted" due to the cropping
    # performed above relative to the unrotated target scan's indices -- for
    # example, padding the rotatedMask to the left requires a correction after
    # rotating. Unfortunately, because the cells are rotated, the padding done in
    # the rotated domain doesn't come out to be nice (dRow,dCol) translations in
    # the unrotated domain. We need to perform some trig to determine what
    # (dRow,dCol) in the rotated domain translates to in the unrotated domain.

    # In the rotated domain:
    #    ------------* <<- the location of target cell after padding in rotatedMask
    #    |    ^dx^
    #    |
    #    | <- dy
    #    |
    #    * <<- where the target cell *should* be relative to the rotated target's indices
    #

    # no consider rotating this whole space, then the dx, dy will be "tilted" by
    # some theta and we will need to calculate via trig what the correct dx', dy'
    # are in the original, unrotated domain. Draw a diagram with a rotated dx, dy
    # by some theta and draw a straight line between the two * above and you
    # should be able to work out the following formulas again

    #TODO: pay attention to the necessary signs (up/down/left/right) of the
    #corrections below
    psi <- atan2(rowPad,colPad)
    hyp <- sqrt(colPad^2 + rowPad^2)
    phi <- pi/2 - (theta*pi/180 + psi)

    newColPad <- sin(phi)*hyp
    newRowPad <- cos(phi)*hyp

  }

  ret <- rotatedMask %>%
    imager::as.cimg() %>%
    as.data.frame() %>%
    mutate(xnew = y,
           ynew = x) %>%
    select(-c(x,y)) %>%
    rename(x=xnew,y=ynew) %>%
    mutate(x = x - min(which(colSums(rotatedMaskCopy,na.rm = TRUE) > 0)),
           # x = x - newColPad,
           y = y - min(which(rowSums(rotatedMaskCopy,na.rm = TRUE) > 0)),
           # y = y - newRowPad,
           # y = max(y) - y
           y = nrow(target$surface.matrix) - y
    ) %>%
    filter(value == 100) %>%
    select(-value) %>%
    group_by(x,y) %>%
    distinct() %>%
    mutate(cellIndex = cellIndex,
           theta = theta,
           cmcClassif = cmcClassif)

  return(ret)

}
