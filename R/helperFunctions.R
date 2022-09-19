# helper function for x3pListPlot. Rotates a surface matrix, but doesn't crop
# back to the original surface matrix's dimensions.
rotateSurfaceMatrix_noCrop <- function(surfaceMat,
                                       theta = 0,
                                       interpolation = 0){
  surfaceMatFake <- (surfaceMat*10^5) + 1 #scale and shift all non-NA pixels up 1 (meter)
  # imFakeRotated <- :bilinearInterpolation(imFake,theta)
  surfaceMatFakeRotated <- surfaceMatFake %>%
    imager::as.cimg() %>%
    imager::imrotate(angle = theta,
                     interpolation = interpolation, #linear interpolation,
                     boundary = 0) %>% #pad boundary with 0s (dirichlet condition)
    as.matrix()

  surfaceMatFakeRotated[surfaceMatFakeRotated == 0] <- NA
  #shift all of the legitimate pixels back down by 1:
  surfaceMatRotated <- (surfaceMatFakeRotated - 1)/(10^5)

  return(surfaceMatRotated)
}

# @name linear_to_matrix
# @param index integer vector of indices, must be between 1 and nrow*ncol
# @param nrow number of rows, integer value defaults to 7
# @param ncol  number of columns, integer value, defaults to number of rows
# @param byrow logical value, is linear index folded into matrix by row (default) or by column (`byrow=FALSE`).
# @examples
# index <- sample(nrow*ncol, 10, replace = TRUE)
# linear_to_matrix(index, nrow=4, ncol = 5, byrow=TRUE)
#
# @keywords internal
# @importFrom rlang .data

linear_to_matrix <- function(index, nrow = 7, ncol = nrow, byrow = TRUE, sep = ", ") {
  index <- as.integer(index)
  stopifnot(all(index <= nrow*ncol), all(index > 0))

  if (byrow) { # column is the fast index
    idx_out_col <- ((index-1) %% ncol) + 1
    idx_out_row <- ((index-1) %/% ncol) + 1
  } else { # row is the fast index
    idx_out_col <- ((index-1) %/% nrow) + 1
    idx_out_row <- ((index-1) %% nrow) + 1
  }
  paste0(idx_out_row, sep, idx_out_col)
}

targetCellCorners <- function(alignedTargetCell,cellIndex,theta,cmcClassif,target){

  targetScanRows <- alignedTargetCell$cmcR.info$regionIndices[c(3)] + alignedTargetCell$cmcR.info$regionRows - 1
  targetScanCols <- alignedTargetCell$cmcR.info$regionIndices[c(1)] + alignedTargetCell$cmcR.info$regionCols - 1

  rotatedMask <- rotateSurfaceMatrix(target$surface.matrix,theta)

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

  rotatedMask <- rotateSurfaceMatrix_noCrop(rotatedMask,theta = -1*theta)
  #make a copy that isn't going to have the target cell indices added so that we
  #know exactly how many rows/cols we need to translate everything to get back to
  #the original scan indices
  rotatedMaskCopy <- rotatedMask#rotateSurfaceMatrix_noCrop(rotatedMaskCopy,theta = 0)#-1*(-30))

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
    dplyr::mutate(xnew = .data$y,
                  ynew = .data$x) %>%
    dplyr::select(-c(.data$x,.data$y)) %>%
    dplyr::rename(x=.data$xnew,y=.data$ynew) %>%
    dplyr::mutate(x = .data$x - min(which(colSums(abs(rotatedMaskCopy),na.rm = TRUE) > 0)),
                  y = .data$y - min(which(rowSums(abs(rotatedMaskCopy),na.rm = TRUE) > 0)),
                  y = nrow(target$surface.matrix) - .data$y
    ) %>%
    dplyr::filter(.data$value == 100) %>%
    dplyr::select(-.data$value) %>%
    dplyr::group_by(.data$x,.data$y) %>%
    dplyr::distinct() %>%
    dplyr::mutate(cellIndex = cellIndex,
                  theta = theta,
                  cmcClassif = cmcClassif)

  return(ret)

}
