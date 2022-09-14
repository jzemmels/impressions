# returns the filter boundary outlines in a data frame
filterBoundaries <- function(averageBinarized){

  suppressWarnings({

    averageMat <- averageBinarized %>%
      dplyr::mutate(x = x+1,
                    y = y+1) %>%
      as.data.frame() %>%
      dplyr::select(x,y,value) %>%
      imager::as.cimg() %>%
      as.matrix()

  })

  averageMat[is.na(averageMat)] <- 0

  # we pad the matrix so that the contours one the edge blobs are properly
  # identified. the padding is removed in the last lines of the creation of
  # the outline object below
  averageMat  <- averageMat %>%
    imager::as.cimg() %>%
    imager::pad(nPix = 10,axes = "xy",val = 0)

  labels <- imager::label(averageMat)

  bounds <- purrr::map(unique(labels[labels > 0]),
                       function(lab){

                         imager::boundary(labels == lab)

                       })

  blobBoundaries <- list(bounds,labels)

  # combine all labeled blobs into one image
  boundaryPx <- Reduce("+",blobBoundaries[[1]] %>%
                         purrr::map(as.matrix)) %>%
    imager::as.cimg()

  # the mask used to dilate the blobs will grow them towards the bottom-right of
  # the matrix
  dilatedPx <- imager::dilate_rect(boundaryPx,sx = 2,sy = 2)
  dilatedPx_labels <- imager::dilate_rect(blobBoundaries[[2]],sx = 2,sy = 2)

  # flip the image and re-apply the dilation to grow the borders to the other
  # corners. flip back after dilation
  dilatedPx_mirrorx <- imager::mirror(imager::dilate_rect(imager::mirror(boundaryPx,axis="x"),sx = 2,sy = 2),axis="x")
  dilatedPx_mirrorx_labels <- imager::mirror(imager::dilate_rect(imager::mirror(blobBoundaries[[2]],axis="x"),sx = 2,sy = 2),axis="x")

  dilatedPx_mirrory <- imager::mirror(imager::dilate_rect(imager::mirror(boundaryPx,axis="y"),sx = 2,sy = 2),"y")
  dilatedPx_mirrory_labels <- imager::mirror(imager::dilate_rect(imager::mirror(blobBoundaries[[2]],axis="y"),sx = 2,sy = 2),"y")

  dilatedPx_mirrorxy <- imager::mirror(imager::dilate_rect(imager::mirror(boundaryPx,axis="xy"),sx = 3,sy = 3),"xy")
  dilatedPx_mirrorxy_labels <- imager::mirror(imager::dilate_rect(imager::mirror(blobBoundaries[[2]],axis="xy"),sx = 3,sy = 3),"xy")

  # combine all of the dilated images together into one image
  dilatedPx_comb <- dilatedPx + dilatedPx_mirrorx + dilatedPx_mirrory + dilatedPx_mirrorxy

  # we just want a binary labeling
  dilatedPx_comb[dilatedPx_comb > 0] <- 1

  # the dilated boundaries will have also grown into the blobs, so we take those
  # pixels out
  dilatedPx_comb[blobBoundaries[[2]] > 0] <- 0

  # from: https://stackoverflow.com/questions/34756755/plot-outline-around-raster-cells
  outline <- dilatedPx_comb %>%
    as.data.frame() %>%
    dplyr::filter(value > 0) %>%
    dplyr::mutate(x = x-1,
                  y = y-1) %>%
    raster::rasterFromXYZ() %>%
    raster::rasterToPolygons(dissolve = TRUE) %>%
    ggplot2::fortify() %>%
    #the boundaries around the filtered blobs all share a common value in the
    #"hole" column of TRUE
    dplyr::filter(hole) %>%
    # remove padding used previously
    dplyr::mutate(lat = lat-5,
                  long = long-5)

  return(outline)

}
