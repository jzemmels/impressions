labelBlobs <- function(reference,target,filterCutoff = 1){

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
    mutate(alpha = ifelse(valueDiff <= filterCutoff,1,0),
           x3pName = "Pixelwise Average") %>%
    select(-valueDiff) %>%
    filter(alpha == 0)

  ret <- matrix(0,nrow = nrow(reference$surface.matrix),ncol = ncol(reference$surface.matrix))

  for(rowInd in 1:nrow(surfaceMat_df)){
    # for some reason the x,y indices here are off by 1 index each from what
    # indices are actually filtered in the original image? Not quite sure why
    # this is the case.
    ret[surfaceMat_df[rowInd,]$y+1,surfaceMat_df[rowInd,]$x+1] <- 1

  }

  contours <- skimage$measure$find_contours(ret,0,fully_connected="high")

  #create a matrix to hold the contours
  dat <- matrix(-3L,nrow = nrow(reference$surface.matrix),ncol = ncol(reference$surface.matrix))

  polys <- map(contours,
               function(cont){

                 rows <- as.integer(cont[,1])
                 cols <- as.integer(cont[,2])

                 ret <- skimage$draw$polygon(rows,cols)

                 return(ret)
               })



  for(polyInd in 1:length(polys)){
    dat[matrix(c(polys[[polyInd]][[1]],polys[[polyInd]][[2]]),ncol = 2) + 1L] <- polyInd
  }

  ret1 <- dat %>%
    imager::as.cimg() %>%
    as.data.frame() %>%
    # cimg labels the x and y image axes incorrectly, so this fixes it
    mutate(xnew = y,
           ynew = x) %>%
    select(-c(x,y)) %>%
    rename(x = xnew,y=ynew) %>%
    arrange(x,y) %>%
    # some edge pixels may not have been flooded with the polygon algorithm
    # above, so we need to label the pixels that need to be associated with a
    # particular blob with -1
    left_join(refTargAverage %>%
                mutate(notFiltered = ifelse(valueDiff <= 1,TRUE,FALSE)) %>%
                mutate(x = x + 1,
                       y = y + 1) %>%
                select(x,y,notFiltered),
              by = c("x","y")) %>%
    mutate(value = ifelse(notFiltered,-1,value)) %>%
    select(-notFiltered) %>%
    rename(blobNum = value) %>%
    mutate(comparisonName = ..4,
           cellIndex = ..3)

  blobNums <- ret1 %>%
    filter(blobNum > 0) %>%
    pull(blobNum) %>%
    unique()

  needABlob <- ret1 %>%
    filter(blobNum == -3)

  # the pixels labeled -3 are either neighbors of a labeled blob, and
  # therefore should be a part of that blob, or should be a blob by
  # themselves. it's common for the -3 labeled pixels to be on the edge of the
  # image. the following code will progressively fill-in the -3 labeled pixels
  # with the neighboring blob label (if there is one) until there are no
  # additional labels added (the if statement at the end of the while loop).

  if(nrow(needABlob) > 0){

    while(any(needABlob$blobNum == -3)){

      # make a copy to check at the end of the while loop if anything updated.
      ret2 <- ret1

      # for each labeled blob, we want to determine if any neighboring pixels
      # (within 1 pixel in the x or y direction) has value -1, in which case it
      # needs to be included
      for(blob in blobNums){

        # filters down to the pixels belonging to a particular blob
        blobPixels <- ret1 %>%
          filter(blobNum == blob) %>%
          select(x,y)

        # for each pixel in the blob...
        for(rowInd in 1:nrow(blobPixels)){

          # if that pixel is within the neighborhood of a pixel needing a blob,
          # then mark that blob-needing pixel for later flood filling
          needABlobNeighbor <- (abs(needABlob$x - blobPixels[rowInd,]$x) <= 1) &
            (abs(needABlob$y - blobPixels[rowInd,]$y) <= 1)

          if(any(needABlobNeighbor)){

            # identify the indices that need updating
            rowsToUpdate <- needABlob$x[needABlobNeighbor]
            colsToUpdate <- needABlob$y[needABlobNeighbor]

            #update the indices directly in ret1
            ret1 <- ret1 %>%
              mutate(blobNum = ifelse(paste0(x,",",y) %in% paste0(rowsToUpdate,",",colsToUpdate),
                                      blob,blobNum))

            # remove the newly labeled indices from the needABlob data frame
            needABlob <- needABlob %>%
              filter(!(paste0(x,",",y) %in% paste0(rowsToUpdate,",",colsToUpdate)))


          }

        }

      }

      # the while loop will end on the loop where nothing is updated.
      if(all(ret2$blobNum == ret1$blobNum,na.rm = TRUE)){
        break
      }

    }
  }

  # there are some edge pixels labeled with -3 that are in their own island
  # and do not have a neighboring blob. This will create a new blob label for
  # those pixels.
  ret1 <- ret1 %>%
    mutate(blobNum = ifelse(blobNum == -3,max(blobNums) + 1,blobNum))

  return(ret1)

}
