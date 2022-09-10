#'Plot Congruent Matching Cells results for a pair of cartridge cases.
#'@name cmcPlot
#'
#'@param reference the scan that is partitioned into a grid of cells
#'@param target the scan to which each reference cell is compared during the
#'  cell-based comparison procedure
#'@param cmcClassifs a data frame containing columns cellHeightValues,
#'  alignedTargetCell, cellIndex, theta, and user-defined cmcCol & corrCol
#'@param type the form of the returned plot object(s). Either "faceted," meaning
#'  the reference and target plot will be shown side-by-side or "list" meaning
#'  each element of the plot (referece, target, and legend) will be returned
#'  separately as elements of a list
#'@param cmcCol name of column containing CMC classifications as returned by the
#'  decision_CMC function. Defaults to "originalMethod"
#'@param corrCol name of column containing correlation values for each cell.
#'  Defaults to "pairwiseCompCor," but "fft_ccf" is a common alternative.
#'@export
#'@importFrom patchwork wrap_plots
#'@importFrom ggplotify as.ggplot
cmcPlot <- function(reference,
                    target,
                    cmcClassifs,
                    type = "faceted",
                    cmcCol = "originalMethod",
                    corrCol = "pairwiseCompCor"){

  #check that the necessary columns are in cmcClassifs

  stopifnot("Make sure that there is a column called 'cellHeightValues' that is the result of the comparison_alignedTargetCell() function." = any(stringr::str_detect(names(cmcClassifs),"cellHeightValues")))

  stopifnot("Make sure that there is a column called 'alignedTargetCell' that is the result of the comparison_alignedTargetCell() function." = any(stringr::str_detect(names(cmcClassifs),"alignedTargetCell")))

  stopifnot("Make sure that there is a column called 'cellIndex'" = any(stringr::str_detect(names(cmcClassifs),"cellIndex")))

  stopifnot("Make sure that there is a column called 'theta'" = any(stringr::str_detect(names(cmcClassifs),"theta")))

  stopifnot(any(stringr::str_detect(names(cmcClassifs),cmcCol)))

  stopifnot(any(stringr::str_detect(names(cmcClassifs),corrCol)))

  # get the indices for the necessary columns
  referenceCellCol <- which(stringr::str_detect(names(cmcClassifs),"cellHeightValues"))

  targetCellCol <- which(stringr::str_detect(names(cmcClassifs),"alignedTargetCell"))

  cellIndexCol <- which(stringr::str_detect(names(cmcClassifs),"cellIndex"))

  thetaCol <- which(stringr::str_detect(names(cmcClassifs),"theta"))

  cmcIndexCol <- which(stringr::str_detect(names(cmcClassifs),cmcCol))

  # cmcClassifs <- cmcClassifs %>%
  #   dplyr::group_by(cellIndex) %>%
  #   dplyr::filter(!!as.name(corrCol) == max(!!as.name(corrCol)))

  targetCellData <- cmcClassifs %>%
    dplyr::select(c(targetCellCol,cellIndexCol,thetaCol,cmcIndexCol)) %>%
    purrr::pmap_dfr(~ targetCellCorners(alignedTargetCell = ..1,
                                        cellIndex = ..2,
                                        theta = ..3,
                                        cmcClassif = ..4,
                                        target = target))

  referenceCells <- cmcClassifs %>%
    dplyr::pull(referenceCellCol)

  cellData <- cmcClassifs %>%
    dplyr::select(c(cellIndexCol,referenceCellCol,cmcIndexCol)) %>%
    purrr::pmap_dfr(~ {

      cellInds <- ..2$cmcR.info$cellRange %>%
        stringr::str_remove("rows: ") %>%
        stringr::str_remove("cols: ") %>%
        stringr::str_split(pattern = ", ")

      cellInds_rows <- stringr::str_split(cellInds[[1]][1]," - ")[[1]]
      cellInds_cols <- stringr::str_split(cellInds[[1]][2]," - ")[[1]]

      return(data.frame(rowStart = as.numeric(cellInds_rows[1]),
                        rowEnd = as.numeric(cellInds_rows[2]),
                        colStart = as.numeric(cellInds_cols[1]),
                        colEnd = as.numeric(cellInds_cols[2])) %>%
               dplyr::mutate(cellIndex = ..1,
                             cmcClassif = ..3))

    }) %>%
    dplyr::mutate(rowStart = max(.data$rowEnd) - .data$rowStart,
                  rowEnd = max(.data$rowEnd) - .data$rowEnd,
                  colMean = purrr::map2_dbl(.data$colStart,.data$colEnd,~ mean(c(.x,.y))),
                  rowMean = purrr::map2_dbl(.data$rowStart,.data$rowEnd,~ mean(c(.x,.y))))

  # ggplot2 complains about the guides
  suppressWarnings({

    plt <- x3pListPlot(list("target" = target),
                       height.colors =
                         c('#1B1B1B','#404040','#7B7B7B','#B0B0B0','#DBDBDB','#F7F7F7','#E4E4E4','#C5C5C5','#999999','#717171','#4E4E4E')) +
      ggplot2::theme(legend.position = "none")

    if(all(targetCellData$cmcClassif == "CMC")){

      refPlt <- x3pListPlot(list("reference" = reference),
                            height.colors =
                              c('#1B1B1B','#404040','#7B7B7B','#B0B0B0','#DBDBDB','#F7F7F7','#E4E4E4','#C5C5C5','#999999','#717171','#4E4E4E')) +
        ggplot2::guides(fill = "none") +
        ggnewscale::new_scale_fill() +
        ggplot2::geom_rect(data = cellData,
                           ggplot2::aes(xmin = .data$colStart,xmax = .data$colEnd,ymin = .data$rowStart,ymax = .data$rowEnd,fill = .data$cmcClassif),
                           alpha = .2,
                           inherit.aes = FALSE) +
        ggplot2::scale_fill_manual(values = c("#313695")) +
        ggplot2::geom_text(data = cellData,
                           ggplot2::aes(x = .data$colMean,y = .data$rowMean,label = .data$cellIndex),inherit.aes = FALSE) +
        ggplot2::guides(fill = ggplot2::guide_legend(order = 1)) +
        ggplot2::theme(
          legend.direction = "horizontal"
        ) +
        ggplot2::labs(fill = "CMC Classif.")

      plt <- plt +
        ggnewscale::new_scale_fill() +
        ggplot2::geom_raster(data = targetCellData,
                             ggplot2::aes(x = .data$x,y = .data$y,fill = .data$cmcClassif),
                             alpha = .2) +
        ggplot2::scale_fill_manual(values = c("#313695")) +
        ggplot2::geom_text(data = targetCellData %>%
                             dplyr::group_by(.data$cellIndex) %>%
                             dplyr::summarise(x = mean(.data$x),
                                              y = mean(.data$y),
                                              theta = unique(.data$theta)),
                           ggplot2::aes(x=.data$x,y=.data$y,label = .data$cellIndex,angle = -1*.data$theta))

    }
    else if(all(targetCellData$cmcClassif == "non-CMC")){

      refPlt <- x3pListPlot(list("reference" = reference),
                            height.colors =
                              c('#1B1B1B','#404040','#7B7B7B','#B0B0B0','#DBDBDB','#F7F7F7','#E4E4E4','#C5C5C5','#999999','#717171','#4E4E4E')) +
        ggplot2::guides(fill = "none") +
        ggnewscale::new_scale_fill() +
        ggplot2::geom_rect(data = cellData,
                           ggplot2::aes(xmin = .data$colStart,xmax = .data$colEnd,ymin = .data$rowStart,ymax = .data$rowEnd,fill = .data$cmcClassif),
                           alpha = .2,
                           inherit.aes = FALSE) +
        ggplot2::scale_fill_manual(values = c("#a50026")) +
        ggplot2::geom_text(data = cellData,
                           ggplot2::aes(x = .data$colMean,y = .data$rowMean,label = .data$cellIndex),inherit.aes = FALSE) +
        ggplot2::guides(fill = ggplot2::guide_legend(order = 1)) +
        ggplot2::theme(
          legend.direction = "horizontal"
        ) +
        ggplot2::labs(fill = "CMC Classif.")

      plt <- plt +
        ggnewscale::new_scale_fill() +
        ggplot2::geom_raster(data = targetCellData,
                             ggplot2::aes(x = .data$x,y = .data$y,fill = .data$cmcClassif),
                             alpha = .2) +
        ggplot2::scale_fill_manual(values = c("#a50026")) +
        ggplot2::geom_text(data = targetCellData %>%
                             dplyr::group_by(.data$cellIndex) %>%
                             dplyr::summarise(x = mean(.data$x),
                                              y = mean(.data$y),
                                              theta = unique(.data$theta)),
                           ggplot2::aes(x=.data$x,y=.data$y,label = .data$cellIndex,angle = -1*.data$theta))

    }
    else{
      refPlt <- x3pListPlot(list("reference" = reference),
                            height.colors =
                              c('#1B1B1B','#404040','#7B7B7B','#B0B0B0','#DBDBDB','#F7F7F7','#E4E4E4','#C5C5C5','#999999','#717171','#4E4E4E')) +
        ggplot2::guides(fill = "none") +
        ggnewscale::new_scale_fill() +
        ggplot2::geom_rect(data = cellData,
                           ggplot2::aes(xmin = .data$colStart,xmax = .data$colEnd,ymin = .data$rowStart,ymax = .data$rowEnd,fill = .data$cmcClassif),
                           alpha = .2,
                           inherit.aes = FALSE) +
        ggplot2::scale_fill_manual(values = c("#313695","#a50026")) +
        ggplot2::geom_text(data = cellData,
                           ggplot2::aes(x = .data$colMean,y = .data$rowMean,label = .data$cellIndex),inherit.aes = FALSE) +
        ggplot2::guides(fill = ggplot2::guide_legend(order = 1)) +
        ggplot2::theme(
          legend.direction = "horizontal"
        ) +
        ggplot2::labs(fill = "CMC Classif.")

      plt <- plt +
        ggnewscale::new_scale_fill() +
        ggplot2::geom_raster(data = targetCellData,
                             ggplot2::aes(x = .data$x,y = .data$y,fill = .data$cmcClassif),
                             alpha = .2) +
        ggplot2::scale_fill_manual(values = c("#313695","#a50026")) +
        ggplot2::geom_text(data = targetCellData %>%
                             dplyr::group_by(.data$cellIndex) %>%
                             dplyr::summarise(x = mean(.data$x),
                                              y = mean(.data$y),
                                              theta = unique(.data$theta)),
                           ggplot2::aes(x=.data$x,y=.data$y,label = .data$cellIndex,angle = -1*.data$theta))

    }

    cmcLegend <- ggplotify::as.ggplot(cowplot::get_legend(refPlt)$grobs[[1]])

    refPlt <- refPlt +
      ggplot2::theme(legend.position = "none")

  })

  # library(patchwork)
  # return((refPlt | plt))
  if(type == "list"){
    return(list("reference" = refPlt,
                "target" = plt,
                "legend" = cmcLegend))
  }

  return(patchwork::wrap_plots(refPlt,plt,cmcLegend,nrow = 2,heights = c(1,.1)))
}
