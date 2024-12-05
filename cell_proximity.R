#' Evaluate Ripley's K function for each pair of cell types.
#'
#' @param DT data.table with cell type and spatial locations for each cell
#' @param R range of radii to evaluate
#' @param cell_type_column which column in DT contains cell type information
#' @param slide_column which column in DT contains slide_number to uniquely
#'   identify fov
#' @param fov_column which column in DT contains fov number
#'
#' @return
#' @export
Kall <- function(DT, R=c(0, 0.1), cell_type_column="cell_type",
                 slide_column="slide_ID_numeric", fov_column="fov") {
    # create window around each FOV's cells
    poly <- lapply(split(DT, by = c(slide_column, fov_column)), create_poly)

    # create point pattern
    cosmx.pp <- spatstat.geom::ppp(DT$sdimx, DT$sdimy,
        poly = poly,
        marks = factor(DT[[cell_type_column]]),
        unitname = c("mm", "mm"))
    
    # run K-function for each combination of cell types
    spatstat.core::alltypes(cosmx.pp,
        spatstat.core::Kcross,
        r = seq(R[1], R[2], length.out=100),
        correction="border")
}

#' Create a heatmap based on the results of the K function.
#'
#' @param K value from Kall function 
#' @param title A title for the plot.
#'
#' @return Heatmap
#' @export
heatmapK <- function(K, title = "Kall") {
    cell_pairs <- K[['which']]
    cell_pairs <- apply(cell_pairs, c(1, 2), function(x) {
        Ki <- K[['fns']][[x]]
        sum(Ki[['border']] - Ki[['theo']])
    })
    col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    # cell_pairs_scaled <- apply(cell_pairs, 2, scale)
    # rownames(cell_pairs_scaled) <- rownames(cell_pairs)
    ComplexHeatmap::Heatmap(cell_pairs, row_names_side = "left",
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_dend = FALSE, show_column_dend = FALSE,
        column_title = title, show_heatmap_legend = FALSE,
        col=col_fun
    )
}

pcf_plot <- function(K, from_cell_type, to_cell_type,
                     title = paste0(from_cell_type, ' to ', to_cell_type)) {
  idx <- K[['which']][from_cell_type, to_cell_type]
  plot(spatstat.core::pcf(K[['fns']][[idx]]),
       main = title, legend = FALSE)
}

#' Create polygon to define window in spatstat
#'
#' @param DT data.table containing spatial locations for each cell
#'
#' @return coordinates for polygon
#'
#' @internal
create_poly <- function(DT) {
  tissue_window <- with(DT, spatstat.geom::convexhull.xy(sdimx, sdimy))
  tissue_window[['bdry']][[1]]
}