
#' Calculate cell adjacencies and neighboring cell expression of genes before DE analysis.
#' 
#' @param counts matrix of counts
#' 
#' @param normalized_data optional, matrix of normalized data.  
#' If not provided, totalCounts normalization is applied to counts.  
#' @param metadata cell-level metadata corresponding to the cells in counts, with one row per cell.
#' Should contain the columns of cell_type_metadata_colname, 
#' cellid_colname, split_neighbors_by_colname, sdimx_colname, and sdimy_colname.
#' @param ref_celltype Reference category (commonly a cell type), used to calculate neighbor expression with respect to.  
#' Can be vector of one or more levels of cell_type_metadata_colname. 
#' Can use keyword "all" for all levels of cell_type_metadata_colname.
#' @param cell_type_metadata_colname column (commonly cell type column) corresponding to categories by which cell expression is 
#' aggregated by to calculate expression among neighboring cells by cell_type_metadata_colname type.
#' @param weight_colname (optional).  If NULL, neighbor expression is unweighted.  
#' If "weight", the column corresponding to "weight" (1/distance) by default is used to weight calculated expression of neighbor cells.
#' @param cellid_colname column name corresponding to cell id.
#' @param sdimx_colname column name corresponding to 'x' axis spatial coordinate of cell.
#' @param sdimy_colname column name corresponding to 'y' axis spatial coordinate of cell.
#' @param split_neighbors_by_colname (optional) If specified, identification of neighboring cells is split by this 
#' column in the meta data.  For example, to avoid identifying cells as neighbors if they were to have close x/y coordinates, but were on different tissue.
#' @param mm_radius maximum euclidean distance in units of sdimx, sdimy by which two cells will be identified as neighbors. 
#' @param verbose print some extra messages during calculation.
#' @param contamination "sum" (default) calculates a total amount of target expressed in neighboring cells.  Anything other values passed will calculate an average.
#' @param nb_celltypes_to_individually_calc vector of cell types for which to individually calculate neighbor expression.
#' default is NULL, which will make a separate list in neighbor_expr_by_ct output with matrix of expression by cell type
#' 
#' @return list object with data.table of cell_adjacencies (one row per pair of adjacent cells),
#'  a list (nblist) with sparse matrices corresponding to neighbor cell expression passed to downstream DE in smi_de or run_de functions.
#'  
#' @examples 
#' library(Giotto)
#' library(data.table)
#' datadir<-system.file("extdata", package="smiDE")
#' gem <- readRDS(paste0(datadir, "/small_nsclc.rds"))
#' metainfo <- data.table::copy(gem@cell_metadata$rna)
#' metainfo <- merge(metainfo, gem@spatial_locs$raw, by="cell_ID")
#'  
#' pre_de_obj <- 
#' pre_de(counts = gem@expression[["rna"]][["raw"]]
#'        ,normalized_data = gem@expression[["rna"]][["normalized"]]
#'        ,metadata = metainfo
#'        ,cell_type_metadata_colname = "cell_type"
#'        ,split_neighbors_by_colname = "tissue"
#'        ,mm_radius = 0.05
#'        ,ref_celltype = "fibroblast"
#'        ,sdimx_colname = "sdimx"
#'        ,sdimy_colname = "sdimy"
#'        ,weight_colname = "weight"
#'        ,contamination = "sum"
#'        ,verbose=TRUE
#'  )
#' names(pre_de_obj)
#' 
#' # Example for measuring neighbor expression with 'all' levels of
#' # ref_celltype = 'all'
#' pre_de_obj <- 
#' pre_de(counts = gem@expression[["rna"]][["raw"]]
#'        ,normalized_data = gem@expression[["rna"]][["normalized"]]
#'        ,metadata = metainfo
#'        ,cell_type_metadata_colname = "cell_type"
#'        ,split_neighbors_by_colname = "tissue"
#'        ,mm_radius = 0.05
#'        ,ref_celltype = c("all")
#'        ,sdimx_colname = "sdimx"
#'        ,sdimy_colname = "sdimy"
#'        ,weight_colname = "weight"
#'        ,contamination = "sum"
#'        ,verbose=TRUE
#' )
#' names(pre_de_obj)
#' dim(pre_de_obj$nblist$neighbor_expr_byct$otherct)  
#' metainfo[, .N]
#'
#'
#' # Example for measuring neighbor expression with a vector of multiple cell types
#' pre_de_obj <- 
#' pre_de(counts = gem@expression[["rna"]][["raw"]]
#'        ,normalized_data = gem@expression[["rna"]][["normalized"]]
#'        ,metadata = metainfo
#'        ,cell_type_metadata_colname = "cell_type"
#'        ,split_neighbors_by_colname = "tissue"
#'        ,mm_radius = 0.05
#'        ,ref_celltype = c("macrophage", "fibroblast")
#'        ,sdimx_colname = "sdimx"
#'        ,sdimy_colname = "sdimy"
#'        ,weight_colname = "weight"
#'        ,contamination = "sum"
#'        ,verbose=TRUE
#' )
#' names(pre_de_obj)
#' dim(pre_de_obj$nblist$neighbor_expr_byct$otherct)  
#' metainfo[cell_type %in% pre_de_obj$nblist$ref_celltype, .N]
#'  
#' @export  
pre_de <- function(counts
                   ,normalized_data = NULL
                   ,metadata
                   ,ref_celltype
                   ,cell_type_metadata_colname
                   ,weight_colname = NULL
                   ,cellid_colname = "cell_ID"
                   ,sdimx_colname = "sdimx"
                   ,sdimy_colname = "sdimy"
                   ,split_neighbors_by_colname = "tissue" 
                   ,mm_radius = 0.05
                   ,verbose=TRUE
                   ,contamination = "sum"
                   ,nb_celltypes_to_individually_calc = NULL
                   ){

  metainfo <- data.table::copy(as.data.table(metadata))
 
  
  stopifnot(!is.null(rownames(metadata)) & (cellid_colname %in% names(metainfo)))
  if(!(cellid_colname) %in% names(metainfo)){
    cellid_colname <- "cell_ID"
    metainfo[[cellid_colname]] <- rownames(metadata)
  }
   
  stopifnot(sdimx_colname %in% names(metainfo)) 
  stopifnot(sdimy_colname %in% names(metainfo)) 
  stopifnot(cell_type_metadata_colname %in% names(metainfo)) 
  
  if("all" %in% ref_celltype) ref_celltype <- unique(metainfo[[cell_type_metadata_colname]])
  stopifnot(all(ref_celltype %in% unique(metainfo[[cell_type_metadata_colname]])))
  stopifnot(all(setdiff(nb_celltypes_to_individually_calc, c("otherct", "allct")) %in% unique(metainfo[[cell_type_metadata_colname]])))
  
  if(length(setdiff(colnames(counts), metainfo[[cellid_colname]])) > 0 ||
     length(setdiff(metainfo[[cellid_colname]], colnames(counts))) > 0 
     ){
    warning("Not identical set of cells between assay matrix and meta.data")
    mcells <- nrow(metainfo)
    acells <- ncol(counts)
    comm <- length(intersect(colnames(counts), metainfo[[cellid_colname]]))
    message(paste0(mcells, " cells in metadata."))
    message(paste0(acells, " cells in assay."))
    message(paste0(comm, " cells common between metadata and assay."))
  }
  
  if(is.null(normalized_data)){
    colsumms <- Matrix::colSums(counts)
    norm_factors <- colsumms / mean(colsumms)
    normalized_data <- Matrix::t(Matrix::t(counts)/(norm_factors))
    ### apply normalization factor by column (cell).
  }
  stopifnot(all.equal(dim(counts), dim(normalized_data)))
   
  if(is.null(split_neighbors_by_colname)){
    split_neighbors_by_colname <- "all_data"
    metainfo[[split_neighbors_by_colname]] <- "all_cells" 
  }
   
  # Identify cell-cell spatial neighbors 
  if(verbose) message(paste0(Sys.time()
                             ,", identifying cell-cell spatial neighbors within "
                             ,mm_radius, " radius."))
  
  tiss <- split(metainfo[,c(cellid_colname, sdimx_colname, sdimy_colname,split_neighbors_by_colname),with=FALSE]
                ,by=split_neighbors_by_colname)
  cell_adjacency_dt <-
    rbindlist(
      lapply(1:length(tiss), function(x){
        nbr <- fast_make_all_neighbors(tiss[[x]]
                                       ,cellid_col = cellid_colname
                                       ,sdimx_col = sdimx_colname
                                       ,sdimy_col = sdimy_colname
                                       ,radius = mm_radius
                                       )
        if(verbose) message(paste0("neighbors calculated for "
                                 ,split_neighbors_by_colname
                                 , ": "
                                 ,tiss[[x]][[split_neighbors_by_colname]][1]
                                 ))
        return(nbr)
      })
    ) 
  rm(tiss); gc()
  
  # Create neighborhood list objects from normalized data. 
  if(verbose) message(paste0(Sys.time()
                             ,". creating neighborhood list object"
                             #, ref_celltype
                             #, " ."
                             ))
 
  #nblistl <- list()
  #nblistl[[paste0(ref_celltype[1])]] <- 
  ict <- 1L 
  if(verbose) message(paste0("Measuring neighbor expression for ", ref_celltype[1]))
  if(verbose) message(paste0("(cluster ",1, " of ", length(ref_celltype), ")"))
  nblist <- 
    measure_neighbor_expression_by_celltype(assay_matrix = normalized_data
                                            ,metadata = metainfo
                                            ,cell_adjacencies = cell_adjacency_dt
                                            ,ref_celltype[1]
                                            ,cell_type_metadata_colname
                                            ,contamination
                                            ,cellid_colname
                                            ,weight_colname
                                            #,Wmat
                                            #,ct_matrix
                                            #,adjacency_counts_by_ct
                                            ,nb_celltypes_to_individually_calc = nb_celltypes_to_individually_calc
    )
   
  if(length(ref_celltype) > 1){
    for(ct in ref_celltype[2:length(ref_celltype)]){
      ict <- ict + 1L
      if(verbose) message(paste0(Sys.time()))
      if(verbose) message(paste0("Measuring neighbor expression for ", ct))
      if(verbose) message(paste0("(cluster ",ict, " of ", length(ref_celltype), ")"))
      #nblistl[[paste0(ct)]] <- 
      nblist_tmp <- 
        measure_neighbor_expression_by_celltype(assay_matrix = normalized_data
                                                ,metadata = metainfo
                                                ,cell_adjacencies = cell_adjacency_dt
                                                ,ct
                                                ,cell_type_metadata_colname
                                                ,contamination
                                                ,cellid_colname
                                                ,weight_colname
                                                ,Wmat = nblist$adjacency_mat 
                                                ,ct_matrix = nblist$ct_matrix
                                                ,adjacency_counts_by_ct = nblist$adjacency_counts_by_ct 
                                                ,nb_celltypes_to_individually_calc = nb_celltypes_to_individually_calc
        )
      
     #for(nm in names(nblist[["neighbor_expr_byct"]])){
     for(nm in names(nblist[["neighbor_expr_byct"]])){
       nblist[["neighbor_expr_byct"]][[nm]] <- 
         Matrix::cbind2(nblist[["neighbor_expr_byct"]][[nm]]
                        ,nblist_tmp[["neighbor_expr_byct"]][[nm]]
       )
     }
     rm(nblist_tmp); gc()
    }
  }
  
  
  #nblist <- nblistl[[1]]
  #if(length(ref_celltype) > 1){
  #  for(ii in 2:length(nblistl)){
  #    for(nm in names(nblist[["neighbor_expr_byct"]])){
  #      nblist[["neighbor_expr_byct"]][[nm]] <- 
  #        Matrix::cbind2(nblist[["neighbor_expr_byct"]][[nm]]
  #                       ,nblistl[[ii]]$neighbor_expr_byct[[nm]]
  #      )
  #    }
  #  } 
  #}
 
  #is_all_zero <- which(unlist(lapply(nblist$neighbor_expr_byct, sum))==0)
  #dims <- dim(nblist[["neighbor_expr_byct"]][[1]])
  #for(nbz in is_all_zero){
  #  outsm <- new("dgCMatrix")
  #  outsm@Dim <- dims
  #  outsm@p <- integer(dims[2] +  1L)
  #  dimnames(outsm) <- dimnames(nblist[["neighbor_expr_byct"]][[1]])
  #  nblist$neighbor_expr_byct[[nbz]] <- outsm
  #} 
  #
  #rslist <- lapply(nblist$neighbor_expr_byct, Matrix::rowSums)
  #nz_rows <-  which(unlist(lapply(rslist, function(x) any(x==0))))
  #nz_rows <- setdiff(nz_rows, is_all_zero)
  #for(nzr in nz_rows){
  #  dgt <- as(nblist$neighbor_expr_byct[[nzr]], "dgTMatrix") 
  #  missing_rows <- setdiff(0:(nrow(dgt)-1L), dgt@i)
  #  dgt@i <- c(dgt@i, missing_rows)
  #  dgt@j <- c(dgt@j, rep(0L, length(missing_rows)))
  #  dgt@x <- c(dgt@x, rep(0L, length(missing_rows)))
  #  nblist$neighbor_expr_byct[[nzr]] <- as(dgt, "dgCMatrix")
  #  #print(nzr)
  #  #print(all.equal(as.matrix(dgt), as.matrix(nblist$neighbor_expr_byct[[nzr]])))
  #  #stopifnot(all.equal(as.matrix(dgt), as.matrix(nblist$neighbor_expr_byct[[nzr]])))
  #} 
  
  nblist$ref_celltype <- ref_celltype 
  
  return(list(
              nblist = nblist
              ,cell_adjacency_dt = cell_adjacency_dt
              ))
}
  