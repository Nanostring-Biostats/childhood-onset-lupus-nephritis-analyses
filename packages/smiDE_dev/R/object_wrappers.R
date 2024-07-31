
#' Calculate cell adjacencies and neighboring cell expression of genes before DE analysis.
#'
#' @param object object of class giotto or seurat.
#' @param config_pre_de A list with the following potential configuration parameters:
#'   * \strong{assay} name of assay to use.
#'   *  \strong{counts} matrix of counts
#'   *  \strong{normalized} optional, matrix of normalized data.  
#'   *  \strong{cellid_colname} cellid_colname column name corresponding to cell id.
#'   *  \strong{cell_type_metadata_colname} column (commonly cell type column) corresponding to categories by which cell expression is aggregated by to calculate expression among neighboring cells by cell_type_metadata_colname type.
#'   *  \strong{split_neighbors_by_colname} (optional) If specified, identification of neighboring cells is split by this column in the meta data.  For example, to avoid identifying cells as neighbors if they were to have close x/y coordinates, but were on different tissue.
#'   *  \strong{mm_radius} maximum euclidean distance in units of sdimx, sdimy by which two cells will be identified as neighbors. 
#'   *  \strong{ref_celltype} Reference category (commonly a cell type), used to calculate neighbor expression with respect to.
#'   *  \strong{weight_colname} (optional).  If NULL, neighbor expression is unweighted. If "weight", the column corresponding to "weight" (1/distance) by default is used to weight calculated expression of neighbor cells.
#'   *  \strong{contamination} "sum" by default. 
#'   *  \strong{sdimx_colname} column name corresponding to 'x' axis spatial coordinate of cell.
#'   *  \strong{sdimy_colname} column name corresponding to 'y' axis spatial coordinate of cell.
#'   *  \strong{verbose} print some extra messages during calculation.
#'   
#' @param ... currently ignored.
#' 
#' @return list of lists, pre_de_object can be passed directly to run_de or smi_de functions. 
#'
#' @examples
#' library(Giotto)
#' datadir<-system.file("extdata", package="smiDE")
#' gem <- readRDS(paste0(datadir, "/small_nsclc.rds"))
#' config_pre_de <- list(assay = "rna"
#'                       ,counts = "raw" 
#'                       ,normalized = "normalized"
#'                       ,cellid_colname = "cell_ID"
#'                       ,cell_type_metadata_colname = "cell_type"
#'                       ,split_neighbors_by_colname = "tissue"
#'                       ,mm_radius = 0.05
#'                       ,ref_celltype = "fibroblast"
#'                       ,sdimx_colname = "sdimx"
#'                       ,sdimy_colname = "sdimy"
#'                       ,weight_colname = "weight"
#'                       ,contamination = "sum"
#'                       ,verbose=TRUE
#'                       )
#' pre_de_obj <- run_pre_de(gem, config_pre_de) 
#' 
#' @export
run_pre_de <- function(object
                       ,config_pre_de
                       ,...){
  UseMethod("run_pre_de")
}

#' @export
run_pre_de.giotto <- function(object
                              ,config_pre_de
                              ,...){
  
  ## required arguments in config
  stopifnot(all(c("assay"
                  ,"counts"
                  ,"normalized"
                  ,"ref_celltype"
                  ,"cell_type_metadata_colname"
                  ) %in% names(config_pre_de)
                ))
  
   
  args <- formals(pre_de)
  for(nm in names(config_pre_de)){
    args[[nm]] <- config_pre_de[[nm]]
  } 
 
  metainfo <- data.table::copy(object@cell_metadata[[args[["assay"]]]])
  metainfo <- merge(metainfo, object@spatial_locs$raw, by="cell_ID")
  
  #use do.call and pass argument list overwriting null defaults 
  pre_de_obj <- 
  pre_de(counts = object@expression[[args[["assay"]]]][[args[["counts"]]]]
         ,normalized_data = object@expression[[args[["assay"]]]][[args[["normalized"]]]]
         ,metadata = metainfo 
         ,cell_type_metadata_colname = args[["cell_type_metadata_colname"]]
         ,split_neighbors_by_colname = args[["split_neighbors_by_colname"]]
         ,mm_radius = args[["mm_radius"]]
         ,ref_celltype = args[["ref_celltype"]]
         ,sdimx_colname = args[["sdimx_colname"]]
         ,sdimy_colname = args[["sdimy_colname"]]
         ,weight_colname =args[["weight_colname"]]
         ,contamination =args[["contamination"]]
         ,verbose= args[["verbose"]]
         ,nb_celltypes_to_individually_calc = args[["nb_celltypes_to_individually_calc"]]
  )
  return(pre_de_obj)
}

#' @export
run_pre_de.Seurat <- function(object
                              ,config_pre_de
                              ,...){
  
  ## required arguments in config
  stopifnot(all(c("assay"
                  ,"counts"
                  ,"normalized"
                  ,"ref_celltype"
                  ,"cell_type_metadata_colname"
                  ) %in% names(config_pre_de)
                ))
  
   
  args <- formals(pre_de)
  args[["sdimx_colname"]] <- "viz_sdimx"
  args[["sdimy_colname"]] <- "viz_sdimy"
  args[["cell_type_metadata_colname"]] <- "cellType"
   
  for(nm in names(config_pre_de)){
    args[[nm]] <- config_pre_de[[nm]]
  } 
 
  metainfo <- data.table(object@meta.data)
  
  if(!("totalcounts" %in% names(metainfo))){
    totalcounts <- data.table(stack(Matrix::colSums(object@assays[[args[["assay"]]]]@counts)))
    totalcounts[,`:=`(cell_ID = ind, totalcounts = values, values=NULL, ind=NULL)]
    metainfo <- merge(metainfo, totalcounts, by="cell_ID")
  }
  
  #use do.call and pass argument list overwriting null defaults 
  pre_de_obj <- 
  pre_de(counts = Seurat::GetAssayData(object, slot=args[["counts"]], assay = args[["assay"]])
         ,normalized_data = Seurat::GetAssayData(object, slot=args[["normalized"]],assay = args[["assay"]])
         ,metadata = metainfo 
         ,cell_type_metadata_colname = args[["cell_type_metadata_colname"]]
         ,split_neighbors_by_colname = args[["split_neighbors_by_colname"]]
         ,mm_radius = args[["mm_radius"]]
         ,ref_celltype = args[["ref_celltype"]]
         ,sdimx_colname = args[["sdimx_colname"]]
         ,sdimy_colname = args[["sdimy_colname"]]
         ,weight_colname =args[["weight_colname"]]
         ,contamination =args[["contamination"]]
         ,verbose= args[["verbose"]]
         ,nb_celltypes_to_individually_calc = args[["nb_celltypes_to_individually_calc"]]
  )
  return(pre_de_obj)
}

#' Run DE analysis 
#'
#' @param object object of class giotto or seurat.
#' @param config_de A list with the following potential configuration parameters:
#'   * \strong{assay} name of assay to use.
#'   *  \strong{data_use} name of data layer in assay used as outcome in DE.  
#'                        Could be 'counts' for count based regression or 'normalized' , for a linear based regression.  
#'   * \strong{formula} right hand side of formula used in DE, possibly containing random effects.  
#'                      Must contain groupVar fixed effect. 
#'                     It's recommended that raw count models with a log-link (i.e., negative binomial,poisson)
#'                     , use a corresponding offset for library size in formula
#'                     (i.e., ~de_variable + offset(log(totalcounts))), where totalcounts is a column in the metadata.
#'   * \strong{groupVar}  DE group variable with two levels,must exist in meta.data slot 
#'   * \strong{groupVar_levels}  optional factor ordering for groupVar, which should have only two-levels (for now) 
#'                        May be used to control for confounding or 'bleed-over' in DE model (see examples)
#'   * \strong{family} a family function, a character string naming a family function, 
#'          or the result of a call to a family function (variance/link function) information. 
#'          See `family` for a generic discussion of families.
#'          Default is nbinom2 (with log link), specifying Negative binomial distribution: quadratic parameterization (Hardin & Hilbe 2007). V=mu*(1+mu/phi) = mu+mu^2/phi.
#'  * \strong{targets} optional vector of targets to run DE for.  If not supplied, run models for all targets in assay.
#'  * \strong{cellid_colname} column name corresponding to cell id in the metadata of giotto / seurat object. 
#'  * \strong{targets} optional vector of targets to run DE for.  If not supplied, run models for all targets in assay.
#'  * \strong{ziformula} optional formula for zero inflation model.  Default formula is ~0, corresponding to no zero inflation.
#' 
#' @param pre_de_obj optional object produced by `pre_de` or `run_pre_de` with cell neighbor expression information. 
#' @param ... optional further arguments to be passed to during model fitting. 
#'            See example below.  
#' 
#' @return list of lists, pre_de_object can be passed directly to run_de or smi_de functions. 
#'
#' @examples
#' library(Giotto)
#' datadir<-system.file("extdata", package="smiDE")
#' gem <- readRDS(paste0(datadir, "/small_nsclc.rds"))
#' config_pre_de <- list(assay = "rna"
#'                       ,counts = "raw" 
#'                       ,normalized = "normalized"
#'                       ,cellid_colname = "cell_ID"
#'                       ,cell_type_metadata_colname = "cell_type"
#'                       ,split_neighbors_by_colname = "tissue"
#'                       ,mm_radius = 0.05
#'                       ,ref_celltype = "fibroblast"
#'                       ,sdimx_colname = "sdimx"
#'                       ,sdimy_colname = "sdimy"
#'                       ,weight_colname = "weight"
#'                       ,contamination = "sum"
#'                       ,verbose=TRUE
#'                       )
#' pre_de_obj <- run_pre_de(gem, config_pre_de) 
#' 
#' config_de <- list(
#'                   assay = "rna"
#'                   ,data_use = "raw"
#'                   ,formula = ~RankNorm(otherct_expr) + niche + tissue + offset(log(totalcounts))
#'                   ,groupVar="niche"
#'                   ,ziformula = ~1
#'                   ,family="nbinom2"
#'                   ,targets=rownames(gem@expression$rna$raw)[1:2]
#'                   )
#' 
#' de_obj <- 
#'   smiDE::run_de(gem, pre_de_obj, config_de)
#'
#' results(de_obj)
#' 
#' @export
run_de <- function(object, pre_de_obj, config_de, ...) UseMethod("run_de")


#' @export
run_de.giotto <- function(object
                          ,pre_de_obj=NULL
                          ,config_de
                          ,...
                          ){

  ## required arguments in config
  stopifnot(all(c("assay", "data_use", "formula", "groupVar") %in% names(config_de)))
  
  args <- formals(smi_de)
  for(nm in names(config_de)){
    args[[nm]] <- config_de[[nm]]
  } 
  
  if(!missing(pre_de_obj)){
    pre_de_obj <- pre_de_obj[["nblist"]]
  } 
  
  metainfo <- data.table::copy(object@cell_metadata[[config_de$assay]]) 
  metainfo <- merge(metainfo, object@spatial_locs$raw, by="cell_ID")
  
  de_obj <-  
  smi_de(assay_matrix = object@expression[[config_de$assay]][[config_de$data_use]]
         ,metadata = metainfo 
         ,formula = args$formula
         ,neighborhood_counts = pre_de_obj
         ,groupVar=args[["groupVar"]]
         ,ziformula = as.formula(args[["ziformula"]])
         ,family=args[["family"]]
         ,targets=args[["targets"]]
         ,nCores = args[["nCores"]]
         ,...
         )
  return(de_obj)
}

#' @export
run_de.Seurat <- function(object
                          ,pre_de_obj=NULL
                          ,config_de
                          ,...
                          ){

  ## required arguments in config
  stopifnot(all(c("assay", "data_use", "formula", "groupVar") %in% names(config_de)))
  
  args <- formals(smi_de)
  for(nm in names(config_de)){
    args[[nm]] <- config_de[[nm]]
  } 
  
  if(!missing(pre_de_obj)){
    pre_de_obj <- pre_de_obj[["nblist"]]
  } 
  
  metainfo <- data.table(object@meta.data)
  
  if(!("totalcounts" %in% names(metainfo))){
    totalcounts <- data.table(stack(Matrix::colSums(object@assays[[args[["assay"]]]]@counts)))
    totalcounts[,`:=`(cell_ID = ind, totalcounts = values, values=NULL, ind=NULL)]
    metainfo <- merge(metainfo, totalcounts, by="cell_ID")
  }
  
  
  de_obj <-  
  smi_de(assay_matrix = Seurat::GetAssayData(object, slot=args[["data_use"]], assay = args[["assay"]])
         ,metadata = metainfo 
         ,formula = args$formula
         ,neighborhood_counts = pre_de_obj
         ,groupVar=args[["groupVar"]]
         ,ziformula = as.formula(args[["ziformula"]])
         ,family=args[["family"]]
         ,targets=args[["targets"]]
         ,nCores = args[["nCores"]]
         ,...
         )
  
 
  return(de_obj)
}


