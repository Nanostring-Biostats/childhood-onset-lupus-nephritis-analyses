#' Rank based inverse normal transformation
#'
#' This function is copied and borrowed from `RNOmni` package. 
#' Applies the rank-based inverse normal transform (INT) to a numeric vector. 
#' The INT can be broken down into a two-step procedure. 
#' In the first, the observations are transformed onto the probability scale using the 
#' empirical cumulative distribution function (ECDF). 
#' In the second, the observations are transformed onto the real line, as Z-scores, using the probit function. 
#'
#' @param u Numeric vector.
#' @param k Offset. Defaults to (3/8), correspond to the Blom transform.
#'
#' @return Numeric vector of rank normalized measurements. 
#'
#' @examples
#'
#' # Draw from chi-1 distribution
#' y <- rchisq(n = 1e3, df = 1)
#' # Rank normalize
#' z <- RankNorm(y)
#' # Plot density of transformed measurement
#' plot(density(z))
#' 
#' @export
#'

RankNorm <- 
function (u, k = 0.375) 
{
  if (!is.vector(u)) {
    stop("A numeric vector is expected for u.")
  }
  if ((k < 0) || (k > 0.5)) {
    stop("Select the offset within the interval (0,0.5).")
  }
  if (sum(is.na(u)) > 0) {
    stop("Please exclude observations with missing measurements.")
  }
  n <- length(u)
  r <- rank(u)
  out <- qnorm((r - k)/(n - 2 * k + 1))
  return(out)
}



#' Run DE analysis models on SMI data.
#'
#' @param assay_matrix counts or normalized assay matrix.  
#'                     Raw counts recommended for count distributions s.a. negative binomial (i.e., family=nbinom2). 
#'                     Normalized is recommended for gaussian models.
#' @param metadata data.table or data.frame of meta.data associated with assay.  
#'                 Must contain groupVar column and any covariates in the `formula` argument used for DE.
#'                 If cell id is not included as a column in metadata, should correspond to rownames in the metadata data.frame object.
#' @param formula right hand side of formula used in DE, possibly containing random effects.  
#'                     Must contain groupVar fixed effect. 
#'                     It's recommended that raw count models with a log-link (i.e., negative binomial,poisson)
#'                     , use a corresponding offset for library size in formula
#'                     (i.e., ~de_variable + offset(log(totalcounts))), where totalcounts is a column in the metadata.
#' @param groupVar  DE group variable with two levels,must exist in meta.data slot 
#' @param groupVar_levels  optional factor ordering for groupVar, which should have only two-levels (for now) 
#' @param neighborhood_counts  optional object created by smiDE::measure_neighbor_expression_by_celltype, which contains lists of 
#'                             the expression counts in neighboring cells by cell type compared to a cell type of reference 
#'                             (see `help(measure_neighbor_expression_by_celltype)`). 
#'                             May be used to control for confounding or 'bleed-over' in DE model (see examples)
#' @param nCores = 1, number of cores to use, set to 1 if running in serial mode
#' @param multiCore = TRUE, set to TRUE to use multiCore, FALSE to run in cluster mode
#' @param family a character string naming a family function, 
#'               See `family` for a generic discussion of families.
#'               Default is nbinom2 (with log link), specifying Negative binomial distribution: quadratic parameterization (Hardin & Hilbe 2007). V=mu*(1+mu/phi) = mu+mu^2/phi.
#' @param targets optional vector of targets to run DE for.  If not supplied, run models for all targets in assay.
#' @param cellid_colname column name in metadata corresponding to cell id.
#' @param ziformula optional formula for model of zero inflation.  Default is ~0 (no zeroinflation).
#' @param ... further arguments to be passed to model fitting function. 
#'            See example below. 
#'
#' @return data.table of DE results, with one row per target, marginal means and their SE's of the DE groups
#'          , estimated fold change and p-values for non-zero difference (identity link) or ratio \eqn{\neq 1} (log link) between groups.
#'
#' @examples
#'
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
#'        ,ref_celltype = c("macrophage", "fibroblast")
#'        ,sdimx_colname = "sdimx"
#'        ,sdimy_colname = "sdimy"
#'        ,weight_colname = "weight"
#'        ,contamination = "sum"
#'        ,verbose=TRUE
#' )
#' 
#' fibroblast_and_macrophage_cells <- gem@cell_metadata$rna[cell_type %in% c(pre_de_obj$nblist$ref_celltype),cell_ID]
#' de_results <-       
#'   smi_de(assay_matrix = gem@expression$rna$raw[,fibroblast_and_macrophage_cells]
#'          ,metadata = metainfo[cell_ID %in% fibroblast_and_macrophage_cells]
#'          ,formula = ~RankNorm(otherct_expr) + niche + tissue + offset(log(totalcounts)) 
#'          ,neighborhood_counts = pre_de_obj$nblist
#'          ,groupVar="niche"
#'          ,family="nbinom2"
#'          ,targets=rownames(gem@expression$rna$raw)[1:2]
#'   ) 
#'
#' results(de_results, "pairwise", variable=de_results$groupVar)
#'
#' fibroblast_and_macrophage_cells <- gem@cell_metadata$rna[cell_type %in% c(pre_de_obj$nblist$ref_celltype),cell_ID]
#' de_results <-       
#'   smi_de(assay_matrix = gem@expression$rna$raw[,fibroblast_and_macrophage_cells]
#'          ,metadata = metainfo[cell_ID %in% fibroblast_and_macrophage_cells]
#'          ,formula = ~RankNorm(otherct_expr) + cell_type + tissue + offset(log(totalcounts)) 
#'          ,neighborhood_counts = pre_de_obj$nblist
#'          ,groupVar="cell_type"
#'          ,family="nbinom2"
#'          ,targets=rownames(gem@expression$rna$raw)[1:2]
#'   ) 
#' results(de_results, "pairwise", variable=de_results$groupVar)
#' 
#' ## fit a mixed model with all cells, cell type as a covariate
#' pre_de_obj <- 
#' pre_de(counts = gem@expression[["rna"]][["raw"]]
#'        ,normalized_data = gem@expression[["rna"]][["normalized"]]
#'        ,metadata = metainfo
#'        ,cell_type_metadata_colname = "cell_type"
#'        ,split_neighbors_by_colname = "tissue"
#'        ,mm_radius = 0.05
#'        ,ref_celltype = "all"
#'        ,sdimx_colname = "sdimx"
#'        ,sdimy_colname = "sdimy"
#'        ,weight_colname = "weight"
#'        ,contamination = "sum"
#'        ,verbose=TRUE
#' )
#' de_results <-       
#'   smi_de(assay_matrix = gem@expression$rna$raw
#'          ,metadata = metainfo
#'          ,formula = ~RankNorm(otherct_expr) + cell_type + (1 | tissue) + offset(log(totalcounts)) 
#'          ,neighborhood_counts = pre_de_obj$nblist
#'          ,groupVar="cell_type"
#'          ,family="nbinom2"
#'          ,targets=rownames(gem@expression$rna$raw)[1:10]
#'   ) 
#' results(de_results, "one.vs.rest", variable=de_results$groupVar)
#'   
#'   
#' @export
#' @import data.table
#' @importFrom stats as.formula coef family formula logLik qnorm terms.formula
#' @importFrom methods as 
#' @importFrom utils stack
smi_de <- function(assay_matrix
                   ,metadata
                   ,formula
                   ,neighborhood_counts=NULL
                   ,groupVar
                   ,groupVar_levels=NULL
                   ,nCores = 1
                   ,multiCore = TRUE
                   ,family="nbinom2"
                   ,targets=NULL
                   ,cellid_colname = "cell_ID"
                   ,ziformula=~0
                   ,...
) {
 
  
  metainfo <- data.table::copy(as.data.table(metadata))
  
  stopifnot(!is.null(rownames(metadata)) & (cellid_colname %in% names(metainfo)))
  if(!(cellid_colname) %in% names(metainfo)){
    cellid_colname <- "cell_ID"
    metainfo[[cellid_colname]] <- rownames(metadata)
  }
  setnames(metainfo, old=cellid_colname, new="cell_ID")
  
  if(length(setdiff(colnames(assay_matrix), metainfo[[cellid_colname]])) > 0 ||
     length(setdiff(metainfo[[cellid_colname]], colnames(assay_matrix))) > 0 
  ){
    warning("Not identical set of cells between assay matrix and meta.data")
    mcells <- nrow(metainfo)
    acells <- ncol(assay_matrix)
    comm <- length(intersect(colnames(assay_matrix), metainfo[[cellid_colname]]))
    message(paste0(mcells, " cells in metadata."))
    message(paste0(acells, " cells in assay_matrix."))
    message(paste0(comm, " cells common between metadata and assay."))
  } 
  
  if(is.null(targets)) targets <- rownames(assay_matrix)
  mTerms <- all.vars(formula)
  mTerms <- setdiff(mTerms, "1")
  
  # check if groupVar is in model formula terms
  if (!groupVar %in% mTerms){
    stop ("Error: groupVar needs to be defined as fixed effect in the model.\n")
  }
 
  if(!missing(neighborhood_counts)) 
    stopifnot("neighborhood_counts must be NULL or an object of neighborexpr class returned by 'measure_neighbor_expr_by_celltype' function" = inherits(neighborhood_counts, "neighborexpr"))
  # check if terms in model are in sData
  if(missing(neighborhood_counts)) neighborhood_counts <- list()
  neighborTerms <- setdiff(mTerms, names(metainfo))
  candidate_neighborTerms <- unlist(lapply(names(neighborhood_counts$neighbor_expr_byct)
                                    ,function(x) paste0(x, c("", "_intensity", "_expr", "_cellcount"))))
  neighborTerms <- intersect(neighborTerms, candidate_neighborTerms)
  mTerms <- setdiff(mTerms, neighborTerms)
  
  missingTerms <- setdiff(c(mTerms, neighborTerms)
                          ,c(names(metainfo)
                             ,names(neighborhood_counts$neighbor_expr_byct)
                             ,candidate_neighborTerms)
                         )
  if (length(missingTerms) > 0){
    stop (paste0("Error: ", paste0(missingTerms, collapse=","), "were not found in in the meta.data slot of the passed seurat object.\n"))
  }
  pDat <- metainfo[,c("cell_ID", mTerms),with=F]
 
  ### Add check to ensure that covariates in model have more than one unique value. 
  uniq_check <- unlist(pDat[,lapply(.SD, uniqueN),.SDcols=c(mTerms)])
  non_uniq <- names(uniq_check)[which(uniq_check==1)]
  if(length(non_uniq) > 0){
    message(paste0("covariates: ", paste0(non_uniq, collapse=", "), " in formula have only one unique value in the dataset.\n"))
    for(ii in 1:length(non_uniq)){
      message(paste0(non_uniq[ii], " unique value: ", pDat[1][[non_uniq[ii]]]))
    }
    message(paste0("(", nrow(pDat), " total cells)"))
    stop("Non-unique covariates will cause error in regression models, such as:\n"
         ,"'contrasts can only be applied to factors with two or more levels'")
  }
  if(!is.numeric(pDat[[groupVar]])){
    if(missing(groupVar_levels)){
      pDat[[groupVar]] <- as.factor(pDat[[groupVar]])
    } else {
      pDat[[groupVar]] <- factor(pDat[[groupVar]],levels=groupVar_levels)
    }
  }
  
  for (i in setdiff(names(pDat), groupVar)){
    if (inherits(pDat[[i]], "character")) {
      pDat[, i] <- as.factor(pDat[[i]])
    }
  }
  updatedFormula <- formula(paste("y", as.character(formula)[2]
                                  , sep = " ~ "))
  
  
  if (nCores > 1) {
    assay_expr <- new.env()
    assay_expr$assay_expr <- assay_matrix 
    #assay_expr$assay_expr <- obj@assays$rna #assayDataElement(object, elt = elt)
    if (multiCore & Sys.info()['sysname'] != "Windows") {
      mixedOut <- parallel::mclapply(targets
                                     ,deFunc
                                     ,groupVar
                                     ,groupVar_levels
                                     ,pDat
                                     ,updatedFormula 
                                     ,family=family
                                     ,assay_expr
                                     ,neighborhood_counts
                                     ,ziformula
                                     ,typ="parallel"
                                     ,mc.cores = nCores
                                     ,...)
    } else {
      cl <- parallel::makeCluster(getOption("cl.cores", nCores))
      mixedOut <- parallel::parLapply(cl
                                      ,targets
                                      ,deFunc
                                      ,groupVar
                                      ,pDat
                                      ,updatedFormula 
                                      ,family
                                      ,assay_expr
                                      ,neighborhood_counts
                                      ,ziformula
                                      ,typ
                                      ,...)
      suppressWarnings(parallel::stopCluster(cl))
    }
    #mOut <- rbindlist(mixedOut)
    #mOut[,target:=targets]
  } else {
   #debugonce(deFunc)
   #debugonce(summarize_model)
   #deFunc(targets[1]
   #     ,groupVar
   #     ,groupVar_levels
   #     ,pDat
   #     ,updatedFormula
   #     ,family
   #     ,obj@assays$rna
   #     ,neighborhood_counts
   #     ,"non-parallel"
   #     )
    mixedOut <- lapply(targets
                       ,deFunc
                       ,groupVar
                       ,groupVar_levels
                       ,pDat
                       ,updatedFormula
                       ,family=family
                       ,assay_matrix
                       ,neighborhood_counts
                       ,ziformula
                       ,typ="non-parallel"
                       ,...)
  }
    names(mixedOut) <- targets
    
    return_obj <- list(results = mixedOut
                       ,groupVar = groupVar
                       ,modelterms = mixedOut[[1]]$terms
                       ,targets = targets
                       )
    
    #mOut <- rbindlist(mixedOut)
    #mOut[,target:=targets]
    #mOut <- assayDataApply(object, 1, deFunc, groupVar, pDat, formula(paste("expr", as.character(formula)[2], sep = " ~ ")), pairwise,  elt = elt, return_model=return_model,...)
  #return(mOut)
  class(return_obj) <- append(class(return_obj), "smide")
  #structure(return_obj, class="smi_de")
  return(return_obj)
}


deFunc <- function(target, groupVar, groupVar_levels, pDat
                   ,formula, family
                   ,assay_expr
                   ,neighbor_counts
                   ,zeroinflation_formula
                   ,typ
                   ,...) {
  message_parallel(paste0("Fitting model to target ",target, "."))
  
  #if(!is.null(neighbor_counts)){
  if(length(neighbor_counts) > 0){
    neighb_dt <- extract_neighborexpr_by_gene(neighbor_counts$neighbor_expr_byct
                                              ,target)
    
    pDat <- merge(pDat,neighb_dt, by="cell_ID")
    
    neighbor_ct_count <- data.table::copy(neighbor_counts$adjacency_counts_by_ct)
    names(neighbor_ct_count) <- gsub("\\-|\\ ","_",names(neighbor_ct_count))
    
    pDat <- merge(pDat, neighbor_ct_count,suffixes=c("_expr", "_cellcount"))
    all_celltypes <- setdiff(names(neighbor_ct_count), "cell_ID")
    ref_celltype <- neighbor_counts$ref_celltype 
    other_ctypes <- setdiff(all_celltypes
                            , c(ref_celltype, "allct", "otherct"))
    #pDat[,otherct_cellcount:=rowSums(.SD),.SDcols=paste0(other_ctypes, "_cellcount")]
    #pDat[,allct_cellcount:=rowSums(.SD),.SDcols=paste0(all_celltypes, "_cellcount")]
    for(cc in intersect(c("allct", "otherct"), names(pDat))){
      setnames(pDat, old=cc, new=paste0(cc, "_expr"))
    }
  }
  
  if (is.character(family)) {
    if (family=="beta") {
      family <- "beta_family"
      warning("please use ",sQuote("beta_family()")," rather than ",
              sQuote("\"beta\"")," to specify a Beta-distributed response")
    }
    family_char <- family
    #family <- get(family, mode = "function", envir = parent.frame())
  } else {
    family_char <- family()$family
  }
  
  if(typ=="parallel"){
    #y <- FetchData(assay_expr$assay_expr, vars=target)
    y <- matrix(assay_expr$assay_expr[target,],ncol=1)
    rownames(y) <- colnames(assay_expr$assay_expr)
  } else {
    #y <- FetchData(assay_expr, vars=target)
    y <- matrix(assay_expr[target,],ncol=1)
    rownames(y) <- colnames(assay_expr)
  }
  if (grepl("(nbinom|pois)",family_char)) {
    if (any(abs(y[,1]-round(y[,1]))>0.001)) {
      warning(sprintf("non-integer counts in a %s model",
                      family_char))
      message_parallel(paste0("converting counts to nearest integer"))
      y[,1] <- round(y[,1])
    }
  }
  colnames(y) <- "y"
  dat <- merge(data.table(cell_ID=rownames(y), y=y[,1])
               ,pDat, by=c("cell_ID"))
  
  re_terms <- lme4::findbars(formula) 
  has_re <- length(re_terms) > 0
  stopifnot("Only 1 random effect grouping level currently supported" = 
              length(re_terms) <=1)
  convergence_error <- FALSE
  err <- NULL
  model_warning_msg <- emmeans_warning_msg <- NULL
  fixed_formula <- NULL ## only needed in nebula::nebula case; remains null for fixed effect only models
  mod_time <- system.time({
  if(!has_re){
    if(family_char=="nbinom2"){
      fittype <- "MASS::glm.nb"
      mod <- 
      withCallingHandlers({
        tryCatch({
          MASS::glm.nb(formula
                              ,data=dat
                              , ...
                              ) 
          
          
        },error = function(e){
            convergence_error <<- TRUE
            err <<- conditionMessage(e)
            "convergence_error"
          })
      }, warning = function(w){
        model_warning_msg <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      })
      
      if(!convergence_error) model_warning_msg <- paste0("converged=",mod$converged, " ",model_warning_msg)
       
    }  else {
      fittype <- "stats::glm"
      
      mod <- 
      withCallingHandlers({
        tryCatch({
          stats::glm(formula
                    ,data=dat
                    ,family=family_char
                    , ...
                    )
          },error = function(e){
            convergence_error <<- TRUE
            err <<- conditionMessage(e)
            "convergence_error"
          })
      }, warning = function(w){
        model_warning_msg <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      })
      
    }
  } else {
    #re_formula <- as.formula(paste0('~', paste0("(",re_terms,")", collapse="+")))
    re_formula <- as.formula(paste0('~', paste0(re_terms, collapse="+")))
    tt <- terms(formula)
    re_vars <- lapply(re_terms, function(x) gsub("1 | ", "", x))
    re_vars <- unlist(lapply(1:length(re_terms), function(ii) gsub("^.*\\|[\\ ]+", "", re_terms[ii])))
    re_frmlas <- unlist(lapply(1:length(re_terms), function(ii) gsub("\\|.*$", "", re_terms[ii])))
    re_frmlas <- lapply(re_frmlas, function(x){
      as.formula(paste0("~",x))
    })
    names(re_frmlas) <- re_vars
    allvar <- as.character(attr(tt, "variables"))[-1]
    response <- allvar[attr(tt,"response")] 
    offsetv <- allvar[attr(tt,"offset")] 
    allterms <- attr(terms.formula(formula),"term.labels")
    fixed_terms <- setdiff(allterms, re_terms)
    if(length(fixed_terms)==0) fixed_terms <- "1"
    if(length(offsetv) > 0) fixed_terms <-  c(fixed_terms, offsetv)
    fixed_formula <-  as.formula(paste0(response, ' ~ ', paste0(fixed_terms, collapse="+")))
    
      if(family_char == "gaussian"){
        fittype <- "lme4::lmer"
        mod <- 
        withCallingHandlers({
          tryCatch({
            lme4::lmer(formula 
                       ,data=dat
                       , ...
                       ) 
            },error = function(e){
              convergence_error <<- TRUE
               err <<- conditionMessage(e)
              "convergence_error"
          })
         }, warning = function(w){
           model_warning_msg <<- conditionMessage(w)
           invokeRestart("muffleWarning")
         })
        
             
      } else if (family_char %in% c("nbinom2","poisson")){
       
       fittype <- "nebula::nebula"
       modelv <- "NBGMM" # neg binomial gamma mixed model
       if(family_char=="poisson") modelv <- "PMM" # poisson gamma mixed model
       data.table::setkeyv(dat, re_vars)
       offsetcol <- NULL
       offsetuse <- NULL
       if(length(offsetv) > 0){
         offsetcol <- all.vars(as.formula(paste0("~",offsetv))) 
         offsetuse <- dat[[offsetcol]]
       }
       mod <- 
       withCallingHandlers({
         tryCatch({
           nebula::nebula(
             count = as(matrix(dat[["y"]],nrow=1), "dgCMatrix")
             ,id=dat[[re_vars]]
             ,pred = model.matrix(fixed_formula,data=dat) 
             ,offset = offsetuse 
             ,covariance=TRUE
             ,ncore=1
             ,model=modelv
           )
           ## attach overdispersions if needed
         },error = function(e){
           convergence_error <<- TRUE
           err <<- conditionMessage(e)
           "convergence_error"
         }) 
       }, warning = function(w){
         model_warning_msg <<- conditionMessage(w)
         invokeRestart("muffleWarning")
       })
      
      }  else {
        fittype <- "MASS::glmmPQL"
        
        mod <- 
        withCallingHandlers({
          tryCatch({
            MASS::glmmPQL(
                                 fixed = fixed_formula
                                 ,random = re_formula
                                 ,family = family
                                 ,data=dat
                                 , ...
                                 ) 
            },error = function(e){
              convergence_error <<- TRUE
              err <<- conditionMessage(e)
              "convergence_error"
         })
       }, warning = function(w){
         model_warning_msg <<- conditionMessage(w)
         invokeRestart("muffleWarning")
       })
      }
  } 
  })

  mod_out <- 
    withCallingHandlers({
                summarize_model(mod
                               ,modtime = mod_time
                               ,groupVar
                               ,dat
                               ,fittype
                               ,original_formula = formula
                               ,fixed_formula = fixed_formula
                               ,error_msg = err
                               )
       }, warning = function(w){
         emmeans_warning_msg <<- conditionMessage(w)
         invokeRestart("muffleWarning")
       })
    
  
  mod_out[["model_summary"]][,target:=target]
  for(dt in mod_out[["pairwise_list"]]) dt[,target:=target]
  for(dt in mod_out[["onevall_list"]]) dt[,target:=target]
  for(dt in mod_out[["onevrest_list"]]) dt[,target:=target]
  for(dt in mod_out[["emm_list"]]) dt[,target:=target]

  if(!is.null(model_warning_msg)){
    mod_out[["model_summary"]][,msg:=model_warning_msg]
  }
  if(!is.null(emmeans_warning_msg)){
    for(dt in mod_out[["pairwise_list"]]) dt[,msg:=emmeans_warning_msg]
    for(dt in mod_out[["onevall_list"]]) dt[,msg:=emmeans_warning_msg]
    for(dt in mod_out[["onevrest_list"]]) dt[,msg:=emmeans_warning_msg]
    for(dt in mod_out[["emm_list"]]) dt[,msg:=emmeans_warning_msg]
  } 
  return(mod_out)
}







