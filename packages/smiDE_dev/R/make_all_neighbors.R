
#' Make table with all cell neighbors within a radius.
#'
#' Make table with all cell neighbors within a radius.
#' Passed to downstream functions  
#' such as `smiDE::measure_neighbor_expression_by_celltype`.
#' 
#' @param dt data.table of metadata containing (at least) columns for cell ids, 
#' x-coordinates, y-coordinates.
#' @param cellid_col name of column with cellid, i.e., "cell_ID". 
#' @param sdimx_col name of column with x-coordinate for cell, i.e., "sdimx". 
#' @param sdimy_col name of column with x-coordinate for cell, i.e., "sdimy". 
#' @param radius radius within which to identify all other neighboring cells, i.e, 0.05.
#' 
#' @examples
#' library(Giotto)
#' library(data.table)
#' datadir<-system.file("extdata", package="smiDE")
#' gem <- readRDS(paste0(datadir, "/small_nsclc.rds"))
#' metainfo <- data.table::copy(gem@cell_metadata$rna)
#' metainfo <- merge(metainfo, gem@spatial_locs$raw, by="cell_ID")
#' 
#' tiss <- split(metainfo, by="tissue")
#' 
#' cell_adjacency_dt <- 
#'   rbindlist(
#'     lapply(1:length(tiss), function(x){
#'       nb <- fast_make_all_neighbors(tiss[[x]][,.(cell_ID, sdimx, sdimy)]
#'                                     ,cellid_col = "cell_ID"
#'                                     ,sdimx_col ="sdimx"
#'                                     ,sdimy_col ="sdimy"
#'                                     ,radius = 0.05) 
#'       return(nb) 
#'     })
#'   )
#'
#' @export
#'

fast_make_all_neighbors <- function(dt
                                    ,cellid_col="cell_ID"
                                    ,sdimx_col="sdimx"
                                    ,sdimy_col="sdimy"
                                    ,radius=0.05){
  dtuse <- copy(dt)
  setnames(dtuse
           ,old=c(cellid_col, sdimx_col, sdimy_col)
           ,new=c("cell_ID", "sdimx", "sdimy"))
  neighbs <- 
    RANN::nn2(dtuse[,c("sdimx","sdimy"),with=F]
              ,dtuse[,c("sdimx","sdimy"),with=F]
              ,k=min(200, nrow(dtuse))
              ,searchtype = 'radius'
              ,radius=radius)$nn.idx
  
  nzcol <- which(Matrix::colSums(neighbs)!=0)
  neighbs <- neighbs[,nzcol]
  dtuse[,from:=1:.N]
  idx <- cbind(dtuse[,c(cellid_col, "from"),with=FALSE], neighbs)
  midx <- melt(idx, id.vars=c(cellid_col, "from")
               ,value.name="to"
  )[from!=to][to!=0]
  
  midx[,`:=`(minid=from,maxid=to)]
  midx[to < from,`:=`(minid=to,maxid=from)]
  midx <- midx[minid==to]
  midx <- merge(midx, dtuse, by.x=c("from",cellid_col), by.y=c("from",cellid_col)) #,suffixes = c("_from", "_to"))
  midx <- merge(midx, dtuse, by.x="to", by.y="from",suffixes = c("_from", "_to"))
  midx[,`:=`(maxid=NULL,minid=NULL,variable=NULL,from=NULL,to=NULL)]
  midx[,distance:=sqrt((sdimx_to-sdimx_from)^2 + (sdimy_to-sdimy_from)^2)]
  midx[,weight:=1/distance] 
  setnames(midx
           ,old=c("cell_ID_to", "sdimx_to", "sdimy_to"
                  ,"cell_ID_from", "sdimx_from", "sdimy_from")
           ,new=c("to", "sdimx_end", "sdimy_end"
                  ,"from","sdimx_begin", "sdimy_begin"))
  setcolorder(midx
              ,c('from','to','sdimx_begin','sdimy_begin'
                 ,'sdimx_end','sdimy_end','distance','weight'))
  return(midx)
}