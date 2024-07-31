rm(list = ls())
library(umap)
library(pheatmap)
library(spatstat)
library(spatstat.geom)
library(InSituCor)
library(smiDE)

#### load data -------------------- 
fixedprofiles = readRDS("processed_data/fixedprofiles.RDS")
cols <- readRDS("processed_data/cols.RDS")
load("processed_data/cleaneddata.RData")
rownames(customlocs) = annot$cell_ID
clusts = clust

#### run pre-DE: -------------------------------

# addendum to standard pre-DE: also screen out low expressers:
meanclusterexpression <- InSituType::getRNAprofiles(x = Matrix::t(raw), 
                                                    neg = annot$negmean, 
                                                    clust = clust)
colMeans(meanclusterexpression < 0.05)

# get contamination ratio metric:
neighbors <- InSituCor:::nearestNeighborGraph(x = customlocs[, 1], y = customlocs[, 2], N = 50, subset = annot$tissuename)
neighborssymm <- 1 * ((neighbors + Matrix::t(neighbors)) != 0)
rownames(neighborssymm) <- colnames(neighborssymm) <- colnames(raw)
contam <- smiDE:::contamination_ratio_metric(assay_matrix = raw,
                                    metadata = data.frame(sdimx = customlocs[, 1],
                                                          sdimy = customlocs[, 2],
                                                          cell_ID = annot$cell_ID,
                                                          tissue = annot$tissuename,
                                                          clust = clust),
                                    adjacency_matrix = neighborssymm,
                                    cluster_col = c("clust"),
                                    cellid_col = "cell_ID",
                                    grouping_col = NULL, # "tissue",
                                    sdimx_col = "sdimx",
                                    sdimy_col = "sdimy",
                                    verbose = TRUE)
contam0 = contam
contam = as.data.frame(contam)
# convert to matrix as expected by later code:
mean_to_neighborhood_ratio <- matrix(NA, length(unique(contam$target)), length(unique(contam$clust)),
                                     dimnames = list(unique(contam$target), unique(contam$clust)))
for (cell in unique(contam$clust)){
  print(cell)
  mean_to_neighborhood_ratio[, cell] <- (contam$ratio[contam$clust == cell])[
    match(rownames(mean_to_neighborhood_ratio), contam$target[contam$clust == cell])] 
}

colMeans(mean_to_neighborhood_ratio[, c("Podocyte", "Mesangial.cell", "macrophage")] > 1.5)
toomuchcontam <- mean_to_neighborhood_ratio > 1.5
saveRDS(toomuchcontam, file = "processed_data/toomuchcontam.RDS")
saveRDS(mean_to_neighborhood_ratio, file = "processed_data/mean_to_neighborhood_ratio.RDS")

hist(colMeans(toomuchcontam))
barplot(colMeans(toomuchcontam), las = 2)
names(which(toomuchcontam[, "Podocyte"]))
names(which(toomuchcontam[, "Mesangial.cell"]))
names(which(toomuchcontam[, "macrophage"]))


# now use pre-DE to get matrices of neighbor expression:

predefile <- "processed_data/prede.RDS"
if (!file.exists(predefile)) {
  prede <- smiDE::pre_de(counts = raw, 
                         ref_celltype = c("Podocyte", "Mesangial.cell", "Glomerular.endothelium", "macrophage", "plasmablast", "PCT"), 
                         cell_type_metadata_colname = "clust",
                         metadata = data.frame(sdimx = customlocs[, 1],
                                               sdimy = customlocs[, 2],
                                               cell_ID = annot$cell_ID,
                                               tissue = annot$tissuename,
                                               clust = clust))
  str(prede)
  
  # save whole object:
  saveRDS(prede, file = predefile)
  # save smaller object:
  neighborcounts.other <- prede$nblist
  neighborcounts.other$neighbor_expr_byct <- neighborcounts.other$neighbor_expr_byct["otherct"]
  saveRDS(neighborcounts.other, file = "processed_data/prede_concise.RDS")
} else {
  prede <- readRDS(predefile)
}



# prede for immune cells: (macs and CD4)
predefile.imm <- "processed_data/prede.imm.RDS"
if (!file.exists(predefile.imm)) {
  prede <- smiDE::pre_de(counts = raw, 
                         ref_celltype = c("macrophage", "T CD4 naive", "T CD4 memory", "Treg"), 
                         cell_type_metadata_colname = "clust",
                         metadata = data.frame(sdimx = customlocs[, 1],
                                               sdimy = customlocs[, 2],
                                               cell_ID = annot$cell_ID,
                                               tissue = annot$tissuename,
                                               clust = clust))
  str(prede)
  
  # save whole object:
  saveRDS(prede, file = predefile.imm)
  # save smaller object:
  neighborcounts.other <- prede$nblist
  neighborcounts.other$neighbor_expr_byct <- neighborcounts.other$neighbor_expr_byct["otherct"]
  saveRDS(neighborcounts.other, file = "processed_data/prede_concise_immune.RDS")
} else {
  prede <- readRDS(predefile.imm)
}


