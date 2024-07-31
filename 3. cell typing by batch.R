
library(InSituCor)
library(RColorBrewer)
library(pheatmap)
library(rlist)
library(InSituType)

# load necessary data:
load("processed_data/data_with_failed_slide.RData")
viz <- readRDS("processed_data/viz.RDS")
fixedprofiles = readRDS("processed_data/fixedprofiles.RDS")
counts <- Matrix::t(raw)
mat = fixedprofiles[intersect(rownames(fixedprofiles), colnames(counts)), ]
rm(raw)
pdf("results/profiles heatmap.pdf")
pheatmap(sweep(mat, 1, pmax(apply(mat,1,max), 200),"/"),
         col = colorRampPalette(c("white", "darkblue"))(101), show_rownames = F)
dev.off()

#### get alternative data types for cell typing: -------------------

## get spatial context:
# get neighbor relationships:
neighborpcsfile <- "processed_data/neighborhoodPCs.RDS"
if (!file.exists(neighborpcsfile)) {
  neighbors <- InSituCor:::nearestNeighborGraph(x = customlocs[, 1], y = customlocs[, 2], N = 50, subset = annot$tissue) 
  # get mean neighborhood expression:
  neighborexpression <- InSituCor:::get_neighborhood_expression(counts = counts, neighbors = neighbors) 
  # save top 20 PCs:
  neighborhoodPCs <- irlba::prcomp_irlba(neighborexpression, n = 10)$x
  saveRDS(neighborhoodPCs, file = neighborpcsfile)
} else {
  neighborhoodPCs <- readRDS(neighborpcsfile)
}
dim(neighborhoodPCs)

# define cohorts:
cohortsfile <- "processed_data/cohorts.RDS"
if (!file.exists(cohortsfile)) {
  # get IF data:
  ifmat = as.matrix(annot[, c("Area", paste0("Mean.", c("CD45", "PanCK", "DAPI")))])
  # cluster to get cohorts:
  cohort = fastCohorting(mat = cbind(ifmat, neighborhoodPCs))
  saveRDS(cohort, file = cohortsfile)
} else {
  cohort <- readRDS(cohortsfile)
}

#### supervised cell typing: ----------------------------------

# genes that misbehave in the controls:
blacklist = c("NPPC", "WIF1", "DUSP5", "AZU1", "COL9A2", "FKBP11", "AZGP1", "KRT18", "MZT2A")
goodgenes = setdiff(colnames(counts), blacklist)


slidelist = list()
slidelist[["newcontrols"]] = c("SMI0016C_SP17SP19", "SMI0016c_SP212319" , "SMI0016c_SP2219")
slidelist[["oldcontrol"]] = c("SP22-1068-A2-S1", "SP22-1068-A2-S2")
slidelist[["sle"]] = setdiff(unique(annot$slidename), unlist(slidelist))

# insitutype for each:
reslist = list()
for (name in names(slidelist)) {
  print(name)
  inds = is.element(annot$slidename, slidelist[[name]])
  reslist[[name]] = insitutype(x = counts[inds, goodgenes],
                               assay_type = "rna",
                               neg = annot$negmean[inds],
                               n_clusts = 0,
                               cohort = cohort[inds],
                               reference_profiles = mat,
                               update_reference_profiles = TRUE,
                               rescale = TRUE, refit = FALSE)
}
saveRDS(reslist, file = "processed_data/cell typing by batch.RDS")
reslist = readRDS(file = "processed_data/cell typing by batch.RDS")


#### refine clusters in the first batch of data -------------------------

name = "sle"

# first delete the failed clusters:
set.seed(0)
refinedsle <- refineClusters(logliks = reslist$sle$logliks, 
                             #to_delete = c("DCT", "Parietal.epithelium"),
                             subcluster = list("Myofibroblast" = 4,      # break into its 3 subtypes
                                               "Distinct.proximal.tubule.1" = 3),   # recover the parietal epi cells lost in this cluster
                             counts = counts[rownames(reslist$sle$logliks), ],
                             neg = annot[reslist$sle$logliks, "negmean"])
inds = is.element(annot$slidename, slidelist[["sle"]])
table(refinedsle$clust)

png("results/PCT xy.png", width = 30, height = 40, res = 500, units = "in")
par(mar = c(0,0,0,0))
plot(viz[inds, ], col = c(NA, "red", "blue", "chartreuse3", "black")[
  1 + 1*(refinedsle$clust == "Distinct.proximal.tubule.1_1") + 2*(refinedsle$clust == "Distinct.proximal.tubule.1_2") + 
    3*(refinedsle$clust == "Distinct.proximal.tubule.1_3") + 4*(refinedsle$clust == "Podocyte")], cex = 0.2, pch = 16,asp= 1)
dev.off()

png("results/myo xy.png", width = 30, height = 40, res = 500, units = "in")
plot(viz[inds, ], col = c(NA, "red", "blue", "chartreuse3", "darkviolet", "black")[
  1 + 1*(refinedsle$clust == "Myofibroblast_1") + 2*(refinedsle$clust == "Myofibroblast_2") + 
    3*(refinedsle$clust == "Myofibroblast_3") + 4*(refinedsle$clust == "Myofibroblast_4") +  
    5*(refinedsle$clust == "Podocyte")], cex = 0.2, pch = 16,asp= 1)
dev.off()


## decisions on re-naming:
newnames = c("Distinct.proximal.tubule.1_3" = "Parietal.epithelium",
             "Distinct.proximal.tubule.1_1" = "PCT",
             "Distinct.proximal.tubule.1_2" = "PCT",
             "Distinct.proximal.tubule.2" = "PCT",
             "Proliferating.Proximal.Tubule" = "PCT",
             "Proximal.tubule" = "PCT",
             "Myofibroblast_1" = "Mesangial.cell",
             "Myofibroblast_2" = "Vascular.pericyte",
             "Myofibroblast_3" = "Myofibroblast",
             "Myofibroblast_4" = "Myofibroblast")

if (FALSE) {
  # rename, and merge all the PCTs:
  refinedsle2 <- refineClusters(logliks = refinedsle$logliks, 
                                merges = c("Distinct.proximal.tubule.1_3" = "Parietal.epithelium",
                                           "Distinct.proximal.tubule.1_1" = "PCT",
                                           "Distinct.proximal.tubule.1_2" = "PCT",
                                           "Distinct.proximal.tubule.2" = "PCT",
                                           "Proliferating.Proximal.Tubule" = "PCT",
                                           "Proximal.tubule" = "PCT",
                                           "Myofibroblast_1" = "Mesangial.cell",
                                           "Myofibroblast_2" = "Vascular.pericyte",
                                           "Myofibroblast_3" = "Myofibroblast",
                                           "Myofibroblast_4" = "Myofibroblast"))
  table(refinedsle2$clust)
}


#### now apply the subcluster profiles from the above to the newer batches: ---------------------------

# global cluster assignments:
bclust = rep(NA, nrow(counts)); names(bclust) = rownames(counts)
inds = is.element(annot$slidename, slidelist[["sle"]])
bclust[inds] = refinedsle$clust

# subtype profiles:
myoprofiles <- refinedsle$profiles[, paste0("Myofibroblast_", 1:4)]
myoprofiles <- myoprofiles[apply(myoprofiles, 1, max) > 0.2, ]
pct1profiles <- refinedsle$profiles[, paste0("Distinct.proximal.tubule.1_", 1:3)]
pct1profiles <- pct1profiles[apply(pct1profiles, 1, max) > 0.2, ]

# apply subtyping logic to new batches:
for (name in c("oldcontrol", "newcontrols")) {
 
  inds = is.element(annot$slidename, slidelist[[name]])
  bclust[inds] = reslist[[name]]$clust
  
  # reclassify myofibroblasts:
  myoinds = which(inds)[which(bclust[inds] == "Myofibroblast")]
  myores <- insitutype(counts[myoinds, ],
                       neg = annot$negmean[myoinds],
                       n_clusts = 0,
                       cohort = cohort[myoinds],
                       reference_profiles = myoprofiles,
                       update_reference_profiles = TRUE,
                       rescale = TRUE, refit = FALSE)
  bclust[myoinds] = myores$clust                    
  
  # reclassify PCT1:
  pctinds = which(inds)[which(bclust[inds] == "Distinct.proximal.tubule.1")]
  pctres <- insitutype(counts[pctinds, ],
                       neg = annot$negmean[pctinds],
                       n_clusts = 0,
                       cohort = cohort[pctinds],
                       reference_profiles = pct1profiles,
                       update_reference_profiles = TRUE,
                       rescale = TRUE, refit = FALSE)
  bclust[pctinds] = pctres$clust                    
}


# rename:
for (oldname in names(newnames)) {
  bclust[bclust == oldname] = newnames[oldname]
}

saveRDS(bclust, file = "processed_data/clusts.RDS")

clusts = bclust



#### full QCs ----------------------------------------------

# abundances
tab = table(bclust, annot$tissuename)
pheatmap(sweep(tab, 2, colSums(tab), "/"))

# spatial QC:
png("results/PCT xy2.png", width = 30, height = 40, res = 500, units = "in")
par(mar = c(0,0,0,0))
plot(viz, col = c(NA, "red", "dodgerblue2", "black")[
  1 + 1*(bclust == "Parietal.epithelium") + 2*(bclust == "PCT") + 
    3*(bclust == "Podocyte")], cex = 0.2, pch = 16,asp= 1)
dev.off()

png("results/myo xy2.png", width = 30, height = 40, res = 500, units = "in")
plot(viz, col = c(NA, "red", "blue", "chartreuse3", "black")[
  1 + 1*(bclust == "Mesangial.cell") + 2*(bclust == "Vascular.pericyte") + 
    3*(bclust == "Myofibroblast") + 4*(bclust == "Podocyte")], cex = 0.2, pch = 16,asp= 1)
dev.off()

#### assign colors: -------------------------------


# bring cell types together:
colnames(fixedprofiles)
uniquecols = c(brewer.pal(12, "Set3")[-2], brewer.pal(8, "Set2"))
cols = iocolors
cols["mast"] = "gold"
cols[c("T CD8 memory","T CD4 memory","Treg","T CD4 naive","T CD8 naive")] = "red"
cols["NK"] = iocolors["NK"]
cols[c("Fibroblast")] = iocolors["fibroblast"]
#cols[c("Proliferating.Proximal.Tubule", "Distinct.proximal.tubule.2", "Distinct.proximal.tubule.1", "Proximal.tubule")] = "peachpuff2"
cols["PCT"] = "peachpuff2"
cols[c("Type.A.intercalated.cell", "Indistinct.intercalated.cell", "Type.B.intercalated.cell")] = uniquecols[4]
cols[c("Thick.ascending.limb.of.Loop.of.Henle")] = uniquecols[5]
cols["Epithelial.progenitor.cell"] = "#FFFF99" #uniquecols[6]
cols["Pelvic.epithelium"] = "#CCCC00" #"chartreuse2"
cols["Parietal.epithelium"] = "yellow" #"olivedrab3" #uniquecols[15]
cols[c("Connecting.tubule")] = uniquecols[7]
cols[c("Podocyte")] = "orange" # "deeppink" # "#FF9966" #uniquecols[8]
cols["Glomerular.endothelium"] = "slateblue1" # "orangered"  #uniquecols[2]
cols["Mesangial.cell"] = "orchid3" #"purple" #uniquecols[9]
cols[c("Myofibroblast")] = uniquecols[10]
cols["Vascular.pericyte"] = uniquecols[14]
cols[c("Peritubular.capillary.endothelium.2", "Ascending.vasa.recta.endothelium",
       "Peritubular.capillary.endothelium.1", "Descending.vasa.recta.endothelium")] = uniquecols[1]
cols[c("Principal.cell")] = uniquecols[3]
cols[c("Transitional.urothelium")] = "brown" #uniquecols[11]
cols["DCT"] = "wheat4"
cols = cols[is.element(names(cols), unique(clusts))]

setdiff(names(cols), unique(clusts))
setdiff(unique(clusts), names(cols))

saveRDS(cols, file = "processed_data/cols.RDS")


### colors w Tcells collapsed:
clusts.T = clusts
clusts.T[clusts.T == "Treg"] = "CD4 T-cell"
clusts.T[clusts.T == "T CD4 naive"] = "CD4 T-cell"
clusts.T[clusts.T == "T CD4 memory"] = "CD4 T-cell"
clusts.T[clusts.T == "T CD8 naive"] = "CD8 T-cell"
clusts.T[clusts.T == "T CD8 memory"] = "CD8 T-cell"
cols.T = cols[unique(clusts.T)]
cols.T["CD8 T-cell"] = "red"
cols.T["CD4 T-cell"] = "firebrick"
cols.T = cols.T[order(names(cols.T))]
cols.T = cols.T[!is.na(cols.T)]

### Shaun colors:


scols = c("Podocyte" = "#ED2526",
          "Mesangial.cell" = "#F47F21",
          "Glomerular.endothelium" = "blue",
          "Parietal.epithelium" = "#FFF100",
          "PCT" = "#C5C7C9",
          "Thick.ascending.limb.of.Loop.of.Henle" = "#2D3962",
          "Connecting.tubule" = "#00A95A",
          "Epithelial.progenitor.cell" = "#8CCDAA",
          "Type.A.intercalated.cell" = "#6CCBDE",
          "Type.B.intercalated.cell" = "#6CCBDE",
          "Indistinct.intercalated.cell" = "#6CCBDE",
          "Principal.cell" = "#760055",
          "Peritubular.capillary.endothelium.1" = "#8ED1C5",
          "Peritubular.capillary.endothelium.2" = "#8ED1C5",
          "Ascending.vasa.recta.endothelium" = "#8ED1C5",
          "Descending.vasa.recta.endothelium" = "#8ED1C5",
          "Vascular.pericyte" = "#555658",
          "Myofibroblast" = "#AED36E",
          
          "Pelvic.epithelium" = "black",
          "Fibroblast" = "tan",
          "Transitional.urothelium" = "pink")
#          "Immune.cell" = "#FFFFFF")  #"#993399"

# immune cells should be left out; give them all the same color:
setdiff(names(scols), unique(clusts.T))
setdiff(unique(clusts.T), names(scols))
scols[setdiff(unique(clusts.T), names(scols))] = "#993399"

saveRDS(scols, file = "processed_data/shaun cols.RDS")

#### summary plots: -------------------

png("results/umap by cell type - Ts condensed.png", res = 400, units = "in", width = 5, height = 5)
par(mar = c(0,0,0,0))
plot(um, pch = 16, cex = 0.1, col = cols.T[clusts.T], asp = 1,
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend("bottomleft", pch = 16, col = cols, legend = names(cols), cex = 0.5)
dev.off()
svg("results/umap by cell type - Ts condensed - color legend.svg")
frame()
legend("center", pch = 16, col = cols.T, legend = names(cols.T), cex = 0.5)
dev.off()



# svg from a single sample:
#name = "control2"
#name = "SP19_1139"
name = "SP20_642"
use = annot$tissue == name
svg(paste0("results/cell type xy - ", name, " - shaun colors.svg"), width = 10, height = 9)
plot(viz[use, ], pch = 16, col = scols[clusts.T[use]], cex = 0.025, asp = 1)
#lines(c(-2,-1), c(-11.25,-11.25))
#text(-1.5,-11.5, "1 mm")
dev.off()
svg(paste0("results/cell type xy - ", name, " - orig colors.svg"), width = 10, height = 9)
plot(viz[use, ], pch = 16, col = cols[clusts.T[use]], cex = 0.025, asp = 1)
dev.off()

png("results/cell type UMAP - Shaun colors.png", width = 8, height = 8, units = "in", res = 600)
par(mar = c(0,0,0,0))
plot(um, pch = 16, cex = 0.1, col = scols[clusts.T], 
     asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
dev.off()

svg("results/Shaun color legend.svg", height = 15, width = 5)
frame()
legend("center", pch = 16, col = scols, legend = names(scols))
dev.off()


png("results/cell type umap - with legend.png", width = 10, height = 10, res = 500, units = "in")
par(mar = c(0,0,0,0))
layout(mat = matrix(c(1,2), nrow = 1), widths = c(8,2))
plot(um, pch = 16, cex = 0.1, col = cols[clusts], 
     asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
frame()
legend("center", pch = 16, col = cols, legend = names(cols), cex = 0.55)
dev.off()
png("results/cell type umap.png", width = 10, height = 10, res = 500, units = "in")
par(mar = c(0,0,0,0))
plot(um, pch = 16, cex = 0.1, col = cols[clusts], 
     asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
dev.off()
png("results/cell type umap - with labels.png", width = 10, height = 10, res = 500, units = "in")
par(mar = c(0,0,0,0))
plot(um, pch = 16, cex = 0.1, col = cols[clusts], 
     asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
for (name in unique(clusts)) {
  text(median(um[clusts == name, 1]), median(um[clusts == name, 2]), name)
}
dev.off()

png("results/cell type xy.png", width = 30, height = 40, res = 600, units = "in")
par(mar = c(0,0,0,0))
plot(customlocs, pch = 16, cex = 0.1, col = cols[clusts], 
     asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
#legend("right", pch = 16, col = cols, legend = names(cols))
dev.off()


png("results/cell type viz.png", width = 14, height = 14, res = 600, units = "in")
par(mar = c(0,0,0,0))
plot(viz, pch = 16, cex = 0.1, col = cols[clusts], 
     asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
#legend("right", pch = 16, col = cols, legend = names(cols))
dev.off()


# color legend
legendcols = cols
legendcols["T-cell"] = legendcols["T CD4 memory"] 
legendcols = legendcols[!is.element(names(legendcols), c("T CD8 memory","T CD8 naive","T CD4 memory","T CD4 naive","Treg"))]
legendcols["Endothelial"] = legendcols["Descending.vasa.recta.endothelium"]
legendcols = legendcols[!grepl("capillary.endothelium", names(legendcols))]
legendcols = legendcols[!grepl("recta.endothelium", names(legendcols))] 
legendcols["Intercalated"] = legendcols["Indistinct.intercalated.cell"]
legendcols = legendcols[!grepl("intercalated.cell", names(legendcols))]
# reorder:
legendcols = legendcols[order(names(legendcols))]

svg("results/cell type color legend - fitted.svg", height = 6, width = 3.5)
par(mar = c(0,0,0,0))
frame()
legend("center", pch = 16, col = legendcols, legend = names(legendcols))
dev.off()


# plots per cell type:
dir.create("results/celltyping")
for (cell in unique(clusts)) {
  png(paste0("results/celltyping/", cell, ".png"), width = 14, height = 14, res = 400, units = "in")
  par(mar = c(0,0,0,0))
  plot(viz, pch = 16, cex = 0.1, col = scales::alpha(cols[clusts], 0.4),
       asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  points(viz[clusts == cell, ], pch = 16, cex = 0.1, col = "black",
         asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  for (tiss in unique(annot$tissuename)) {
    text(median(range(viz[annot$tissuename == tiss, 1])), max(viz[annot$tissuename == tiss, 2]) + 0.1, tiss, col = "black")
    rect(min(viz[annot$tissuename == tiss, 1]),
         min(viz[annot$tissuename == tiss, 2]),
         max(viz[annot$tissuename == tiss, 1]),
         max(viz[annot$tissuename == tiss, 2]))
  }
  
  dev.off()
  
}

#### summarize abundances: --------------------------

tab = table(clusts, annot$tissuename)
write.csv(tab, file = "results/cell type abundances all data including flagged tissues.csv")

badslide = "SMI0016c_SP212319"
isgood = annot$slidename != badslide
tab = table(clusts[isgood], annot$tissuename[isgood])
write.csv(tab, file = "results/cell type abundances without flagged slide.csv")
