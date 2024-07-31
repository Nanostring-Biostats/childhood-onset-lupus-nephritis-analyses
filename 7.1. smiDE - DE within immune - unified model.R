# Goal: for each glomerular cell type, run DE vs. various annotations
  

rm(list = ls())
library(lmerTest)
library(smiDE)
library(spatstat.geom)
library(spatstat)
library(RColorBrewer)

# load data:
gloms = readRDS("processed_data/gloms with cell annotations.RDS")
load("processed_data/cannot with distances from gloms.RData") # loads "cannot"
cols <- readRDS("processed_data/cols.RDS")
load("processed_data/cleaneddata.RData")
clusts = clust
rownames(customlocs) = annot$cell_ID
fixedprofiles <- readRDS("processed_data/fixedprofiles.RDS")

mean_to_neighborhood_ratio = readRDS("processed_data/mean_to_neighborhood_ratio.RDS")
neighborhood_counts <- readRDS("processed_data/prede_concise.RDS")
toomuchcontam <- readRDS("processed_data/toomuchcontam.RDS")

load("processed_data/polygons.RData")
#### extract targeted environment variables: -------------------------

# define spatial context in terms of number of immune cell neighbors:
#context = readRDS("processed_data/context.RDS")
pp <- spatstat.geom::ppp(customlocs[, 1], customlocs[, 2], xrange = range(customlocs[, 1]), yrange = range(customlocs[, 2]))
marks(pp) <- clust
spatstat.geom::marks(pp) <- as.factor(spatstat.geom::marks(pp))
# count neighbors of each db cluster:
context <- spatstat.core::marktable(X = pp, R = NULL, N = 50, exclude=TRUE, collapse=FALSE)
context <- as.matrix(context)
rownames(context) <- rownames(customlocs)
head(context)
context <- cbind(context, rowSums(context[, c("T CD4 naive", "T CD4 memory", "T CD8 naive", "T CD8 memory", "Treg")]))
colnames(context)[ncol(context)] <- "T"


# initialize deannot with immune context:
deannot <- cbind(rowSums(context[, c("B-cell", "plasmablast")]),
             context[, "T"],
             rowSums(context[, c("macrophage", "mDC", "monocyte")]))
# log-transform immune variables:
deannot = log2(1 + deannot)
colnames(deannot) <- c("B-plasma", "T", "MoMacDC")
# add position vs. glom:
deannot = as.data.frame(deannot)
#deannot$in.glom = 1 * (cannot$position.vs.glom == "inside glomerulous")
#deannot$near.glom = 1 * (cannot$position.vs.glom == "bordering glomerulous")
deannot$totalcounts <- annot$totalcounts
deannot$tissue = annot$tissuename
deannot$cell_ID = annot$cell_ID
deannot$path = cannot$position.vs.glom

#### run DE over all variables and cell types: -----------------------


delist = list()
if (!file.exists("processed_data/DE results - immcells.RDS")) {
  celltypes = c("macrophage")
  deres <- list()
  for (name in celltypes) {
    print(name)
    tempgenes <- names(which(!toomuchcontam[, name]))
    tempcells <- annot$cell_ID[(clust == name) & !grepl("Control", annot$tissuename)]
    # call smi_de:
    deres[[name]] <- smiDE::smi_de(
      assay_matrix = raw[tempgenes, tempcells],
      metadata = deannot[tempcells, ],
      formula = ~offset(log(totalcounts)) + RankNorm(otherct_expr) + (1|tissue) +
        path + `B-plasma` + T + MoMacDC,
      neighborhood_counts = neighborhood_counts,
      groupVar = "path",
      family="nbinom2"
    )
    
    saveRDS(deres, file = "processed_data/DE results - immcells.RDS")
  }
} else {
  deres <- readRDS("processed_data/DE results - immcells.RDS")
}


# check out volcplots for SLE vs. control:
par(mfrow = c(2,3))
for (cell in names(deres)) {
  res = results(deres[[cell]])
  for (contrastname in c("bordering glomerulous / inside glomerulous",  "bordering glomerulous / tubulointerstitium",
                         "inside glomerulous / tubulointerstitium", "B-plasma 1.992 / B-plasma 0.7662",
                         "T 0.7936 / T 0.245", "MoMacDC 3.314 / MoMacDC 2.06")) {
    temp <- res$pairwise[contrast == contrastname, ]
    plot(log2(temp$fold_change), -log10(temp$p.value), main = contrastname)
  }
}


## export volc plot values:
for (cell in names(deres)) {
  res = results(deres[[cell]])
  temp <- res$pairwise[contrast == "inside glomerulous / tubulointerstitium", ]
  temp$fdr <- p.adjust(temp$p.value, "BH")
  temp$log2foldchange = log2(temp$fold_change)
  write.csv(temp[, c("target", "log2foldchange", "p.value", "fdr", "contrast")],
            file = paste0("results/DE immune cells - glom vs tubulointerstitium - ", cell, ".csv"), row.names = F)
}


plot(temp$log2foldchange, -log10(temp$p.value), col = 0, xlab = "log2(fold-change): in glom / tubulointerstitium")
text(temp$log2foldchange, -log10(temp$p.value), temp$target, cex = 0.6)



#### spatial plots of macrophage expression --------------------------------------
temp <- res$pairwise[contrast == "inside glomerulous / tubulointerstitium"]
temp$fdr <- p.adjust(temp$p.value, "BH")
temp[order((fold_change > 1.5) * (fdr < 0.05), decreasing = TRUE)][1:20, ]
macgenes <- as.vector(temp$target)[(temp$fold_change > 1.5) & (temp$fdr < 0.05)]
macgenes

score <- Matrix::colMeans(raw[macgenes, ]) / Matrix::colMeans(raw) * mean(Matrix::colMeans(raw))
ismac = clust == "macrophage"

png("results/macrophage expression.png", width = 30, height = 30, res = 600, units = "in")
par(mar = c(0,0,0,0))
plot(customlocs, col = "grey80", pch = 16, cex = 0.05, asp = 1)  # col = cols[clust], 
points(customlocs[is.element(clust, c("Podocyte", "Mesangial.cell", "Glomerular.endothelium")), ],
       cex = 0.05, pch = 16, col = "blue")
points(customlocs[ismac, ], pch = 16, cex = 0.15,
       col = colorRampPalette(brewer.pal(8, "YlOrRd")[-1])(101)[
       #col = colorRampPalette(brewer.pal(11, "BrBg")[-(5:7)])(71)[
         1 + pmin(score[ismac] * 100, 100)
       ])
for (tiss in unique(annot$tissuename)) {
  inds = annot$tissuename == tiss
  rect(min(customlocs[inds, 1]), min(customlocs[inds, 2]), max(customlocs[inds, 1]), max(customlocs[inds, 2]))
  text(median(range(customlocs[inds,1])), max(customlocs[inds,2]), tiss)
}
dev.off()

par(mar = c(15,2,1,1))
boxplot(score[ismac] ~ deannot$path[ismac] + annot$tissuename[ismac], las = 2) 

# svgs from selected tissues:
for (tiss in c("SLE3", "SLE6", "Control2")) {
  use = annot$tissuename == tiss
  svg(paste0("results/macrophage metagene - ", tiss, ".svg"), width = 12, height = 12)
  par(mar = c(0,0,0,0))
  plot(customlocs[use, ], col = "grey80", pch = 16, cex = 0.05, asp = 1)  # col = cols[clust], 
  points(customlocs[use & is.element(clust, c("Podocyte", "Mesangial.cell", "Glomerular.endothelium")), ],
         cex = 0.05, pch = 16, col = "cornflowerblue")
  points(customlocs[use & ismac, ], pch = 16, cex = 0.15,
         col = colorRampPalette(brewer.pal(8, "YlOrRd")[-(1:2)])(101)[
           1 + pmin(score[ismac] * 100, 100)
         ])
  dev.off()

  use = annot$tissuename == tiss
  svg(paste0("results/macrophage metagene - ", tiss, " - glom polygons.svg"), width = 12, height = 12)
  par(mar = c(0,0,0,0))
  plot(customlocs[use, ], col = "grey80", pch = 16, cex = 0.05, asp = 1)  # col = cols[clust], 
  #points(customlocs[use & is.element(clust, c("Podocyte", "Mesangial.cell", "Glomerular.endothelium")), ],
  #       cex = 0.05, pch = 16, col = "cornflowerblue")
  for (name in names(gpoly.customlocs)) {
    polygon(gpoly.customlocs[[name]], border = "black", lwd = 0.5)
  }
  points(customlocs[use & ismac, ], pch = 16, cex = 0.15,
         col = colorRampPalette(brewer.pal(8, "YlOrRd")[-(1:2)])(101)[
           1 + pmin(score[ismac] * 100, 100)
         ])
  dev.off()
}

svg("results/macrophage metagene legend.svg")
frame()
legend("center", pch = 16,
       col = c("cornflowerblue", colorRampPalette(brewer.pal(8, "YlOrRd")[-(1:2)])(101)[c(1,51,101)]),
       legend = c("resident glomerular cells", paste0("macrophage: glomerulus metagene score = ", c(0, 0.5, ">1"))))
dev.off()

par(mar = c(15,2,1,1))
boxplot(score[ismac] ~ deannot$path[ismac] + annot$tissuename[ismac], las = 2) 


#### immune cell DE: in hotspots vs. in tubulointerstitium ---------------------------

# initialize deannot with immune context:
deannot <- data.frame(hotspot = c("other", "hotspot")[1 + (annot$context == "hotspot")])
deannot$hotspot <- factor(deannot$hotspot, levels = c("other", "hotspot"))
deannot$totalcounts <- annot$totalcounts
deannot$tissue = annot$tissuename
deannot$cell_ID = annot$cell_ID
rownames(deannot) = deannot$cell_ID


#### run DE over all variables and cell types: -----------------------
neighborhood_counts <- readRDS("processed_data/prede_concise_immune.RDS")


delist = list()
if (!file.exists("processed_data/DE results - immcells vs hotspots.RDS")) {
  celltypes = c("macrophage", "CD4")
  deres <- list()
  for (name in celltypes) {
    print(name)
    if (name == "macrophage") {
      tempgenes <- names(which(!toomuchcontam[, name]))
      tempcells <- annot$cell_ID[(clust == name) & !grepl("Control", annot$tissuename)]
    }
    if (name == "CD4") {
      celltypes = c("T CD4 naive", "T CD4 memory", "Treg")
      tempgenes <- names(which(rowSums(!toomuchcontam[, celltypes]) >= 3))
      tempcells <- annot$cell_ID[is.element(clust, celltypes) & !grepl("Control", annot$tissuename)]
    }
    # call smi_de:
    deres[[name]] <- smiDE::smi_de(
      assay_matrix = raw[tempgenes, tempcells],
      metadata = deannot[tempcells, ],
      formula = ~offset(log(totalcounts)) + RankNorm(otherct_expr) + #(1|tissue) +
        hotspot,
      neighborhood_counts = neighborhood_counts,
      groupVar = "hotspot",
      family="nbinom2"
    )
    
    saveRDS(deres, file = "processed_data/DE results - immcells vs hotspots.RDS")
  }
} else {
  deres <- readRDS("processed_data/DE results - immcells vs hotspots.RDS")
}



## export volc plot values:
for (cell in names(deres)) {
  res = results(deres[[cell]])
  temp <- res$pairwise[contrast == "other / hotspot", ]
  temp$fdr <- p.adjust(temp$p.value, "BH")
  temp$log2foldchange = log2(temp$fold_change)
  write.csv(temp[, c("target", "log2foldchange", "p.value", "fdr", "contrast")],
            file = paste0("results/DE immune cells - hotspot vs elsewhere - ", cell, ".csv"), row.names = F)
}

