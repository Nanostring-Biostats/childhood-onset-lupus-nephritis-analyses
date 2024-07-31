if (FALSE) {
  devtools::install("cellPoly")
}
rm(list = ls())
library(umap)
library(scales)
library(pheatmap)
library(spatstat)
library(spatstat.geom)
source("cell_proximity.R")
library(CellPoly)

#### load data -------------------- 
gloms = readRDS("processed_data/gloms with cell annotations.RDS")
cols <- readRDS("processed_data/cols.RDS")
load("processed_data/cleaneddata.RData")
rownames(customlocs) = annot$cell_ID
load("processed_data/cannot with distances from gloms.RData") # loads "cannot"
load("processed_data/polygons.RData")

#### categorize spatial context: --------------------------

# positition vs glom:
annot$context = cannot$position.vs.glom
# try dbscan on imm cells:
immcells = c("B-cell","macrophage","mast","mDC","monocyte"  ,"neutrophil",
             "NK","pDC","plasmablast","T CD4 memory" ,"T CD4 naive" , "T CD8 memory",
             "T CD8 naive" , "Treg")
isimmune = is.element(clust, immcells)
db = dbscan::dbscan(customlocs[isimmune, ], eps = 0.05, minPts = 5)   #<----- need to check behavior of this - getting weird result: no clusters if minPts = 20
bigclusters = setdiff(names(table(db$cluster))[table(db$cluster) > 50], "0")
annot$inhotspot = 0
annot$inhotspot[isimmune][is.element(db$cluster, bigclusters)] = 1

tempcols = alpha(cols, 0.4)
tempcols[immcells] = cols[immcells]
png(paste0("results/immune hostspot viz.png"), width = 14, height = 14, res = 300, units = "in")
par(mar = c(0,0,0,0))
plot(viz, pch = 16, cex = 0.1, col = tempcols[clust], 
     asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")

points(viz[annot$inhotspot == 1, ], pch = 1, cex = .2, col = scales::alpha("red", 0.3))
dev.off()


annot$context[(annot$inhotspot == 1) & (annot$context !="inside glomerulous")] = "hotspot"
saveRDS(annot$context, file = "processed_data/context4.RDS")

contextcols = c("inside glomerulous" = "dodgerblue4", 
                "bordering glomerulous" = "dodgerblue2", 
                "tubulointerstitium" = "peachpuff3", 
                "hotspot"= "red")

tempcols = alpha(cols, 0.3)
#tempcols[immcells] = contextcols[annot$context[immcells]]
png(paste0("results/immune viz - color by spatial context.png"), width = 14, height = 14, res = 300, units = "in")
par(mar = c(0,0,0,0))
plot(viz, pch = 16, cex = 0.1, col = tempcols[clust], 
     asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
for (i in 1:length(gpoly.viz)) {
  polygon(gpoly.viz[[i]], border = "black", lwd = 0.3)
}
points(viz[isimmune, ], cex = 0.1, pch =16, col = contextcols[annot$context[isimmune]])
dev.off()


## plot spatial context for manuscript figure:
use = (annot$tissuename == "SLE4") & (annot$fov == "SP21_213_R1080_S15")
#png("results/spatial context example.png", width = 6.3, height = 2, units = "in", res = 600)
svg("results/spatial context example.svg", width = 6.3, height = 2)
par(mfrow = c(1,2))
par(mar = c(0,.5,0,.5))
plot(customlocs[use, ], cex = 0.3,
     col = contextcols[annot$context[use]], 
     pch = 16, asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
for (i in 1:length(gpoly.customlocs)) {
  polygon(gpoly.customlocs[[i]], border = "black", lwd = 1)
}
legend("topleft", pch = 16, col = contextcols, legend = c(names(contextcols[1:3]), "immune cell in hotspot"), cex = 0.5, bty = "n")
plot(customlocs[use, ], cex = 0.3,
     col = cols[clust[use]], 
     pch = 16, asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
for (i in 1:length(gpoly.customlocs)) {
  polygon(gpoly.customlocs[[i]], border = "black", lwd = 1)
}
dev.off()



#### summarize cell type abundance by tissue ----------------------------


### high-level abundances:
# prop of immune cells within SLE only:
sleimmune = isimmune & grepl("SLE", annot$tissuename)
table(clust[sleimmune]) / sum(sleimmune) * 100
# each cell type's spread by spatial context:
tab = table(clust[sleimmune], annot$context[sleimmune])
sweep(tab, 1, rowSums(tab), "/") * 100
write.csv(tab, file = "results/immune cells per spatial context - SLE only.csv")

# immune composition of controls:
controlimmune = isimmune & !grepl("SLE", annot$tissuename)
table(clust[controlimmune]) / sum(controlimmune) * 100
# each cell type's spread by spatial context:
tab = table(clust[controlimmune], annot$context[controlimmune])
sweep(tab, 1, rowSums(tab), "/") * 100
write.csv(tab, file = "results/immune cells per spatial context - controls only.csv")


#T-cells:
tcellnames = c("T CD4 memory", "T CD4 naive" ,"T CD8 memory", "T CD8 naive", "Treg")
sum(table(clust[sleimmune])[tcellnames] / sum(sleimmune) * 100)
colSums(tab[tcellnames, ]) / sum(colSums(tab[tcellnames, ])) * 100 


# immune counts per context * tissue:
tab = table(clust[isimmune], paste0(annot$tissuename, "_", annot$context)[isimmune])
write.csv(tab, file = "results/immune cells per tissue x spatial context.csv")
tab = table(clust, paste0(annot$tissuename, "_", annot$context))
write.csv(tab, file = "results/all cells per tissue x spatial context.csv")


#### immune neighbor details --------------
use = TRUE
pp <- spatstat.geom::ppp(customlocs[use, 1], customlocs[use, 2],
                         xrange = range(customlocs[use, 1]), yrange = range(customlocs[use, 2]))
marks(pp) <- clust[use]
marks(pp) = as.factor(marks(pp))
# count neighbors of each db cluster:
mt50 <- marktable(X = pp, R = NULL, N = 50, exclude=TRUE, collapse=FALSE)

# plasmablast neighbor details:
mean(rowSums(mt50[sleimmune & (clust == "plasmablast"), immcells]))
mean(rowSums(mt50[sleimmune & (clust == "plasmablast"), immcells]) > 0)
#mean(rowSums(mt50[sleimmune & (clust == "plasmablast"), tcellnames]))
#mean(rowSums(mt50[sleimmune & (clust == "plasmablast"), tcellnames]) > 0)


# T-cell neighbor details:
is.tcell = is.element(clust, tcellnames)
mean(rowSums(mt50[sleimmune & is.tcell, tcellnames]) == 0)
mean(rowSums(mt50[sleimmune & is.tcell, tcellnames]) > 2)

# b-cells neighbor details:
tab["B-cell", ] / sum(tab["B-cell", ]) * 100
mean(mt50[sleimmune & (clust == "B-cell"), "B-cell"] > 0)
mean(rowSums(mt50[sleimmune & (clust == "B-cell"), immcells]) > 0)

# CD4 neighbor details:
thesecelltypes = c("T CD4 naive", "T CD4 memory", "Treg")
mean(rowSums(mt50[sleimmune & is.element(clust, thesecelltypes), thesecelltypes]) > 0)
mean(rowSums(mt50[sleimmune & is.element(clust, thesecelltypes), immcells]) > 0)

# CD8 neighbor details:
thesecelltypes = c("T CD8 naive", "T CD8 memory")
mean(rowSums(mt50[sleimmune & is.element(clust, thesecelltypes), thesecelltypes]) > 0)
mean(rowSums(mt50[sleimmune & is.element(clust, thesecelltypes), immcells]) > 0)

# hotspot summary:
table(db$cluster, annot$tissuename[isimmune])[bigclusters, ]
mean(annot$context[grepl("SLE", annot$tissuename)] == "hotspot")
mean(annot$context[sleimmune] == "hotspot")

tab = table(is.element(clust, c("T CD4 naive", "T CD4 memory", "Treg")), annot$context == "hotspot")
tab[2,2] / sum(tab[2, ])
tab = table(is.element(clust, c("T CD8 naive", "T CD8 memory")), annot$context == "hotspot")
tab[2,2] / sum(tab[2, ])

#### figures --------------------------------------------

## 2a: cell type abundance by patient
colsT = cols
#colsT["Treg"] = "pink"
colsT["CD8 T-cell"] = "red"
colsT["CD4 T-cell"] = "orange"
colsT["B-cell"] = "blue"
clustsT = clust
clustsT[grepl("CD8", clust)] = "CD8 T-cell"
clustsT[grepl("CD4", clust)] = "CD4 T-cell"
clustsT[clust == "Treg"] = "CD4 T-cell"
immcellsT = c(immcells[!grepl("CD", immcells)], "CD4 T-cell", "CD8 T-cell")
immcellsT = immcellsT[c(2,4,5,6,8,9,1,12,11,10,7,3)]
immcellsT = intersect(immcellsT, clustsT)
mat0 = table(clustsT, annot$tissuename)
mat = sweep(mat0, 2, colSums(mat0), "/")[immcellsT, ]

svg("results/fig_immune a - abundance barplots.svg", width = 2.5, height = 3.5)
layout(mat = c(1,2), height = c(.7, 2.5))
par(mar = c(.7,3,0,0))
# abundances:
barplot(mat, col = colsT[immcellsT],
        border = NA,
        names.arg = rep("", ncol(mat)),
        yaxt = "n", 
        ylab = "Total cells")
axis(2, at = c(0, 0.03, 0.06), las = 2, cex.axis = 0.7)
# proportions:
par(mar = c(5,3,0,0))
barplot(sweep(mat, 2, colSums(mat), "/"), col = colsT[immcellsT],
        border = NA,
        names.arg = colnames(mat),
        cex.names = 0.75,
        las = 2,
        yaxt = "n",
        ylab = "Proportion")
axis(2, at = c(0, 0.25, 0.5, 0.75, 1), las = 2, cex.axis = 0.7)
dev.off()

svg("results/fig_immune a - legend.svg", height = 2.5, width = 1.25)
par(mar = c(0,0,0,0))
frame()
legend("center", fill = colsT[immcellsT], legend = names(colsT[immcellsT]), bty = "n", cex = 0.8)
dev.off()


## barplots of imm cells per context

contextcols2 = c("inside glomerulous" = "red", 
                "bordering glomerulous" = "orange", 
                "tubulointerstitium" = "grey60", 
                "hotspot"= "blue")

svg("results/fig_immune b - context barplots.svg", width = 2.5, height = 4.5)
par(mar = c(1,3.5,1,0))
layout(mat = c(1:5), height = c(1,1,1,1,.9))
for (cell in c("macrophage", "plasmablast", "B-cell", "T-cell")) { # "CD8 T-cell", "CD4 T-cell")) {
  if (cell == "T-cell") {
    use = is.element(clustsT, c("Treg", "CD8 T-cell", "CD4 T-cell"))
  } else {
    use = clustsT == cell
  }
  
  mat = table(annot$context[use], annot$tissuename[use])
  mat = sweep(mat, 2, table(annot$tissuename)[colnames(mat)], "/")
  mat = mat[c(3,1,4,2), ]
  bp=barplot(mat, 
          main = paste0(cell, "s"),
          col = contextcols[rownames(mat)], 
          las = 2,
          names.arg = rep("", ncol(mat)),
          cex.axis = 0.7, cex.names = 0.75, cex.main = 1)
}
par(mar = c(6,3.5,0,0))
plot(0,0,xaxt = "n", yaxt = "n", col = 0,xlab = "", ylab = "", bty = "n", xlim = c(0,ncol(mat)+2.2))
axis(1, at = bp-0.2, las = 2, cex.axis = 0.8, labels = colnames(mat))
dev.off()

svg("results/fig_immune b - legend.svg", height = 1, width = 1.75)
par(mar = c(0,0,0,0))
frame()
legend("center", fill = rev(contextcols), 
       legend = rev(c("immune foci", names(contextcols)[2:4])), 
       bty = "n", cex = 0.8)
dev.off()


## 2b: boxplot of #macs per glom x tissue
macsperglom = table(cannot$closest.glom[(clust == "macrophage")])
macsperglom = macsperglom[match(rownames(gloms), names(macsperglom))]
macsperglom[is.na(macsperglom)] = 0

svg("results/immune fig c - macs per glom.svg", height = 3, width = 3.5)
par(mar = c(5,5,0,0))
boxplot(macsperglom ~ unique(annot$tissuename)[match(gloms$tissue, unique(annot$tissuename))], 
        ylab = "# macrophages\nin/bordering glomerulus", xlab = "", las = 2,
        outline = F, ylim = c(0, max(macsperglom)),
        cex.axis = 0.7, cex.lab = 0.75)
points(jitter(as.numeric(as.factor(unique(annot$tissuename)[match(gloms$tissue, unique(annot$tissuename))]))), macsperglom,
       pch = 16, col = scales::alpha("dodgerblue4", 0.75), cex = .75)
dev.off()



#### plot of immune cells in xy for a single sample --------------------

#tissname = "SLE8.1"
tissname = "SLE4"

use = annot$tissuename == tissname

tempcols = colsT
tempcols[!is.element(names(tempcols), immcellsT)] = scales::alpha(tempcols[!is.element(names(tempcols), immcellsT)], 0.5)

svg(paste0("results/immune viz - ", tissname, ".svg"), width = 9.5, height = 4)
par(mar = c(0,0,0,0))
plot(viz[use, ], pch = 16, cex = 0.25 + 0.2 * is.element(clustsT[use], immcellsT), col = tempcols[clustsT[use]],
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", asp = 1)
#xlims = range(tdf$x[is.element(tdf$cell_id, tempids)])
#ylims = range(tdf$y[is.element(tdf$cell_id, tempids)])
lines(c(9,10), rep(-14.55, 2))
text(9.5,-14.35, "1 mm")
dev.off()
