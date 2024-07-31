# This script: annotate cells for whether they're in a glom spatial context or not. 
# No need to distinguish between individual gloms - the glom IDs used for glom-level analysis are in script 4.1.
# This is just to have a more sensitive way to annotate cells as in/out of gloms. 



rm(list = ls())
library(RColorBrewer)
library(scales)
library(dbscan)

load("processed_data/cleaneddata.RData")
cols = readRDS("processed_data/cols.RDS")
db <- readRDS(file = "processed_data/db.rds")
load("processed_data/polygons.RData")
load("processed_data/cannot with distances from gloms.RData")

reslist = readRDS(file = "processed_data/cell typing by batch.RDS")

#### explore glom cell landscape: ---------------------------------------------

png("results/glom cells and old glom annot.png", width = 30, height = 30, res = 600, units = "in")
par(mar = c(0,0,0,0))
plot(customlocs, col = alpha(cols[clust], 0.3), pch = 16, cex = 0.1, asp = 1)
points(customlocs[clust == "Podocyte", ], col = "red", pch = 16, cex = 0.2)
points(customlocs[clust == "Mesangial.cell", ], col = "blue", pch = 16, cex = 0.15)
points(customlocs[clust == "Glomerular.endothelium", ], col = "forestgreen", pch = 16, cex = 0.15)
for (cl in names(gpoly.customlocs)) {
  polygon(gpoly.customlocs[[cl]])
}
for (tiss in unique(annot$tissuename)) {
  inds = annot$tissuename == tiss
  rect(min(customlocs[inds, 1]), min(customlocs[inds, 2]), max(customlocs[inds, 1]), max(customlocs[inds, 2]))
  text(median(range(customlocs[inds,1])), max(customlocs[inds,2]), tiss)
}
dev.off()

# conclude: old glom calling was strict and missed some small gloms. 

#### more sensitive glom definition:  ---------------------------------------------


# cell type makeup of nearest neighbors:
pp <- spatstat.geom::ppp(customlocs[, 1], customlocs[, 2], xrange = range(customlocs[, 1]), yrange = range(customlocs[, 2]))
marks(pp) <- clust
spatstat.geom::marks(pp) <- as.factor(spatstat.geom::marks(pp))
# count neighbors of each db cluster:
neighbormarks <- spatstat.core::marktable(X = pp, R = NULL, N = 50, exclude=TRUE, collapse=FALSE)
neighbormarks <- as.matrix(neighbormarks)
rownames(neighbormarks) <- rownames(customlocs)
head(neighbormarks)

# number of nearby glom cells:
cannot$nglomcellneighbors <- rowSums(neighbormarks[, c("Podocyte", "Mesangial.cell", "Glomerular.endothelium")])

boxplot(cannot$nglomcellneighbors ~ clust, las = 2)

tab = table(clust, cannot$nglomcellneighbors >= 5)
barplot(tab[, "TRUE"] / rowSums(tab), las = 2)


cannot$inglomcontext <- cannot$nglomcellneighbors >= 3
tab = table(clust, cannot$inglomcontext)
round(tab[, "TRUE"] / rowSums(tab), 2)
barplot(tab[, "TRUE"] / rowSums(tab), las = 2)



inds = clust == "Podocyte"
tab = table(annot$tissuename[inds], cannot$inglomcontext[inds])
tab
tab[,2]/rowSums(tab)
inds = clust == "Mesangial.cell"
tab = table(annot$tissuename[inds], cannot$inglomcontext[inds])
tab[,2]/rowSums(tab)

# what if we only look at confident calls?
conf <- rep(NA, length(clust))
names(conf) = names(clust)
for (name in names(reslist)) {
  tempprobs = reslist[[name]]$prob
  tempprobs = tempprobs[is.element(names(tempprobs), names(conf))]
  conf[match(names(tempprobs), names(conf))] = tempprobs
}

inds = (clust == "Podocyte") & (conf > 0.95)
tab = table(annot$tissuename[inds], cannot$inglomcontext[inds])
tab
tab[,2]/rowSums(tab)



#### next attempt: same logic as before, just let dbscan find smaller gloms:
glom.cell.types = c("Podocyte", "Mesangial.cell", "Glomerular.endothelium")
endo.cell.types = unique(clust)[grepl("endo", unique(clust))]
glomcells = is.element(clust, c(glom.cell.types, endo.cell.types))

set.seed(0)
gl = dbscan(customlocs[glomcells, ],
            eps = 0.04, minPts = 10)$cluster  # was 0.04
names(gl) = names(clust)[glomcells]
gl[gl == 0] = NA

# clean up gloms without actual glom cells: (i.e clusters of only endo):
for (cl in setdiff(unique(gl), NA)) {
  if (mean(is.element(clust[names(which((gl == cl) & !is.na(gl)))], glom.cell.types)) < 0.1) {
    gl[(gl == cl) & !is.na(gl)] <- NA
  }
}


png("results/glom def - new vs. old.png", width = 30, height = 30, res = 600, units = "in")
par(mar = c(0,0,0,0))
plot(customlocs, col = cols[clust], pch = 16, cex = 0.05, asp = 1)
for (cl in setdiff(unique(db), NA)) {
  inds = which(glomcells)[(db == cl) & !is.na(db)]
  hull = chull(customlocs[inds], )
  polygon(customlocs[inds[c(hull, hull[1])], ], lwd = 0.3, border = "blue")
}
for (cl in setdiff(unique(gl), NA)) {
  inds = which(glomcells)[(gl == cl) & !is.na(gl)]
  hull = chull(customlocs[inds], )
  polygon(customlocs[inds[c(hull, hull[1])], ], lwd = 0.3, border = "red")
}
points(customlocs[clust == "Podocyte", ], col = "red", pch = 16, cex = 0.2)
points(customlocs[clust == "Mesangial.cell", ], col = "blue", pch = 16, cex = 0.15)
points(customlocs[clust == "Glomerular.endothelium", ], col = "forestgreen", pch = 16, cex = 0.15)
for (tiss in unique(annot$tissuename)) {
  inds = annot$tissuename == tiss
  rect(min(customlocs[inds, 1]), min(customlocs[inds, 2]), max(customlocs[inds, 1]), max(customlocs[inds, 2]))
  text(median(range(customlocs[inds,1])), max(customlocs[inds,2]), tiss)
}
dev.off()


#### summarize frequency in gloms under more generous definition: -------------------------

