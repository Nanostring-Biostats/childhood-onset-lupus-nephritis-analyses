# annotate cells by glom characteristics
# output: "cannot", a data frame storing 1. which glom each cell is closest to, and 2. where it is within the glom (within, bordering, near)

library(spatstat.geom)
library(spatstat)

rm(list = ls())
gloms = readRDS("processed_data/gloms.rds")
db = readRDS("processed_data/db.rds")
cols <- readRDS("processed_data/cols.RDS")
load("processed_data/cleaneddata.RData")
rownames(customlocs) = annot$cell_ID
rownames(viz) = annot$cell_ID

#### new attempt 5-12: annotate cells by position vis-a-vis glom --------------------------------------------

### for each glom, add neighboring Parietal.epithelium cells to the db clust:
## for each parietal epi cell, record nearest neighbor from !is.na(db) cells, and map to it
# point pattern of pe cells and db!=NA cells:
#use = is.element(names(clust), names(db)[!is.na(db)]) | (clust == "Parietal.epithelium")
use = rep(TRUE, nrow(customlocs)); names(use) = rownames(customlocs)
pp <- spatstat.geom::ppp(customlocs[use, 1], customlocs[use, 2],
                         xrange = range(customlocs[use, 1]), yrange = range(customlocs[use, 2]))
marks(pp) <- db[match(names(which(use)), names(db))]
#marks(pp)[is.na(marks(pp))] = "na"
marks(pp) = as.factor(marks(pp))
# count neighbors of each db cluster:
mt05 <- marktable(X = pp, R = 0.05, N = NULL, exclude=TRUE, collapse=FALSE)
mt01 <- marktable(X = pp, R = 0.01, N = NULL, exclude=TRUE, collapse=FALSE)
rownames(mt01) <- rownames(mt05) <- names(which(use))

## assign cells to their most neighboring glom:
# which glom wins the vote for each cell?
closestglom <- colnames(mt05)[apply(mt05,1,which.max)]
closestglom[rowSums(mt05) == 0] = NA
names(closestglom) = rownames(mt05)


## assign parietal epi cells to gloms with a more stringent radius:
closestglom.pe <- colnames(mt01)[apply(mt01,1,which.max)]
closestglom.pe[rowSums(mt01) == 0] = NA

is.glom.or.pe = is.element(names(clust), names(db)[!is.na(db)]) | (clust == "Parietal.epithelium")
closestglom.pe[!is.glom.or.pe] = NA
names(closestglom.pe) = rownames(mt01)

# made db object that includes pe cells along with the previous glom cells:
db.pe = db
db.pe[names(which(!is.glom.or.pe))] = closestglom.pe[names(which(!is.glom.or.pe))]
#db.pe[names(closestglom.pe)] = closestglom.pe
# db.pe is now a list of glom membership that INCLUDES parietal epi cells


### start new cannot:
cannot = data.frame(id = annot$cell_ID)
rownames(cannot) = cannot$id

# variable for members of a glom (glomerular cells only):
cannot$in.glom = NA
cannot[names(db.pe), "in.glom"] = db.pe
# remove glom ids not present in "gloms" object:
removedgloms = setdiff(db.pe, c(rownames(gloms), NA))
print("removing these gloms from db.pe:")
print(removedgloms)
cannot$in.glom[is.element(cannot$in.glom, removedgloms)] = NA
table(clust[!is.na(cannot$in.glom)])

### for each glom, define a polygon:
gpolys = list()
for (gid in setdiff(unique(db.pe), NA)) {
  tempcells = names(db.pe)[(db.pe == gid) & !is.na(db.pe)]
  tempinds = chull(customlocs[tempcells, ])
  tempinds = c(tempinds, tempinds[1])
  gpolys[[gid]] = customlocs[tempcells[tempinds], ]
}

# test the polys:
png(paste0("results/glom polys test xy.png"), width = 17, height = 17, res = 600, units = "in")
par(mar = c(0,0,0,0))
plot(customlocs, pch = 16, 
     cex = 0.01 + 0.05 * is.element(clust, c("Podocyte", "Mesangial.cell","Glomerular.endothelium","Parietal.epithelium")),
     col = cols[clust], 
     asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
for (i in 1:length(gpolys)) {
  polygon(gpolys[[i]], border = "black", lwd = 0.3)
}
dev.off()


### record all cells within each polygon:
cannot$inside.glom = NA
for (i in 1:length(gpolys)) {
  if (i%%20==0){print(i)}
  insidethispoly = sp::point.in.polygon(point.x = customlocs[, 1], 
                                        point.y = customlocs[, 2], 
                                        pol.x = gpolys[[i]][, 1], 
                                        pol.y = gpolys[[i]][, 2], mode.checked=FALSE)
  cannot$inside.glom[insidethispoly > 0] = names(gpolys)[i]
}
#table(cannot$inside.glom)
#table(is.na(cannot$in.glom))
#table(is.na(cannot$inside.glom)) # great

### for each cell: what glom is it near?
cannot$closest.glom = closestglom

cannot$position.vs.glom = "tubulointerstitium"  
cannot$position.vs.glom[!is.na(cannot$closest.glom)] = "bordering glomerulous"
cannot$position.vs.glom[!is.na(cannot$inside.glom)] = "inside glomerulous"

png(paste0("results/glom position.vs.glom test xy.png"), width = 17, height = 17, res = 600, units = "in")
par(mar = c(0,0,0,0))
plot(customlocs, 
     pch = 16, #c(1, 16, 2)[(cannot$position.vs.glom == "inside glomerulous") + 2*(cannot$position.vs.glom == "bordering glomerulous") + 3*(cannot$position.vs.glom == "tubulointerstitium")], 
     cex = 0.02,
     col = c("red", "blue", "grey80")[
       (cannot$position.vs.glom == "inside glomerulous") + 2*(cannot$position.vs.glom == "bordering glomerulous") + 3*(cannot$position.vs.glom == "tubulointerstitium")
     ],
     #col = cols[clust], 
     asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
for (i in 1:length(gpolys)) {
  polygon(gpolys[[i]], border = "black", lwd = 0.3)
}
dev.off()

# save
save(cannot, file = "processed_data/cannot with distances from gloms.RData")
write.csv(cannot, file = "processed_data/cell positions wrt gloms.csv", row.names = F)


#### annotate gloms further based on their cellular context --------------------------------

# count immune cells inside/nearby: (split by T&NK/B/plasma/myl)
gloms$T.in = gloms$T.near = gloms$T.out = gloms$T.tot = 
  gloms$B.in = gloms$B.near = gloms$B.out = gloms$B.tot = 
  gloms$myeloid.in = gloms$myeloid.near = gloms$myeloid.out = gloms$myeloid.tot =
  gloms$plasma.in = gloms$plasma.near = gloms$plasma.out = gloms$plasma.tot = NA

celllists = list("T" = c("Treg", "T CD4 memory", "T CD4 naive", "T CD8 memory", "T CD8 naive"),
                 "B" = "B-cell",
                 "myeloid" = c("mDC", "monocyte", "macrophage"),
                 "plasma" = "plasmablast")
for (cell in names(celllists)) {
  print(cell)
  for (gid in rownames(gloms)) {
    is.in = ((cannot$closest.glom == gid) & (cannot$position.vs.glom == "inside glomerulous"))
    is.near = ((cannot$closest.glom == gid) & (cannot$position.vs.glom == "bordering glomerulous"))
    is.out = ((cannot$closest.glom == gid) & (cannot$position.vs.glom == "tubulointerstitium"))
    gloms[gid, paste0(cell, ".in")] = sum(is.element(clust, celllists[[cell]]) & is.in, na.rm = T)
    gloms[gid, paste0(cell, ".near")] = sum(is.element(clust, celllists[[cell]]) & is.near, na.rm = T)
    gloms[gid, paste0(cell, ".out")] = sum(is.element(clust, celllists[[cell]]) & is.out, na.rm = T)
  }
  gloms[[paste0(cell, ".tot")]] = rowSums(gloms[, paste0(cell, ".", c("in", "near", "out"))])
}

# record proportion of glom cell types:
glom.cell.types = c("Podocyte", "Mesangial.cell", "Glomerular.endothelium")
other.endo.cell.types = setdiff(unique(clust)[grepl("endo", unique(clust))], "Glomerular.endothelium")
gloms$prop.Podocyte = gloms$prop.Mesangial.cell = 
  gloms$prop.Glomerular.endothelium = gloms$other.endo = NA

for (cell in glom.cell.types) {
  print(cell)
  for (gid in rownames(gloms)) {
    thisglom = (cannot$in.glom == gid) & !is.na(cannot$in.glom)
    gloms[gid, paste0("prop.", cell)] = mean(clust[thisglom] == cell)
  }
}
for (gid in rownames(gloms)) {
  thisglom = (cannot$in.glom == gid) & !is.na(cannot$in.glom)
  gloms[gid, "other.endo"] = mean(is.element(clust[thisglom], other.endo.cell.types))
}

saveRDS(gloms, file = "processed_data/gloms with cell annotations.RDS")

# sanity check:
if (FALSE) {
  gloms0 = readRDS("processed_data/gloms with cell annotations - old.RDS")
  identical(gloms$prop.Glomerular.endothelium, gloms0$prop.Glomerular.endothelium)
  identical(gloms$prop.Mesangial.cell, gloms0$prop.Mesangial.cell)
  identical(gloms$other.endo, gloms0$other.endo)
  identical(gloms$tissue, gloms0$tissue)
  identical(gloms$disease, gloms0$disease)
} #OK.

#### define glom bounding polygons -----------------------
gpoly.viz = gpoly.customlocs = list()
for (gid in gloms$glomnames) {
  tempcells = rownames(cannot)[!is.na(cannot$in.glom) & (cannot$in.glom == gid)]
  # what's the majority FOV?
  mainfov = names(which.max(table(annot$fov[match(tempcells, annot$cell_ID)])))
  tempcells = tempcells[annot$fov[match(tempcells, annot$cell_ID)] == mainfov]
  
  tempxy = viz[tempcells, ][chull(viz[tempcells, ]), ]
  gpoly.viz[[gid]] = tempxy #polygon(tempxy)
  
  tempxy = customlocs[tempcells, ][chull(customlocs[tempcells, ]), ]
  gpoly.customlocs[[gid]] = tempxy #polygon(tempxy)
  
  #if (!is.na(vizgrid[tempcells, 1])) {
  #  tempxy = vizgrid[tempcells, ][chull(vizgrid[tempcells, ]), ]
  #  gpoly.vizgrid[[gid]] = tempxy #polygon(tempxy)
  #}
}
save(gpoly.viz, gpoly.customlocs, file = "processed_data/polygons.RData")

  
length(gpoly.viz)
length(gpoly.customlocs)
c(sapply(gpoly.viz, function(x){(x[1,1])}))


