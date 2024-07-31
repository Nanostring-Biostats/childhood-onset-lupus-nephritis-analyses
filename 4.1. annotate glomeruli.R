# create initial labelling of glomeruli, to be reviewed by pathologist then updated with a subsequent script

rm(list = ls())
library(RColorBrewer)
library(dbscan)
library(mclust)

load("processed_data/cleaneddata.RData")
cols = readRDS("processed_data/cols.RDS")

# run dbscan on glomeruli cells: 
glom.cell.types = c("Podocyte", "Mesangial.cell", "Glomerular.endothelium")
endo.cell.types = unique(clust)[grepl("endo", unique(clust))]
glomcells = is.element(clust, c(glom.cell.types, endo.cell.types))
table(glomcells)
table(is.na(glomcells))

set.seed(0)
db = dbscan(customlocs[glomcells, ],
            eps = 0.04, minPts = 20)$cluster  # was 0.04
names(db) = names(clust)[glomcells]
db[db == 0] = NA
  


png("results/glom clusters pre-split.png", width = 30, height = 30, res = 600, units = "in")
par(mar = c(0,0,0,0))
#plot(customlocs[glomcells, ], col = brewer_pal(12, "Set3")(12)[1+(db %% 12)], pch = 16, cex = 0.2)
plot(customlocs, col = cols[clust], pch = 16, cex = 0.05, asp = 1)
for (cl in setdiff(unique(db), NA)) {
  inds = which(glomcells)[(db == cl) & !is.na(db)]
  hull = chull(customlocs[inds], )
  polygon(customlocs[inds[c(hull, hull[1])], ], lwd = 0.3)
  text(mean(customlocs[inds, 1]), mean(customlocs[inds, 2]), cl, cex = 0.2)
}
dev.off()
 

png("results/glom clusters black-grey.png", width = 30, height = 30, res = 600, units = "in")
par(mar = c(0,0,0,0))
#plot(customlocs[glomcells, ], col = brewer_pal(12, "Set3")(12)[1+(db %% 12)], pch = 16, cex = 0.2)
plot(customlocs, col = c("grey80", "grey30", "blue")[
  1 + is.element(rownames(customlocs), names(db)) + !is.na(db[rownames(customlocs)])], 
  pch = 16, cex = 0.1, asp = 1)
for (cl in setdiff(unique(db), NA)) {
  inds = which(glomcells)[(db == cl) & !is.na(db)]
  hull = chull(customlocs[inds], )
  polygon(customlocs[inds[c(hull, hull[1])], ], lwd = 0.3)
  text(mean(customlocs[inds, 1]), mean(customlocs[inds, 2]), cl, cex = 0.2)
}
dev.off()


# hand-flag clusters with multiple gloms, with number of gloms inside:
# (referring to the file "results/glom clusters pre-split.png)
tosplit = c(
  "137" = 2,
  "133" = 2,
  "233" = 2,
  "256" = 2,
  "240" = 2,
  "335" = 2,
  "349" = 2,
  "381" = 2,
  "581" = 2,
  "580" = 2,
  "494" = 2,
  "404" = 2,
  "442" = 2,
  "449" = 2,
  "528" = 2,
  "406" = 3,
  "407" = 2,
  "493" = 2,
  "427" = 2,
  "424" = 2,
  "469" = 2,
  "472" = 2,
  "31" = 2,
  "118" = 3,
  "79" = 2,
  "105" = 2,
  "102" = 3,
  "85" = 2,
  "93" = 2,
  "11" = 2,
  "191" = 2
)

  
# look like dups but just singles per old annots:
#"34" = 2,  "93" = 2,  "129" = 2,   "239" = 2,  "258" = 2,






db0 = db
# split cluster:
for (cl in names(tosplit)) {
  cl = as.numeric(cl)
  inds = which(glomcells)[(db == cl) & !is.na(db)]
  mc = Mclust(customlocs[inds, ], G = tosplit[as.character(cl)])
  db[(db == cl) & !is.na(db)] = paste0(cl, ".", mc$classification)
  
}


pdf("results/glomsplitting.pdf")
par(mfrow = c(2,2))
for (cl in names(tosplit)) {
  cl = as.numeric(cl)
  
  inds = which(glomcells)[(db0 == cl) & !is.na(db0)]
  plot(customlocs[inds, ], pch = 16, main = cl, asp = 1,
       col = as.numeric(as.factor(db[(substr(db,1,nchar(cl)) == cl) & !is.na(db)])))
}
dev.off()

# to merge:  (decisions made based on pathology review)
to_merge = list(c("13", "14"),
                c("31.1", "31.2"),
                c("39", "41"),
                c("58", "60"),
                c("82", "83"),
                c("93.1", "93.2"),
                c("142", "143"),
                c("145", "146"),
                c("151", "152"),
                c("210", "211"),
                c("218", "215"),
                c("240.1", "240.2"),
                c("191.1", "191.2")
)
for (i in 1:length(to_merge)) {
  tempids = to_merge[[i]]
  db[is.element(db, tempids)] = paste0(tempids, collapse = "_")
}


png("results/glom clusters post-split.png", width = 30, height = 30, res = 600, units = "in")
par(mar = c(0,0,0,0))
plot(customlocs, col = cols[clust], pch = 16, cex = 0.1, asp = 1)
for (cl in setdiff(unique(db), NA)) {
  inds = which(glomcells)[(db == cl) & !is.na(db)]
  hull = chull(customlocs[inds], )
  polygon(customlocs[inds[c(hull, hull[1])], ], lwd = 0.3)
  text(mean(customlocs[inds, 1]), mean(customlocs[inds, 2]), cl, cex = 0.2)
}
for (tiss in unique(annot$tissue)) {
  inds = which(annot$tissue == tiss)
  rect(min(customlocs[inds, 1]), min(customlocs[inds, 2]), max(customlocs[inds, 1]), max(customlocs[inds, 2]), border = "black")
  text(max(customlocs[inds, 1]), max(customlocs[inds, 2]), tiss, cex = 2)
}
dev.off()

#### annotate gloms: ------------------------

# center xy coords
# cell composition
gannot = data.frame(new.id = as.character(setdiff(unique(db), NA)))
rownames(gannot) = gannot$new.id
gannot$y = gannot$x = gannot$tissue = NA
for (cl in gannot$new.id) {
  cells = names(db)[(db == cl) & !is.na(db)]
  gannot[cl, "x"] = mean(customlocs[match(cells, annot$cell_ID), 1])
  gannot[cl, "y"] = mean(customlocs[match(cells, annot$cell_ID), 2])
  gannot[cl, "tissue"] = annot$tissue[match(cells[1], annot$cell_ID)]
}
gannot$rename_to = ""

### align to ids from the previous version:
#load("~/smiqumulo/01 SMI TAP project/SMI-0016_ShaunJackson_SeattleCH/6.99 custom analysis/processed_data/cannot.RData")
load("data/original_cannot.RData")
head(cannot)
old.glom.ids = setdiff(unique(cannot$in.glom), NA)
gannot$old.id = NA
for (i in 1:nrow(gannot)) {
  cellids = names(db)[(db == gannot$new.id[i]) & !is.na(db)]
  cannot.gloms = cannot$in.glom[match(cellids, rownames(cannot))]
  tab = table(cannot.gloms)
  tab = tab[order(tab, decreasing = T)]
  if (length(tab) > 0) {
    gannot$old.id[i] = names(tab)[1]
    if ((length(tab) > 1) & (tab[2] / tab[1] > 0.1)) {
      print(paste0("row ", i, " might be off"))
      print(tab)
    }
  } 
  if (length(tab) == 0) {
    gannot$old.id[i] = paste0("new", gannot$new.id[i])
  }
}

# look at duplicates:
dups =  names(which(table(gannot$old.id) > 1))
gannot[is.element(gannot$old.id, dups), ]

if (length(dups) > 0) {
  pdf("results/glom that might need merging.pdf")
  par(mfrow = c(2,2))
  for (oldid in dups) {
    newids = gannot$new.id[gannot$old.id == oldid]
    
    inds = which(glomcells)[is.element(db, newids) & !is.na(db)]
    plot(customlocs[inds, ], pch = 16, asp = 1, main = oldid) #,  col = as.numeric(as.factor(db[(substr(db,1,nchar(cl)) == cl) & !is.na(db)])))
  }
  dev.off()
}

# decision: let the pathologist's annotation overrule the clustering results: merge all gloms that map to a single old annotation


write.csv(gannot, file = "processed_data/glomeruli annotations.csv")
saveRDS(db, file = "processed_data/db.rds")

