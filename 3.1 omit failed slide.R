
rm(list = ls())

# load necessary data:
load("processed_data/data_with_failed_slide.RData")
viz <- readRDS("processed_data/viz.RDS")
clust <- readRDS("processed_data/clusts.RDS")

ls()

badslide = "SMI0016c_SP212319"
isgood = annot$slidename != badslide

raw = raw[, isgood]
annot = annot[isgood, ]
viz = viz[isgood, ]
customlocs = customlocs[isgood, ]
clust = clust[isgood]
um = um[isgood, ]

# adjust viz:
inds = (viz[,2] > 15) & (viz[, 1] > 5)
plot(viz, col = 1 + inds)

viz[inds, 2] <- viz[inds, 2] - 3.5

save(raw, annot, viz, customlocs, clust, um, file = "processed_data/cleaneddata.RData")
