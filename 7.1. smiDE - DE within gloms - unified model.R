# Goal: for each glomerular cell type, run DE vs. various annotations

rm(list = ls())
library(lmerTest)
library(glmmTMB)
library(smiDE)

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



#### define and explore predictors: ----------------------------------

# define a disease variable: "control", "normal-like", "involved
gloms$path = NA
# define what makes a glom normal-like:
gloms$nopathology = ((gloms$involved == 0) & !is.na(gloms$involved)) *
  ((gloms$segmental.sclerosis == 0) & !is.na(gloms$segmental.sclerosis)) *
  ((gloms$endocapillary10 == 0) & !is.na(gloms$endocapillary10)) *
  ((gloms$wire10 == 0) & !is.na(gloms$wire10)) *
  ((gloms$wire10 == 0) & !is.na(gloms$wire10)) *
  ((gloms$`circulating inflamm cells` == 0) & !is.na(gloms$`circulating inflamm cells`)) *
  is.na(gloms$other) *
  is.na(gloms$...10) # (the "comments" field - all gloms with pathologist comments aren't normal)
gloms$path[gloms$nopathology == 1] = "non-pathological"
gloms$path[gloms$nopathology == 0] = "pathological"
gloms$path[grepl("Control", gloms$tissue)] = "control"

# total lymphoid:
gloms$lymphoid.tot = gloms$B.tot + gloms$T.tot + gloms$plasma.tot

gloms$lymphoid.transformed = log2(1 + gloms$lymphoid.tot)
gloms$myeloid.transformed = log2(1 + gloms$myeloid.tot)

#### run DE over all variables and cell types: -----------------------

# build a data frame of glomerular annots at the cell level:
deannot = cannot
deannot$tissue = annot$tissuename
deannot$totalcounts = annot$totalcounts
deannot$cell_ID = annot$cell_ID
glomvars = c("path", "lymphoid.transformed", "myeloid.transformed")

for (varname in glomvars) {
  deannot[[varname]] = gloms[[varname]][match(cannot$closest.glom, rownames(gloms))]
}


#### run main DE model (all data) --------------------------

if (!file.exists("processed_data/DE results - glomcells.RDS")) {
  celltypes = c("Podocyte", "Mesangial.cell", "Glomerular.endothelium") 
  deres <- list()
  for (name in celltypes) {
    print(name)
    tempgenes <- names(which(!toomuchcontam[, name]))
    tempcells <- annot$cell_ID[(clust == name) & !is.na(cannot$in.glom)]
    # call smi_de:
    deres[[name]] <- smiDE::smi_de(
      assay_matrix = raw[tempgenes, tempcells],
      metadata = deannot[tempcells, ],
      formula = ~offset(log(totalcounts)) + RankNorm(otherct_expr) + (1|in.glom) + #(1|tissue) +
        path + lymphoid.transformed + myeloid.transformed,
      neighborhood_counts = neighborhood_counts,
      groupVar = "path",
      family="nbinom2"
    )
    saveRDS(deres, file = "processed_data/DE results - glomcells.RDS")
  }
} else {
  deres <- readRDS("processed_data/DE results - glomcells.RDS")
}

results(deres[[1]])


#### summarize results -------------------------

# get data for volcano plots 
reslist <- deres

# check out volcplots for SLE vs. control:
par(mfrow = c(1,3))
for (cell in names(reslist)) {
  res = results(reslist[[cell]])
  temp <- res$one.vs.rest[res$one.vs.rest$contrast == "control vs. avg.rest", ]
  plot(log2(temp$fold_change), -log10(temp$p.value), col = 0, main = cell)
  text(log2(temp$fold_change), -log10(temp$p.value), temp$target, cex = 0.7)
}
#par(mfrow = c(3,2))
#for (cell in names(deres.newcontrols)) {
#  resnew = results(deres.newcontrols[[cell]])
#  resold = results(deres.oldcontrol[[cell]])
#  tempnew <- resnew$one.vs.rest[resnew$one.vs.rest$contrast == "control vs. avg.rest", ]
#  tempold <- resold$one.vs.rest[resold$one.vs.rest$contrast == "control vs. avg.rest", ]
#  plot(log2(tempnew$fold_change), -log10(tempnew$p.value))
#  plot(log2(tempold$fold_change), -log10(tempold$p.value))
#}


#### run DE model using just subsets of controls --------------------------

oldcontrols = c("Control2")
newcontrols = paste0("Control", c(1,3,4))

# old control only:
if (!file.exists("processed_data/DE results - glomcells - oldcontrolonly.RDS")) {
  deres.oldcontrol <- list()
  for (name in celltypes) {
    print(name)
    tempgenes <- names(which(!toomuchcontam[, name]))
    tempcells <- annot$cell_ID[((clust == name) & !is.na(cannot$in.glom)) &
                                 (!is.element(annot$tissuename, newcontrols))]
    # call smi_de:
    deres.oldcontrol[[name]] <- smiDE::smi_de(
      assay_matrix = raw[tempgenes,tempcells],
      metadata = deannot[tempcells, ],
      formula = ~offset(log(totalcounts)) + RankNorm(otherct_expr) + (1|in.glom) + #(1|tissue) +
        path + lymphoid.transformed + myeloid.transformed,
      neighborhood_counts = neighborhood_counts,
      groupVar = "path",
      family="nbinom2"
    )
    saveRDS(deres.oldcontrol, file = paste0("processed_data/DE results - glomcells - oldcontrolonly.RDS"))
  }
} else {
  deres.oldcontrol <- readRDS("processed_data/DE results - glomcells - oldcontrolonly.RDS")
}

# new controls only:
if (!file.exists("processed_data/DE results - glomcells - newontrolsonly.RDS")) {
  deres.newcontrols <- list()
  for (name in celltypes) {
    print(name)
    tempgenes <- names(which(!toomuchcontam[, name]))
    tempcells <- annot$cell_ID[((clust == name) & !is.na(cannot$in.glom)) &
                                 (!is.element(annot$tissuename, oldcontrols))]
    # call smi_de:
    deres.newcontrols[[name]] <- smiDE::smi_de(
      assay_matrix = raw[tempgenes, tempcells],
      metadata = deannot[tempcells, ],
      formula = ~offset(log(totalcounts)) + RankNorm(otherct_expr) + (1|in.glom) + #(1|tissue) +
        path + lymphoid.transformed + myeloid.transformed,
      neighborhood_counts = neighborhood_counts,
      groupVar = "path",
      family="nbinom2"
    )
    saveRDS(deres.newcontrols, file = paste0("processed_data/DE results - glomcells - newontrolsonly.RDS"))
  }
} else {
  deres.newcontrols <- readRDS("processed_data/DE results - glomcells - newontrolsonly.RDS")
}


#### As a sanity check against batch effects, compare the analysis one of batch of controls at a time: ----------------------
par(mfrow = c(3,3))
for (cell in names(reslist)) {
  res = results(reslist[[cell]])
  resnew = results(deres.newcontrols[[cell]])
  resold = results(deres.oldcontrol[[cell]])
  temp <- res$one.vs.rest[res$one.vs.rest$contrast == "control vs. avg.rest", ]
  tempnew <- resnew$one.vs.rest[resnew$one.vs.rest$contrast == "control vs. avg.rest", ]
  tempold <- resold$one.vs.rest[resold$one.vs.rest$contrast == "control vs. avg.rest", ]
  
  plot(-log10(temp$p.value) * sign(log2(temp$fold_change)),
       -log10(tempnew$p.value) * sign(log2(tempnew$fold_change)), cex = 0.3)
  points(-log10(temp$p.value) * sign(log2(temp$fold_change)),
         -log10(tempold$p.value) * sign(log2(tempold$fold_change)), cex = 0.3, col = 2)
  
  abline(0,1)
  
  topgenes <- temp$target[p.adjust(temp$p.value, "BH") < 0.05]
  topgenes = topgenes[order(temp$fold_change[match(topgenes, temp$target)])]
  plot(temp[match(topgenes, target), log2(fold_change)], 1:length(topgenes), pch = 16, main = cell)
  points(tempnew[match(topgenes, target), log2(fold_change)], 1:length(topgenes), pch = 16, col = "red")
  points(tempold[match(topgenes, target), log2(fold_change)], 1:length(topgenes), pch = 16, col = "blue")
  abline(v = 0)
  
  plot(tempnew[, log2(fold_change)], tempold[, log2(fold_change)])
  abline(0,1)
  abline(h = 0)
  abline(v = 0)
}

par(mfrow = c(1,3))
for (cell in names(deres)) {
  res = results(deres[[cell]])
  blacklist = c("AZU1", "COL9A2", "DUSP5", "MZT2A", "SPP1", "AZGP1", "EGF")
  
  ## get SLE vs. control results:
  temp <- res$one.vs.rest[res$one.vs.rest$contrast == "control vs. avg.rest", ]
  # remove blacklist genes
  temp <- temp[!is.element(temp$target, blacklist), ]
  
  temp$fdr <- p.adjust(temp$p.value, method = "BH")
  temp$log2foldchange = log2(temp$fold_change)
  
  # get results from subset models:
  resnew = results(deres.newcontrols[[cell]])
  resold = results(deres.oldcontrol[[cell]])
  tempnew <- resnew$one.vs.rest[resnew$one.vs.rest$contrast == "control vs. avg.rest", ]
  tempold <- resold$one.vs.rest[resold$one.vs.rest$contrast == "control vs. avg.rest", ]
  tempnew <- tempnew[!is.element(tempnew$target, blacklist), ]
  tempold <- tempold[!is.element(tempold$target, blacklist), ]
  
  # do the subsets confirm the full data's result?
  temp$confirmed.samedirection = ((tempnew$fold_change > 1) == (tempold$fold_change > 1)) 
  temp$confirmed.samedirectionANDalwayssignificant = (((tempnew$fold_change > 1) == (tempold$fold_change > 1))) & 
    ((tempold$p.value < 0.05) & (tempnew$p.value < 0.05))
  
  
  write.csv(temp[, c("target", "log2foldchange", "p.value", "fdr", "contrast", "confirmed.samedirection", "confirmed.samedirectionANDalwayssignificant" )],
            file = paste0("results/DE glom cells - control vs. SLE - ", cell, ".csv"))
  
  # quick volc plot:
  plot(log2(temp$fold_change), -log10(temp$p.value), col = 0)
  text(log2(temp$fold_change), -log10(temp$p.value), temp$target, cex = 0.6)
  
  ## pairwise contrasts:
  for (contrastval in c("control / non-pathological", "control / pathological", "non-pathological / pathological")) {
    temp <- res$pairwise
    temp <- temp[contrast == eval(contrastval), ]
    temp$fdr <- p.adjust(temp$p.value, method = "BH")
    temp$log2foldchange = log2(temp$fold_change)
    write.csv(temp[, c("target", "log2foldchange", "p.value", "fdr", "contrast")],
              file = paste0("results/DE glom cells - ", make.names(contrastval), " - ", cell, ".csv"), row.names = F)
  }
  
}

