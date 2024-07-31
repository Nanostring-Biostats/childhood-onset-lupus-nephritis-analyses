# get each gloms total expression from the major cell types, and explore

library(ggplot2)
#library(ppcor)
library(harmony)
rm(list = ls())

#### load data: ------------------------

gloms = readRDS("processed_data/gloms with cell annotations.RDS")
cols <- readRDS("processed_data/cols.RDS")
load("processed_data/cleaneddata.RData")
load("processed_data/cannot with distances from gloms.RData") # loads "cannot"
load("processed_data/polygons.RData")
toomuchcontam <- readRDS("processed_data/toomuchcontam.RDS")
tissuecols = readRDS(file = "processed_data/tissuecols.RDS")
blacklist = c("AZU1", "COL9A2", "DUSP5", "MZT2A", "SPP1", "AZGP1", "EGF")

#### run Harmony for batch correction: ------------------------

if (!file.exists("processed_data/harmony embeddings.RDS")) {
  pcs <- irlba::prcomp_irlba(norm[, !is.element(colnames(norm), blacklist)], n = 20)$x
  print("havepcs")
  harmony_embeddings <- RunHarmony(data_mat = pcs, meta_data = annot$slidename) #, "dataset")
  saveRDS(harmony_embeddings, file = "processed_data/harmony embeddings.RDS")
} else {
  harmony_embeddings <- readRDS("processed_data/harmony embeddings.RDS")
}

x <- harmony_embeddings
rm(harmony_embeddings)

#### get each glom's total profile from each cell type:

if (!file.exists("processed_data/glomwidedata.RData")) {
  glomneg = glomraw = glomx = list()
  
  celllists = list("Podocyte" = "Podocyte", 
                   "Mesangial.cell" = "Mesangial.cell", 
                   "Endothelial" = unique(clust)[grepl("ndothel", unique(clust))])
  celllists
  # get total counts and negmean per glom
  for (name in names(celllists)) {
    glomraw[[name]] = matrix(NA, nrow(gloms), nrow(raw),
                             dimnames = list(rownames(gloms), rownames(raw)))
    glomneg[[name]] = rep(NA, nrow(gloms))
    names(glomneg[[name]]) = rownames(gloms)
    
    for (i in 1:nrow(gloms)) {
      cell.inds = which(((cannot$inside.glom == gloms$glomnames[i]) & !is.na(cannot$inside.glom)) & is.element(clust, celllists[[name]]))
      if (length(cell.inds) >= 5) { # requiring at least 5 cells to report data
        glomraw[[name]][i, ] = rowMeans(raw[, cell.inds, drop = F])
        glomneg[[name]][i] = mean(annot$negmean[cell.inds])
      }
    }
  }
  
  
  # get normalized bg-sub data per glom:
  glomnorm = list()
  for (name in names(glomraw)) {
    # bg-subtract:
    temp = pmax(sweep(glomraw[[name]], 1, glomneg[[name]], "-"), 0)
    # normalize (for gloms where the cell type is present):
    #hasdata = !is.na(temp[,1])
    #temp9 = rowSums(temp[hasdata, ])
    #temp[hasdata,] = sweep(temp[hasdata, ], 1, tempq9, "/") * mean(tempq9)
    tempq9 = rowSums(temp)
    temp = sweep(temp, 1, tempq9, "/") * mean(tempq9, na.rm = T)
    glomnorm[[name]] = temp
  }
  
  # get harmony embeddings per glom:
  for (name in names(celllists)) {
    glomx[[name]] = matrix(NA, nrow(gloms), ncol(x),
                           dimnames = list(rownames(gloms), colnames(x)))
    for (i in 1:nrow(gloms)) {
      cell.inds = which(((cannot$inside.glom == gloms$glomnames[i]) & !is.na(cannot$inside.glom)) & is.element(clust, celllists[[name]]))
      if (length(cell.inds) >= 5) { # requiring at least 5 cells to report data
        glomx[[name]][i, ] = colMeans(x[cell.inds, , drop = F])
      }
    }
  }
  save(glomraw, glomnorm, glomx, file = "processed_data/glomwidedata.RData")
} else {
  load("processed_data/glomwidedata.RData")
}

## list of genes to use: -----------------------------
uselist = list("Podocyte" = rownames(toomuchcontam)[!toomuchcontam[, "Podocyte"]],
               "Mesangial.cell" = rownames(toomuchcontam)[!toomuchcontam[, "Mesangial.cell"]])
uselist[["Endothelial"]] = c()
for (name in  unique(clust)[grepl("ndothel", unique(clust))]) {
  uselist[["Endothelial"]] = unique(c(uselist[["Endothelial"]], rownames(toomuchcontam)[!toomuchcontam[, name]]))
}
# also remove too-low genes:
for (name in names(glomraw)) {
  toolow = colMeans(glomraw[[name]], na.rm = T) < 0.1
  print(paste0(name, " ", sum(toolow)))
  uselist[[name]] = setdiff(uselist[[name]], c(names(which(toolow)), blacklist))
}


#### make wide glom-level data: ---------------------------

## get triple-wide data for each glom:
# normalized data:
wglom = c()
for (name in names(glomnorm)) {
  temp = glomnorm[[name]][, uselist[[name]]]
  colnames(temp) = paste0(colnames(temp), "_", name)
  wglom = cbind(wglom, temp)
}
# harmony embeddings:
hglom = c()
for (name in names(glomx)) {
  temp = glomx[[name]]
  colnames(temp) = paste0(colnames(temp), "_", name)
  hglom = cbind(hglom, temp)
}


#### cluster only SLE gloms ------------------------------------

slegloms = rownames(gloms)[grepl("SLE", gloms$tissue)]
temp <- wglom[is.element(rownames(wglom), slegloms), ]
temp <- temp[rowSums(is.na(temp)) == 0, ]

sleum <- uwot::umap(temp, n_neighbors = 50)

## umap plots:
svg("results/glomeruli umaps - sle only.svg", widt = 9, height = 3.5)
par(mfrow = c(1,3))
par(mar = c(0,0.1,2,0.1))
# by tissue:
plot(sleum, pch = 16, main = "By tissue", 
     xaxt = "n", yaxt = "n", xlab = "",ylab = "",
     col = tissuecols[gloms[rownames(sleum), "tissue"]])
legend("topleft", pch = 16, col = tissuecols[paste0("SLE", c(1:7, 8.1, 8.2, 8.3))],
       legend = names(tissuecols[paste0("SLE", c(1:7, 8.1, 8.2, 8.3))]))
# by slide:
annot$slideid = substr(annot$slidename, nchar(annot$slidename) - 5, nchar(annot$slidename))
glomslide = as.factor(annot$slideid[match(gloms$tissue, annot$tissuename)])
plot(sleum, pch = 16, main = "By slide", 
     xaxt = "n", yaxt = "n", xlab = "",ylab = "",
     col = as.numeric(glomslide))
legend("topleft", pch = 16, col = 1:length(unique(glomslide)),
       legend = levels(glomslide))
# by involved:
umapinvolved <- gloms[rownames(sleum), "involved"]
plot(sleum, pch = 16, main = "By involvement", 
     xaxt = "n", yaxt = "n", xlab = "",ylab = "",
     col = c("grey50", "red")[1 + (umapinvolved == "1")])
legend("topleft", pch = 16, col = c("grey50", "red"),
       legend = c("Not involved", "Involved"))
dev.off()

#### cluster gloms based on their wide data -------------------
temp = sweep(wglom, 2, pmax(colMeans(wglom, na.rm = T), 0.3), "/")
temp = sqrt(temp)
keep = rowSums(is.na(temp)) == 0
wtemp = temp[keep, ]
htemp = hglom[keep, ]

# umaps:
wum = uwot::umap(wtemp, n_neighbors = 50)
hum = uwot::umap(htemp)

pdf("results/glom umaps.pdf", width = 8, height = 8)
par(mfrow = c(2,2))
# color by tissue:
plot(wum, pch = 16, col = tissuecols[gloms[rownames(wum), "tissue"]],
     main = "normalized counts")
plot(hum, pch = 16, col = tissuecols[gloms[rownames(hum), "tissue"]], 
     main = "batch-corrected PCs")
# color by involvement:
gloms$involved2 = gloms$involved; gloms$involved2[grepl("Control", gloms$tissue)] = "control"
umapinvolved <- gloms[rownames(wum), "involved2"]
plot(wum, pch = 16, 
     col = c("grey80", "cornflowerblue", "orange")[1 + (umapinvolved == "0") + 2 * (umapinvolved == "1")], 
  main = "normalized counts") 
legend("topleft", pch = 16, col = c("grey80", "cornflowerblue", "orange"),
       legend = c("control", "SLE not involved", "SLE involved"))
plot(hum, pch = 16, 
     col = c("grey80", "cornflowerblue", "orange")[1 + (umapinvolved == "0") + 2 * (umapinvolved == "1")], 
     main = "batch-corrected PCs")
dev.off()


## umaps of gloms by various variables
# tissue ID:
svg(paste0("results/glom umap by tissue - ", subname, ".svg"))
plot(ums[[subname]], pch = 16, col = mycols$tissue[gloms$tissue], asp = 1, main = subname)
legend("topright", pch = 16, col = mycols$tissue, legend =  names(mycols$tissue), cex = 0.6)
dev.off()

#  involved (gloms$involved = 1 * (((gloms$`wire-loop` > 0) | (gloms$double.contours > 0)) | (gloms$endocapillary > 0)) )
svg(paste0("results/glom umap by involved - ", subname, ".svg"))
gloms$involved2 = gloms$involved; gloms$involved2[grepl("control", gloms$tissue)] = "control"
plot(ums[[subname]], pch = 16, col = c("grey40", "orange", "cornflowerblue")[1+(gloms$involved2=="1")+2*(gloms$involved2=="control")], main = "involved")
legend("topright", pch = 16, col = c("grey40", "orange", "cornflowerblue"), legend =  c("not involved", "involved", "control"), cex = 0.6)
dev.off()

# disease state
gloms$class = annot$class[match(gloms$tissue, annot$tissue)]
classcols = c("control" = "chartreuse3", "III" = "cornflowerblue", "IV" = "red", "IV+V" = "grey10")
svg(paste0("results/glom umap by disease stats - ", subname, ".svg"))
plot(ums[[subname]], pch = 16, col = classcols[gloms$class], main = "segmental.sclerosis")
legend("topright", pch = 16, col = classcols, legend =  names(classcols), cex = 0.6)
dev.off()

# histopathologically normal (all 0s per Robyn):
gloms$histopath0 = c("abnormal", "normal")[1 + ((1-gloms$mesang10) * (1-gloms$endocapillary10) * (1-gloms$wire10) * (1-gloms$`circulating inflamm cells`) * 
  (1-gloms$double.contours) * (1-gloms$protein.droplets) * (1-gloms$podocyte.hyperplasia) * (1-gloms$segmental.sclerosis))]
histopathcols = c("normal" = "black", "abnormal" = "orange")
svg(paste0("results/glom umap by histopath normal vs. not - ", subname, ".svg"))
plot(ums[[subname]], pch = 16, col = histopathcols[(gloms$histopath0)], main = "segmental.sclerosis")
legend("topright", pch = 16, col = histopathcols, legend =  names(histopathcols), cex = 0.6)
dev.off()






