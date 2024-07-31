
#### try combining ioprofiles with cellprofilelibrary's kidney scRNAseq:

## get ioprofiles:
library(InSituType)
library(pheatmap)
data("ioprofiles")

## get HCA kidney profiles:
load("data/Kidney_HCA.RData")
kid = as.matrix(profile_matrix)
colnames(kid)
immnames = c("CD4.T.cell" , "NK.cell", "B.cell", "CD8.T.cell", "Neutrophil", "NKT.cell", 
             "MNP.d.Tissue.macrophage", "Mast.cell", "Plasmacytoid.dendritic.cell",
             "MNP.c.dendritic.cell","MNP.a.classical.monocyte.derived","MNP.b.non.classical.monocyte.derived")

## merge them all
sharedgenes = intersect(rownames(ioprofiles), rownames(kid))
fixedprofiles = cbind(kid[sharedgenes, setdiff(colnames(kid), immnames)], 
                      ioprofiles[sharedgenes, setdiff(colnames(ioprofiles), c("fibroblast", "endothelial"))])
fixedprofiles = sweep(fixedprofiles, 2, colSums(fixedprofiles), "/") * 1e6
saveRDS(fixedprofiles, file = "processed_data/fixedprofiles.RDS")

if (FALSE) {
  # quick check:
  temp <- fixedprofiles[is.element(rownames(fixedprofiles), colnames(counts)), ]
  pheatmap(sweep(temp, 1, pmax(apply(temp, 1, max), 0.2), "/"), 
           col = colorRampPalette(c("white", "darkblue"))(100))
  
  # what genes separate DCT from PCT?
  plot(rowMeans(temp[, grepl("roximal", colnames(temp))]), temp[, "DCT"],
       col = 0, log = "xy")
  text(rowMeans(temp[, grepl("roximal", colnames(temp))]), temp[, "DCT"], rownames(temp), cex = 0.7)
  
  DCTgenes = c("MECOM", "EGF", "FGF13", "CASR", "KITLG", "TOX", "EFNA5", "MAML2", "CALB1")
  pheatmap(temp[DCTgenes, grepl("roximal", colnames(temp))])
  
  intersect(c("NCC", "CALB1", "SLC12A3", "CLDN14"), rownames(temp))
}
