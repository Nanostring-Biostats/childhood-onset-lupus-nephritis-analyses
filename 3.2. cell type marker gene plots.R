
rm(list = ls())
library(RColorBrewer)
library(ggplot2)
library(InSituType)
library(pheatmap)
library(viridis)

#### load data -----------------
load("processed_data/cleaneddata.RData")
cols = readRDS("processed_data/cols.RDS")
reslist = readRDS(file = "processed_data/cell typing by batch.RDS")

# subset on only high-confidence cell type calls:
prob = rep(NA, length(clust)); names(prob) = names(clust)
for (name in names(reslist)) {
  tempprobs = reslist[[name]]$prob
  ids = intersect(names(prob), names(tempprobs))
  prob[ids] = tempprobs[ids]
}

minprob = rep(0.95, length(unique(clust))); names(minprob) = unique(clust)
minprob["T CD4 naive"] = minprob["T CD8 naive"] = minprob["T CD4 memory"] = minprob["T CD8 memory"] = minprob["Treg"] = 0.75
minprob[grepl("ndothel", names(minprob))] = 0.9

propabove = nabove = c()
highprob = prob*NA
for (cell in names(minprob)) {
  propabove[cell] = mean(prob[clust == cell] > minprob[cell])
  nabove[cell] = sum(prob[clust == cell] > minprob[cell])
  highprob[clust == cell] = prob[clust == cell] > minprob[cell]
}
highprob = as.logical(highprob)
names(highprob) = names(prob)

#### define marker genes: -------------------------

# bigger clusters:
bigclust <- clust
#bigclust <- replace(clust, grepl("ndothel", clust), "endothelial")
bigclust <- replace(bigclust, grepl("CD4", bigclust), "CD4")
bigclust <- replace(bigclust, grepl("CD8", bigclust), "CD8")
bigclust <- replace(bigclust, grepl("Treg", bigclust), "CD4")
bigclust <- replace(bigclust, grepl("intercalated", bigclust), "intercalated")
bigclust <- replace(bigclust, grepl("ithelial", bigclust), "epithelial")
bigclust <- replace(bigclust, grepl("ithelium", bigclust), "epithelial")
unique(bigclust)

profiles <- InSituType:::Estep(counts = Matrix::t(raw)[highprob, ], clust[highprob], annot$negmean[highprob])$profiles
bigprofiles <- InSituType:::Estep(counts = Matrix::t(raw)[highprob, ], bigclust[highprob], annot$negmean[highprob])$profiles

profiles = profiles[, (colnames(profiles) != "Transitional.urothelium") ]
bigprofiles = bigprofiles[, (colnames(bigprofiles) != "Transitional.urothelium") ]

meanexpressing <- bigprofiles * NA
for (cell in colnames(meanexpressing)) {
  tempmat <- raw[rownames(meanexpressing), (bigclust == cell) & highprob]
  meanexpressing[, cell] <- Matrix::rowMeans(tempmat > 0)
}


# which genes act as markers?
markers <- (apply(bigprofiles, 1, max) > 0.1) & 
  (apply(bigprofiles, 1, max) > 2 * apply(bigprofiles, 1, function(x){
    x[order(x, decreasing = T)[2]]}))
# list of each cell type's markers:
celltypemarkers <- list()
for (cell in colnames(bigprofiles)) {
  celltypemarkers[[cell]] <- names(which(
    (bigprofiles[, cell] > 0.1) & 
      (bigprofiles[, cell] > 1.75 * apply(bigprofiles, 1, function(x){
        x[order(x, decreasing = T)[2]]}))
  ))
  
  celltypemarkers[[cell]] <- celltypemarkers[[cell]][
    order(bigprofiles[celltypemarkers[[cell]], cell] / 
            apply(bigprofiles[celltypemarkers[[cell]], setdiff(colnames(bigprofiles), cell), drop = F], 1, max),
          decreasing = T)[1:pmin(5, length(celltypemarkers[[cell]]))]]
}
celltypemarkers
celltypemarkers$endothelial = c("CDH5","FLT1","TIE1","TGFBR2" , "CLEC14A")

pheatmap(sweep(bigprofiles[setdiff(unlist(celltypemarkers), NA), ], 1, 
               apply(bigprofiles[setdiff(unlist(celltypemarkers), NA), ],1,max), "/"))

pheatmap(pmin(meanexpressing[setdiff(unlist(celltypemarkers), NA), ], 0.5))

# for each gene x cell type, get % expressing and mean bg-sub expression
df = data.frame(gene = rep(rownames(bigprofiles), ncol(bigprofiles)),
                clust = rep(colnames(bigprofiles), each = nrow(bigprofiles)),
                mean = as.vector(bigprofiles),
                prop.expressing = as.vector(meanexpressing))
df <- df[is.element(df$gene, unlist(celltypemarkers)), ]
df$scaled.mean = df$mean
for (gene in unique(df$gene)) {
  df$scaled.mean[df$gene == gene] = df$scaled.mean[df$gene == gene] / max(df$scaled.mean[df$gene == gene])
}
# order genes:
genetree = hclust(dist(sweep(bigprofiles[setdiff(unlist(celltypemarkers), NA), ], 1, 
                             apply(bigprofiles[setdiff(unlist(celltypemarkers), NA), ],1,max), "/")))
geneorder = genetree$labels[genetree$order]
df$gene = factor(df$gene, levels = geneorder)
cellorder = c("CD4", "CD8", 
              "Thick.ascending.limb.of.Loop.of.Henle",
              "Connecting.tubule", "Principal.cell", 
              "mDC", "macrophage", "Fibroblast", "Myofibroblast", 
              "Mesangial.cell", "Vascular.pericyte",
              "neutrophil", "intercalated", "Podocyte",
              "monocyte", "PCT", "mast", "NK",
              "pDC", "plasmablast", "B-cell",
              "Peritubular.capillary.endothelium.1",
              "Peritubular.capillary.endothelium.2",
              "Ascending.vasa.recta.endothelium",
              "Descending.vasa.recta.endothelium",
              "Glomerular.endothelium", "epithelial")
df$clust = factor(df$clust, levels = cellorder)
pdf("results/markers heatmap - supplemental figure.pdf", height = 15, width = 12)
ggplot(df, aes(x=clust, y = gene, color = scaled.mean, size = prop.expressing)) + 
  geom_point() + scale_color_viridis(option = "A", direction = -1) + ggthemes::theme_few() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
dev.off()
