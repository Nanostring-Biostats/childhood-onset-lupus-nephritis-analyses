rm(list = ls())
library(ComplexHeatmap)
library(umap)
library(pheatmap)
library(spatstat)
library(spatstat.geom)
library(data.table)
library(readxl)
library(InSituCor)
library(ggplot2)
library(ggthemes)

#### load data -------------------- 
rm(list = ls())
cols <- readRDS("processed_data/cols.RDS")
load("processed_data/cleaneddata.RData")
rownames(customlocs) = annot$cell_ID
load( "processed_data/polygons.RData")
tissuecols = readRDS(file = "processed_data/tissuecols.RDS")
gloms <- readRDS("processed_data/gloms.rds")
norm <- Matrix::t(sweep(raw, 1, Matrix::colMeans(raw), "/") * mean(Matrix::colMeans(raw)))

#### run insitucor -------------------------------

insitucorres <- "processed_data/insitucor results.RDS"
if (!file.exists(insitucorres)) {
  set.seed(0)
  res <- insitucor(counts = norm, 
                   conditionon = cbind(clust, annot$tissuename, annot$negmean, annot$totalcounts), 
                   celltype = clust, xy = customlocs, 
                   k = 50, roundcortozero = 1e-2, max_cells = 20000, max_module_size = 20)
  saveRDS(res, file = insitucorres)
} else {
  res <- readRDS(insitucorres)
}


#### explore results ---------------------------


#### correlation networks ---------------------
# plot correlation network:
pdf("results/insitucor/correlation networks.pdf")
set.seed(0)
InSituCor::plotCorrelationNetwork(res$condcor * (abs(res$condcor) > 0.1), 
                                  modules = res$modules, show_gene_names = F)
dev.off()
pdf("results/insitucor/correlation networks with gene names.pdf")
set.seed(0)
InSituCor::plotCorrelationNetwork(res$condcor * (abs(res$condcor) > 0.1), 
                                  modules = res$modules, show_gene_names = F)
dev.off()

pdf("results/insitucor/IFIT correlation network.pdf", width = 6, height = 6)
par(mar = c(0,0,0,0))
tempgenes = res$modules$gene[res$modules$module == "IFITM3_IFI27_25"]
x = (res$condcor * (abs(res$condcor) > 0.1))
xum <- uwot::umap(as.matrix(x[tempgenes, tempgenes]), spread = 25, 
                  min_dist = 0.1, n_neighbors = max(min(length(tempgenes) - 2, 15), 1))
gr0 <- InSituCor:::get_coregulation_network(cormat = x[tempgenes, tempgenes], 
                                corthresh = 0.1)
set.seed(0)
igraph::plot.igraph(gr0, layout = xum, 
                    vertex.label = unlist(list(NA, names(igraph::V(gr0)))[1 + TRUE]), vertex.size = 0, 
                    vertex.label.color = "black")
dev.off()

pdf("results/insitucor/COL correlation network.pdf", width = 5.5, height = 5)
par(mar = c(0,1,0,1))
tempgenes = res$modules$gene[res$modules$module == "COL3A1_COL1A1_10"]
x = (res$condcor * (abs(res$condcor) > 0.1))
xum <- uwot::umap(as.matrix(x[tempgenes, tempgenes]), spread = 25, 
                  min_dist = 0.1, n_neighbors = max(min(length(tempgenes) - 2, 15), 1))
gr0 <- InSituCor:::get_coregulation_network(cormat = x[tempgenes, tempgenes], 
                                            corthresh = 0.1)
set.seed(0)
igraph::plot.igraph(gr0, layout = xum, 
                    vertex.label = unlist(list(NA, names(igraph::V(gr0)))[1 + TRUE]), vertex.size = 0, 
                    vertex.label.color = "black")
dev.off()


pdf("results/insitucor/macrophage correlation network.pdf", width = 6.5, height = 6)
par(mar = c(0,1,0,1))
tempgenes = res$modules$gene[is.element(res$modules$module, c("C1QA_CD163_11", "CD74_HLA.DRA_9"))]
tempgenes2 = res$modules$gene[is.element(res$modules$module, c("CD74_HLA.DRA_9"))]

x = (res$condcor * (abs(res$condcor) > 0.1))
xum <- uwot::umap(as.matrix(x[tempgenes, tempgenes]), spread = 25, 
                  min_dist = 0.1, n_neighbors = max(min(length(tempgenes) - 2, 15), 1))
gr0 <- InSituCor:::get_coregulation_network(cormat = x[tempgenes, tempgenes], 
                                            corthresh = 0.1)
set.seed(0)
igraph::plot.igraph(gr0, layout = xum, 
                    vertex.label = unlist(list(NA, names(igraph::V(gr0)))[1 + TRUE]), vertex.size = 0, 
                    vertex.label.color = c("red", "blue")[1 + is.element(tempgenes, tempgenes2)])
dev.off()


#### attribution results ---------------------


pdf("results/insitucor/insitucor attribution matrices.pdf")
for (name in names(res$attributionmats)) {
  print(pheatmap(res$attributionmats[[name]], main = name, fontsize_row = 6,
           col = colorRampPalette(c("white", "darkblue"))(100),
           breaks = seq(0,1,length.out=101)))
}
dev.off()

# fibrosis module heatmap:
# tidy cell types:
# - combine T-cell subtypes into CD4, CD8
# - no urothelium
# - combine peritubularcapillary endo
# - 

name <- "COL3A1_COL1A1_10"
mat <- res$attributionmats[[name]]
newmat <- rbind(
  colMeans(mat[c("T CD8 memory", "T CD8 naive"), ]),
  colMeans(mat[c("T CD4 memory", "T CD4 naive", "Treg"), ]),
  colMeans(mat[c("Peritubular.capillary.endothelium.1", "Peritubular.capillary.endothelium.2"), ])
)
rownames(newmat) <- c("CD8 T cell", "CD4 T cell", "Peritubular capillary endothelium")
mat <- rbind(newmat, mat)
mat <- mat[!is.element(rownames(mat), c("Transitional.urothelium", "T CD8 memory", "T CD8 naive",
                                        "T CD4 memory", "T CD4 naive", "Treg",
                                        "Peritubular.capillary.endothelium.1", "Peritubular.capillary.endothelium.2")), ]

write.csv(mat, file = "results/insitucor/fibrosis module attribution scores.csv")

#### spatial plots ---------------------

#for (name in colnames(res$scores_env)) {
for (name in c("COL3A1_COL1A1_10", "CD74_HLA.DRA_9", "C1QA_CD163_11", "IFITM3_IFI27_25")) {
  png(paste0("results/insitucor/", name, ".png"), width = 14, height = 14, res = 400, units = "in")
  par(mar = c(0,0,0,0))
  plot(viz, pch = 16, cex = 0.2,
       col = viridis_pal(option = "B")(101)[
         pmin(101, 1 + round(100 * res$scores_env[, name] / quantile(res$scores_env[, name], 0.995)))
       ],
       asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  legend("topright", legend = name)
  for (tiss in unique(annot$tissuename)) {
    text(median(range(viz[annot$tissuename == tiss, 1])),
         max(viz[annot$tissuename == tiss, 2] + 0.2),
         tiss, col = "darkred")
  }
  for (g in names(gpoly.viz)) {
    polygon(gpoly.viz[[g]], border = "cyan")
  }
  dev.off()
  
}

#### custom spatial plot of COL network ---------------------

name = "COL3A1_COL1A1_10"
png(paste0("results/insitucor/", name, " - with fibroblasts.png"), width = 14, height = 14, res = 400, units = "in")
par(mar = c(0,0,0,0))
plot(viz, pch = 16, cex = 0.2,
     col = viridis_pal(option = "B")(101)[
       pmin(101, 1 + round(100 * res$scores_env[, name] / quantile(res$scores_env[, name], 0.995)))
     ],
     asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
legend("topright", legend = name)
for (tiss in unique(annot$tissuename)) {
  text(median(range(viz[annot$tissuename == tiss, 1])),
       max(viz[annot$tissuename == tiss, 2] + 0.2),
       tiss, col = "darkred")
}
for (g in names(gpoly.viz)) {
  polygon(gpoly.viz[[g]], border = "cyan")
}
points(viz[clust == "Fibroblast", ], cex = 0.2, pch = 16, col = "blue")
dev.off()

####  macropahge scores scatterplot ------------------------------

png("results/insitucor/macrophage scores.png", width = 5, height = 5, units = "in", res = 500)
o <- sample(1:nrow(annot))
plot(res$scores_env[o, c("C1QA_CD163_11", "CD74_HLA.DRA_9")],
     pch = 16, cex = 0.2, col = tissuecols[annot$tissuename][o],
     xlab = "Compliments score", ylab = "Antigen presentation score")
dev.off()
pdf("results/tissue color legend.pdf", width = 1.4, height = 3.5)
par(mar = c(0,0,0,0))
frame()
legend('center', col = tissuecols, pch = 16, legend = names(tissuecols))
dev.off()


#### subset spatial plots: max and median FOVs per tissue:

for (name in c("COL3A1_COL1A1_10", "IFITM3_IFI27_25", "C1QA_CD163_11", "CD74_HLA.DRA_9")) {
  # get median and max scores per tissue:
  keepcells <- c()
  fovpairs <- c()
  for (tiss in unique(annot$tissuename)) {
    inds <- annot$tissuename == tiss
    means <- by(res$scores_env[inds, name], annot$fov[inds], mean)
    keepfovs <- names(means)[order(means, decreasing = T)[c(1, floor(length(means) / 2))]]
    keepcells <- c(keepcells, which(inds & (is.element(annot$fov, keepfovs))))
    fovpairs <- rbind(fovpairs, c(tiss, keepfovs[1], "max"))
    fovpairs <- rbind(fovpairs, c(tiss, keepfovs[2], "mean"))
  }
  
  # make df for ggplot:
  df <- data.frame(tissue = annot$tissuename[keepcells],
                   score = res$scores_env[keepcells, name],
                   score99 = pmin(res$scores_env[keepcells, name], quantile(res$scores_env[keepcells, name], 0.99)),
                   fov = annot$fov[keepcells],
                   x = customlocs[keepcells, 1],
                   y = customlocs[keepcells, 2])
  df$fovtype <- fovpairs[,3][match(df$fov, fovpairs[, 2])]
  
  # df of polygon data:
  polydf <- data.frame(
    x = unlist(sapply(gpoly.customlocs, function(x){x[,1]})),
    y = unlist(sapply(gpoly.customlocs, function(x){x[,2]})),
    name = unlist(sapply(1:length(gpoly.customlocs), 
                         function(i){rep(names(gpoly.customlocs)[i], 
                                         nrow(gpoly.customlocs[[i]]))})))
  # associate with tissues:
  polydf$tissue <- annot$tissuename[match(
    substr(rownames(polydf), unlist(gregexpr("c", rownames(polydf))), nchar(rownames(polydf))),
    annot$cell_ID
  )]
  polydf$fov <- annot$fov[match(
    substr(rownames(polydf), 
           sapply(rownames(polydf), function(str){gregexpr("c", str)[[1]][1]}),
           #unlist(gregexpr("c", rownames(polydf))), 
           nchar(rownames(polydf))),
    annot$cell_ID
  )]
  polydf$fov[is.na(polydf$fov)] = ""
  polydf = polydf[is.element(polydf$fov, df$fov), ]
  polydf = polydf[!is.na(polydf$tissue), ]
  
  # center the xy data within each FOV:
  for (fov in unique(df$fov)) {
    polydf$x[polydf$fov == fov] <- polydf$x[polydf$fov == fov] - min(df$x[df$fov == fov])
    polydf$y[polydf$fov == fov] <- polydf$y[polydf$fov == fov] - min(df$y[df$fov == fov])
    df$x[df$fov == fov] <- df$x[df$fov == fov] - min(df$x[df$fov == fov])
    df$y[df$fov == fov] <- df$y[df$fov == fov] - min(df$y[df$fov == fov])
   }
  # offset the mean fovs:
  meanfovs <- df$fov[df$fovtype == "mean"]
  df$y[is.element(df$fov, meanfovs)] <- df$y[is.element(df$fov, meanfovs)] - 0.8
  polydf$y[is.element(polydf$fov, meanfovs)] <- polydf$y[is.element(polydf$fov, meanfovs)] - 0.8
  
  # plot:
  g <- ggplot(df, aes(x = x, y = y, color = score99, size = I(0.05))) +
    geom_point() + 
    geom_polygon(data = polydf, aes(x = x, y = y, group = name), 
                 fill = NA, color = "cyan", linewidth = .25) +
    scale_color_viridis_c(option = "B") +
    labs(title = "",
         x = "",
         y = "",
         color = "") +  theme_few() + coord_fixed(ratio = 1) + 
    facet_wrap(~tissue, nrow = 2) +
    theme(axis.text = element_text(size = 4),
          legend.position = "none") 
  png(paste0("results/insitucor/fovsubset ", name, ".png"), width = 6, height = 3.2, units = "in", res = 600)
  print(g)
  dev.off()
  
}


### custom plot for paired mac scores:
# get median and max scores per tissue:
macscore = res$scores_env[, "C1QA_CD163_11"] / mean(res$scores_env[, "C1QA_CD163_11"]) +
  res$scores_env[, "CD74_HLA.DRA_9"] / mean(res$scores_env[, "CD74_HLA.DRA_9"])
  
keepcells <- c()
fovpairs <- c()
for (tiss in unique(annot$tissuename)) {
  inds <- annot$tissuename == tiss
  means <- by(macscore[inds], annot$fov[inds], mean)
  keepfovs <- names(means)[order(means, decreasing = T)[c(1, floor(length(means) / 2))]]
  keepcells <- c(keepcells, which(inds & (is.element(annot$fov, keepfovs))))
  fovpairs <- rbind(fovpairs, c(tiss, keepfovs[1], "max"))
  fovpairs <- rbind(fovpairs, c(tiss, keepfovs[2], "mean"))
}

# make df for ggplot:
df <- data.frame(tissue = annot$tissuename[keepcells],
                 score99.h = pmin(res$scores_env[keepcells, "CD74_HLA.DRA_9"] / quantile(res$scores_env[keepcells, "CD74_HLA.DRA_9"], 0.99), 1),
                 score99.c = pmin(res$scores_env[keepcells, "C1QA_CD163_11"] / quantile(res$scores_env[keepcells, "C1QA_CD163_11"], 0.99), 1),
                 fov = annot$fov[keepcells],
                 x = customlocs[keepcells, 1],
                 y = customlocs[keepcells, 2])
df$fovtype <- fovpairs[,3][match(df$fov, fovpairs[, 2])]
# rgb color:
df$col = rgb(df$score99.c, 0, df$score99.h, 1)
# df of polygon data:
polydf <- data.frame(
  x = unlist(sapply(gpoly.customlocs, function(x){x[,1]})),
  y = unlist(sapply(gpoly.customlocs, function(x){x[,2]})),
  name = unlist(sapply(1:length(gpoly.customlocs), 
                       function(i){rep(names(gpoly.customlocs)[i], 
                                       nrow(gpoly.customlocs[[i]]))})))
# associate with tissues:
polydf$tissue <- annot$tissuename[match(
  substr(rownames(polydf), unlist(gregexpr("c", rownames(polydf))), nchar(rownames(polydf))),
  annot$cell_ID
)]
polydf$fov <- annot$fov[match(
  substr(rownames(polydf), 
         sapply(rownames(polydf), function(str){gregexpr("c", str)[[1]][1]}),
         #unlist(gregexpr("c", rownames(polydf))), 
         nchar(rownames(polydf))),
  annot$cell_ID
)]
polydf$fov[is.na(polydf$fov)] = ""
polydf = polydf[is.element(polydf$fov, df$fov), ]
polydf = polydf[!is.na(polydf$tissue), ]

# center the xy data within each FOV:
for (fov in unique(df$fov)) {
  polydf$x[polydf$fov == fov] <- polydf$x[polydf$fov == fov] - min(df$x[df$fov == fov])
  polydf$y[polydf$fov == fov] <- polydf$y[polydf$fov == fov] - min(df$y[df$fov == fov])
  df$x[df$fov == fov] <- df$x[df$fov == fov] - min(df$x[df$fov == fov])
  df$y[df$fov == fov] <- df$y[df$fov == fov] - min(df$y[df$fov == fov])
}
# offset the mean fovs:
meanfovs <- df$fov[df$fovtype == "mean"]
df$y[is.element(df$fov, meanfovs)] <- df$y[is.element(df$fov, meanfovs)] - 0.8
polydf$y[is.element(polydf$fov, meanfovs)] <- polydf$y[is.element(polydf$fov, meanfovs)] - 0.8

# plot:
g <- ggplot(df, aes(x = x, y = y, color = I(col), size = I(0.05))) +
  geom_point() + 
  geom_polygon(data = polydf, aes(x = x, y = y, group = name), 
               fill = NA, color = "yellow", linewidth = .25) +
  scale_color_viridis_c(option = "B") +
  labs(title = "",
       x = "",
       y = "",
       color = "") +  theme_few() + coord_fixed(ratio = 1) + 
  facet_wrap(~tissue, nrow = 2) +
  theme(axis.text = element_text(size = 4),
        legend.position = "none") 
png("results/insitucor/fovsubset - macrophages - .png", width = 6, height = 3.2, units = "in", res = 600)
print(g)
dev.off()

svg("results/insitucor/fovsubset - macrophages - legend.svg", width = 3, height = 2)
par(bg = "black")
par(mar = c(0,0,0,0))
frame()
legend("center", pch = 16,
       col = c(rgb(0.1, 0, 0.9, 1),
               rgb(0.9, 0, 0.1, 1),
               rgb(0.9, 0, 0.9, 1)),
       legend = c("high compliments",
                  "high HLA class II",
                  "high both"),
       text.col = "white")
dev.off()

