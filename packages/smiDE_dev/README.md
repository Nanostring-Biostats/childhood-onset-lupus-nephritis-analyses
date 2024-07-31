Example DE usage
================
Dan McGuire
2023-05-11

``` r
library(data.table); setDTthreads(1)
library(smiDE)
```

### 0. Load example dataset

A small subset of giotto object from NSCLC data. This is FOV 1 from
Lung6, Lung12, and Lung13.

``` r
datadir <- system.file("extdata", package="smiDE")
gem <- readRDS(paste0(datadir, "/small_nsclc.rds"))
```

### 1. pre-de step: Calculate all cell-cell neighbors within a specified radius of each other. For a cell-type of interest, calculate expression of genes in neighbors of those cells.

This can be done in one step. For example, if one wanted to find DE
genes in fibroblast cells, comparing spatial regions specified by the
column `niche` in the meta data.

``` r
config_pre_de <- list(assay = "rna"
                       ,counts = "raw" 
                       ,normalized = "normalized"
                       ,cellid_colname = "cell_ID"
                       ,cell_type_metadata_colname = "cell_type"
                       ,split_neighbors_by_colname = "tissue"
                       ,mm_radius = 0.05
                       ,ref_celltype = "fibroblast"
                       ,sdimx_colname = "sdimx"
                       ,sdimy_colname = "sdimy"
                       ,weight_colname = "weight"
                       ,contamination = "sum"
                       ,verbose=TRUE
                       )

pre_de_obj <- run_pre_de(gem, config_pre_de) 
#> Loading required package: Giotto
#> Loading required package: Matrix
#> 2023-05-11 23:38:43, identifying cell-cell spatial neighbors within 0.05 radius.
#> neighbors calculated for tissue: Lung6
#> neighbors calculated for tissue: Lung12
#> neighbors calculated for tissue: Lung13
#> 2023-05-11 23:38:44. creating neighborhood list object
#> Measuring neighbor expression for fibroblast
#> (cluster 1 of 1)
#> 'as(<dgCMatrix>, "dgTMatrix")' is deprecated.
#> Use 'as(., "TsparseMatrix")' instead.
#> See help("Deprecated") and help("Matrix-deprecated").
names(pre_de_obj)
#> [1] "nblist"            "cell_adjacency_dt"
```

The `config_de` below specifies the model formula, for a negative
binomial regression, where for each gene, we include the expression of
neighboring cells which are not fibroblasts as a covariate to adjust
for. `niche` is a primary variable of interest, and tissue is a fixed
effect below which accounts for some gene expression attributable to
differences in the tissue / biological sample rather than the spatial
region (`niche`).

If we remove the `targets` parameter or set to `NULL`, then a model
would be fit for all genes.

``` r
config_de <- list(
                   assay = "rna"
                  ,data_use = "raw"
                  ,formula = ~RankNorm(otherct_expr) + niche + tissue + offset(log(totalcounts))
                  ,groupVar="niche"
                  ,family="nbinom2"
                  ,targets=rownames(gem@expression$rna$raw)[1:20]
                  )

de_obj <- 
  run_de(gem, pre_de_obj, config_de)
```

Fit negative binomial mixed models, with random effect for tissue
variable instead of a fixed effects.

``` r
config_de <- list(
                   assay = "rna"
                  ,data_use = "raw"
                  ,formula = ~RankNorm(otherct_expr) + niche + (1 | tissue) + offset(log(totalcounts))
                  ,groupVar="niche"
                  ,family="nbinom2"
                  ,targets=rownames(gem@expression$rna$raw)[1:20]
                  )

de_obj_mixedmod <- 
  run_de(gem, pre_de_obj, config_de)
#> Remove  0  genes having low expression.
#> Analyzing  1  genes with  3  subjects and  1097  cells.
#> Loading required package: foreach
#> Loading required package: future
#> Loading required package: rngtools
#> Remove  0  genes having low expression.
#> Analyzing  1  genes with  3  subjects and  1097  cells.
#> Remove  0  genes having low expression.
#> Analyzing  1  genes with  3  subjects and  1097  cells.
#> Remove  0  genes having low expression.
#> Analyzing  1  genes with  3  subjects and  1097  cells.
#> Remove  0  genes having low expression.
#> Analyzing  1  genes with  3  subjects and  1097  cells.
#> Remove  0  genes having low expression.
#> Analyzing  1  genes with  3  subjects and  1097  cells.
#> Remove  0  genes having low expression.
#> Analyzing  1  genes with  3  subjects and  1097  cells.
#> Remove  0  genes having low expression.
#> Analyzing  1  genes with  3  subjects and  1097  cells.
#> Remove  0  genes having low expression.
#> Analyzing  1  genes with  3  subjects and  1097  cells.
#> Remove  0  genes having low expression.
#> Analyzing  1  genes with  3  subjects and  1097  cells.
#> Remove  0  genes having low expression.
#> Analyzing  1  genes with  3  subjects and  1097  cells.
#> Remove  0  genes having low expression.
#> Analyzing  1  genes with  3  subjects and  1097  cells.
#> Remove  0  genes having low expression.
#> Analyzing  1  genes with  3  subjects and  1097  cells.
#> Remove  0  genes having low expression.
#> Analyzing  1  genes with  3  subjects and  1097  cells.
#> Remove  0  genes having low expression.
#> Analyzing  1  genes with  3  subjects and  1097  cells.
#> Remove  0  genes having low expression.
#> Analyzing  1  genes with  3  subjects and  1097  cells.
#> Remove  0  genes having low expression.
#> Analyzing  1  genes with  3  subjects and  1097  cells.
#> Remove  0  genes having low expression.
#> Analyzing  1  genes with  3  subjects and  1097  cells.
#> Remove  0  genes having low expression.
#> Analyzing  1  genes with  3  subjects and  1097  cells.
#> Remove  0  genes having low expression.
#> Analyzing  1  genes with  3  subjects and  1097  cells.
```

Note here that we are fitting models with a very small (1,097) number of
cells, and with all categories of niche (), some of which have few cells
per niche, which can cause stability issues in model fitting.

``` r
gem@cell_metadata$rna[,.N,by=.(cell_type,niche)][order(-N)][cell_type=="fibroblast"]
#>     cell_type                       niche   N
#> 1: fibroblast                      immune 399
#> 2: fibroblast       tumor-stroma boundary 226
#> 3: fibroblast              tumor interior 193
#> 4: fibroblast                 macrophages 143
#> 5: fibroblast                 neutrophils  48
#> 6: fibroblast plasmablast-enriched stroma  40
#> 7: fibroblast                      stroma  38
#> 8: fibroblast          lymphoid structure   5
#> 9: fibroblast     myeloid-enriched stroma   5
gem@cell_metadata$rna[cell_type=="fibroblast",.N,by=.(cell_type,niche)][,sum(N)]
#> [1] 1097
```

To check out the results, we can use the `results` function.

Some documentation on options here:

``` r
help(results)
```

For example, to look at comparisons of each spatial `niche` compared to
others:

``` r
results(de_obj, comparisons = "one.vs.rest", variable = "niche")
#> $one.vs.rest
#>                                      contrast        ratio           SE  df
#>   1:                      immune vs. avg.rest 1.160359e+00 2.512442e+03 Inf
#>   2:          lymphoid structure vs. avg.rest 3.582174e-11 7.483889e-06 Inf
#>   3:                 macrophages vs. avg.rest 9.324604e-01 1.477206e+03 Inf
#>   4:     myeloid-enriched stroma vs. avg.rest 4.654640e-11 1.016789e-05 Inf
#>   5:                 neutrophils vs. avg.rest 1.294860e+00 1.865547e+03 Inf
#>  ---                                                                       
#> 176:                 neutrophils vs. avg.rest 1.750875e+00 8.418212e+04 Inf
#> 177: plasmablast-enriched stroma vs. avg.rest 1.833648e+00 8.749456e+04 Inf
#> 178:                      stroma vs. avg.rest 7.131330e+00 3.396368e+05 Inf
#> 179:              tumor interior vs. avg.rest 2.360016e+00 1.316699e+05 Inf
#> 180:       tumor-stroma boundary vs. avg.rest 1.790789e+00 1.036971e+05 Inf
#>      null       z.ratio   p.value  fold_change  term target
#>   1:    1  6.869016e-05 0.9999452 1.160359e+00 niche   AATK
#>   2:    1 -1.151275e-04 0.9999081 3.582174e-11 niche   AATK
#>   3:    1 -4.414123e-05 0.9999648 9.324604e-01 niche   AATK
#>   4:    1 -1.089081e-04 0.9999131 4.654640e-11 niche   AATK
#>   5:    1  1.793548e-04 0.9998569 1.294860e+00 niche   AATK
#>  ---                                                       
#> 176:    1  1.164966e-05 0.9999907 1.750875e+00 niche ADGRE1
#> 177:    1  1.270655e-05 0.9999899 1.833648e+00 niche ADGRE1
#> 178:    1  4.124843e-05 0.9999671 7.131330e+00 niche ADGRE1
#> 179:    1  1.539054e-05 0.9999877 2.360016e+00 niche ADGRE1
#> 180:    1  1.006215e-05 0.9999920 1.790789e+00 niche ADGRE1
```

``` r
results(de_obj_mixedmod, comparisons = "one.vs.rest", variable = "niche")
#> $one.vs.rest
#>                                      contrast        ratio           SE  df
#>   1:                      immune vs. avg.rest 7.529357e-01  9.401866581 Inf
#>   2:          lymphoid structure vs. avg.rest 1.679008e-06  0.001912346 Inf
#>   3:                 macrophages vs. avg.rest 6.753535e-01  6.181835397 Inf
#>   4:     myeloid-enriched stroma vs. avg.rest 1.791263e-06  0.002361590 Inf
#>   5:                 neutrophils vs. avg.rest 1.835044e+00 15.292438414 Inf
#>  ---                                                                       
#> 176:                 neutrophils vs. avg.rest 2.472171e+00 24.996431175 Inf
#> 177: plasmablast-enriched stroma vs. avg.rest 2.566615e+00 25.764928490 Inf
#> 178:                      stroma vs. avg.rest 5.633489e+00 56.364306716 Inf
#> 179:              tumor interior vs. avg.rest 1.236094e+00 14.484663076 Inf
#> 180:       tumor-stroma boundary vs. avg.rest 2.189612e+00 26.592458214 Inf
#>      null     z.ratio   p.value  fold_change  term target
#>   1:    1 -0.02272577 0.9818690 7.529357e-01 niche   AATK
#>   2:    1 -0.01167481 0.9906851 1.679008e-06 niche   AATK
#>   3:    1 -0.04288194 0.9657956 6.753535e-01 niche   AATK
#>   4:    1 -0.01003690 0.9919918 1.791263e-06 niche   AATK
#>   5:    1  0.07284628 0.9419284 1.835044e+00 niche   AATK
#>  ---                                                     
#> 176:    1  0.08951497 0.9286727 2.472171e+00 niche ADGRE1
#> 177:    1  0.09389740 0.9251907 2.566615e+00 niche ADGRE1
#> 178:    1  0.17278266 0.8628223 5.633489e+00 niche ADGRE1
#> 179:    1  0.01808793 0.9855687 1.236094e+00 niche ADGRE1
#> 180:    1  0.06453153 0.9485470 2.189612e+00 niche ADGRE1
```

To get pairwise comparisons for fold change comparing tissue.

``` r
results(de_obj, comparisons="pairwise", variable="tissue")[[1]][1:10]
#>            contrast     ratio        SE  df null    z.ratio     p.value
#>  1: Lung12 / Lung13 1.9678351 1.4162603 Inf    1  0.9405718 0.346924334
#>  2:  Lung12 / Lung6 3.1628070 3.8891143 Inf    1  0.9364203 0.349056829
#>  3:  Lung13 / Lung6 1.6072521 2.2279235 Inf    1  0.3423290 0.732103324
#>  4: Lung12 / Lung13 0.7070481 0.2706988 Inf    1 -0.9054450 0.365229708
#>  5:  Lung12 / Lung6 0.3276122 0.2750379 Inf    1 -1.3292369 0.183769836
#>  6:  Lung13 / Lung6 0.4633521 0.4107256 Inf    1 -0.8678348 0.385484790
#>  7: Lung12 / Lung13 0.3241189 0.1377688 Inf    1 -2.6505780 0.008035417
#>  8:  Lung12 / Lung6 1.4042977 1.2913603 Inf    1  0.3692320 0.711954817
#>  9:  Lung13 / Lung6 4.3326620 4.1847687 Inf    1  1.5179983 0.129014835
#> 10: Lung12 / Lung13 0.5015197 0.2514298 Inf    1 -1.3765468 0.168652420
#>     fold_change   term target
#>  1:   1.9678351 tissue   AATK
#>  2:   3.1628070 tissue   AATK
#>  3:   1.6072521 tissue   AATK
#>  4:   0.7070481 tissue   ABL1
#>  5:   0.3276122 tissue   ABL1
#>  6:   0.4633521 tissue   ABL1
#>  7:   0.3241189 tissue   ABL2
#>  8:   1.4042977 tissue   ABL2
#>  9:   4.3326620 tissue   ABL2
#> 10:   0.5015197 tissue    ACE
```

To check out the model summary, say for the “ADGRB2” gene which gave a
warning for ‘bad fit’

``` r
results(de_obj, comparisons="model_summary", targets="ADGRB2")$model_summary
#>                                 term          est           se             z
#>  1:                      (Intercept)   -6.9428979 1.115571e+00 -6.223630e+00
#>  2:           RankNorm(otherct_expr)    0.4167882 2.126407e-01  1.960059e+00
#>  3:          nichelymphoid structure  -31.2360925 2.605508e+06 -1.198848e-05
#>  4:                 nichemacrophages   -0.1956178 7.077558e-01 -2.763917e-01
#>  5:     nichemyeloid-enriched stroma  -30.1344737 2.673545e+06 -1.127135e-05
#>  6:                 nicheneutrophils   -1.2820174 1.271037e+00 -1.008639e+00
#>  7: nicheplasmablast-enriched stroma   -2.0623221 1.380254e+00 -1.494161e+00
#>  8:                      nichestroma  -28.2634106 9.582602e+05 -2.949451e-05
#>  9:              nichetumor interior   -2.1676230 1.440364e+00 -1.504914e+00
#> 10:       nichetumor-stroma boundary   -0.7706849 1.076428e+00 -7.159652e-01
#> 11:                     tissueLung13   -2.2981408 1.064439e+00 -2.159015e+00
#> 12:                      tissueLung6    0.4466936 1.129780e+00  3.953810e-01
#> 13:                            theta    0.2847745 1.696979e-01            NA
#> 14:                           logLik -151.2148528           NA            NA
#>             pval    termtype user.self sys.self elapsed target             msg
#>  1: 4.857825e-10       fixed     10.15   15.152   0.468 ADGRB2 converged=TRUE 
#>  2: 4.998893e-02       fixed     10.15   15.152   0.468 ADGRB2 converged=TRUE 
#>  3: 9.999904e-01       fixed     10.15   15.152   0.468 ADGRB2 converged=TRUE 
#>  4: 7.822472e-01       fixed     10.15   15.152   0.468 ADGRB2 converged=TRUE 
#>  5: 9.999910e-01       fixed     10.15   15.152   0.468 ADGRB2 converged=TRUE 
#>  6: 3.131478e-01       fixed     10.15   15.152   0.468 ADGRB2 converged=TRUE 
#>  7: 1.351335e-01       fixed     10.15   15.152   0.468 ADGRB2 converged=TRUE 
#>  8: 9.999765e-01       fixed     10.15   15.152   0.468 ADGRB2 converged=TRUE 
#>  9: 1.323463e-01       fixed     10.15   15.152   0.468 ADGRB2 converged=TRUE 
#> 10: 4.740129e-01       fixed     10.15   15.152   0.468 ADGRB2 converged=TRUE 
#> 11: 3.084898e-02       fixed     10.15   15.152   0.468 ADGRB2 converged=TRUE 
#> 12: 6.925618e-01       fixed     10.15   15.152   0.468 ADGRB2 converged=TRUE 
#> 13:           NA       fixed     10.15   15.152   0.468 ADGRB2 converged=TRUE 
#> 14:           NA model_stats     10.15   15.152   0.468 ADGRB2 converged=TRUE
```

Compare with the ‘AATK’ gene, which didn’t give a warning.

``` r
results(de_obj, comparisons="model_summary", targets="AATK")$model_summary
#>                                 term           est           se             z
#>  1:                      (Intercept)   -8.51678983 7.931928e-01 -10.737350954
#>  2:           RankNorm(otherct_expr)    0.10369471 2.226824e-01   0.465661914
#>  3:          nichelymphoid structure  -24.03747163 2.089179e+05  -0.000115057
#>  4:                 nichemacrophages   -0.15544696 6.870626e-01  -0.226248620
#>  5:     nichemyeloid-enriched stroma  -23.77677063 2.184442e+05  -0.000108846
#>  6:                 nicheneutrophils    0.15246187 1.041453e+00   0.146393417
#>  7: nicheplasmablast-enriched stroma    0.20353068 1.053066e+00   0.193274319
#>  8:                      nichestroma   -0.02814041 1.101373e+00  -0.025550305
#>  9:              nichetumor interior    0.40504499 9.888351e-01   0.409618351
#> 10:       nichetumor-stroma boundary    0.28726984 7.173280e-01   0.400472108
#> 11:                     tissueLung13   -0.67693400 7.197048e-01  -0.940571810
#> 12:                      tissueLung6   -1.15145993 1.229640e+00  -0.936420301
#> 13:                            theta    0.29709920 2.010725e-01            NA
#> 14:                           logLik -145.68567492           NA            NA
#>             pval    termtype user.self sys.self elapsed target             msg
#>  1: 6.796663e-27       fixed     6.028     8.79   0.288   AATK converged=TRUE 
#>  2: 6.414575e-01       fixed     6.028     8.79   0.288   AATK converged=TRUE 
#>  3: 9.999082e-01       fixed     6.028     8.79   0.288   AATK converged=TRUE 
#>  4: 8.210081e-01       fixed     6.028     8.79   0.288   AATK converged=TRUE 
#>  5: 9.999132e-01       fixed     6.028     8.79   0.288   AATK converged=TRUE 
#>  6: 8.836108e-01       fixed     6.028     8.79   0.288   AATK converged=TRUE 
#>  7: 8.467441e-01       fixed     6.028     8.79   0.288   AATK converged=TRUE 
#>  8: 9.796160e-01       fixed     6.028     8.79   0.288   AATK converged=TRUE 
#>  9: 6.820859e-01       fixed     6.028     8.79   0.288   AATK converged=TRUE 
#> 10: 6.888088e-01       fixed     6.028     8.79   0.288   AATK converged=TRUE 
#> 11: 3.469243e-01       fixed     6.028     8.79   0.288   AATK converged=TRUE 
#> 12: 3.490568e-01       fixed     6.028     8.79   0.288   AATK converged=TRUE 
#> 13:           NA       fixed     6.028     8.79   0.288   AATK converged=TRUE 
#> 14:           NA model_stats     6.028     8.79   0.288   AATK converged=TRUE
```

## Deeper explanations…

##### The pre_de step:

The object `pre_de_obj$cell_adjacency_dt` contains a row for each pair
of cells within the 0.05 distance of each other specified by
`mm_radius=0.05`. The `split_neighbors_by_colname=tissue` specifies that
neighbors are calculated separately by by the column in meta data named
“tissue”.

``` r

config_pre_de <- list(assay = "rna"
                       ,counts = "raw" 
                       ,normalized = "normalized"
                       ,cellid_colname = "cell_ID"
                       ,cell_type_metadata_colname = "cell_type"
                       ,split_neighbors_by_colname = "tissue"
                       ,mm_radius = 0.05
                       ,ref_celltype = "fibroblast"
                       ,sdimx_colname = "sdimx"
                       ,sdimy_colname = "sdimy"
                       ,weight_colname = "weight"
                       ,contamination = "sum"
                       ,verbose=TRUE
                       )


show <-copy(pre_de_obj$cell_adjacency_dt)
show <- merge(show, gem@cell_metadata$rna[,.(from=cell_ID,cell_type_from = cell_type)], by="from", sort=FALSE)
show <- merge(show, gem@cell_metadata$rna[,.(to=cell_ID,cell_type_to = cell_type)], by="to", sort=FALSE)
show[]
#>                to       from sdimx_begin sdimy_begin sdimx_end sdimy_end
#>      1:  c_4_1_10   c_4_1_11    29.13752   -18.97240  29.17568 -18.97204
#>      2:  c_4_1_10 c_4_1_1458    29.12618   -18.97492  29.17568 -18.97204
#>      3:  c_4_1_10 c_4_1_1459    29.15030   -18.97996  29.17568 -18.97204
#>      4:  c_4_1_10 c_4_1_1460    29.15966   -18.97186  29.17568 -18.97204
#>      5:  c_4_1_10 c_4_1_1461    29.16812   -18.97582  29.17568 -18.97204
#>     ---                                                                 
#> 196198: c_8_1_993  c_8_1_995    15.11668    -2.87298  15.13630  -2.87082
#> 196199: c_8_1_993  c_8_1_998    15.16204    -2.87334  15.13630  -2.87082
#> 196200: c_8_1_994  c_8_1_999    15.19444    -2.87496  15.21532  -2.87262
#> 196201: c_8_1_995  c_8_1_998    15.16204    -2.87334  15.11668  -2.87298
#> 196202: c_8_1_998  c_8_1_999    15.19444    -2.87496  15.16204  -2.87334
#>            distance    weight tissue_from tissue_to cell_type_from cell_type_to
#>      1: 0.038161698  26.20428       Lung6     Lung6        tumor 6      tumor 6
#>      2: 0.049583711  20.16791       Lung6     Lung6        tumor 6      tumor 6
#>      3: 0.026587042  37.61231       Lung6     Lung6        tumor 6      tumor 6
#>      4: 0.016021011  62.41803       Lung6     Lung6        tumor 6      tumor 6
#>      5: 0.008452337 118.31047       Lung6     Lung6        tumor 6      tumor 6
#>     ---                                                                        
#> 196198: 0.019738541  50.66231      Lung13    Lung13       tumor 13   fibroblast
#> 196199: 0.025863062  38.66518      Lung13    Lung13           mast   fibroblast
#> 196200: 0.021010712  47.59477      Lung13    Lung13   T CD4 memory     tumor 13
#> 196201: 0.045361429  22.04516      Lung13    Lung13           mast     tumor 13
#> 196202: 0.032440475  30.82569      Lung13    Lung13   T CD4 memory         mast
```

Because the ref_celltype is specified as fibroblast
(`ref_celltype=fibroblast`), the total normalized, distance weighted
expression of each gene is calculated for the neighbors of each of the
1097 fibroblast cells in this small sample dataset.

``` r
lapply(pre_de_obj$nblist$neighbor_expr_byct, dim)
#> $allct
#> [1]  960 1097
#> 
#> $otherct
#> [1]  960 1097
#> 
#> $tumor_6
#> [1]  960 1097
#> 
#> $T_CD4_memory
#> [1]  960 1097
#> 
#> $fibroblast
#> [1]  960 1097
#> 
#> $tumor_12
#> [1]  960 1097
#> 
#> $T_CD4_naive
#> [1]  960 1097
#> 
#> $macrophage
#> [1]  960 1097
#> 
#> $Treg
#> [1]  960 1097
#> 
#> $neutrophil
#> [1]  960 1097
#> 
#> $tumor_5
#> [1]  960 1097
#> 
#> $pDC
#> [1]  960 1097
#> 
#> $T_CD8_naive
#> [1]  960 1097
#> 
#> $epithelial
#> [1]  960 1097
#> 
#> $tumor_9
#> [1]  960 1097
#> 
#> $NK
#> [1]  960 1097
#> 
#> $tumor_13
#> [1]  960 1097
#> 
#> $endothelial
#> [1]  960 1097
#> 
#> $plasmablast
#> [1]  960 1097
#> 
#> $mDC
#> [1]  960 1097
#> 
#> $mast
#> [1]  960 1097
#> 
#> $T_CD8_memory
#> [1]  960 1097
#> 
#> $B_cell
#> [1]  960 1097
#> 
#> $monocyte
#> [1]  960 1097
```

This is calculated separately for all of the cell types. For example,
the total expression of endothelial cells which are the neighbors of the
fibroblast cells can be accessed like this.

``` r
pre_de_obj$nblist$neighbor_expr_byct$endothelial[1:20,1:20]
#> 20 x 20 sparse Matrix of class "dgCMatrix"
#>   [[ suppressing 20 column names 'c_4_1_101', 'c_4_1_1021', 'c_4_1_1033' ... ]]
#>                                                                     
#> AATK     .      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 . 0   .        .      
#> ABL1     .      . . . . . . . . . . . . . . . . .   .       59.98763
#> ABL2     .      . . . . . . . . . . . . . . . . .   .        .      
#> ACE      .      . . . . . . . . . . . . . . . . . 155.16464  .      
#> ACE2   208.1959 . . . . . . . . . . . . . . . . .   .        .      
#> ACKR1    .      . . . . . . . . . . . . . . . . .   .        .      
#> ACKR3    .      . . . . . . . . . . . . . . . . . 118.62663  .      
#> ACKR4    .      . . . . . . . . . . . . . . . . . 118.62663 59.98763
#> ACTA2    .      . . . . . . . . . . . . . . . . .   .        .      
#> ACTG2  333.1853 . . . . . . . . . . . . . . . . .   .        .      
#> ACVR1    .      . . . . . . . . . . . . . . . . .   .        .      
#> ACVR1B   .      . . . . . . . . . . . . . . . . .   .        .      
#> ACVR2A   .      . . . . . . . . . . . . . . . . .   .        .      
#> ACVRL1   .      . . . . . . . . . . . . . . . . .   .        .      
#> ADGRA2   .      . . . . . . . . . . . . . . . . . 118.62663  .      
#> ADGRA3   .      . . . . . . . . . . . . . . . . .   .        .      
#> ADGRB2   .      . . . . . . . . . . . . . . . . .   .        .      
#> ADGRB3   .      . . . . . . . . . . . . . . . . .   .        .      
#> ADGRD1   .      . . . . . . . . . . . . . . . . .  40.78229  .      
#> ADGRE1   .      . . . . . . . . . . . . . . . . .   .        .
```

There’s also a list for expression from any cell type *other* than the
fibroblast cells (*otherct*). There’s also a list for expression from
all cell types including other fibroblast cells (*allct*).

``` r
pre_de_obj$nblist$neighbor_expr_byct$otherct[1:10,1:10]
#> 10 x 10 sparse Matrix of class "dgCMatrix"
#>   [[ suppressing 10 column names 'c_4_1_101', 'c_4_1_1021', 'c_4_1_1033' ... ]]
#>                                                                           
#> AATK   30.17248   .         .         .         .         .         .     
#> ABL1  222.37865 248.29146 176.90788 136.65606 311.92902 611.28069   .     
#> ABL2  124.36308 143.02026 160.27138  92.54601   .        30.44504 185.6214
#> ACE   174.67647   .         .         .       122.25106 313.12621   .     
#> ACE2  208.19591   .        10.44656  17.91831   .         .       163.7922
#> ACKR1  30.17248   .         .         .        28.83073  40.30382   .     
#> ACKR3   .         .        10.44656   .         .         .         .     
#> ACKR4   .        44.71815  49.65259  67.37350  14.56298   .         .     
#> ACTA2   .         .         .         .         .         .         .     
#> ACTG2 333.18530   .         .         .         .        21.07011   .     
#>                                    
#> AATK    .         .         .      
#> ABL1  150.09725 421.54171 298.17699
#> ABL2   81.53878  79.18371 125.67410
#> ACE   167.88788 137.31107  66.67863
#> ACE2    .         .         .      
#> ACKR1   .        16.45510  31.59523
#> ACKR3   .         .         .      
#> ACKR4   .        30.80876  30.71691
#> ACTA2   .        22.82957  23.48007
#> ACTG2   .        49.83408  41.68309
pre_de_obj$nblist$neighbor_expr_byct$allct[1:10,1:10]
#> 10 x 10 sparse Matrix of class "dgCMatrix"
#>   [[ suppressing 10 column names 'c_4_1_101', 'c_4_1_1021', 'c_4_1_1033' ... ]]
#>                                                                           
#> AATK   30.17248   .         .         .         .         .         .     
#> ABL1  222.37865 535.73222 176.90788 390.58938 311.92902 611.28069   .     
#> ABL2  124.36308 143.02026 160.27138  92.54601   .        30.44504 406.3436
#> ACE   473.04048   .         .         .       122.25106 313.12621   .     
#> ACE2  208.19591   .        10.44656  17.91831   .         .       163.7922
#> ACKR1  30.17248   .         .         .        28.83073  40.30382   .     
#> ACKR3   .         .        10.44656   .         .         .         .     
#> ACKR4   .        44.71815  49.65259  67.37350  14.56298   .         .     
#> ACTA2   .         .         .         .         .         .         .     
#> ACTG2 333.18530   .         .         .         .        21.07011   .     
#>                                    
#> AATK    .         .         .      
#> ABL1  150.09725 421.54171 298.17699
#> ABL2   81.53878  79.18371 125.67410
#> ACE   723.52038 137.31107  66.67863
#> ACE2    .         .         .      
#> ACKR1   .        16.45510  31.59523
#> ACKR3   .         .         .      
#> ACKR4   .        30.80876  30.71691
#> ACTA2   .       341.32640  23.48007
#> ACTG2   .        49.83408  41.68309
```

Using the *normalized* counts and weighting the expression by it’s
distance to the cell of interest (`weight_colname="weight"`) was found
to work well in some real applications and is set as a default. In some
applications, this helped best for controlling known false positive
genes that were initially indicated to be DE due to cell-typing or
cell-segmentation error.

Under the hood , calculations are done using the function
`measure_neighbor_expression_by_celltype`. It can be easier to visualize
if we look at raw counts which are not weighted by distance.

``` r

metainfo <- copy(gem@cell_metadata$rna)
metainfo <- merge(metainfo, gem@spatial_locs$raw, by="cell_ID")
nblist_raw_unweighted <- 
  measure_neighbor_expression_by_celltype(
                                          assay_matrix = gem@expression$rna$raw
                                          ,metadata = metainfo
                                          ,cell_adjacencies = pre_de_obj$cell_adjacency_dt
                                          ,contamination = "sum"
                                          ,weight_colname = NULL
                                          ,ref_celltype = "fibroblast"
                                          ,cellid_colname = "cell_ID"
                                          ,cell_type_metadata_colname = "cell_type"
                                          )
```

Using the neighbor cell expression, we can decompose how much a gene is
expressed in neighboring cells, and by which cell types.  
For example, looking at the tumor gene “KRT17”, we can see that most of
the expressionof KRT17 in cells which are neighbors to fibroblast cells
are coming from tumor cells.

``` r
krt17_neighbor_expr <- 
  extract_neighborexpr_by_gene(nblist_raw_unweighted$neighbor_expr_byct, gene= "KRT17")

krt17_neighbor_expr[1:5]
#>    allct otherct tumor_6 T_CD4_memory fibroblast tumor_12 T_CD4_naive
#> 1:    51      51      49            2          0        0           0
#> 2:   161     161     160            0          0        1           0
#> 3:   225     225     225            0          0        0           0
#> 4:   142     139     138            0          3        1           0
#> 5:    81      77      72            0          4        0           0
#>    macrophage Treg neutrophil tumor_5 pDC T_CD8_naive epithelial tumor_9 NK
#> 1:          0    0          0       0   0           0          0       0  0
#> 2:          0    0          0       0   0           0          0       0  0
#> 3:          0    0          0       0   0           0          0       0  0
#> 4:          0    0          0       0   0           0          0       0  0
#> 5:          2    0          0       0   0           0          0       0  2
#>    tumor_13 endothelial plasmablast mDC mast T_CD8_memory B_cell monocyte
#> 1:        0           0           0   0    0            0      0        0
#> 2:        0           0           0   0    0            0      0        0
#> 3:        0           0           0   0    0            0      0        0
#> 4:        0           0           0   0    0            0      0        0
#> 5:        0           0           0   0    0            1      0        0
#>       cell_ID
#> 1:  c_4_1_101
#> 2: c_4_1_1021
#> 3: c_4_1_1033
#> 4: c_4_1_1065
#> 5: c_4_1_1074
```

A few other notes on the ‘nblist_raw_unweighted’ object.

The object ‘adjacency_mat’ is a sparse matrix indicating which cells are
neighbors/within 0.05 distance of each other.

``` r
nblist_raw_unweighted$adjacency_mat[1:10,1:10]
#> 10 x 10 sparse Matrix of class "dgCMatrix"
#>   [[ suppressing 10 column names 'c_4_1_10', 'c_4_1_100', 'c_4_1_1001' ... ]]
#>                               
#> c_4_1_10   . . . . . . . . . .
#> c_4_1_100  . . . . . . . . . .
#> c_4_1_1001 . . . . . . . . . .
#> c_4_1_1002 . . . . 1 . 1 1 1 .
#> c_4_1_1004 . . . 1 . . 1 1 1 .
#> c_4_1_1005 . . . . . . . . . .
#> c_4_1_1006 . . . 1 1 . . 1 1 .
#> c_4_1_1007 . . . 1 1 . 1 . 1 .
#> c_4_1_1008 . . . 1 1 . 1 1 . .
#> c_4_1_1009 . . . . . . . . . .
```

The object ‘adjacency_counts_by_ct’ indicates how many cells of each
cell type are ‘neighbors’ of each fibroblast cell.

``` r
nblist_raw_unweighted$adjacency_counts_by_ct
#>       B-cell NK T CD4 memory T CD4 naive T CD8 memory T CD8 naive Treg
#>    1:      0  0            0           0            0           0    0
#>    2:      0  0            1           0            0           1    0
#>    3:      0  0            0           0            0           0    0
#>    4:      0  0            0           0            0           0    0
#>    5:      0  0            0           0            0           0    0
#>   ---                                                                 
#> 9346:      0  0           13           0            5           0    3
#> 9347:      0  0            2           0            1           0    0
#> 9348:      0  1            0           0            2           0    0
#> 9349:      0  0           18           0            3           0    1
#> 9350:      0  0           12           0            2           0    1
#>       endothelial epithelial fibroblast mDC macrophage mast monocyte neutrophil
#>    1:           0          0          0   0          0    0        0          0
#>    2:           0          3          6   0          0    0        0          0
#>    3:           0          0          5   0          2    0        0          0
#>    4:           0          0          1   0          0    0        0          0
#>    5:           0          0          5   0          2    0        0          0
#>   ---                                                                          
#> 9346:           0          0          7   0         22    2        0          0
#> 9347:           0          0          4   0          9    0        0          1
#> 9348:           0          0          0   0          2    0        0          1
#> 9349:           2          0         13   0         11    2        0          0
#> 9350:           2          0         11   0          4    4        0          0
#>       pDC plasmablast tumor 12 tumor 13 tumor 5 tumor 6 tumor 9    cell_ID
#>    1:   0           0        0        0       0      25       0   c_4_1_10
#>    2:   0           0        0        0       0      16       0  c_4_1_100
#>    3:   0           0        0        0       0      33       0 c_4_1_1001
#>    4:   0           0        0        0       0      54       0 c_4_1_1002
#>    5:   0           0        0        0       0      45       0 c_4_1_1004
#>   ---                                                                     
#> 9346:  11           2        0        4       0       0       0  c_8_1_995
#> 9347:   0           0        0       13       0       0       0  c_8_1_996
#> 9348:   0           0        0       41       0       0       1  c_8_1_997
#> 9349:   2           2        0        8       0       0       0  c_8_1_998
#> 9350:   2           1        0       19       0       0       0  c_8_1_999
nblist_raw_unweighted$adjacency_counts_by_ct[,lapply(.SD,sum),.SDcols=1:22]
#>    B-cell   NK T CD4 memory T CD4 naive T CD8 memory T CD8 naive Treg
#> 1:   3905 1145        20841        4238         3831        1649 3803
#>    endothelial epithelial fibroblast  mDC macrophage mast monocyte neutrophil
#> 1:       13402      17543      49149 1473      42188 8156      282      30027
#>     pDC plasmablast tumor 12 tumor 13 tumor 5 tumor 6 tumor 9
#> 1: 6054       22945    39016    46095     224   75933     505
```

#### Modeling options

It’s possible to specify a number of different DE models.

In addition to trying different likelihoods (`family`) argument, any of
the cell types in the ‘nblist’ with ’\_expr’ suffix are available to use
as arguments in modeling formula.

``` r

## Neg binomial regression on raw counts adjusting specifically for neighboring expression in certain class of *tumor* cells
nbmod <- 
smi_de(assay_matrix = gem@expression$rna$raw
       ,metadata=metainfo
       ,neighborhood_counts = pre_de_obj$nblist
       ,groupVar="niche"
       ,targets="KRT17"
       ,formula = ~niche + tumor_12_expr  + offset(log(totalcounts))
       ,family="nbinom2"
       )

results(nbmod, comparisons="model_summary")
#> $model_summary
#>                                 term           est           se             z
#>  1:                      (Intercept) -8.019933e+00 1.919162e-01 -4.178873e+01
#>  2:          nichelymphoid structure  9.920198e-01 1.265765e+00  7.837317e-01
#>  3:                 nichemacrophages  1.750051e-02 3.731816e-01  4.689544e-02
#>  4:     nichemyeloid-enriched stroma -3.395112e+01 2.859729e+07 -1.187214e-06
#>  5:                 nicheneutrophils  1.159321e+00 4.176131e-01  2.776065e+00
#>  6: nicheplasmablast-enriched stroma -6.289158e-01 7.860827e-01 -8.000631e-01
#>  7:                      nichestroma -5.438803e-01 7.728622e-01 -7.037222e-01
#>  8:              nichetumor interior  2.159640e+00 2.703850e-01  7.987278e+00
#>  9:       nichetumor-stroma boundary  2.957090e-01 3.188968e-01  9.272874e-01
#> 10:                    tumor_12_expr  2.785111e-04 3.318112e-04  8.393663e-01
#> 11:                            theta  5.615629e-01 1.711268e-01            NA
#> 12:                           logLik -4.146546e+02           NA            NA
#>             pval    termtype user.self sys.self elapsed target             msg
#>  1: 0.000000e+00       fixed    11.757   17.187   0.546  KRT17 converged=TRUE 
#>  2: 4.331976e-01       fixed    11.757   17.187   0.546  KRT17 converged=TRUE 
#>  3: 9.625966e-01       fixed    11.757   17.187   0.546  KRT17 converged=TRUE 
#>  4: 9.999991e-01       fixed    11.757   17.187   0.546  KRT17 converged=TRUE 
#>  5: 5.502126e-03       fixed    11.757   17.187   0.546  KRT17 converged=TRUE 
#>  6: 4.236742e-01       fixed    11.757   17.187   0.546  KRT17 converged=TRUE 
#>  7: 4.816058e-01       fixed    11.757   17.187   0.546  KRT17 converged=TRUE 
#>  8: 1.379508e-15       fixed    11.757   17.187   0.546  KRT17 converged=TRUE 
#>  9: 3.537773e-01       fixed    11.757   17.187   0.546  KRT17 converged=TRUE 
#> 10: 4.012638e-01       fixed    11.757   17.187   0.546  KRT17 converged=TRUE 
#> 11:           NA       fixed    11.757   17.187   0.546  KRT17 converged=TRUE 
#> 12:           NA model_stats    11.757   17.187   0.546  KRT17 converged=TRUE

## Linear regression on normalized counts adjusting specifically for neighboring expression in certain class of *tumor* cells.  And a random effect for tissue.

linearmod <- 
smi_de(assay_matrix = gem@expression$rna$raw
       ,metadata=metainfo
       ,neighborhood_counts = pre_de_obj$nblist
       ,groupVar="niche"
       ,targets="KRT17"
       ,formula = ~niche + tumor_12_expr + (1|tissue) 
       ,family="gaussian"
       )

results(linearmod, comparisons="model_summary")
#> $model_summary
#>                                 term           est           se           t
#>  1:                      (Intercept)  1.620675e-01 1.268332e-01  1.27780060
#>  2:          nichelymphoid structure  1.794967e-01 2.018694e-01  0.88917218
#>  3:                 nichemacrophages  3.715583e-03 4.070916e-02  0.09127142
#>  4:     nichemyeloid-enriched stroma -1.893267e-01 1.959670e-01 -0.96611535
#>  5:                 nicheneutrophils  1.358642e-01 8.676098e-02  1.56596002
#>  6: nicheplasmablast-enriched stroma -2.989479e-02 9.488809e-02 -0.31505311
#>  7:                      nichestroma -2.756892e-02 7.091032e-02 -0.38878577
#>  8:              nichetumor interior -6.492975e-02 7.719303e-02 -0.84113485
#>  9:       nichetumor-stroma boundary -5.112555e-02 4.995546e-02 -1.02342263
#> 10:                    tumor_12_expr  3.494502e-04 8.571859e-05  4.07671379
#> 11:                         Residual  4.176800e-01           NA          NA
#> 12:                       1 | tissue  2.037000e-01           NA          NA
#> 13:                           logLik -6.229322e+02           NA          NA
#>                             termtype user.self sys.self elapsed target
#>  1:                            fixed     0.038        0   0.039  KRT17
#>  2:                            fixed     0.038        0   0.039  KRT17
#>  3:                            fixed     0.038        0   0.039  KRT17
#>  4:                            fixed     0.038        0   0.039  KRT17
#>  5:                            fixed     0.038        0   0.039  KRT17
#>  6:                            fixed     0.038        0   0.039  KRT17
#>  7:                            fixed     0.038        0   0.039  KRT17
#>  8:                            fixed     0.038        0   0.039  KRT17
#>  9:                            fixed     0.038        0   0.039  KRT17
#> 10:                            fixed     0.038        0   0.039  KRT17
#> 11: random,Residual,Residual,std.dev     0.038        0   0.039  KRT17
#> 12:     random,tissue,tissue,std.dev     0.038        0   0.039  KRT17
#> 13:                      model_stats     0.038        0   0.039  KRT17
```
