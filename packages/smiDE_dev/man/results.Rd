% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarize_model.R
\name{results}
\alias{results}
\title{Extract results from DE object by contrast comparison, model variable, gene}
\usage{
results(
  smide_results,
  comparisons = c("pairwise", "one.vs.rest", "one.vs.all", "emmeans", "model_summary",
    "all"),
  variable = c(smide_results$groupVar, smide_results$modelterms, "all"),
  targets = c("all")
)
}
\arguments{
\item{smide_results}{object of class "smide" returned by call to run_de or smi_de functions.}

\item{comparisons}{one or more of c("pairwise", "one.vs.rest", "one.vs.all", "emmeans", "model_summary").
Corresponding results from emmeans package are returned for variables which were included
as fixed effects in the DE model formula.
For continuous variables, only the 'emmeans' returns a result,
corresponding to that of the mean of the untransformed covariate and the corresponding model estimated marginal response of the outcome at that mean.
\itemize{
\item pairwise : pairwise contrasts between levels of a categorical fixed effect variable in DE model. See pairwise.emmc in ?emmeans::\code{contrast-methods} .
\item one.vs.rest : one category vs. average of all other category contrasts, for categorical fixed effect variables.

\if{html}{\out{<div class="sourceCode">}}\preformatted{          This contrast will be empty for categorical variables with only 2 levels, 
          as the pairwise contrast is equivalent. The weights argument 'wts' is set to correspond to # of cells in a category, 
          so that these contrasts don't emphasize niche or categorical groups with very few cells. See del.eff.emmc in ?emmeans::`contrast-methods`.  
          (also see here for more discussion on using weights : \\href\{https://github.com/rvlenth/emmeans/issues/346\}\{emmeans github issue\}).
}\if{html}{\out{</div>}}
\item one.vs.all : one category vs. average of all categories (including the one being compared to), for categorical fixed effect variables.
See eff.emmc in ?emmeans::\code{contrast-methods}.
This contrast will be empty for categorical variables with only 2 levels. The weights argument 'wts' is set to correspond to # of cells in a category,
so that these contrasts don't emphasize niche or categorical groups with very few cells.
\item emmeans : Model estimated marginal means of the response level, at different levels of the predictor.
Ratios of these can be used to estimate fold change. Lots of vignettes and examples available online.
See here, for example (\href{https://cran.r-project.org/web/packages/emmeans/vignettes/}{emmeans-vignettes})
\item model_summary : summary tables for fixed effects and std deviations of variance components.
(i.e., this is equal to sqrt(mse) in linear regression models, or the dispersion parameter in a negative binomial model.)
\item all : if 'all' is contained anywhere in comparisons argument, a list of all available results are returned.
}}

\item{variable}{vector of one or more of the fixed effect variables in DE model.  Can be "all".}

\item{targets}{vector of one or more of the targets analyzed in smide_results object. Can be "all".}
}
\value{
list of data.tables corresponding to targets, comparisons, and fixed effect variables from DE models.
}
\description{
Extract results from DE object by contrast comparison, model variable, gene
}
\examples{
library(Giotto)
library(data.table)
datadir<-system.file("extdata", package="smiDE")
gem <- readRDS(paste0(datadir, "/small_nsclc.rds"))
metainfo <- data.table::copy(gem@cell_metadata$rna)
metainfo <- merge(metainfo, gem@spatial_locs$raw, by="cell_ID")
pre_de_obj <- 
  pre_de(counts = gem@expression[["rna"]][["raw"]]
         ,normalized_data = gem@expression[["rna"]][["normalized"]]
         ,metadata = metainfo
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

fibroblast_cells <- gem@cell_metadata$rna[cell_type=="fibroblast",cell_ID]
de_results <-       
  smi_de(assay_matrix = gem@expression$rna$raw[,fibroblast_cells]
         ,metadata = metainfo[cell_ID \%in\% fibroblast_cells]
         ,formula = ~RankNorm(otherct_expr) + niche + tissue + offset(log(totalcounts)) 
         ,neighborhood_counts = pre_de_obj$nblist
         ,groupVar="niche"
         ,family="nbinom2"
         ,targets=rownames(gem@expression$rna$raw)[1:2]
  )

results(de_results, comparisons="pairwise",  variable="niche")[[1]]
results(de_results, comparisons="model_summary", targets="ABL1")[[1]]
results(de_results, comparisons="one.vs.rest", variable="tissue")[[1]] 
results(de_results, comparisons="emmeans", variable="otherct_expr")
results(de_results, comparisons="emmeans")
results(de_results)

}
