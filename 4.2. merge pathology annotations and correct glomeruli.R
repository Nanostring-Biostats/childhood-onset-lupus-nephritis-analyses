# inputs: 
# - glom annotations from 4.1. annotate glomeruli ('gannot' - includes gloms not annotated by pathologist)
# - pathology annotations

# outputs:
# - updated cell->glomeruli assignments ("db")
# - merged glomerui annotations ("gloms")



library(readxl)

# load data:
load("processed_data/cleaneddata.RData")

# load original glom annot:
gannot = read.csv("processed_data/glomeruli annotations.csv", header = T, row.names = 1)
db = readRDS("processed_data/db.rds")

# load path annots:
pathannot = read_xlsx("data/Scoring glomeruli.xlsx", sheet = gsub("_", "-", unique(gannot$tissue)[1]))
pathannot = as.data.frame(pathannot)
for (tiss in unique(gannot$tissue)[-1]) {
  if (annot$class[annot$tissue == tiss][1] != "control") {
    temp = read_xlsx("data/Scoring glomeruli.xlsx", sheet = gsub("_", "-", tiss))
    if (identical(colnames(pathannot), colnames(temp))) {
      pathannot = rbind(pathannot, temp)
    } else {
      stop(paste0(tiss, " misalignment"))
    }
  }
}
head(pathannot)

# correct IDs with errors introduced by excel:
pathannot$`annotated #` = gsub("19999999999999", "2", pathannot$`annotated #`)
pathannot$`annotated #` = gsub("10000000000002", "1", pathannot$`annotated #`)

# extract more variables:
pathannot$double.contours = grepl("double contours", pathannot$other) * 1
pathannot$protein.droplets = grepl("protein droplets", pathannot$other) * 1
pathannot$podocyte.hyperplasia = grepl("podocyte hyperplasia", pathannot$other) * 1
pathannot$segmental.sclerosis = grepl("segmental sclerosis", pathannot$other) * 1


sort(setdiff(pathannot$`annotated #`, as.character(gannot$old.id)))
sort(setdiff(as.character(gannot$old.id), pathannot$`annotated #`))


#### merge gannot with pathannot:
sharedids = intersect(gannot$old.id, pathannot$`annotated #`)
print(sharedids)
print(setdiff(gannot$old.id, pathannot$`annotated #`))
print(setdiff(pathannot$`annotated #`, gannot$old.id))

# build "gloms" data frame of path-annotated gloms:
gloms = cbind(gannot[match(sharedids, gannot$old.id), ], pathannot[match(sharedids, pathannot$`annotated #`), ])
gloms$originalid = gloms$id
gloms$id = gloms$new.id
rownames(gloms) = gloms$id
gloms$pathology.id = gloms$glom
gloms$glom = NULL

# add in gloms missing from path:
glomnames = setdiff(unique(db), NA)
setdiff(glomnames, gloms$new.id)
setdiff(gloms$new.id, glomnames)
gloms = cbind(glomnames, gloms[match(glomnames, gloms$new.id), ])
rownames(gloms) = glomnames
gloms$new.id = glomnames

# map to clin variables:
#gloms$tissue = gannot[rownames(gloms), "tissue"]
gloms$tissue = NA
for (i in 1:nrow(gloms)) {
  onecellinglom = names(which(db == gloms$new.id[i]))[1]
  gloms$tissue[i] = annot$tissuename[match(onecellinglom, annot$cell_ID)]
}
gloms$disease = annot$class[match(gloms$tissue, annot$tissue)]
gloms$x = gloms$y = NA
for (i in 1:nrow(gloms)) {
  tempcells = names(db)[(db == gloms$new.id[i]) & !is.na(db)]
  gloms$x[i] = mean(customlocs[match(tempcells, annot$cell_ID), 1])
  gloms$y[i] = mean(customlocs[match(tempcells, annot$cell_ID), 2])
}

# look at available variables:
table(gloms$tissue, gloms$mesang)  # dichotomize?
table(gloms$tissue, gloms$endocapillary) # dichotomize?
table(gloms$tissue, gloms$double.contours) # only in one sample, where it's 50% of gloms
table(gloms$tissue, gloms$protein.droplets) # very rare
table(gloms$tissue, gloms$podocyte.hyperplasia) # too few
table(gloms$tissue, gloms$segmental.sclerosis) # only 11 gloms in 1 sample.
table(gloms$tissue, gloms$`wire-loop`) # probably too few
table(gloms$tissue, gloms$`circulating inflamm cells`)  # too few to use


# define new glom variables: (added here 6-15, moved from the DE script. )
#gloms$hasTcells = 1 * (gloms$T.tot > 0)
#gloms$hasBcells = 1 * (gloms$B.tot > 0)
#gloms$hasplasmacells = 1 * (gloms$plasma.tot > 0)
gloms$wire10 = 1 * (gloms$`wire-loop` > 0)
gloms$endocapillary10 = 1 * (gloms$endocapillary > 0)
gloms$mesang10 = 1 * (gloms$mesang != "0")

#gloms$hypercellular = 1 * ((gloms$mesang > 0) | (gloms$endocapillary > 0))
gloms$involved = 1 * (((gloms$`wire-loop` > 0) | (gloms$double.contours > 0)) | (gloms$endocapillary > 0)) 

rownames(gloms) = gloms$new.id
saveRDS(gloms, file = "processed_data/gloms.rds")
#saveRDS(db, file = "processed_data/db - updated glom membership.rds")
