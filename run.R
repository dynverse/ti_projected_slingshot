#!/usr/local/bin/Rscript

library(readr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(purrr, warn.conflicts = FALSE)
library(dyndimred, warn.conflicts = FALSE)
library(dynwrap, warn.conflicts = FALSE)
library(dyncli, warn.conflicts = FALSE)

library(princurve, warn.conflicts = FALSE)
library(cluster, warn.conflicts = FALSE)
suppressWarnings(library(slingshot, warn.conflicts = FALSE))

#####################################
###           LOAD DATA           ###
#####################################

task <- dyncli::main()

#' @example
#' task <- dyncli::main(
#'   args = "--dataset ~/example/test.loom --dimred landmark_mds --output ~/example/output.h5" %>% strsplit(" ") %>% first(),
#'   definition_location = "~/Workspace/dynverse/methods/ti_angle/definition.yml"
#' )

params <- task$params
counts <- task$counts
start_id <- task$priors$start_id
end_id <- task$priors$end_id

#   ____________________________________________________________________________
#   Preprocessing                                                           ####

start_cell <- if (!is.null(start_id)) { sample(start_id, 1) }  else { NULL }

# normalization & preprocessing
# from the vignette of slingshot
FQnorm <- function(counts){
  rk <- apply(counts, 2, rank, ties.method = "min")
  counts.sort <- apply(counts, 2, sort)
  refdist <- apply(counts.sort, 1, median)
  norm <- apply(rk, 2, function(r) refdist[r])
  rownames(norm) <- rownames(counts)
  return(norm)
}

expr <- t(log1p(FQnorm(t(counts))))

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

#   ____________________________________________________________________________
#   Dimensionality reduction                                                ####
pca <- prcomp(expr)

# this code is adapted from the expermclust() function in TSCAN
# the only difference is in how PCA is performed
# (they specify scale. = TRUE and we leave it as FALSE)
x <- 1:20
optpoint1 <- which.min(sapply(2:10, function(i) {
  x2 <- pmax(0, x - i)
  sum(lm(pca$sdev[1:20] ~ x + x2)$residuals^2 * rep(1:2,each = 10))
}))

# this is a simple method for finding the "elbow" of a curve, from
# https://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve
x <- cbind(1:20, pca$sdev[1:20])
line <- x[c(1, nrow(x)),]
proj <- project_to_curve(x, line)
optpoint2 <- which.max(proj$dist_ind)-1

# we will take more than 3 PCs only if both methods recommend it
optpoint <- max(c(min(c(optpoint1, optpoint2)), 3))
dimred <- pca$x[, seq_len(optpoint)]

#   ____________________________________________________________________________
#   Clustering                                                              ####
clusterings <- lapply(3:10, function(K){
  pam(dimred, K) # we generally prefer PAM as a more robust alternative to k-means
})

# take one more than the optimal number of clusters based on average silhouette width
# (max of 10; the extra cluster improves flexibility when learning the topology,
# silhouette width tends to pick too few clusters, otherwise)
wh.cl <- which.max(sapply(clusterings, function(x){ x$silinfo$avg.width })) + 1
grouping <- clusterings[[min(c(wh.cl, 8))]]$clustering %>% {set_names(paste0("M", .), names(.))}

start.clus <-
  if(!is.null(start_cell)) {
    grouping[[start_cell]]
  } else {
    NULL
  }
end.clus <-
  if(!is.null(end_id)) {
    unique(grouping[end_id])
  } else {
    NULL
  }

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####
sds <- slingshot(
  dimred,
  grouping,
  start.clus = start.clus,
  end.clus = end.clus,
  shrink = params$shrink,
  reweight = params$reweight,
  reassign = params$reassign,
  thresh = params$thresh,
  maxit = params$maxit,
  stretch = params$stretch,
  smoother = params$smoother,
  shrink.method = params$shrink.method
)

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

#   ____________________________________________________________________________
#   Create output                                                           ####

# extract information on clusters
lineages <- slingshot::slingLineages(sds)
lineage_ctrl <- slingshot::slingParams(sds)

# collect milestone network
milestone_network <- lineages %>%
  map_df(~ data_frame(from = .[-length(.)], to = .[-1])) %>%
  unique() %>%
  mutate(
    length = lineage_ctrl$dist[cbind(from, to)],
    directed = TRUE
  )

# calculate cluster dimred_milestones
milestone_ids <- unique(grouping)
dimred_milestones <- t(sapply(milestone_ids, function(cli){
  colMeans(dimred[names(which(cli == grouping)), , drop = FALSE])
}))

#   ____________________________________________________________________________
#   Save output                                                             ####

output <-
  wrap_data(
    cell_ids = rownames(counts)
  ) %>%
  dynwrap::add_dimred_projection(
    milestone_network = milestone_network,
    dimred = dimred,
    dimred_milestones = dimred_milestones,
    grouping = grouping
  ) %>%
  add_timings(checkpoints)

dyncli::write_output(output, task$output)
