#!/usr/local/bin/Rscript

requireNamespace("dyncli", quietly = TRUE)
task <- dyncli::main()

library(dplyr, warn.conflicts = FALSE)
library(purrr, warn.conflicts = FALSE)
library(dynwrap, warn.conflicts = FALSE)

requireNamespace("princurve", quietly = TRUE)
requireNamespace("cluster", quietly = TRUE)
requireNamespace("irlba", quietly = TRUE)

suppressWarnings(library(slingshot, warn.conflicts = FALSE))

#####################################
###           LOAD DATA           ###
#####################################

parameters <- task$parameters
expression <- task$expression
start_id <- task$priors$start_id
end_id <- task$priors$end_id

#####################################
###        INFER TRAJECTORY       ###
#####################################
#   ____________________________________________________________________________
#   Preprocessing                                                           ####

start_cell <- if (!is.null(start_id)) { sample(start_id, 1) }  else { NULL }

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

#   ____________________________________________________________________________
#   Dimensionality reduction                                                ####
pca <- irlba::prcomp_irlba(expression, n = 20)

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
  shrink = parameters$shrink,
  reweight = parameters$reweight,
  reassign = parameters$reassign,
  thresh = parameters$thresh,
  maxit = parameters$maxit,
  stretch = parameters$stretch,
  smoother = parameters$smoother,
  shrink.method = parameters$shrink.method
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
    cell_ids = rownames(expression)
  ) %>%
  dynwrap::add_dimred_projection(
    milestone_network = milestone_network,
    dimred = dimred,
    dimred_milestones = dimred_milestones,
    grouping = grouping
  ) %>%
  add_timings(checkpoints)

dyncli::write_output(output, task$output)
