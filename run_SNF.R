library(cluster)
library(CancerSubtypes)

data_path <- "./data"
result_path <- "./results"
code_path <- "./code"

source(sprintf("%s/SNF.R", code_path))

pdist <- function(X1, X2) {
  if (identical(X1, X2) == TRUE) {
    D <- as.matrix(dist(X1))
  }
  else {
    D <- as.matrix(dist(rbind(X1, X2)))
    D <- D[1:nrow(X1), (nrow(X1) + 1):(nrow(X1) + nrow(X2))]
  }
  return(D)
}

args <- commandArgs(trailingOnly = TRUE)
cluster_count <- as.numeric(args[[1]]) %% 10 + 1
pathway <- args[[2]]
cohort <- args[[3]]

if (file.exists(sprintf("%s/state_cohort_%s_pathway_%s_SNF_%s.RData", result_path, cohort, pathway, cluster_count)) == FALSE) {
  print(sprintf("running %s cohort with %d clusters ...", cohort, cluster_count))
  load(sprintf("%s/TCGA-%s.RData", data_path, cohort))
  load(sprintf("./msigdb/%s.RData", pathway))
  N_pathway <- length(pathways)
  
  X_train <- log2(TCGA$mrna[unique(rownames(TCGA$mrna)),] + 1)
  X_train <- scale(X_train)
  X_train <- X_train[,colMeans(is.na(X_train)) == 0]
  
  pathway_names <- sapply(1:length(pathways), FUN = function(m) pathways[[m]]$name)
  N_train <- nrow(X_train)
  
  K_train <- list()
  for (m in 1:N_pathway) {
    hit_features <- which(colnames(X_train) %in% pathways[[m]]$symbols)
    D_train <- pdist(X_train[, hit_features, drop = FALSE], X_train[, hit_features, drop = FALSE])
    sigma <- mean(D_train)
    K_train[[m]] <- exp(-D_train^2 / (2 * sigma ^2))
  }
  
  parameters <- list()
  parameters$cluster_count <- cluster_count
  state <- SNF(K_train, parameters)
  save("state", file = sprintf("%s/state_cohort_%s_pathway_%s_SNF_%s.RData", result_path, cohort, pathway, cluster_count))
}
