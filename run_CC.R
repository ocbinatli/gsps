library(cluster)
library(CancerSubtypes)

data_path <- "./data"
result_path <- "./results"
code_path <- "./code"

source(sprintf("%s/CC.R", code_path))

args <- commandArgs(trailingOnly = TRUE)
cluster_count <- as.numeric(args[[1]]) %% 10 + 1
pathway <- args[[2]]
cohort <- args[[3]]

if (file.exists(sprintf("%s/state_cohort_%s_pathway_%s_CC_%s.RData", result_path, cohort, pathway, cluster_count)) == FALSE) {
  print(sprintf("running %s cohort with %d clusters with method CC...", cohort, cluster_count))
  load(sprintf("%s/TCGA-%s.RData", data_path, cohort))
  load(sprintf("./msigdb/%s.RData", pathway))
  N_pathway <- length(pathways)
  
  X_train <- log2(TCGA$mrna[unique(rownames(TCGA$mrna)),] + 1)
  X_train <- scale(X_train)
  X_train <- X_train[,colMeans(is.na(X_train)) == 0]
  
  K_train <- list()
  for (m in 1:N_pathway) {
    hit_features <- which(colnames(X_train) %in% pathways[[m]]$symbols)
    K_train[[m]] <- t(X_train[, hit_features, drop = FALSE])
  }
  parameters <- list()
  parameters$cluster_count <- cluster_count
  state <- CC(K_train, parameters)
  save("state", file = sprintf("%s/state_cohort_%s_pathway_%s_CC_%s.RData", result_path, cohort, pathway, cluster_count))  
}
