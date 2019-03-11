CC <- function(Km, parameters) {
  state <- list()
  state$time <- system.time({
    result = ExecuteCC(clusterNum=parameters$cluster_count, d=Km, maxK=parameters$cluster_count, clusterAlg="hc", distance="pearson")
    state$clustering <- result$group
    state$distance <- result$distanceMatrix
    state$parameters <- parameters
  })
  return(state)  
}
