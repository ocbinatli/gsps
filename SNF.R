SNF <- function(Km, parameters) {
  state <- list()
  state$time <- system.time({
    result = ExecuteSNF(datasets=Km, clusterNum=parameters$cluster_count)
    state$clustering <- result$group
    state$distance <- result$distanceMatrix
    state$parameters <- parameters
  })
  return(state)  
}
