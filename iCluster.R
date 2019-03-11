iCluster <- function(Km, parameters) {
  state <- list()
  state$time <- system.time({
    result = ExecuteiCluster(datasets=Km, k=parameters$cluster_count, max.iter = 50)
    state$clustering <- result$group
    state$parameters <- parameters
  })
  return(state)  
}
