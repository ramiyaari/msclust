source('utils.R')
library('sigclust2')
library('mclust')

# running monte-carlo tests of the sigclust2 algorithm


if(!exists('runSeriesTests') || !runSeriesTests) {
  alpha <- 0.05
  M <- 40
  m <- 8
  n <- 40
  delta <- 0.05
  sigma <- 100
  mc <- 1000 
}

res <- generate_clusters(M,m)
actual_clusters <- res$clusters
actual_sep <- res$sep
source('run_sir.R')

obs_list1 <- generate_obs_incidence_using_sir_model(mc,M,m,actual_clusters,delta,n,noise_type='gaussian',sigma,seed=1)
obs_list2 <- generate_obs_incidence_using_sir_model(mc,M,m,actual_clusters,delta,n,noise_type='gaussian',sigma,seed=2)

results <- data.frame()
for(rr in 1:mc) {

  # if(rr%%200==0)
    # print(rr)

  obs <- t(obs_list1[[rr]]+obs_list2[[rr]])

  shc_result <- shc(obs, alpha=alpha, n_min=3)
  
  p_sig <- which(shc_result$p_emp<=alpha)
  clusters <- rep(0,M)
  cl_ind <- 0
  for(i in p_sig) {
    if(!any(diff(sort(shc_result$idx_hc[i,1][[1]]))>1))
      clusters[shc_result$idx_hc[i,1][[1]]] <- cl_ind+1
    if(!any(diff(sort(shc_result$idx_hc[i,1][[1]]))>1))
      clusters[shc_result$idx_hc[i,2][[1]]] <- cl_ind+2
    cl_ind <- cl_ind+2
  }
  
  ARI <- adjustedRandIndex(clusters,actual_clusters)
  sep <- which(clusters[1:(M-1)]!=clusters[2:M])
  sep <- paste(sep,collapse=',')

  row <- data.frame(sim=rr,sep=sep,ARI=ARI)
  results <- rbind(results,row)
}

test_summary <- analyze_clustering_results(results,M,actual_sep)
# print(test_summary)


