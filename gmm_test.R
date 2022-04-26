source('utils.R')
library('mclust')

# running monte-carlo tests of the gmm algorithm


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
  res <- Mclust(obs,verbose=F)
  clusters <- apply(res$z,1,which.max)
  uclusters <- unique(clusters)
  for(i in uclusters) {
    x <- which(clusters==i)
    if(length(x) > 0) {
      clusters[x[1]:x[length(x)]] <- i
    }
  }
  ARI <- adjustedRandIndex(clusters,actual_clusters)
  sep <- which(clusters[1:(M-1)]!=clusters[2:M])
  sep <- paste(sep,collapse=',')
  
  row <- data.frame(sim=rr,sep=sep,ARI=ARI)
  results <- rbind(results,row)
}

test_summary <- analyze_clustering_results(results,M,actual_sep)
# print(test_summary)

