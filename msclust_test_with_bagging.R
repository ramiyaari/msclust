source('msclust.R')
source('utils.R')
library('mclust')

# running monte-carlo tests of the msclust+ algorithm


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

if(!exists('noise_type')) {
  obs_list1 <- generate_obs_incidence_using_sir_model(mc,M,m,actual_clusters,delta,n,noise_type='gaussian',sigma,seed=1)
  obs_list2 <- generate_obs_incidence_using_sir_model(mc,M,m,actual_clusters,delta,n,noise_type='gaussian',sigma,seed=2)
} else {
  obs_list1 <- generate_obs_incidence_using_sir_model(mc,M,m,actual_clusters,delta,n,noise_type=noise_type,sigma,seed=1)
  obs_list2 <- generate_obs_incidence_using_sir_model(mc,M,m,actual_clusters,delta,n,noise_type=noise_type,sigma,seed=2)
}
results <- data.frame()
for(rr in 1:mc) {
  
  # if(rr%%200==0)
    # print(rr)
  
  obs1 <- obs_list1[[rr]]
  obs2 <- obs_list2[[rr]]
  
  V <- c()
  for(k in 1:M) {
    obsk <- cbind(obs1[,k],obs2[,k])
    obs_k_mean <- rowMeans(obsk)
    V[k] <- mean(apply(obsk-obs_k_mean,1,var))
  }

  seps <- list()
  for(b in 1:B) {
    sel <- rbinom(M,1,0.5)
    sel1 <- which(sel==1)
    sel2 <- which(sel==0)
    obs1b <- matrix(0,n,M)
    obs2b <- matrix(0,n,M)
    obs1b[,sel1] <- obs1[,sel1]
    obs2b[,sel1] <- obs2[,sel1]
    obs1b[,sel2] <- obs2[,sel2]
    obs2b[,sel2] <- obs1[,sel2]
    res <- msclust(obs1b,obs2b,V,alpha)
    seps[[b]] <- res$sep 
  }
  res <- select_clusters_from_bagging_results(seps,M)
  ARI <- adjustedRandIndex(res$clusters,actual_clusters)
  sep <- paste(res$sep,collapse=',')
  row <- data.frame(sim=rr,sep=sep,ARI=ARI)
  results <- rbind(results,row)
}

test_summary <- analyze_clustering_results(results,M,actual_sep)
# print(test_summary)

