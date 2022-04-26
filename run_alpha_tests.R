
# running tests to verify type-I error - assuming Gaussian noise


runSeriesTests <- TRUE

M_vals <- 40 #c(20,40)
m_vals <- c(1,2,4,6,8)

alpha <- 0.05
n <- 40
delta <- 0.05
sigma <- 10
mc <- 1000

loadTestResults <- TRUE

if(loadTestResults) {
  
  alpha_test_results <- readRDS('alpha_test_results.RDS')
  print(alpha_test_results)
  
} else {

  alpha_test_results <- data.frame(delta=NA,M=NA,m=NA,msclust_known_sig=NA,msclust_unknown_sig=NA,gmm=NA,sigclust2=NA)
  
  test <- 1
  for(M in M_vals) {
    for(m in m_vals) {
  
      alpha_test_results[test,1:3] <- c(delta,M,m)
  
      fixedSigma <- TRUE
      source('msclust_test.R')
      alpha_test_results[test,c(4)] <- 1-test_summary$bullseye/mc
  
      fixedSigma <- FALSE
      source('msclust_test.R')
      alpha_test_results[test,c(5)] <- 1-test_summary$bullseye/mc
  
      source('gmm_test.R')
      alpha_test_results[test,c(6)] <- 1-test_summary$bullseye/mc
  
      source('sigclust2_test.R')
      alpha_test_results[test,c(7)] <- 1-test_summary$bullseye/mc
      
      print(alpha_test_results[test,])
      
      test <- test+1
    }
  }
  alpha_test_results$delta <- delta
  saveRDS(alpha_test_results,'alpha_test_results.RDS')
}


########################################################################

# running tests to verify type-I error - assuming Poisson noise

if(loadTestResults) {
  
  alpha_test_results2 <- readRDS('alpha_test_results_poisson.RDS')
  print(alpha_test_results2)
  
} else {
  
  noise_type <- 'poisson'
  
  alpha_test_results2 <- data.frame(delta=NA,M=NA,m=NA,msclust=NA,msclust100=NA)
  
  test <- 1
  for(M in M_vals) {
    for(m in m_vals) {
      
      alpha_test_results2[test,1:3] <- c(delta,M,m)
      
      fixedSigma <- FALSE
      source('msclust_test.R')
      alpha_test_results2[test,c(4)] <- 1-test_summary$bullseye/mc
      
      B <- 100
      source('msclust_test_with_bagging.R')
      alpha_test_results2[test,c(5)] <- 1-test_summary$bullseye/mc
      
      
      test <- test+1
    }
  }
  
  alpha_test_results2$delta <- delta
  saveRDS(alpha_test_results2,'alpha_test_results_poisson.RDS')
}


