
# running tests to establish power of the algorithms - assuming Gaussian noise


runSeriesTests <- TRUE

M_vals <- 40 #c(20,40)
m_vals <- c(2,4,6,8)
delta_vals <- seq(0.01,0.04,by=0.005)

alpha <- 0.05
n <- 40
sigma <- 100 
mc <- 1000

loadTestResults <- TRUE

if(loadTestResults) {
  power_test_results <- readRDS('power_test_results.RDS')
  print(power_test_results)
  
} else {
  
  power_test_results <- data.frame(delta=NA,M=NA,m=NA,
                                   power_msclust=NA,power_msclust100=NA,power_gmm=NA,power_sigclust2=NA,
                                   mARI_msclust=NA,mARI_msclust100=NA,mARI_gmm=NA,mARI_sigclust2=NA)
  
  test <- 1
  for(M in M_vals) {
    for(m in m_vals) {
      for(delta in delta_vals) {
        
        power_test_results[test,1:3] <- c(delta,M,m)
        
        source('msclust_test.R')
        power_test_results[test,c(4,8)] <- c(test_summary$bullseye/mc,test_summary$mARI)
        
        B <- 100
        source('msclust_test_with_bagging.R')
        power_test_results[test,c(5,9)] <- c(test_summary$bullseye/mc,test_summary$mARI)
        
        source('gmm_test.R')
        power_test_results[test,c(6,11)] <- c(test_summary$bullseye/mc,test_summary$mARI)
        
        source('sigclust2_test.R')
        power_test_results[test,c(7,12)] <- c(test_summary$bullseye/mc,test_summary$mARI)
        
        # print(power_test_results[test,])
        
        test <- test+1
      }
    }
  }
  saveRDS(power_test_results,'power_test_results.RDS')
  
} 


################################################################################

# running tests to establish power of the algorithms - assuming Poisson noise

if(loadTestResults) {
  power_test_results2 <- readRDS('power_test_results_poisson.RDS')
  print(power_test_results2)
  
} else {
  
  noise_type <- 'poisson'
  power_test_results2 <- data.frame(delta=NA,M=NA,m=NA,
                                    power_msclust=NA,power_msclust100=NA,
                                    mARI_msclust=NA,mARI_msclust100=NA)
  
  test <- 1
  for(M in M_vals) {
    for(m in m_vals) {
      for(delta in delta_vals) {
        
        power_test_results2[test,1:3] <- c(delta,M,m)
        
        source('msclust_test.R')
        power_test_results2[test,c(4,6)] <- c(test_summary$bullseye/mc,test_summary$mARI)
        
        B <- 100
        source('msclust_test_with_bagging.R')
        power_test_results2[test,c(5,7)] <- c(test_summary$bullseye/mc,test_summary$mARI)
        
        print(power_test_results2[test,])
        
        test <- test+1
      }
    }
  }
  saveRDS(power_test_results2,'power_test_results_poisson.RDS')
}


################################################################################

# code to compare effect of B on the power of the algorithm - generating figure 3 in the manuscript

if(loadTestResults) {
  power_test_results3 <- readRDS('power_test_results_effect_of_B.RDS')
} else {
  
  power_test_results3 <- power_test_results[22:28,1:5]
  power_test_results3$power_msclust10 <- 0
  M <- 40
  m <- 8
  for(test in 1:length(delta_vals)) {
    delta <- delta_vals[test]
    B <- 10
    source('msclust_test_with_bagging.R')
    power_test_results3$power_msclust10[test] <- test_summary$bullseye/mc
  }
  saveRDS(power_test_results3,'power_test_results_effect_of_B.RDS')
}

x11();
matplot( power_test_results3$delta,power_test_results3[,c(4,6,5)],
        type='l',lty=1,lwd=3,col=2:4,xlab=expression(delta),ylab='power')
matplot( power_test_results3$delta,power_test_results3[,c(4,6,5)],type='p',pch=0:2,col=2:4,add=T,cex=1.5)
legend('topleft',c('B=1','B=10','B=100'),lty=c(1,1,1),lwd=c(3,3,3),pch=0:2,col=2:4)


#######################################################################

# code to show detailed results of power simulations - generating figure 2 in the manuscript

M <- 40
m <- 8 
delta_vals <- seq(0.01,0.04,by=0.01)

if(loadTestResults) {
  power_test_results4 <- readRDS('power_test_results_detailed.RDS')
  
} else {
  
  power_test_results4 <- list()
  test <- 1
  for(test in 1:length(delta_vals)) {
    
    delta <- delta_vals[test]
    power_test_results4[[test]] <- list()
    
    source('msclust_test.R')
    power_test_results4[[test]][[1]] <- test_summary
    
    B <- 100
    source('msclust_test_with_bagging.R')
    power_test_results4[[test]][[2]] <- test_summary
    
    source('gmm_test.R')
    power_test_results4[[test]][[3]] <- test_summary
    
    source('sigclust2_test.R')
    power_test_results4[[test]][[4]] <- test_summary
    
  }
  saveRDS(power_test_results4,'power_test_results_detailed.RDS')
}

source('utils.R')
source('run_sir.R')
res <- generate_clusters(M,m)
actual_clusters <- res$clusters
actual_sep <- res$sep

delta <- 0.01
obs_list1 <- generate_obs_incidence_using_sir_model(mc,M,m,actual_clusters,delta,n,noise_type='gaussian',sigma,seed=1)
delta <- 0.04
obs_list4 <- generate_obs_incidence_using_sir_model(mc,M,m,actual_clusters,delta,n,noise_type='gaussian',sigma,seed=1)

Tmax <- 40

x11(width=25,heigh=16); par(mfrow=c(2,2),mar=c(5,5,3,3),mgp=c(3,0.5,0))
matplot(seq(1,Tmax,length.out=n),obs_list1[[1]],type='l',lty=2,col=actual_clusters,
        xlab='time',ylab='incidence',main=bquote(delta ~ '=' ~ 0.01),cex.axis=1.2,cex.lab=1.2,yaxt='n',ylim=c(0,4200))
matplot(seq(1,Tmax,length.out=n),obs_list1[[1]],pch=actual_clusters,col=actual_clusters,add=T)
axis(2,at=seq(0,4000,1000),labels=seq(0,4000,1000),las=2,cex.axis=1.2)
legend('topright',paste0('cluster',1:m),col=1:m,pch=1:m,cex=1.2)
text(1,4000,'A',cex=2)
matplot(seq(1,Tmax,length.out=n),obs_list4[[1]],type='l',lty=2,col=actual_clusters,
        xlab='time',ylab='incidence',main=bquote(delta ~ '=' ~ 0.04),cex.axis=1.2,cex.lab=1.2,yaxt='n',ylim=c(0,4200))
matplot(seq(1,Tmax,length.out=n),obs_list4[[1]],pch=actual_clusters,col=actual_clusters,add=T)
axis(2,at=seq(0,4000,1000),labels=seq(0,4000,1000),las=2,cex.axis=1.2)
legend('topright',paste0('cluster',1:m),col=1:m,pch=1:m,cex=1.2)
text(1,4000,'B',cex=2)

cols1 <- c('red3','blue3','green3','yellow2')
barplot(t(sapply(1:4, function(i) power_test_results4[[1]][[i]]$size_count[1:16])),
        xlab='number of clusters',ylab='frequency',main=bquote(delta ~ '=' ~ 0.01),beside=T,col=cols1,ylim=c(0,1000),cex.axis=1.2,cex.lab=1.2,yaxt='n')
axis(2,at=seq(0,1000,200),labels=seq(0,1000,200),las=2,cex.axis=1.2)
legend('topright',c('msclust','msclust+','gmm','sigclust2'),fill=cols1,cex=1.2)
text(1,900,'C',cex=2)
box()
barplot(t(sapply(1:4, function(i) power_test_results4[[4]][[i]]$size_count[1:16])),
        xlab='number of clusters',ylab='frequency',main=bquote(delta ~ '=' ~ 0.04),beside=T,col=cols1,ylim=c(0,1000),cex.axis=1.2,cex.lab=1.2,yaxt='n')
axis(2,at=seq(0,1000,200),labels=seq(0,1000,200),las=2,cex.axis=1.2)
legend('topright',c('msclust','msclust+','gmm','sigclust2'),fill=cols1,cex=1.2)
text(1,900,'D',cex=2)
box()
