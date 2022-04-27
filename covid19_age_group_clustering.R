source('msclust.R')
source('utils.R')

# running the msclust+ algorithm on real incidence data of Covid-19 from Israel, 
# to obtain age-group clustering per epidemic wave and sector of the population. 

############################################################

# aggregate daily incidence to weekly incidence

daily_obs_to_weekly_obs <- function(obs) {
  d <- 7
  m <- ncol(obs)
  n <- nrow(obs)
  if(n%%d<=3) 
    obs <- obs[1:(n-n%%d),]
  else
    obs <- rbind(obs,matrix(0,d-n%%d,m))
  n <- nrow(obs)/d
  weekly_obs <- matrix(NA,n,m) 
  for(t in 1:n) {
    weekly_obs[t,] <- apply(as.matrix(obs[((t-1)*d+1):(t*d),]),2,sum)
  }
  return (weekly_obs)
}

############################################################

# plot results of clustering

plot_clustered_incidence <- function(M, sep_vals, wave, sector) {
  
  part <- as.numeric(unlist(strsplit(as.character(sep_vals),',')))
  part <- c(0,part,M)
  clusters <- list()
  cluster_labels <- c()
  for(i in 1:(length(part)-1)) {
    clusters[[i]] <- (part[i]):((part[i+1])-1)+1
    if(i<(length(part)-1))
      cluster_labels[i] <- paste0(part[i],'-',part[i+1]-1)
    else
      cluster_labels[i] <- paste0(part[i],'+')
  }
  
  clusters_num <- length(clusters)
  cols <- rainbow(clusters_num) 
  
  dt1 <- as.Date('2020-03-01') + waves[[wave]][1]-1
  dates <- dt1+seq(0,n*7,by=7)
  title <- paste0('wave: ',wave,', sector: ',sector)
  rates <- matrix(0,nrow=n,ncol=clusters_num)
  for(i in 1:clusters_num) {
    if(length(clusters[[i]])==1)
      rates[,i] <- 1e4*cases_per_age[,clusters[[i]]]/N[clusters[[i]]]
    else
      rates[,i] <- 1e4*rowSums(cases_per_age[,clusters[[i]]])/sum(N[clusters[[i]]])
  }
  matplot(rates,type='l',lty=1,lwd=3,col=cols,xaxt='n',xlab='',main=title,
          ylim=c(0,max(rates)), ylab='incidence per 10,000',cex.lab=1.2,cex.axis=1.2)
  axis(1, at=1:n, labels=format(dates[1:n],'%d/%m/%y'),las=2,cex.lab=1.2,cex.axis=1.2)
  grid(nx=0,ny=NULL)
  legend('topleft',legend=cluster_labels,lty=rep(1,clusters_num),lwd=rep(3,clusters_num),col=cols,cex=1)
}

############################################################

alpha <- 0.05

sectors <- list()
sectors[[1]] <- 'all'
sectors[[2]] <- 'general'
sectors[[3]] <- 'haredim'
sectors[[4]] <- 'arabs'


for(sec in 1:length(sectors)) {

  sector <- sectors[[sec]]
  if(sector=='all') {
    cases_per_age0 <- read.csv('detected_cases_upto_all.csv',header=F) 
    N <- read.csv("Israel population general 2020.csv",header=F)[,1]+
         read.csv("Israel population haredim 2020.csv",header=F)[,1]+
         read.csv("Israel population arabs 2020.csv",header=F)[,1]
  } else if(sector=='general') {
    cases_per_age0 <- read.csv('detected_cases_upto_general.csv',header=F)
    N <- read.csv("Israel population general 2020.csv",header=F)[,1]
  } else if(sector=='haredim') {
    cases_per_age0 <- read.csv('detected_cases_upto_haredim.csv',header=F)
    N <- read.csv("Israel population haredim 2020.csv",header=F)[,1]
  } else if(sector=='arabs') {
    cases_per_age0 <- read.csv('detected_cases_upto_arabs.csv',header=F) 
    N <- read.csv("Israel population arabs 2020.csv",header=F)[,1]
  } else {
    stop('unknown sector')
  }
  
  cases_per_age0 <- as.matrix(cases_per_age0)
  Tdat <- nrow(cases_per_age0)
  
  if(sector=='all') {
    #generate figure 1 in manuscript
    
    dt1 <- as.Date('2020-03-01') + 1:Tdat-1
    months1 <- c(1,32,62,93,123,154,185,215,246,276,307,338,366,397,427,458,488,519,550,580,611,641,672,703,731)
    x11(height=20,width=25); par(mar=c(8,6,4,4))
    matplot(1:nrow(cases_per_age0),rowSums(cases_per_age0),type='l',lwd=3,col='blue',
            xaxt='n',yaxt='n',xlab='',ylab='')
    axis(1, at=months1,labels=format(dt1[months1],'%B %Y'), las=2)
    axis(2, at=seq(0,80000,by=20000),labels=c('0','20,000','40,000','60,000','80,000'), las=2)
    title(ylab='confirmed Covid-19 cases',line=5)
    abline(v=93,col='black',lty=2)
    abline(v=246,col='black',lty=2)
    abline(v=458,col='black',lty=2)
    abline(v=641,col='black',lty=2)
  }

  
  #############################################################################
  
  waves <- list()
  waves[[1]] <- 1:92 
  waves[[2]] <- 93:245 
  waves[[3]] <- 246:457
  waves[[4]] <- 458:640
  waves[[5]] <- 641:Tdat
  
  #combine age 0 and age 1 incidence and population 
  cases_per_age0[,2] <- cases_per_age0[,2]+cases_per_age0[,1]
  N[2] <- N[2]+N[1]
  cases_per_age0 <- cases_per_age0[,-1]
  N <- N[-1]
  
  # defined basic age groups
  groups <- list()
  groups[[1]] <- 1:2
  groups[[2]] <- 3:5
  groups[[3]] <- 6:8
  groups[[4]] <- 9:11
  groups[[5]] <- 12:14
  groups[[6]] <- 15:17
  groups[[7]] <- 18:20
  groups[[8]] <- 21:24
  groups[[9]] <- 25:29
  groups[[10]] <- 30:34
  groups[[11]] <- 35:39
  groups[[12]] <- 40:44
  groups[[13]] <- 45:49
  groups[[14]] <- 50:54
  groups[[15]] <- 55:59
  groups[[16]] <- 60:64
  groups[[17]] <- 65:69
  groups[[18]] <- 70:74
  groups[[19]] <- 75:79
  groups[[20]] <- 80:85

  M0 <- length(groups)
  
  #############################################################################
  
  #select wave
  for(w in 1:length(waves)) {
    
    selected_wave <- w
    cases_per_age <- cases_per_age0[waves[[selected_wave]],]
    
    
    #accumulate daily cases to weekly cases
    cases_per_age <- daily_obs_to_weekly_obs(cases_per_age)
    
    M <- ncol(cases_per_age)
    n <- nrow(cases_per_age)
    
    #############################################################################
    
    #calculate variance
    V <- c()
    for(i in 1:M0) {
      obsi <- c()
      L <- length(groups[[i]])
      for(l in 1:L) {
        j <- groups[[i]][l] 
        obsi <- cbind(obsi,cases_per_age[,j]/N[j]) 
      }
      obs_i_mean <- rowMeans(obsi)
      V[i] <- mean(apply(obsi-obs_i_mean,1,var))
    }
    
    #run clustering algorithm with bagging
    set.seed(100)
    sep_vals <- list()
    for(rr in 1:1000)
    {
      j <- 0
      obs1 <- matrix(0,nrow=n,ncol=M0)
      obs2 <- matrix(0,nrow=n,ncol=M0)
      for(i in 1:M0) {
        L <- length(groups[[i]])
        s1 <- sample(1:L,ifelse(pracma::rand()<0.5,floor(L/2),ceiling(L/2)))
        s2 <- setdiff(1:L,s1)
        j1 <- j+s1
        j2 <- j+s2
        if(length(j1) > 1)
          obs1[,i] <- rowSums(cases_per_age[,j1])/sum(N[j1])
        else
          obs1[,i] <- (cases_per_age[,j1])/sum(N[j1])
        if(length(j2) > 1)
          obs2[,i] <- rowSums(cases_per_age[,j2])/sum(N[j2])
        else
          obs2[,i] <- (cases_per_age[,j2])/sum(N[j2])
        j <- j+L
      }

      res <- msclust(obs1,obs2,V,alpha)
      
      if(!is.null(res$sep)) {
        sep_vals[[rr]] <-  res$sep 
      }
      else {
        sep_vals[[rr]] <- list()
      }
    }
    
    #summarize results of 1000 runs with bagging
    
    sep_vals <- sapply(sep_vals, function(s) { 
      sep <- c()
      if(length(s)>0) {
        for(j in 1:length(s)) { 
          sep[j] <- groups[[s[j]+1]][1] 
        } 
      }
      return (sep)
    })
    res_final <- select_clusters_from_bagging_results(sep_vals,M)
    
    #generate plot for figure 4-5 in manuscript
    plot_clustered_incidence(M, res_final$sep, selected_wave, sector)
  }
}
