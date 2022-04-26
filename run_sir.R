library("simode")


# running SIR model to produce simulated age group incidence

# mc - number of monte-carlo simulations
# M - number of groups
# m - number of clusters
# clusters - the clusters
# delta - difference in R between adjacent clusters
# n - number of observed time points
# sigma - sd for gaussian noise

generate_obs_incidence_using_sir_model <- function(mc,M,m,clusters,delta,n,noise_type=c('gaussian','poisson','none'),sigma,seed=1) {
  
  gamma <- 1/4              #recovery rate
  N  <- rep(1e5,M)          #group sizes
  S0 <- rep(0.6,M)*N        #initial susceptibles
  I0 <- 1e-2*S0             #initial infected
  Tmax <- 40                #max simulation run time

  noise_type <- match.arg(noise_type)
  stopifnot(clusters[M]==m)
  
  gamma_name <- 'gamma'
  beta_names <- matrix(apply(expand.grid(1:M,1:M),1,function(row) paste0('beta',row[1],'_',row[2])),nrow=M)
  N_names <- paste0('N',1:M)
  S_names <- paste0('S',1:M)
  I_names <- paste0('I',1:M)
  vars <- c(S_names,I_names)
  x0 <- c(S0,I0)
  
  # -----------------------------------------------------------------------------
  # generate equations ----------------------------------------------------------
  
  S_eq <- list()
  I_eq <- list()
  for (i in 1:M) {
    lambda_str <- ''
    for (j in 1:M) {
      lambda_str <- paste0(lambda_str,beta_names[i,j],"*",I_names[j],"/",N_names[j])
      if(j<M)
        lambda_str <- paste0(lambda_str,"+")
    }
    lambda_str <- paste0("(",lambda_str,")")
    S_eq <- c(S_eq,paste0("-",S_names[i],"*",lambda_str))
    I_eq <- c(I_eq,paste0(S_names[i],"*",lambda_str,"-",gamma_name,"*",I_names[i]))
  }
  equations <- unlist(c(S_eq,I_eq))
  names(equations) <- vars
  
  # -----------------------------------------------------------------------------
  # set transmission matrix -----------------------------------------------------
  
  beta <- matrix(0,M,M) 
  diag(beta) <- rep(1,M) - delta*(clusters-1)
  # diag(beta) <- rep(1,M) - delta*((clusters-1)%%2)
  

  # -----------------------------------------------------------------------------
  # run deterministic SIR model  ------------------------------------------------
  
  names(N) <- N_names
  names(beta) <- beta_names
  names(gamma) <- gamma_name
  names(S0) <- S_names
  names(I0) <- I_names
  names(x0) <- vars
  
  time <- seq(1,Tmax,length.out=n+1) #so that we have n points after diff 
  model_out <- solve_ode(equations,c(beta,gamma,N),x0,time)
  S_det <- model_out[,S_names]
  I_det <- model_out[,I_names]
  
  # -----------------------------------------------------------------------------
  # generate obs incidence rates ------------------------------------------------
  
  set.seed(seed)
  obs_list <- list()
  for(rr in 1:mc) {
    obs <- matrix(0,nrow=n,ncol=M)
    for(i in 1:M) {
      inc <- -diff(S_det[,i])
      # inc <- I_det[-Tmax,i] #
      if(noise_type=='poisson')
        obs[,i] <- rpois(n,inc) 
      else if(noise_type=='gaussian')
        obs[,i] <- inc+rnorm(n,0,sigma) #inc+rnorm(n,0,N[i]*sigma) #round(inc+rnorm(n,0,N[i]*sigma)) 
      else
        obs[,i] <- inc #round(inc)
    }
    colnames(obs) <- I_names
    obs_list[[rr]] <- obs
  }

  return (obs_list)
}

# -----------------------------------------------------------------------------


# source('utils.R')
# M <- 40
# m <- 8
# delta <- 0.04
# sigma <- 100
# mc <- 1000
# n <- 40
# Tmax <- 40
# actual_clusters <- generate_clusters(M,m)$clusters
# 
# obs_list <- generate_obs_incidence_using_sir_model(mc,M,m,actual_clusters,delta,n,noise_type='gaussian',sigma,seed=1)
# # obs_list <- generate_obs_incidence_using_sir_model(mc,M,m,actual_clusters,delta,n,noise_type='poisson',seed=1)
# # obs_list <- generate_obs_incidence_using_sir_model(mc,M,m,actual_clusters,delta,n,noise_type='none',seed=1)
# 
# matplot(seq(1,Tmax,length.out=n),obs_list[[1]],type='l',lty=2,col=actual_clusters,xlab='time',ylab='incidence',main=bquote(delta ~ '=' ~ 0.04))
# matplot(seq(1,Tmax,length.out=n),obs_list[[1]],pch=actual_clusters,col=actual_clusters,add=T)
# legend('topright',paste0('cluster',1:m),col=1:m,pch=1:m)