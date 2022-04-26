
# implementation of the msclust algorithm

msclust <- function(X1,X2,V,alpha)
{
  stopifnot(all(dim(X1)==dim(X2)))
  M0 <- ncol(X1) 
  stopifnot(length(V)==M0)

  parts <- msclustR(X1,X2,V,alpha,M0,0)
   
  parts$clusters <- rep(1,M0)
  if(is.null(parts$sep))
    parts$sep <- list()
  else {
    sep <- c(parts$sep,M0)
    for(i in 2:length(sep)) {  
      parts$clusters[(sep[i-1]+1):sep[i]] <- i
    }
  }
  return (parts)
}

msclustR <- function(X1,X2,V,alpha,M0,part_start) {
  
  n <- nrow(X1)
  M <- ncol(X1)
  
  if(M==1)
    return (list())
  
  q_vals <- c()
  for(k in 1:(M-1)) {
    q_vals[k] <- calcQ(X1,V,M,k)
  }
    
  k_sep <- which.max(q_vals)
  q.val <- calcQ(X2,V,M,k_sep) 
  p.val <- 1-pchisq(q.val,n)
  p.adj <- p.val*M0/M 
  H0_rej <- (p.adj <= alpha)
  
  # if(is.nan(p.adj)) {
  #   H0_rej <- F
  # }
  
  parts <- list()
  if(H0_rej==T) {
    parts1 <- list()
    parts2 <- list()
    edges <- c()
    part_mid <- part_start+k_sep
    part_end <- part_start+M-1
    if(part_end == M0-1)
      part_end <- paste0(part_end,"+")
    node <- paste0(part_start,'-',part_end)
    nodeL <- paste0(part_start,'-',part_mid-1)
    nodeR <- paste0(part_mid,'-',part_end)
    edges <- c(edges,node,nodeL,node,nodeR)
    if(M > 2) {
      idx1 <- 1:k_sep
      idx2 <- (k_sep+1):M
      X11 <- as.matrix(X1[,idx1])
      X21 <- as.matrix(X2[,idx1])
      X12 <- as.matrix(X1[,idx2])
      X22 <- as.matrix(X2[,idx2])
      V1 <- V[idx1]
      V2 <- V[idx2]
      if(length(idx1) > 1) {
        parts1 <- msclustR(X11,X21,V1,alpha,M0,part_start)
      }
      if(length(idx2) > 1) {
        parts2 <- msclustR(X12,X22,V2,alpha,M0,part_mid)
      }
    }
    parts$sep <- c(parts1$sep,k_sep,parts2$sep+k_sep)
    parts$p.adj <- c(pmax(parts1$p.adj,p.adj),p.adj,pmax(parts2$p.adj,p.adj))
    parts$q.val <- c(parts1$q.val,q.val,parts2$q.val)
    parts$nodes <- c(parts1$nodes,node,parts2$nodes)
    parts$edges <- c(edges,parts1$edges,parts2$edges)
  }
  return (parts)
}

calcQ <- function(X,V,M,k) {
  idx1 <- 1:k
  idx2 <- (k+1):M
  x1 <- rowSums(as.matrix(X[,idx1])) 
  x2 <- rowSums(as.matrix(X[,idx2])) 
  v1 <- sum(V[idx1])
  v2 <- sum(V[idx2])
  C <- sum(x2*x1)/sum(x1^2)
  d <- C*x1-x2
  v <- C^2*v1+v2
  q <- sum(d^2/v)
  if(is.nan(q)) {
    return (Inf)
  }
  return (q)
}

