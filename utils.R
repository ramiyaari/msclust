
# utility functions

library('mclust')

generate_clusters <- function(M,m) {
  sep <- c()
  clusters <- c()
  for(i in 1:m) {
    cluster_size <- round((M-length(clusters))/(m-i+1))
    clusters <- c(clusters,rep(i,cluster_size))
    if(i < m) {
      if(i == 1)
        sep[i] <- cluster_size
      else
        sep[i] <- sep[i-1]+cluster_size
    }
  }
  return (list(clusters=clusters,sep=sep))
}


analyze_clustering_results <- function(results,M,actual_sep) {
  
  row.names(results) <- NULL
  partition_summary <- table(results$sep)
  partition_summary <- partition_summary[order(partition_summary,decreasing = T)]
  # print(partition_summary)
  
  partition_size_count <- rep(0,M)
  names(partition_size_count) <- 1:M
  for(i in 1:length(partition_summary)){
    part_name <- names(partition_summary[i])
    j <- length(unlist(strsplit(part_name,',')))
    partition_size_count[j+1] <- partition_size_count[j+1] + partition_summary[i]
  }

  mARI <- round(mean(results$ARI),3)
  bullseye <- 0
  contained <- 0
  not_contained <- 0
  seps <- list()
  for(i in 1:nrow(results)) {
    res_sep <- as.numeric(strsplit(as.character(results[i,]$sep),split=',')[[1]])
    seps[[i]] <- res_sep
    if(setequal(res_sep,actual_sep))
      bullseye <- bullseye+1
    else if(length(setdiff(actual_sep,res_sep))==0 && length(setdiff(res_sep,actual_sep))>0)
      contained <- contained+1
    else
      not_contained <- not_contained+1
  }
  
  return (list(summary=partition_summary,size_count=partition_size_count,
               mARI=mARI,bullseye=bullseye,contained=contained,not_contained=not_contained))
}

sep_to_clusters <- function(sep,M) {
  clusters <- rep(1,M)
  if(length(sep) > 0 && !is.na(sep[1])) { 
    if(sep[length(sep)]!=M)
      sep <- c(sep,M)
    for(i in 2:length(sep)) {
      clusters[(sep[i-1]+1):sep[i]] <- i
    }
  }
  return (clusters)
}


select_clusters_from_bagging_results <- function(seps,M,method=c('max_count','max_agreement'),mARI_thresh=0.01) {
  
  method <- match.arg(method)
  B <- length(seps)
  seps_str <- unlist(lapply(seps,function(x) paste(x,collapse=',')))
  names(seps) <- seps_str
  unique_seps <- unique(seps)
  names(unique_seps) <- unique(seps_str)
  num_unique_parts <- length(unique_seps)
  sum_part <- summary(as.factor(seps_str),maxsum=num_unique_parts)
  sum_part <- sum_part[order(sum_part,decreasing = T)]
  unique_seps <- unique_seps[names(sum_part)]
  
  if(method=='max_count') {
    max_count <- max(sum_part)
    max_part <- names(which(sum_part==max_count))
    selected_sep <- seps[max_part]
  } else {
    unique_clusters <- lapply(unique_seps, function(sep) sep_to_clusters(sep,M))
    U <- length(unique_seps)
    agree_mat <- matrix(apply(expand.grid(1:U,1:U),1,function(r) {
      ifelse(r[1]==r[2],1,adjustedRandIndex(unique_clusters[[r[1]]],unique_clusters[[r[2]]]))}),U,U)
    agree_means <- sapply(1:U,function(i) sum(agree_mat[i,]*sum_part)/B)
    names(agree_means) <- names(unique_seps)
    max_agree <- max(agree_means)
    selected_sep <- unique_seps[which(max_agree-agree_means<=mARI_thresh)]
  }
  
  if(length(selected_sep)==0) {
    warning('length(selected_sep)==0')
    max_count <- max(sum_part)
    max_part <- names(which(sum_part==max_count))
    selected_sep <- seps[max_part]
  }
  if(length(selected_sep)==1) {
    selected_sep <- selected_sep[[1]]
  } else {
    part_size <- unlist(lapply(selected_sep,length))+1
    selected_sep <- selected_sep[which(part_size==min(part_size))]
    if(length(selected_sep) > 1) {
      part_count <- c()
      part_name <- c()
      for(i in 1:length(selected_sep)) {
        part_name[i] <- paste(selected_sep[[i]],collapse=',')
        part_count[i] <- sum_part[which(names(sum_part)==part_name[i])]
      }
      names(part_count) <- part_name
      selected_sep <- selected_sep[which.max(part_count)]
      if(length(selected_sep) > 1) {
        warning(paste0('more than 1 max. partition: ',selected_sep),immediate.=T)
      }
    }
    selected_sep <- selected_sep[[1]]
  }
  
  res <- list()
  res$sum_part <- sum_part
  res$sep <- selected_sep
  res$clusters <- sep_to_clusters(selected_sep,M)
  if(method=='max_agreement') {
    res$max_agree <- max_agree
    res$agree_means <- sort(agree_means,decreasing=T)
  }
  return (res)
}