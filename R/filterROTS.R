# Filtering for ROTS
ROTS.filtered <- function(data, groups, B, K, paired=FALSE, progress=FALSE) {
  # List rows having less than two non-missing values
  data1 <- data[,groups==unique(groups)[1]]
  data2 <- data[,groups==unique(groups)[2]]
  filter1 <- which(rowSums(is.na(data1))>=ncol(data1)-1)
  filter2 <- which(rowSums(is.na(data2))>=ncol(data2)-1)
  # Store indexes and names
  filter <- sort(union(filter1,filter2))
  names <- rownames(data)[filter]
  if (length(filter)>0) {
   # Remove rows
    data <- data[-filter,]
    # Run ROTS
    rots.out <- ROTS(data=data, groups=groups, B=B, K=K, paired=paired, progress=progress)
    # Fill back filtered rows with NA
    for(i in 1:length(filter)) {
      rots.out$d <- append(rots.out$d, NA, after=filter[i]-1)
      names(rots.out$d)[filter[i]] <- names[i]
      rots.out$p <- append(rots.out$p, NA, after=filter[i]-1)
      rots.out$FDR <- append(rots.out$FDR, NA, after=filter[i]-1)
    }
    # Return results
    return(rots.out)
  } else {
    # Run ROTS
    rots.out <- ROTS(data=data, groups=groups, B=B, K=K, paired=paired, progress=progress)
    # Return results
    return(rots.out)
  }
}
