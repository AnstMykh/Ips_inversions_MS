install.packages("devtools")
library(devtools)
devtools::install_github("petrelharp/local_pca/lostruct")
library(lostruct)

install.packages("ggplot2")
library(ggplot2)
install.packages("tidyverse")
library(tidyverse)
install.packages("gridExtra")
library(gridExtra)
install.packages("grid")
library(grid)


setwd("path/to/your/data")



k_kept <- 40
data_name <- "whole_genome.bcf"
window_size <- 1e4 
min_windows <- 4
n_permutations <- 1000
max_distance_between_outliers <- 20
min_cor <- 0.8


sites <- vcf_positions(data_name)
win.fn.snp <- vcf_windower(data_name, size=window_size, type='bp', sites=sites) # opening .bcf file



system.time( snp.pca <- eigen_windows(win.fn.snp ,k=2, mc.cores=10) )
system.time( pcdist <- pc_dist( snp.pca ) )


#pcdist_na <- which(is.na(pcdist), TRUE)
col_num <- 1
na.inds <- is.na(pcdist[,col_num])
while (sum(na.inds) == length(na.inds)){
  col_num <- col_num+1
  na.inds <- is.na(pcdist[,col_num])
}

mds <- cmdscale( pcdist[!na.inds,!na.inds], eig=TRUE, k=k_kept )

mds.coords <- mds$points
colnames(mds.coords) <- paste("MDS coordinate", 1:ncol(mds.coords))
win.regions <- region(win.fn.snp)()
win.regions$n <- 1:nrow(win.regions)
win.regions <- win.regions[!na.inds,]
win.regions %>% mutate(mid = (start + end) / 2) ->  win.regions


for (k in 1:k_kept){
  str_pad(k, 2, pad = "0")
  
  name = paste("mds",str_pad(k, 2, pad = "0"),sep="")
  win.regions$tmp <- "NA"
  win.regions <- win.regions %>% rename(!!name := tmp)
}

#Add the MDS coordinates to each window.
for (i in 1:k_kept){
  j = i + 5
  win.regions[,j] <- mds.coords[,i]
}

write.csv(win.regions , 'path/to/intermediary/table.txt', row.names = TRUE)





# Look for "chromosomally clustered" windows along each MDS coordinate (in both positive and negetive directions)
win.regions <-  read.table('path/to/intermediary/table.txt', header=TRUE, skip=0, sep=",") [,2:46]



names_df <- colnames(win.regions)[6:45]
win.regions <- select(win.regions, -X, -mds01)
colnames(win.regions)[5:44] <- names_df
win.regions <- relocate(win.regions, mid, .after = 4)    # you can maybe use it while calling mutate() 


mds_pcs <- colnames(win.regions)[6:(ncol(win.regions))]

mds_clustering <- tibble(mds_coord=character(), direction=character(), clust_pvalue=numeric(), outliers=numeric(), n1_outliers=numeric(), 
                         high_cutoff=numeric(), lower_cutoff=numeric(), chr=character())


for (i in mds_pcs){
  print(paste("Processing",i))
  win.regions %>%
  mutate_(the_mds=i) %>%    
  summarize(sd_mds=sd(the_mds)) %>% pull() -> sd_mds
  mds_high_cutoff <- sd_mds*2
  mds_low_cutoff <- sd_mds
  
  win.regions %>%
    mutate_(the_mds=i) %>%   
    mutate(sd_mds=sd(the_mds)) %>%
    filter(the_mds > sd_mds) -> pos_windows
  if (nrow(pos_windows) >= min_windows){
    permutations <- matrix(nrow=n_permutations, ncol=1)
    for (i in 1:n_permutations){
      max_1 <- win.regions %>% sample_n(nrow(pos_windows)) %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=1) %>% 
        ungroup() %>% summarize(sum=sum(count)) %>% pull(sum)
      permutations[i,1] <- max_1
    }
    pos_windows %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=2) %>% 
      ungroup() %>% summarize(sum=sum(count)) %>% pull(sum) -> sampled_max_1
    pos_windows %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=2) %>%
      pull(chrom) -> clustered_chr
    x <- as.tibble(permutations) %>% filter(V1 >= sampled_max_1) %>% nrow() 
    pvalue <- (x+1)/(n_permutations+1) #this is how p-value is calculated 
    tmp <- tibble(mds_coord=as.character(i), direction=as.character("pos"), clust_pvalue=as.numeric(pvalue), outliers=as.numeric(nrow(pos_windows)),  
                  n1_outliers=as.numeric(sampled_max_1), high_cutoff=as.numeric(mds_high_cutoff), lower_cutoff=as.numeric(mds_low_cutoff), chr=as.character(clustered_chr))
    mds_clustering <- rbind(mds_clustering, tmp)
  }else{
    tmp <- tibble(mds_coord=as.character(i), direction=as.character("pos"), clust_pvalue=as.numeric(NA), outliers=as.numeric(nrow(pos_windows)),    
                  n1_outliers=as.numeric(NA), high_cutoff=as.numeric(NA), lower_cutoff=as.numeric(NA), chr=as.numeric(NA))
    mds_clustering <- rbind(mds_clustering, tmp)
  }
  
  win.regions %>%
    mutate_(the_mds=i) %>%     
    mutate(sd_mds=sd(the_mds)) %>%
    filter(the_mds < -(sd_mds)) -> neg_windows
  if (nrow(neg_windows) >= min_windows){
    permutations <- matrix( nrow=n_permutations, ncol=1)
    for (i in 1:n_permutations){
      max_1 <- win.regions %>% sample_n(nrow(neg_windows)) %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=1) %>% 
        ungroup() %>% summarize(sum=sum(count)) %>% pull(sum)
      permutations[i,1] <- max_1
    }
    neg_windows %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=2) %>% 
      ungroup() %>% summarize(sum=sum(count)) %>% pull(sum) -> sampled_max_1
    neg_windows %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=2) %>%
      pull(chrom) -> clustered_chr
    x <- as.tibble(permutations) %>%  filter(V1 >= sampled_max_1) %>% nrow()
    pvalue <- (x+1)/(n_permutations+1)
    tmp <- tibble(mds_coord = as.character(i),direction = as.character("neg"),clust_pvalue = as.numeric(pvalue), outliers=as.numeric(nrow(neg_windows)), 
                  n1_outliers=as.numeric(sampled_max_1), high_cutoff = as.numeric(-mds_high_cutoff),lower_cutoff=as.numeric(-mds_low_cutoff),
                  chr=as.character(clustered_chr))
    mds_clustering <- rbind(mds_clustering, tmp)
  }else{
    tmp <- tibble(mds_coord = as.character(i),direction = as.character("neg"),clust_pvalue = as.numeric(NA), outliers=as.numeric(nrow(neg_windows)),     
                  n1_outliers=as.numeric(NA), high_cutoff = as.numeric(NA),lower_cutoff=as.numeric(NA),
                  chr=as.numeric(NA))
    mds_clustering <- rbind(mds_clustering, tmp)
  }
}


mds_clustering %>% filter(clust_pvalue < 0.01) -> sig_mds_clusters


write.csv(sig_mds_clusters, 'path/to/sig_mds_clusters_25_e4.csv', row.names = TRUE)   #### This file is used in the downstream analysis and visualisation of data 


# Generate datasets of outlier windows and genotypes
outlier_windows <- tibble(chrom=character(), start=numeric(), end=numeric(), mid=numeric(), 
                          the_mds=numeric(), mds_coord=character(), outlier=character(), 
                          n=numeric())
cluster_genotypes <- tibble(mds_coord=character(), name=character(), PC1=numeric(), genotype=character())
for (i in 1:nrow(sig_mds_clusters)){
  coord <- pull(sig_mds_clusters[i,1])    # i
  direction <- pull(sig_mds_clusters[i,2])  # i
  high_cutoff <- pull(sig_mds_clusters[i,6])  # i
  low_cutoff <- pull(sig_mds_clusters[i,7])  # i
  cluster_chr <- pull(sig_mds_clusters[i,8])  # i
  coord_direction <- paste(coord, "-", direction, sep="") #i 
  print(paste("Testing",coord_direction))
  
  if (direction == "pos"){
    current_windows <- win.regions %>%
      mutate_(the_mds=coord ) %>% 
      mutate(outlier=case_when((the_mds > high_cutoff) & (chrom == cluster_chr) ~ "Outlier", 
                               TRUE ~ "Non-outlier")) %>%
      filter(outlier != "Non-outlier") %>%
      select(chrom,start,end,mid,the_mds,outlier,n) %>%
      mutate(mds_coord=coord_direction) %>%
      mutate(ahead_n=n-lag(n), behind_n=abs(n-lead(n))) %>%
      mutate(min_dist=pmin(ahead_n, behind_n, na.rm=T)) %>%
      filter(min_dist < max_distance_between_outliers) %>%
      select(-ahead_n,-behind_n,-min_dist)
    windows <- current_windows %>% pull(n)
    outlier_windows <- rbind(outlier_windows, current_windows)
    if(nrow(current_windows) < 1){
      next;
    }
  }else{
    current_windows <- win.regions %>%
      mutate_(the_mds=coord) %>% 
      mutate(outlier=case_when((the_mds < high_cutoff) & (chrom == cluster_chr) ~ "Outlier", 
                               TRUE ~ "Non-outlier")) %>%
      filter(outlier != "Non-outlier") %>%
      select(chrom,start,end,mid,the_mds,outlier,n) %>%
      mutate(mds_coord=coord_direction) %>%
      mutate(ahead_n=n-lag(n), behind_n=abs(n-lead(n))) %>%
      mutate(min_dist=pmin(ahead_n, behind_n, na.rm=T)) %>%
      filter(min_dist < max_distance_between_outliers) %>%
      select(-ahead_n,-behind_n,-min_dist)
    windows <- current_windows %>% pull(n)
    outlier_windows <- rbind(outlier_windows, current_windows)
    if(nrow(current_windows) < 1){
      next;
    }
  }
  
  out <- cov_pca(win.fn.snp(windows), k=2)
  matrix.out <- t(matrix(out[4:length(out)], ncol=nrow(samples), byrow=T))
  out <- as_tibble(cbind(samples, matrix.out)) %>% 
    rename(name=sample, PC1="1", PC2="2") %>%       
    mutate(PC1=as.double(PC1), PC2=as.double(PC2))
  try_3_clusters <-try(kmeans(matrix.out[,1], 3, centers=c(min(matrix.out[,1]), (min(matrix.out[,1])+max(matrix.out[,1]))/2, max(matrix.out[,1]))))
  if("try-error" %in% class(try_3_clusters)){
    kmeans_cluster <-kmeans(matrix.out[,1], 2, centers=c(min(matrix.out[,1]), max(matrix.out[,1])))
  }else{
    kmeans_cluster <- kmeans(matrix.out[,1], 3, centers=c(min(matrix.out[,1]), (min(matrix.out[,1])+max(matrix.out[,1]))/2, max(matrix.out[,1])))
  }
  out$cluster <- kmeans_cluster$cluster - 1
  out$cluster <- as.character(out$cluster)
  out$mds_coord <- paste(coord, direction, sep="-")
  genotype.out <- out %>% select(mds_coord,name,PC1,cluster) %>% rename(genotype=cluster)
  cluster_genotypes <- rbind(cluster_genotypes, genotype.out)
}


# Counting
mds_counts <- tibble(mds_coord=character(), n_outliers=numeric())
for (mds in unique(outlier_windows$mds_coord)){
  count_outliers <- outlier_windows %>%
    filter(mds_coord == mds) %>% nrow()
  tmp_tibble <- tibble(mds_coord=as.character(mds), n_outliers=as.numeric(count_outliers))
  mds_counts <- rbind(mds_counts, tmp_tibble)
}


# Collapsing correlated MDS
correlated_mds <- tibble(mds1=character(), mds2=character(), correlation=numeric())
total_mds_coords <- unique(outlier_windows$mds_coord)
for (i in 1:(length(total_mds_coords)-1)){
  for (j in (i+1):length(total_mds_coords)){
    mds1 <- total_mds_coords[i]
    mds2 <- total_mds_coords[j]
    chr1 <- outlier_windows %>% filter(mds_coord == mds1) %>% select(chrom) %>% unique() %>% pull()
    chr2 <- outlier_windows %>% filter(mds_coord == mds2) %>% select(chrom) %>% unique() %>% pull()
    if (chr1 != chr2){next;}
    cluster_genotypes %>% mutate(mds_coord=gsub("_", "-", mds_coord)) %>% filter(mds_coord == mds1 | mds_coord == mds2) %>%
      select(-PC1) %>%
      spread(mds_coord, genotype) %>% select(-name) -> tmp
    x <- tmp %>% pull(1) %>% as.numeric()
    y <- tmp %>% pull(2) %>% as.numeric()
    test_result <- cor.test(x, y, na.rm=T)
    tmp_tibble <- tibble(mds1=as.character(mds1), mds2=as.character(mds2), correlation=as.numeric(abs(test_result$estimate)))
    correlated_mds <- rbind(correlated_mds, tmp_tibble)
  }
}
for (i in 1:nrow(correlated_mds)){
  if (correlated_mds[i,3] >= min_cor){
    count1 <- mds_counts %>% filter(mds_coord == as.character(correlated_mds[i,1])) %>% pull(n_outliers)
    count2 <- mds_counts %>% filter(mds_coord == as.character(correlated_mds[i,2])) %>% pull(n_outliers)
    if (count1 < count2){
      total_mds_coords[which(total_mds_coords != as.character(correlated_mds[i,1]))] -> total_mds_coords
    }else{
      total_mds_coords[which(total_mds_coords != as.character(correlated_mds[i,2]))] -> total_mds_coords
    }
  }
}

