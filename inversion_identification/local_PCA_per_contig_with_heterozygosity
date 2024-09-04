install.packages("devtools")
library(devtools)
devtools::install_github("petrelharp/local_pca/lostruct")
library(lostruct)
library(tidyverse)
library(ggplot2)

library(zoo)
library(templater)

setwd("path/to/your/data")

label=read.table("path/to/populations_meta.txt", header=TRUE, skip=0, sep="\t")
samples <- read_tsv("path/to/sample_info.tsv") [,1]
colnames(samples) <- c("sample")


k_kept <- 10
data_name <- ".bcf"
window_size <- 1e5
chrom <- ""

sites <- vcf_positions(data_name)
win.fn.snp <- vcf_windower(data_name, size=window_size, type='bp', sites=sites) # opening .bcf file

system.time( snp.pca <- eigen_windows(win.fn.snp ,k=2, mc.cores=10) )
system.time( pcdist <- pc_dist( snp.pca ) )

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



#for contig9
# label %>% filter(SEX == "F") -> samples_females_only
# samples_females_only <- as_tibble(samples_females_only[ ,1]) 
# colnames(samples_females_only)[1] <- c("sample")
# 
# label_females_only <- label %>% filter(SEX == "F") 


# for northern populations
# label %>% filter(COUNTRY == "FIN" | COUNTRY == "NOR" | COUNTRY == "SWE") -> label_northern
# 
# label %>% filter(COUNTRY == "FIN" | COUNTRY == "NOR" | COUNTRY == "SWE") -> samples_northern
# samples_northern <- as_tibble(samples_northern[ ,1]) 
# colnames(samples_northern)[1] <- c("sample")




#clustering with k-means

win.regions.clustered <- kmeans(win.regions[, 6:35], 6,iter.max = 10000000, nstart = 50, trace=FALSE)  # munber of clusters depends on the preliminary run with the arbitrary number
main_groups <- data.frame(win.regions.clustered[["cluster"]])
centers <- data.frame(win.regions.clustered[["centers"]])
grouping_by_MDS <- data.frame(win.regions, main_groups)
grouping_by_MDS <- rename(grouping_by_MDS, clusters = win.regions.clustered...cluster...)

# CLUSTER coloring 
print(
  grouping_by_MDS %>%
    gather(., mds, value, colnames(win.regions)[6:7]) %>% 
    ggplot(.,aes(x=mid,y=value, colour= as.factor(clusters))) + geom_point() + facet_grid(mds~.,scales = "free") +
    theme_bw()
  
)


#plotting PCA with the outlying windows after k-means
windows <- grouping_by_MDS %>% filter(clusters == 4) %>% pull(n) 

pca.test <- cov_pca(win.fn.snp(windows), k=2) 
out <- pca.test
matrix.out <- t(matrix(out[4:length(out)], ncol=nrow(samples) , byrow=T))


PCA_of_a_cluster <- data.frame(matrix.out)
PCA_of_a_cluster <- rename(PCA_of_a_cluster, PC1 = X1, PC2 = X2)

names(label)[1] <- "name"
df <- cbind(PCA_of_a_cluster, label)
df$POP=as.factor(df$POP)

PCA_plot <- ggplot(df,aes(PC1, PC2, colour=COUNTRY))+geom_point()+theme_bw()+
  xlab("PC1") + ylab("PC2") + theme(text = element_text(size=15))
summary(df$POP)
print(PCA_plot)


#heterozigosity check
win.regions %>% filter(n %in% windows) %>% summarize(start = min(start), end = max(end)) -> tmp.region
tmp.region$start
tmp.region$end

out <- pca.test
out <- matrix(out[4:length(out)],ncol=nrow(samples),byrow=T) %>% as_tibble() 
colnames(out) <- pull(samples)


out <- as_tibble(cbind(name = names(out), t(out))) %>%
  rename(PC1=V2, PC2=V3) %>% mutate(PC1 = as.double(PC1), PC2 = as.double(PC2))


system(paste("bcftools query -H -f '%END [ %GT]\n' -r ", chrom, ":", tmp.region$start,"-", tmp.region$end," ",data_name,
             '| sed s/END/pos/g > tmp.geno.txt',sep=""))


read_delim("tmp.geno.txt",delim=" ",col_names = c("pos", as.vector(samples$sample)),skip=1, trim_ws = T) %>%
  mutate_if(., 
            is.character, 
            str_replace_all, pattern = '0\\|1', replacement = "1") %>%
  mutate_if(.,
            is.character, 
            str_replace_all, pattern = '1\\|1', replacement = "2") %>%
  mutate_if(.,
            is.character, 
            str_replace_all, pattern = '0\\|0', replacement = "0") %>%
  mutate_if(., 
            is.character, 
            str_replace_all, pattern = '0/0', replacement = "0") %>%
  mutate_if(.,
            is.character, 
            str_replace_all, pattern = '1/1', replacement = "2") %>%
  mutate_if(., 
            is.character, 
            str_replace_all, pattern = '0/1', replacement = "1") %>%
  mutate_if(., 
            is.character, 
            str_replace_all, pattern = './.', replacement = "NA") %>%
  mutate_if(., 
            is.character, 
            str_replace_all, pattern = '\\.', replacement = "NA") -> snps



snps %>%  group_by(pos) %>% gather("name","genotype",2:ncol(snps)) %>%
  group_by(name, genotype) %>%
  summarize(count=n()) %>%
  spread(genotype, count) %>%
  summarize(het=`1`/(`0` + `1` + `2`)) -> heterozygosity

out %>% 
  inner_join(.,heterozygosity) %>%
  ggplot(.,aes(x=PC1,y=PC2)) + geom_point(aes(color=het)) + theme_bw() + scale_colour_viridis_c()


