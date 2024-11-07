library(ggplot2)
library("tidyverse")
library("dplyr")



major_minor <- read.table("path/to/majour_min_al_rec.txt", header=T, skip=0, sep="\t")
names <- read.table("path/to/names.txt", header=F, skip=0, sep="\t")
names$region <- paste(names$V2,names$V3, names$V4, names$V5, sep = '_')


label <-read.table("path/to/Labels_Ips_240_HC_data.txt", header=TRUE, skip=0, sep="\t")
path_to_file <- ""
list_of_files <- dir(path_to_file, pattern="*.genotyped", all.files=FALSE, full.names=FALSE)



allele_fr_count <- tibble(A=numeric(), B=numeric(), C=numeric(), file_name_no_ext = character(), name= character())

file_name <- 'IpsContig22_1917434_2117433_I.genotyped'
for (file_name in list_of_files){
  
  file_name_no_ext <- tools::file_path_sans_ext(file_name)
  genotyping <- read.table(paste(path_to_file, file_name, sep = '/'), header=F, skip=0)
  name <- as_tibble(names) %>% filter(region == file_name_no_ext) %>% unlist(., use.names=FALSE)
  name <- name[1]
  
  names(genotyping)[1] <- "SampleID"
  df <- left_join(genotyping, label)
  #tmp <- df %>% filter(POP == pops[i]) %>% count(V2)
  df %>%count(V2) %>% pivot_wider(., names_from = V2, values_from = n, values_fill = 0) -> genotype_count_s
  #write.table(genotype_count_g, paste(' file_name_no_ext, '_gen_fr_g', '.txt', sep = ''), sep ="\t", row.names = F, col.names = T)
  
  contig_ref <- str_split(file_name_no_ext, "_")[[1]] [1]
  start_ref <- str_split(file_name_no_ext, "_")[[1]] [2]
  end_ref <- str_split(file_name_no_ext, "_")[[1]] [3]
  inversion_ref <- str_split(file_name_no_ext, "_")[[1]] [4]
  
  for (i in 1:nrow(genotype_count_s)){
  if (ncol(genotype_count_s) == 3){
    
    N_AA <-pull(genotype_count_s[i,1])
    N_AB <-pull(genotype_count_s[i,2])
    N_BB <- pull(genotype_count_s[i,3])
    
    
    N_ind <- N_AA + N_AB + N_BB
    fr_A <- (2*N_AA + N_AB)/(2*N_ind)
    fr_B <- (2*N_BB + N_AB)/(2*N_ind)
    
    #Gen_freq <- c(N_AA, N_AB, N_BB)
    
    tmp <- tibble(A=fr_A, B=fr_B, C= 0, file_name_no_ext = file_name_no_ext, name= name)
    
    allele_fr_count <- as_tibble(rbind(allele_fr_count, tmp))
  }
  
  if (ncol(genotype_count_s) == 6){
    
    N_AA <-pull(genotype_count_s[i,1])
    N_AB <-pull(genotype_count_s[i,2])
    N_AC <- pull(genotype_count_s[i,3])
    
    N_BB <-pull(genotype_count_s[i,4])
    N_BC <-pull(genotype_count_s[i,5])
    N_CC <- pull(genotype_count_s[i,6])
    
    
    N_ind <- N_AA + N_AB + N_BB + N_CC + N_AC + N_BC
    fr_A <- (2*N_AA + N_AB + N_AC)/(2*N_ind)
    fr_B <- (2*N_BB + N_AB + N_BC)/(2*N_ind)
    fr_C <- (2*N_CC + N_AC + N_BC)/(2*N_ind)
    
    tmp <- tibble(A=fr_A, B=fr_B, C= fr_C, file_name_no_ext = file_name_no_ext, name= name)
    
    allele_fr_count <- as_tibble(rbind(allele_fr_count, tmp))
  }
  
  if (ncol(genotype_count_s) == 5){
    

    N_AA <-pull(genotype_count_s[i,1])
    N_AB <-pull(genotype_count_s[i,2])
    N_AC <- pull(genotype_count_s[i,3])
    
    N_BB <-pull(genotype_count_s[i,4])
    N_BC <-pull(genotype_count_s[i,5])
    N_CC <- 0
    
    
    N_ind <- N_AA + N_AB + N_BB + N_CC + N_AC + N_BC
    fr_A <- (2*N_AA + N_AB + N_AC)/(2*N_ind)
    fr_B <- (2*N_BB + N_AB + N_BC)/(2*N_ind)
    fr_C <- (2*N_CC + N_AC + N_BC)/(2*N_ind)
    
    tmp <- tibble(A=fr_A, B=fr_B, C= fr_C, file_name_no_ext = file_name_no_ext, name= name)
    
    allele_fr_count <- as_tibble(rbind(allele_fr_count, tmp))
  }
}
} 

long <- allele_fr_count %>% pivot_longer(cols = c('A','B','C'), names_to = 'allele', values_to = 'allele_fr')
major_min_no_fr <- major_minor %>% select(-allele_fr)
allele_fr_count_per_inv <- left_join(long, major_min_no_fr, by = c("file_name_no_ext" = 'Inversion' , 'allele' = 'allele' )) %>% filter(is.na(T) == F)

write.table(allele_fr_count_per_inv,  'allele_fr_count_per_inv.txt', sep ="\t", row.names = F, col.names = T)
