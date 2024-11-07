library("tidyverse")
library("dplyr")
library('ggplot2')

path_to_file <- ""
path_to_output <- ""
list_of_files <- dir(path_to_file, pattern="*.genotyped", all.files=FALSE, full.names=FALSE)


allele_fr_count <- tibble(Inversion=character(), A=numeric(), B=numeric(), C=numeric())

for (filename in list_of_files) {

  file_name_no_ext <- tools::file_path_sans_ext(filename)
  genotypes <- as_tibble(read.table(paste(path_to_file, filename, sep = ""), header = F, skip=0, sep="\t"))


  genotypes %>% count(V2) %>% pivot_wider(., names_from = V2, values_from = n, values_fill = 0, ) -> genotype_count


  if (ncol(genotype_count) == 3){
    
    N_AA <-pull(genotype_count[1,1])
    N_AB <-pull(genotype_count[1,2])
    N_BB <- pull(genotype_count[1,3])
    
    
    N_ind <- N_AA + N_AB + N_BB
    fr_A <- (2*N_AA + N_AB)/(2*N_ind)
    fr_B <- (2*N_BB + N_AB)/(2*N_ind)
    
    tmp <- tibble(Inversion=file_name_no_ext, A=fr_A, B=fr_B, C= 0)
    
    allele_fr_count <- as_tibble(rbind(allele_fr_count, tmp))
  }
  
  if (ncol(genotype_count) == 6){
    
    N_AA <-pull(genotype_count[1,1])
    N_AB <-pull(genotype_count[1,2])
    N_AC <- pull(genotype_count[1,3])
    
    N_BB <-pull(genotype_count[1,4])
    N_BC <-pull(genotype_count[1,5])
    N_CC <- pull(genotype_count[1,6])
    
    
    N_ind <- N_AA + N_AB + N_BB + N_CC + N_AC + N_BC
    fr_A <- (2*N_AA + N_AB + N_AC)/(2*N_ind)
    fr_B <- (2*N_BB + N_AB + N_BC)/(2*N_ind)
    fr_C <- (2*N_CC + N_AC + N_BC)/(2*N_ind)
    
    tmp <- tibble(Inversion=file_name_no_ext, A=fr_A, B=fr_B, C= fr_C)
    
    allele_fr_count <- as_tibble(rbind(allele_fr_count, tmp))
  }
  if (ncol(genotype_count) == 5){
    
    N_AA <-pull(genotype_count[1,1])
    N_AB <-pull(genotype_count[1,2])
    N_AC <- pull(genotype_count[1,3])
    
    N_CC <- 0
    N_BB <-pull(genotype_count[1,4])
    N_BC <- pull(genotype_count[1,5])
    
    
    N_ind <- N_AA + N_AB + N_BB + N_CC + N_AC + N_BC
    fr_A <- (2*N_AA + N_AB + N_AC)/(2*N_ind)
    fr_B <- (2*N_BB + N_AB + N_BC)/(2*N_ind)
    fr_C <- (2*N_CC + N_AC + N_BC)/(2*N_ind)
    
    tmp <- tibble(Inversion=file_name_no_ext, A=fr_A, B=fr_B, C= fr_C)
    
    allele_fr_count <- as_tibble(rbind(allele_fr_count, tmp))
  }
}

#allele_fr_count %>% pivot_longer(cols = c('A','B','C'), names_to = 'allele', values_to = 'allele_fr')  #### wery useful!!!!

long <-allele_fr_count %>% pivot_longer(cols = c('A','B','C'), names_to = 'allele', values_to = 'allele_fr')

long %>% group_by(Inversion) %>% mutate('T' = case_when(
  (sum(allele_fr==0)!= 0) & (allele_fr == max(allele_fr)) ~ "MJ",
  (sum(allele_fr==0)!= 0) & (!(allele_fr %in% c(max(allele_fr),0))) ~ "MN",
  (sum(allele_fr==0)== 0) & (allele_fr == max(allele_fr)) ~ "MJ",
  (sum(allele_fr==0)== 0) & (allele_fr == min(allele_fr)) ~ "MN",
  (sum(allele_fr==0)== 0) & (!(allele_fr %in% c(max(allele_fr),min(allele_fr)))) ~ "INT",
  T ~ "N"
  
  )) -> majour_min


majour_min <- majour_min %>% filter(T != 'N')
write.table(majour_min, paste(path_to_output, "majour_min_al.txt"), row.names = F,  col.names = T, sep = '\t')

long_no_C <- long %>% filter(allele != 'C')
long_no_C %>% group_by(Inversion) %>% mutate('T' = case_when(
  (allele_fr == max(allele_fr)) ~ "MJ",
  (allele_fr == min(allele_fr)) ~ "MN"
  
)) -> majour_min_no_C

write.table(majour_min_no_C, paste(path_to_output, "majour_min_al_no_C_5north.txt", sep = ''), row.names = F,  col.names = T, sep = '\t')


long_C <- long %>% filter(allele == 'C') %>% filter(allele_fr != 0) %>% mutate('T' = 'R')
maj_min_rec <- rbind(majour_min_no_C, long_C)
write.table(maj_min_rec, paste(path_to_output, "majour_min_al_rec.txt", sep = ''), row.names = F,  col.names = T, sep = '\t')

