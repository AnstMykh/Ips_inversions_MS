library(ggplot2)
library("ggfortify")
library("tidyverse")
library("dplyr")
library("HardyWeinberg")


major_minor <- read.table("path/to/majour_min_gen_rec.txt", header=T, skip=0, sep="\t")
names <- read.table("path/to/names.txt", header=F, skip=0, sep="\t")


label <- read.table("path/to/metadata/Labels_Ips_240_HC_data.txt", header=T, skip=0)
path_to_file <- ""
list_of_files <- dir(path_to_file, pattern="*.genotyped", all.files=FALSE, full.names=FALSE)


#### for geographical region no FDR
#AA =numeric(), AB= numeric(), AC = numeric(),  BB = numeric(), BC = numeric(), CC= numeric()
WH_geography <- tibble(contig = character(), start=numeric(), end=numeric(), inversion = character(), G = character(),
                       AA =numeric(), AB= numeric(), AC = numeric(),  BB = numeric(), BC = numeric(), CC= numeric(),
                       exAA =numeric(), exAB= numeric(), exAC = numeric(),  exBB = numeric(), exBC = numeric(), exCC= numeric(),HW_pval = numeric()) # do it now 

HW_sample <- tibble(contig = character(), start=numeric(), end=numeric(), inversion = character(),
                    AA =numeric(), AB= numeric(), AC = numeric(),  BB = numeric(), BC = numeric(), CC= numeric(), HW_pval = numeric()) 


for (file_name in list_of_files){
  
  file_name_no_ext <- tools::file_path_sans_ext(file_name)
  genotyping <- read.table(file_name, header=F, skip=0)
  names(genotyping)[1] <- "SampleID"
  df <- left_join(genotyping, label)
  #tmp <- df %>% filter(POP == pops[i]) %>% count(V2)
  df %>% group_by(GROUP) %>% count(V2)%>% pivot_wider(., names_from = V2, values_from = n, values_fill = 0, names_sort = T) -> genotype_count_g
  df %>%count(V2) %>% pivot_wider(., names_from = V2, values_from = n, values_fill = 0) -> genotype_count_s
  #write.table(genotype_count_g, paste('/home/anstmykh/Desktop/allele_fr/results/', file_name_no_ext, '_gen_fr_g', '.txt', sep = ''), sep ="\t", row.names = F, col.names = T)
  
  contig_ref <- str_split(file_name_no_ext, "_")[[1]] [1]
  start_ref <- str_split(file_name_no_ext, "_")[[1]] [2]
  end_ref <- str_split(file_name_no_ext, "_")[[1]] [3]
  inversion_ref <- str_split(file_name_no_ext, "_")[[1]] [4]
  
  
  ##### HW for the whole sample 
  
  #### contig 9 
  if (contig_ref == "IpsContig9"){
    N_AA <-pull(genotype_count_s[1,1])
    N_AB <-pull(genotype_count_s[1,2])
    N_BB <- pull(genotype_count_s[1,3])
    
    Nind <- N_AA + N_AB + N_BB
    fr_A <- (2*N_AA + N_AB)/ (2*Nind)
    fr_B <- (2*N_BB + N_AB) / (2*Nind)
    
    Exp_AA <- (fr_A^2)*Nind
    Exp_AB <- (2*fr_A*fr_B)*Nind
    Exp_BB <- (fr_B^2)*Nind
    
    F_count <- df %>% filter(SEX=='F') %>% count(V2) %>% pivot_wider(., names_from = V2, values_from = n, values_fill = 0)
    
    F_AA <-pull(F_count[1,1])
    F_AB <-pull(F_count[1,2])
    F_BB <- pull(F_count[1,3])
    
    M_count <- df %>% filter(SEX=='M') %>% count(V2) %>% pivot_wider(., names_from = V2, values_from = n, values_fill = 0)
    #M_num <- df %>% filter(SEX=='M') %>% group_by(V2) %>% nrow()
    #M_fraq <- round((M_num*100/nrow(df))/100, digits = 2)
    N_A <-pull(M_count[1,1])
    N_B <-pull(M_count[1,2])
    
    
    Gen_freq <- c(A=N_A, B=N_B,AA= F_AA, AB=F_AB,BB= F_BB)
    
    
    HWE_X_linked <- HWExact(Gen_freq, x.linked = T)#, phifixed = M_fraq)   
    
    tmp <- tibble(contig = contig_ref, start=start_ref, end=end_ref, inversion = inversion_ref,
                  AA = N_AA, AB= N_AB, AC = 'N',  BB = N_BB, BC = 'N', CC= 'N' ,
                  exAA =Exp_AA, exAB= Exp_AB, exAC = 'N',  exBB = Exp_BB, exBC = 'N', exCC= 'N',HW_pval = HWE_X_linked$pval)
    HW_sample <- rbind(HW_sample, tmp)   
    
  } else {
    
    if (ncol(genotype_count_s) == 3){
      N_AA <-pull(genotype_count_s[1,1])
      N_AB <-pull(genotype_count_s[1,2])
      N_BB <- pull(genotype_count_s[1,3])
      
      Nind <- N_AA + N_AB + N_BB
      fr_A <- (2*N_AA + N_AB)/ (2*Nind)
      fr_B <- (2*N_BB + N_AB) / (2*Nind)
      
      Exp_AA <- (fr_A^2)*Nind
      Exp_AB <- (2*fr_A*fr_B)*Nind
      Exp_BB <- (fr_B^2)*Nind
      
      Gen_freq <- c(N_AA, N_AB, N_BB)
      
      Exc_HWE <- HWExact(Gen_freq)   ## Hardy_Weinberg exact when biallelic 
      
      tmp <- tibble(contig = contig_ref, start=start_ref, end=end_ref, inversion = inversion_ref,
                    AA = N_AA, AB= N_AB, AC = 'N',  BB = N_BB, BC = 'N', CC= 'N' , 
                    exAA =Exp_AA, exAB= Exp_AB, exAC = 'N',  exBB = Exp_BB, exBC = 'N', exCC= 'N',HW_pval = Exc_HWE$pval)
      HW_sample <- rbind(HW_sample, tmp)
      
    }
    if (ncol(genotype_count_s) == 6){
      
      N_AA <-pull(genotype_count_s[1,1])
      N_AB <-pull(genotype_count_s[1,2])
      N_AC <-pull(genotype_count_s[1,3])
      
      N_BB <- pull(genotype_count_s[1,4])
      N_BC <-pull(genotype_count_s[1,5])
      N_CC <- pull(genotype_count_s[1,6])
      
      Nind <- N_AA + N_AB + N_BB + N_AC + N_BC + N_CC
      fr_A <- (2*N_AA + N_AB + N_AC)/ (2*Nind)
      fr_B <- (2*N_BB + N_AB + N_BC) / (2*Nind)
      fr_C <- (2*N_CC + N_AC + N_BC) / (2*Nind)
      
      
      Exp_AA <- (fr_A^2)*Nind
      Exp_AB <- (2*fr_A*fr_B)*Nind
      Exp_BB <- (fr_B^2)*Nind
      
      Exp_AC <- (2*fr_A*fr_C)*Nind
      Exp_BC <- (2*fr_B*fr_C)*Nind
      Exp_CC <- (fr_C^2)*Nind
      
      Gen_freq <- c(N_AA, N_AB, N_BB)
      
      Gen_freq <- c(AA=N_AA,AB=N_AB,AC=N_AC,BB=N_BB,BC=N_BC,CC=N_CC)
      Gen_freq <- toTriangular(Gen_freq)
      
      HWE_perm_mult <- HWPerm.mult(Gen_freq)   #  approximates exact test probabilities for joint tests for HWE and equality of allele frequencies for variants with multiple alleles
      
      tmp <- tibble(contig = contig_ref, start=start_ref, end=end_ref, inversion = inversion_ref, 
                    AA = N_AA, AB= N_AB, AC = N_AC,  BB = N_BB, BC = N_BC, CC= N_CC ,
                    exAA =Exp_AA, exAB= Exp_AB, exAC = Exp_AC,  exBB = Exp_BB, exBC = Exp_BC, exCC= Exp_CC, HW_pval = HWE_perm_mult$pval)
      HW_sample <- rbind(HW_sample, tmp)   
      
    }  
    
    if (ncol(genotype_count_s) == 5){
      
      N_AA <-pull(genotype_count_s[1,1])
      N_AB <-pull(genotype_count_s[1,2])
      N_BB <- pull(genotype_count_s[1,3])
      
      N_AC <-pull(genotype_count_s[1,4])
      N_BC <-pull(genotype_count_s[1,5])
      N_CC <- 0
      
      
      Nind <- N_AA + N_AB + N_BB + N_AC + N_BC + N_CC
      fr_A <- (2*N_AA + N_AB + N_AC)/ (2*Nind)
      fr_B <- (2*N_BB + N_AB + N_BC) / (2*Nind)
      fr_C <- (2*N_CC + N_AC + N_BC) / (2*Nind)
      
      Exp_AA <- (fr_A^2)*Nind
      Exp_AB <- (2*fr_A*fr_B)*Nind
      Exp_BB <- (fr_B^2)*Nind
      
      Exp_AC <- (2*fr_A*fr_C)*Nind
      Exp_BC <- (2*fr_B*fr_C)*Nind
      Exp_CC <- (fr_C^2)*Nind
      
      
      Gen_freq <- c(AA=N_AA,AB=N_AB,AC=N_AC,BB=N_BB,BC=N_BC,CC=0)
      Gen_freq <- toTriangular(Gen_freq)
      
      HWE_perm_mult <- HWPerm.mult(Gen_freq)   #  approximates exact test probabilities for joint tests for HWE and equality of allele frequencies for variants with multiple alleles
      
      tmp <- tibble(contig = contig_ref, start=start_ref, end=end_ref, inversion = inversion_ref, 
                    AA = N_AA, AB= N_AB, AC = N_AC,  BB = N_BB, BC = N_BC, CC= 0 ,
                    exAA =Exp_AA, exAB= Exp_AB, exAC = Exp_AC,  exBB = Exp_BB, exBC = Exp_BC, exCC= Exp_CC,HW_pval = HWE_perm_mult$pval)
      HW_sample <- rbind(HW_sample, tmp)   
      
    }  
  }
  
  
  
  ##### HW for S/N/P
  i <- 1
  for (i in 1:nrow(genotype_count_g)){
    
    #### contig 9 
    if (contig_ref == "IpsContig9"){
      pop <- pull(genotype_count_g[i,1])
      N_AA <-pull(genotype_count_g[i,2])
      N_AB <-pull(genotype_count_g[i,3])
      N_BB <- pull(genotype_count_g[i,4])
      
      
      Nind <- df %>% filter(GROUP == pop) %>% nrow()
      fr_A <- (2*N_AA + N_AB)/ (2*Nind)
      fr_B <- (2*N_BB + N_AB) / (2*Nind)
      
      Exp_AA <- (fr_A^2)*Nind
      Exp_AB <- (2*fr_A*fr_B)*Nind
      Exp_BB <- (fr_B^2)*Nind
      
      
      F_count <- df %>% filter(SEX=='F') %>% group_by(GROUP) %>% count(V2) %>% pivot_wider(., names_from = V2, values_from = n, values_fill = 0)
      
      
      F_AA <-pull(F_count[i,2])
      F_AB <-pull(F_count[i,3])
      F_BB <- pull(F_count[i,4])
      
      M_count <- df %>% filter(SEX=='M') %>% group_by(GROUP) %>% count(V2) %>% pivot_wider(., names_from = V2, values_from = n, values_fill = 0)
      #M_num <- df %>% filter(SEX=='M') %>% group_by(V2) %>% nrow()
      #M_fraq <- round((M_num*100/nrow(df))/100, digits = 2)
      N_A <-pull(M_count[i,2])
      N_B <-pull(M_count[i,3])
      
      
      Gen_freq <- c(A=N_A, B=N_B, AA=F_AA, AB=F_AB, BB=F_BB)
      
      
      HWE_X_linked <- HWExact(Gen_freq, x.linked = T)#, phifixed = M_fraq)   
      
      tmp <- tibble(contig = contig_ref, start=start_ref, end=end_ref, inversion = inversion_ref, G = pop, 
                    AA = N_AA, AB= N_AB, AC = 'N',  BB = N_BB, BC = 'N', CC= 'N' ,
                    exAA =Exp_AA, exAB= Exp_AB, exAC = 'N',  exBB = Exp_BB, exBC = 'N', exCC= 'N',HW_pval = HWE_X_linked$pval)
      WH_geography <- rbind(WH_geography, tmp)   
      
    } else {
      
      if (ncol(genotype_count_g) == 4){
        pop <- pull(genotype_count_g[i,1])
        N_AA <-pull(genotype_count_g[i,2])
        N_AB <-pull(genotype_count_g[i,3])
        N_BB <- pull(genotype_count_g[i,4])
        
        Nind <- df %>% filter(GROUP == pop) %>% nrow()
        fr_A <- (2*N_AA + N_AB)/ (2*Nind)
        fr_B <- (2*N_BB + N_AB) / (2*Nind)
        
        Exp_AA <- (fr_A^2)*Nind
        Exp_AB <- (2*fr_A*fr_B)*Nind
        Exp_BB <- (fr_B^2)*Nind
        
        
        Gen_freq <- c(N_AA, N_AB, N_BB)
        
        Exc_HWE <- HWExact(Gen_freq)   ## Hardy_Weinberg exact when biallelic 
        
        tmp <- tibble(contig = contig_ref, start=start_ref, end=end_ref, inversion = inversion_ref, G = pop, 
                      AA = N_AA, AB= N_AB, AC = 'N',  BB = N_BB, BC = 'N', CC= 'N' ,
                      exAA =Exp_AA, exAB= Exp_AB, exAC = 'N',  exBB = Exp_BB, exBC = 'N', exCC= 'N',HW_pval = Exc_HWE$pval)
        WH_geography <- rbind(WH_geography, tmp)
        
      }
      if (ncol(genotype_count_g) == 7){
        pop <- pull(genotype_count_g[i,1])
        
        N_AA <-pull(genotype_count_g[i,2])
        N_AB <-pull(genotype_count_g[i,3])
        N_AC <-pull(genotype_count_g[i,4])
        
        N_BB <- pull(genotype_count_g[i,5])
        N_BC <-pull(genotype_count_g[i,6])
        N_CC <- pull(genotype_count_g[i,7])
        
        Nind <- df %>% filter(GROUP == pop) %>% nrow()
        fr_A <- (2*N_AA + N_AB + N_AC)/ (2*Nind)
        fr_B <- (2*N_BB + N_AB + N_BC) / (2*Nind)
        fr_C <- (2*N_CC + N_AC + N_BC) / (2*Nind)
        
        
        Exp_AA <- (fr_A^2)*Nind
        Exp_AB <- (2*fr_A*fr_B)*Nind
        Exp_BB <- (fr_B^2)*Nind
        
        Exp_AC <- (2*fr_A*fr_C)*Nind
        Exp_BC <- (2*fr_B*fr_C)*Nind
        Exp_CC <- (fr_C^2)*Nind
        
        Gen_freq <- c(AA=N_AA,AB=N_AB,AC=N_AC,BB=N_BB,BC=N_BC,CC=N_CC)
        Gen_freq <- toTriangular(Gen_freq)
        
        HWE_perm_mult <- HWPerm.mult(Gen_freq)   #  approximates exact test probabilities for joint tests for HWE and equality of allele frequencies for variants with multiple alleles
        
        tmp <- tibble(contig = contig_ref, start=start_ref, end=end_ref, inversion = inversion_ref, G = pop, 
                      AA = N_AA, AB= N_AB, AC = N_AC,  BB = N_BB, BC = N_BC, CC= N_CC ,
                      exAA =Exp_AA, exAB= Exp_AB, exAC = Exp_AC,  exBB = Exp_BB, exBC = Exp_BC, exCC= Exp_CC, HW_pval = HWE_perm_mult$pval)
        WH_geography <- rbind(WH_geography, tmp)   
        
      }  
      
      if (ncol(genotype_count_g) == 6){
        pop <- pull(genotype_count_g[i,1])
        
        N_AA <-pull(genotype_count_g[i,2])
        N_AB <-pull(genotype_count_g[i,3])
        N_BB <- pull(genotype_count_g[i,4])
        
        N_AC <-pull(genotype_count_g[i,5])
        N_BC <-pull(genotype_count_g[i,6])
        N_CC <- 0
        
        
        Nind <- df %>% filter(GROUP == pop) %>% nrow()
        fr_A <- (2*N_AA + N_AB + N_AC)/ (2*Nind)
        fr_B <- (2*N_BB + N_AB + N_BC) / (2*Nind)
        fr_C <- (2*N_CC + N_AC + N_BC) / (2*Nind)
        
        Exp_AA <- (fr_A^2)*Nind
        Exp_AB <- (2*fr_A*fr_B)*Nind
        Exp_BB <- (fr_B^2)*Nind
        
        Exp_AC <- (2*fr_A*fr_C)*Nind
        Exp_BC <- (2*fr_B*fr_C)*Nind
        Exp_CC <- (fr_C^2)*Nind
        
        
        Gen_freq <- c(AA=N_AA,AB=N_AB,AC=N_AC,BB=N_BB,BC=N_BC,CC=0)
        Gen_freq <- toTriangular(Gen_freq)
        
        HWE_perm_mult <- HWPerm.mult(Gen_freq)   #  approximates exact test probabilities for joint tests for HWE and equality of allele frequencies for variants with multiple alleles
        
        tmp <- tibble(contig = contig_ref, start=start_ref, end=end_ref, inversion = inversion_ref, G = pop, 
                      AA = N_AA, AB= N_AB, AC = N_AC,  BB = N_BB, BC = N_BC, CC= 0 ,
                      exAA =Exp_AA, exAB= Exp_AB, exAC = Exp_AC,  exBB = Exp_BB, exBC = Exp_BC, exCC= Exp_CC, HW_pval = HWE_perm_mult$pval)
        WH_geography <- rbind(WH_geography, tmp)   
        
      }  
    }
  }
  
  
}

singificant_WH_g <- WH_geography %>% filter(HW_pval < 0.05)

write.csv(WH_geography, 'path/to/HW_geography.csv', row.names = TRUE)
write.csv(HW_sample, 'path/to/HW_sample.csv', row.names = TRUE)
