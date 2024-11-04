library("stringr")
library("tidyverse")
library("dplyr")
library('ggplot2')

win_size <- 2e5


path_to_file <- ""
path_to_output <- ""
list_of_files <- dir(path_to_file, pattern="*.txt", all.files=FALSE, full.names=FALSE)


major_minor <- read.table("path/to/table/defining/names/based/on/haplotype/frequencies/.txt", header=T, skip=0, sep="\t")
names <- read.table("path/to/names.txt", header=F, skip=0, sep="\t")
collinear_coordinates <- read.table("path/to/collinear_coordinates.txt", header=F, skip=0, sep="\t")


piN_piS_table <- tibble(contig = character(),start = numeric(), end = numeric(), inversion = character(), piN_piS = numeric(), 
                        MJ_MN_COL = character(), gen_fr=numeric(), gene_density=numeric(),
                        names = character(), win_start = numeric(), win_end = numeric())


filename <- "mut_load_sorted_IpsContig10_collinear.txt"
for (filename in list_of_files) {
  
  codon_results <- as_tibble(read.table(paste(path_to_file,filename, sep = "/"), header=T, skip=0, sep="\t"))
  
  no_n_test <- codon_results %>% mutate(is_N = str_detect(codon, 'N'))
  codon_results <- no_n_test %>% filter(is_N != "TRUE")
  

  file_name <- codon_results$file[1] 
  region <- str_split(file_name, "/")[[1]][length(str_split(file_name, "/")[[1]])]
  region <- str_split(region, "\\.")[[1]][1]
  genotyping <- str_split(region, "_")[[1]][length(str_split(region, "_")[[1]])]
  
  
  
  if (genotyping == "collinear" ) {      # else if
    contig <- str_split(region, "_")[[1]] [1]
    
    col_coord_for_contig <- collinear_coordinates %>% filter (V1 == contig)
    
    if(nrow(col_coord_for_contig) == 1) {
      
     
      first_pos <- col_coord_for_contig$V2
      last_pos <- col_coord_for_contig$V3
      
      region_with_data <- codon_results %>% filter(site > first_pos) %>% filter(site < last_pos)
      
      win_n <- ceiling(last_pos/win_size)
      
      site_min <- first_pos
      site_max <- site_min
      
      for(n in 1:win_n){
        
        site_min <- site_max +1  
        site_max <- site_min + win_size 
        
        data_in_window <- region_with_data %>% filter(site > site_min) %>% filter(site < site_max)
        gene_density <- unique(data_in_window$product) %>% length
        
        if (nrow(data_in_window) == 0) {
          next; 
        } else {
        N_diffs <- data_in_window  %>% pull(N_diffs) %>% sum()
        S_diffs <- data_in_window  %>% pull(S_diffs) %>% sum()
        
        N_sites <- data_in_window  %>% pull(N_sites) %>% sum()
        S_sites <- data_in_window  %>% pull(S_sites) %>% sum()
        
        piN_piS <- (N_diffs/N_sites)/(S_diffs/S_sites)
        tmp <- tibble(contig = contig, start = genotyping, end = genotyping, inversion = genotyping, piN_piS = piN_piS, MJ_MN_COL = genotyping,
                      gen_fr= NA, gene_density=gene_density, names = genotyping, win_start = site_min, win_end = site_max)
        piN_piS_table <- rbind(piN_piS_table, tmp)
        }
      }
    }
    
    if(nrow(col_coord_for_contig) == 2) {
      
      #my shitty windowing          ### rewrite!!!! 

      first_pos1 <- col_coord_for_contig$V2[1]
      last_pos1 <- col_coord_for_contig$V3[1] - 1
      
      first_pos2 <- col_coord_for_contig$V2[2] +1
      last_pos2 <- col_coord_for_contig$V3[2]
      
      region_with_data_1 <- codon_results %>% filter(site > first_pos1) %>% filter(site < last_pos1)
      region_with_data_2 <- codon_results %>% filter(site > first_pos2) %>% filter(site < last_pos2)
      
      win_n1 <- ceiling(last_pos1/win_size)
      win_n2 <- ceiling(last_pos2/win_size)
      
      site_min <- first_pos1
      site_max <- site_min
      
      
      for(n in 1:win_n1){
        
        site_min <- site_max +1  
        site_max <- site_min + win_size 
        
        data_in_window <- region_with_data_1 %>% filter(site > site_min) %>% filter(site < site_max)
        gene_density <- unique(data_in_window$product) %>% length 
        
        if (nrow(data_in_window) == 0) {
          next 
        } else {
        
        N_diffs <- data_in_window  %>% pull(N_diffs) %>% sum()
        S_diffs <- data_in_window  %>% pull(S_diffs) %>% sum()
        
        N_sites <- data_in_window  %>% pull(N_sites) %>% sum()
        S_sites <- data_in_window  %>% pull(S_sites) %>% sum()
        
        piN_piS <- (N_diffs/N_sites)/(S_diffs/S_sites)
        pN_pS <- N_sites/S_sites
        
        tmp <- tibble(contig = contig, start = genotyping, end = genotyping, inversion = genotyping, piN_piS = piN_piS,MJ_MN_COL = genotyping,
                      gen_fr= NA, gene_density=gene_density, names = genotyping, win_start = site_min, win_end = site_max)
        piN_piS_table <- rbind(piN_piS_table, tmp)
        }
      }
      
      
      site_min <- first_pos2
      site_max <- site_min
      
      for(n in 1:win_n2){
        
        site_min <- site_max +1  
        site_max <- site_min + win_size 
        
        data_in_window <- region_with_data_2 %>% filter(site > site_min) %>% filter(site < site_max)
        gene_density <- unique(data_in_window$product) %>% length 
        
        if (nrow(data_in_window) == 0) {
          next 
        } else {
        
        N_diffs <- data_in_window  %>% pull(N_diffs) %>% sum()
        S_diffs <- data_in_window  %>% pull(S_diffs) %>% sum()
        
        N_sites <- data_in_window  %>% pull(N_sites) %>% sum()
        S_sites <- data_in_window  %>% pull(S_sites) %>% sum()
        
        piN_piS <- (N_diffs/N_sites)/(S_diffs/S_sites)

        
        tmp <- tibble(contig = contig, start = genotyping, end = genotyping, inversion = genotyping, piN_piS = piN_piS, MJ_MN_COL = genotyping,
                      gen_fr= NA, gene_density=gene_density, names = genotyping, win_start = site_min, win_end = site_max)
        piN_piS_table <- rbind(piN_piS_table, tmp)
        }
      }
    }
    
  } else {
    
    genotyping <- str_split(region, "_")[[1]] [5]     
    contig <- str_split(region, "_")[[1]] [1]
    start <- str_split(region, "_")[[1]] [2]
    end <- str_split(region, "_")[[1]] [3]
    inversion <- str_split(region, "_")[[1]] [4]
    
    actual_region <- str_split(region, "_")[[1]] [1:4]
    actual_region <- paste(actual_region, collapse = '_')
    
    first_pos <- as.numeric(start)
    last_pos <-as.numeric(end)
    
    region_with_data <- codon_results %>% filter(site > first_pos) %>% filter(site < last_pos)
    
    win_n <- ceiling(last_pos/win_size)
    
    site_min <- first_pos
    site_max <- site_min

    for(n in 1:win_n){
      
      site_min <- site_max +1  
      site_max <- site_min + win_size 
      
      data_in_window <- region_with_data %>% filter(site > site_min) %>% filter(site < site_max)
      gene_density <- unique(data_in_window$product) %>% length 
      
      if (nrow(data_in_window) == 0) {
        next 
      } else {
        
      N_diffs <- data_in_window  %>% pull(N_diffs) %>% sum()
      S_diffs <- data_in_window  %>% pull(S_diffs) %>% sum()
      
      N_sites <- data_in_window  %>% pull(N_sites) %>% sum()
      S_sites <- data_in_window  %>% pull(S_sites) %>% sum()
      
      piN_piS <- (N_diffs/N_sites)/(S_diffs/S_sites)
      
      name <- as_tibble(names) %>% filter(V2 == contig & V3 == start & V4 == end & V5 == inversion) %>% unlist(., use.names=FALSE)
      name <- name[1]
      
      gen_fr <- major_minor %>% filter(Inversion == actual_region, T == genotyping) %>% unlist(., use.names=FALSE)
      gen_fr <- as.numeric(gen_fr[3])
      
      
      tmp <- tibble(contig = contig, start = start, end = end, inversion = inversion, piN_piS = piN_piS, MJ_MN_COL = genotyping, 
                    gen_fr=gen_fr, gene_density=gene_density, names = name, win_start = site_min, win_end = site_max)
      piN_piS_table <- rbind(piN_piS_table, tmp)
      }
    }
  }
  
  
}



piN_piS_table <- piN_piS_table  %>%  filter(is.nan(piN_piS) == F) %>% filter(piN_piS != Inf)

write.table(piN_piS_table, '/home/anstmykh/Desktop/mut_load/piN_piS_table.txt', row.names = TRUE)

