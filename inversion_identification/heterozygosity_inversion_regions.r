#here we calculate the levels of heterozygosity per individual per inverison region and synchronize them with the PCA coordinates
# for the following plotting

library(lostruct)
library(ggplot2)
library(tidyverse)




setwd("")
path_to_file <- ""
path_to_output <- ""
list_of_files <- dir(path_to_file, pattern="*.eigenvec", all.files=FALSE, full.names=FALSE)


list_of_files <- dir(path_to_file, pattern="*.eigenvec", all.files=FALSE, full.names=FALSE)



for (filename in list_of_files) {
  file_name_no_ext <- tools::file_path_sans_ext(filename)
  
  
  
  label <-read.table("", header=TRUE, skip=0, sep="\t")
  samples <- label[1]
  eigenvec_filename <- file.path(path_to_file, paste(file_name_no_ext, "eigenvec", sep="."))
  eigenvec <- read.table(eigenvec_filename, header=TRUE, skip=0, sep=" ")
  names(eigenvec)[1] <- "name"
  out <- eigenvec[, c(1, 3:4)]
  
  
  system(paste("bcftools query -H -f '%END [ %GT]\n' ", " ",file_name_no_ext,".bcf",
               '| sed s/END/pos/g > tmp.geno.txt',sep=""))
  
  system(paste("sed 's/  / /' tmp.geno.txt"))
  
  read_delim("tmp.geno.txt",delim=" ",col_names = c("pos", as.vector(samples$SampleID)),skip=1, trim_ws = T)  %>%
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
  
  
  
  write.table(heterozygosity, paste(path_to_output,'/',file_name_no_ext,'_het.txt', sep = ''), row.names = F, col.names = T)
}
