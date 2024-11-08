BiocManager::install("rtracklayer")
library(rtracklayer)

install.packages("data.table")
library(data.table)

library(tidyverse)

path_to_file <- ""
path_to_output <- ""


inversion_coords <- read.table('path/to/Inversion_coords.txt', header=F, skip=0, sep="\t")


list_of_files <- list.files(path_to_file, pattern = '.gtf') 

filename <- list_of_files[1]

matched_all <- tibble()
non_matched_all <- tibble()

for (filename in list_of_files) {
  
  file_name_no_ext <- tools::file_path_sans_ext(filename)
  cont <- str_split(file_name_no_ext, "_")[[1]][2]
  
  Inv_coords <- inversion_coords %>% filter(str_detect(V2, paste0("^", cont, "\\b")))
  
  
  gtf <- read.table(paste(path_to_file, filename, sep = ""), header = F, skip = 1, sep = "\t")
  gtf <- rtracklayer::import(paste(path_to_file, filename, sep = ""))
  gtf=as.data.frame(gtf)
  
  

  Inv_coords <- as.data.table(Inv_coords)
  gtf <- as.data.table(gtf)
  
  # Add a key column to `gtf`
  gtf[, key := .I]
  
  matches <- gtf[Inv_coords, on = .(
                                    start >= V3, 
                                    end <= V4), 
                 nomatch = 0L, 
                 .(key, x.start, x.end, i.V3, i.V4), allow.cartesian=TRUE]
  print(nrow(matches))
  print(head(matches))
  
  # Extract matched rows
  matched <- gtf[key %in% matches$key]
  non_matched <- gtf[!key %in% matches$key]
  
  # Display results to confirm correct subsetting
  print(nrow(matched))
  print(nrow(non_matched))
  
  # Clean up, remove the 'key' column
  gtf[, key := NULL]
  matched[, key := NULL]
  non_matched[, key := NULL]
  
  matched_all <- rbind(matched_all, matched)
  non_matched_all <- rbind(non_matched_all, non_matched)
  
  
  
}

inversion_genes <- unique(matched_all$gene_id)
collinear_genes <- unique(non_matched_all$gene_id)

length(inversion_genes) / 47880549
length(collinear_genes) / 169134853


### plotting

i <- 1

gene_dens_per_inv <- tibble(Inv_n=character(), gene_density = numeric())

for (i in 1:nrow(inversion_coords)) {
  
  
  inv_n <- inversion_coords[i,1]
  cont <- inversion_coords[i,2]
  
  
  start_i <- inversion_coords[i,3]
  end_i <- inversion_coords[i,4]
  
  tmp_inv <- matched_all %>% filter(start_i < start & end_i > end & cont == seqnames )
  
  inversion_genes <- unique(tmp_inv$gene_id)
  
  g_d <- length(inversion_genes) / (end_i - start_i)
  
  gene_dens_per_inv_tmp <- tibble(Inv_n=inv_n, gene_density = g_d)
  
  gene_dens_per_inv <- rbind(gene_dens_per_inv, gene_dens_per_inv_tmp)
  
  
}
#let's check now!

mean_gd <- mean(gene_dens_per_inv$gene_density)
  
write.table(gene_dens_per_inv, paste(path_to_output,'gene_density_per_inversion.txt', sep = '/'), , row.names = F, col.names = T) # should be finished?



