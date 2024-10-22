#After the initial LD calculations we then summarize LD per 10000 pb window for each of the analyzed contigs 
#It is done so the upper triangle has the LD calculated only in the individuals that are homozygous for one orientation with the major
#frequency. Lower triangle shows the LD calculated for all individuals
#It is better to run this code on the backround concerning that it is considerably slow 


### loading packages


if (!require("stringr")) {
  install.packages("stringr", dependencies = TRUE)
  library(stringr)
}


if (!require("Matrix")) {
  install.packages("Matrix", dependencies = TRUE)
  library(Matrix)
}

# library(stringr)
# library(Matrix)

path_to_file <- "~/" 
path_to_output <- "~/"
path_to_contig_ld <- '~/'

list_of_data_files <- dir(path_to_file, pattern="*\\.ld", all.files=FALSE, full.names=FALSE)
list_of_pos_files <- dir(path_to_file, pattern="*\\.pos.txt", all.files=FALSE, full.names=FALSE)


for (filename in list_of_data_files) {
  file_name_no_ext <- tools::file_path_sans_ext(filename)
  contig <- str_split(file_name_no_ext, "_")[[1]] [1]
  
  #data from the homozygous  sample 
  data_filename_inv <- file.path(path_to_file, paste(file_name_no_ext, "ld", sep="."))
  ld_inv <- read.table(data_filename_inv, strip.white = TRUE, header = FALSE, stringsAsFactors=FALSE)
  #data for the whole inversion 
  data_filename_contig <- file.path(path_to_contig_ld, paste(contig, "ld", sep="."))
  ld_contig <- read.table(data_filename_contig, strip.white = TRUE, header = FALSE, stringsAsFactors=FALSE)
  
  
  #SNP positions for the homozygous sample
  positions_filename_inv <- file.path(path_to_file, paste(file_name_no_ext, "pos","txt", sep="."))
  SNP_positions_inv <- read.table(positions_filename_inv, header = TRUE)
  SNP_positions_inv <- as.vector(SNP_positions_inv[, 1])
  #SNP positions for the whole inversion 
  positions_filename_contig <- file.path(path_to_contig_ld, paste(contig, "pos","txt", sep="."))
  SNP_positions_contig <- read.table(positions_filename_contig, header = TRUE)
  SNP_positions_contig<- as.vector(SNP_positions_contig[, 1]) 
  
  
  
  #matrix of windows of 10 kb 
  max_kb <- max(SNP_positions_inv[length(SNP_positions_inv)],SNP_positions_contig[length(SNP_positions_contig)]) 
  window_size <- 10000
  window_max <- ceiling(max_kb/window_size) 
  
  
  ###inversion 
  column_names_inv <- c()
  row_names_inv <- c()
  ld_matrix_inv <- matrix( 0, nrow = window_max, ncol = window_max)
  
  
  window_x <- 1
  window_y <- 1
  for (window_x in 1:window_max) {
    for (window_y in window_x:window_max) {
      col_x_start <- (window_x - 1) * window_size + 1 
      col_x_end <- (window_x) * window_size
      
      col_y_start <- (window_y - 1) * window_size + 1
      col_y_end <- (window_y) * window_size
      
      #print(which(SNP_positions >= col_x_start & SNP_positions < col_x_end))
      #print(which(SNP_positions >= col_y_start & SNP_positions < col_y_end))
      x_indeces <- which(SNP_positions_inv >= col_x_start & SNP_positions_inv < col_x_end)
      y_indeces <- which(SNP_positions_inv >= col_y_start & SNP_positions_inv < col_y_end)
      ld_inv_cut <- ld_inv[y_indeces , x_indeces]
      
      if (length(x_indeces) == 0 | length(y_indeces) == 0) {
        R2_2nd <- 0
      } else {
        R2_2nd <- sapply(ld_inv_cut, max)
        R2_2nd <- sort(as.vector(R2_2nd), decreasing = T)
        # print(c(window_x, window_y, R2_2nd))
        if (length(R2_2nd) == 1) {
          R2_2nd <- R2_2nd[1]
        } else {
          R2_2nd <- R2_2nd[2]
        }
      }
      row_names_inv <- append(row_names_inv, paste(col_y_end))
      row_names_inv <- unique(row_names_inv)
      ld_matrix_inv[window_y, window_x] <- R2_2nd
      ld_matrix_inv[window_x, window_y] <- R2_2nd
      
    }
    column_names_inv <- append(column_names_inv, paste(col_x_end))
  }
  
  colnames(ld_matrix_inv) <- c(column_names_inv)
  rownames(ld_matrix_inv) <-c(row_names_inv)
      
      
  ###whole contig
  column_names_contig <- c()
  row_names_contig <- c()
  ld_matrix_contig <- matrix( 0, nrow = window_max, ncol = window_max)
  
  
  window_x <- 1
  window_y <- 1
  for (window_x in 1:window_max) {
    for (window_y in window_x:window_max) {
      col_x_start <- (window_x - 1) * window_size + 1 
      col_x_end <- (window_x) * window_size
      
      col_y_start <- (window_y - 1) * window_size + 1
      col_y_end <- (window_y) * window_size
      
      #print(which(SNP_positions >= col_x_start & SNP_positions < col_x_end))
      #print(which(SNP_positions >= col_y_start & SNP_positions < col_y_end))
      x_indeces <- which(SNP_positions_contig >= col_x_start & SNP_positions_contig < col_x_end)
      y_indeces <- which(SNP_positions_contig >= col_y_start & SNP_positions_contig < col_y_end)
      ld_contig_cut <- ld_contig[y_indeces , x_indeces]
      
      if (length(x_indeces) == 0 | length(y_indeces) == 0) {
        R2_2nd <- 0
      } else {
        R2_2nd <- sapply(ld_contig_cut, max)
        R2_2nd <- sort(as.vector(R2_2nd), decreasing = T)
        # print(c(window_x, window_y, R2_2nd))
        if (length(R2_2nd) == 1) {
          R2_2nd <- R2_2nd[1]
        } else {
          R2_2nd <- R2_2nd[2]
        }
      }
      row_names_contig <- append(row_names_contig, paste(col_y_end))
      row_names_contig <- unique(row_names_contig)
      ld_matrix_contig[window_y, window_x] <- R2_2nd
      ld_matrix_contig[window_x, window_y] <- R2_2nd
      
    }
    column_names_contig <- append(column_names_contig, paste(col_x_end))
  }
  
  colnames(ld_matrix_contig) <- c(column_names_contig)
  rownames(ld_matrix_contig) <-c(row_names_contig)
  
  

  
  ld_matrix_final <- matrix( 0, nrow = window_max, ncol = window_max)
  ld_matrix_final[upper.tri(ld_matrix_final)] <- ld_matrix_inv[upper.tri(ld_matrix_inv)]
  ld_matrix_final[lower.tri(ld_matrix_final)] <-  ld_matrix_contig[lower.tri(ld_matrix_contig)]
  
  colnames(ld_matrix_final) <- c(column_names_inv)
  rownames(ld_matrix_final) <-c(row_names_contig)
  
  write.table(ld_matrix_final, paste(path_to_output,'/',file_name_no_ext,'.txt', sep = ''), row.names = T, col.names = T)
  

}


