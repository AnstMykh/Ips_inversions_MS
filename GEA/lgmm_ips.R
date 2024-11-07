library(LEA)
library(tidyverse)
library(patchwork)

lfdir <- "LFMM/"


corename <- "merged_maf10_misSNP20_misind30_inversions_as_single_loci"
#corename <- "merged_maf10_misSNP20_misind30"
#corename <- "merged_maf10_misSNP20_misind30_LDpruned_inversions_as_single_loci"

#convert .ped to .geno ----
p <- ped2geno(paste0(corename, ".ped"), output.file = paste0(lfdir, corename, ".geno"))

#convert .geno to .lfmm ----
geno2lfmm(input.file = paste0(lfdir, corename, ".geno"), output.file = paste0(lfdir, corename, ".lfmm"))

#get snmf object ----
#assuming three genetic clusters
#potentially multiple runs are needed
#computationally intensive for the entire dataset
sn <- snmf(p, K = 3, project = "new", CPU = 3)
saveRDS(sn, file=paste0(lfdir, corename, "_snmf_object.rds"))

#impute missing values ----
#Note, that there are only 5% missing data, so the fect of imputation should be minor
impute(sn, paste0(lfdir, corename, ".lfmm"), method = "mode", K = 3)
lfmm_f <- paste0(lfdir, corename, ".lfmm_imputed.lfmm")

#estimate the number of latent factors ----
#looks like there's large effect of imputation on the relative magnitude of the first few eigenvalues
d <- read.lfmm(lfmm_f)
pc <- prcomp(d)
plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")
#it looks like K should be 5 for all SNPs
#and K should be 2 for inversions encoded as haplotypes, but for consistency we'll keep 5 
#as there's a second drop after  K= 5

#generate env data ----
#ids in the same order as in .ped
ids <- read_tsv("ids_from_ped.txt")

#assignment of ids to pops
popinfo <- read_tsv("Labels_Ips_240_HC_data.txt")
n_pop <- popinfo %>% pull(POP) %>% unique %>% length
n_loc <- popinfo %>% pull(LOCALITY) %>% unique %>% length

#environmental data
envinfo <- read_tsv("ips_pop_env_data.txt") %>% rename(LOCALITY = Site_ID) %>% select(LOCALITY, PC1, PC2, PC3)

#adding pop to samples ordered as in .ped
ids_pop <- ids %>% left_join(select(popinfo, SampleID, LOCALITY), by = "SampleID")

#df with individual ids, pops and values of PC1-PC3
ids_env <- ids_pop %>% left_join(envinfo, by = "LOCALITY")

#run lfmm2 ----
kval <- 5
#loop to have separate analysis analysis for each PC
for (pc in c(1:3)){
  #env has to be matrix/vector, so extract from dataframe, +2 as two first columns are descriptors
  e <- ids_env[[pc + 2]]
  l <- lfmm2(lfmm_f, e, K = kval)
  saveRDS(l, file=paste0(lfdir, corename, "_lfmm2_res_object_K_", kval, "_envPC_", pc, ".rds"))
  lt <- lfmm2.test(l, lfmm_f, e)
  saveRDS(lt, file=paste0(lfdir, corename, "_lfmm2_test_object_K_", kval, "_envPC_", pc, ".rds"))
}

#extract p values and prepare for plotting ----

coord <- read_tsv(paste0(corename, ".map"), col_select = c(1, 4), col_names = FALSE) %>% 
  mutate(chr = factor(str_extract(X1, "[0-9]+$"), c(1:36))) %>% select(chr, X4) %>% rename(pos = X4)

inversion_coord <- read_tsv("inversions.map", col_select = c(1, 4), col_names = FALSE) %>% 
  mutate(chr = factor(str_extract(X1, "[0-9]+$"), c(1:36))) %>% select(chr, X4) %>% rename(pos = X4)

# get scaffold lengths
chroms <- read_tsv("ips_contig_lengths.txt") %>% select(-contig) %>% 
  mutate(chr = as.factor(chr),
         cum_end = cumsum(length),
         cum_start = ifelse(chr == "1", 0, lag(cum_end))) 

#ticks in the middle of contigs
tick_labels <- chroms %>% mutate(tick_pos = (cum_start + cum_end)/2)

chr_pos <- coord %>% left_join(chroms, by = "chr")

res_df <- NULL
for (pc in c(1:3)) {
  lfmmt <- readRDS(paste0(lfdir, corename, "_lfmm2_test_object_K_", kval, "_envPC_", pc, ".rds"))
  pval_df <- data.frame(p = lfmmt$pvalues, pfdr = p.adjust(lfmmt$pvalues, method = "fdr")) %>% 
    cbind(chr_pos) %>% 
    mutate(SNP_nr = row_number(),
           cum_pos = pos + cum_start,
           PC = pc) %>% 
    filter(!is.na(p)) %>% 
    select(PC, chr, pos, cum_pos, SNP_nr, p, pfdr)
  res_df <- bind_rows(res_df, pval_df)
}


#save p values for all SNPs
#saveRDS(res_df, paste0("IpsLFMM_results_", corename, "_PC1-3.rds"))

######################################

#Plotting ----
#First load two dataframes with results for all snps and inversions as single loci

res_all <- readRDS("IpsLFMM_results_merged_maf10_misSNP20_misind30_PC1-3.rds") 
res_sl <- readRDS("IpsLFMM_results_merged_maf10_misSNP20_misind30_inversions_as_single_loci_PC1-3.rds")
res_sl_inv_only <- inversion_coord %>% left_join(res_sl, by = c("chr", "pos"))


#to see where's fdr threshold
aaa_all <- res_all %>% filter(pfdr < 0.05) %>% arrange(PC, pfdr)
aaa_sl <- res_sl %>% filter(pfdr < 0.05) %>% arrange(PC, pfdr)

# #df to plot lines at fdr 0.05 threshold, PC3 very much approximate
thr_all <- data.frame(PC = c(1, 2, 3),
                      y = c(-log10(2e-6), -log10(3e-6), -log10(3e-6)))

thr_sl <- data.frame(PC = c(1, 2, 3),
                      y = c(-log10(3.5e-6), -log10(1.6e-6), -log10(5e-7)))

tick_labels_updated <- as.character(c(1:20, "", 22, "", 24, "", 26, "", 28, "", 30, "", 32, "", 34, "", 36))

pall <- ggplot(res_all, aes(x = cum_pos, y = -log10(p), colour = chr)) + geom_point(size = 1) + 
  scale_x_continuous(name = "contig", breaks = tick_labels$tick_pos, labels = tick_labels_updated) +
  geom_hline(data = thr_all, aes(yintercept = y), linetype = "dashed", colour = "red", size = 0.5) +
  ggtitle("All SNPs with MAF > 0.1 ") +
  theme_classic(base_size = 16) + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill="grey80", colour = "grey80"),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 12)) + 
  facet_wrap(vars(PC), nrow = 3,
             labeller = labeller(PC = c("1" = "PC1", "2" = "PC2", "3" = "PC3")))

psl <- ggplot(res_sl, aes(x = cum_pos, y = -log10(p), colour = chr)) + geom_point(size = 1) + 
  geom_point(data = res_sl_inv_only, aes(x = cum_pos, y = -log10(p)), pch = 1, size = 1.5, colour = "blue") +
  scale_x_continuous(name = "contig", breaks = tick_labels$tick_pos, labels = tick_labels_updated) +
  geom_hline(data = thr_sl, aes(yintercept = y), linetype = "dashed", colour = "red", size = 0.5) +
  ggtitle("Inversions coded as single loci") +
  theme_classic(base_size = 16) + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill="grey80", colour = "grey80"),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 12)) + 
  facet_wrap(vars(PC), nrow = 3,
             labeller = labeller(PC = c("1" = "PC1", "2" = "PC2", "3" = "PC3")))
p <- pall / psl



ggsave("LFMM_ips_for_pub.png", p, width = 15, height = 10)
