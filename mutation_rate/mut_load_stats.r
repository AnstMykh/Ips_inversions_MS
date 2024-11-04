library("tidyverse")
library('ggplot2')



piN_piS_table <- read.table('path/to/piN_piS_table.txt')
#piN_piS_table_no_north <- piN_piS_table %>% filter(names != 'Inv5' 
 #                                                  & names != 'Inv9' )

piN_piS_table_filtered_only_col <- piN_piS_table  %>% group_by(MJ_MN_COL, names) %>% mutate(win_n = n()) %>% ungroup() %>% filter(MJ_MN_COL == "collinear")
piN_piS_table_filtered_only_col_strict <- piN_piS_table  %>% group_by(MJ_MN_COL, names) %>% mutate(win_n = n()) %>% ungroup() %>% filter(MJ_MN_COL == "collinear") %>%
  filter(contig == 'IpsContig1'
           | contig == 'IpsContig3' 
         | contig == 'IpsContig4' 
          | contig == 'IpsContig8' 
         | contig == 'IpsContig11')


# First step: filtering, changing filtering criteria and rewrinning the ttest for all
piN_piS_cutoff <- 0.3
gene_density_cutoff <- 10
min_hapl_count <- 10
min_win_n <- 4

#!!!!!
cutoff <- '10/10/4'

piN_piS_table_filtered_only_col_filtered <- piN_piS_table_filtered_only_col_strict  %>% group_by(MJ_MN_COL, names) %>% mutate(win_n = n()) %>% ungroup() %>% 
  filter(MJ_MN_COL == "collinear" & gene_density >= gene_density_cutoff)  # & win_n >= min_win_n 

#!!!!!
strict <- 'yes'

# 
 piN_piS_table_filtered_ftest <- piN_piS_table %>% group_by(MJ_MN_COL, names) %>% mutate(win_n = n()) %>% ungroup() %>% 
   filter(MJ_MN_COL != "collinear" & gene_density >= gene_density_cutoff & gen_fr >= min_hapl_count & win_n >= min_win_n) %>% 
   filter(names != 'Inv5' & names != 'Inv9' )
 


piN_piS_table_filtered_ftest_arranged <- piN_piS_table_filtered_ftest  %>% arrange(win_n)
names_vec <-  unique(piN_piS_table_filtered_ftest_arranged$names)


 
piN_piS_table_filtered_for_col_dist <- piN_piS_table %>% group_by(MJ_MN_COL, names) %>% mutate(win_n = n()) %>% ungroup() %>% 
  filter(MJ_MN_COL != "collinear" & (gen_fr < min_hapl_count | win_n < min_win_n))  

piN_piS_table_filtered_for_col_dist %>% filter(gene_density > gene_density_cutoff) -> piN_piS_table_filtered_for_col_dist


# t-test: col to all MIN haplotypes:



#### stats ####
stats <- tibble(strict = factor(), cutoff = numeric(), mdn_COL = numeric(), group = factor(), 
                mdn_group=numeric(), p_val_ttest= numeric(), p_val_M_W= numeric())
stats <- stats[1:6,]

pinpis_COL <- piN_piS_table_filtered_only_col_filtered %>% pull(piN_piS)
piN_piS_COL_median <- median(pinpis_COL, na.rm = T)

group <- 'INV'
pinpis_inv <- piN_piS_table_filtered_ftest  %>% pull(piN_piS) 
piN_piS_inv_median <- median(pinpis_inv, na.rm = T)
ttest_INV_COL <- t.test(pinpis_inv, pinpis_COL, alternative = "less")$p.value
Wilcx_INV_COL <- wilcox.test(pinpis_inv, pinpis_COL, alternative = 'less')$p.value
tmp1 <- tibble(strict = strict, cutoff = cutoff, mdn_COL = piN_piS_COL_median, group = group, 
               mdn_group=piN_piS_inv_median, p_val_ttest= ttest_INV_COL, p_val_M_W= Wilcx_INV_COL)

group <- 'MN'
pinpis_MN <- piN_piS_table_filtered_ftest %>% filter(MJ_MN_COL == "MN") %>% pull(piN_piS)
piN_piS_MN_median <- median(pinpis_MN, na.rm = T) 
ttest_MN_COL <- t.test(pinpis_MN, pinpis_COL, alternative = "less")$p.value
Wilcx_MN_COL <- wilcox.test(pinpis_MN, pinpis_COL, alternative = 'less')$p.value
tmp2 <- tibble(strict = strict, cutoff = cutoff, mdn_COL = piN_piS_COL_median, group = group, 
               mdn_group=piN_piS_MN_median, p_val_ttest= ttest_MN_COL, p_val_M_W= Wilcx_MN_COL) 

group <- 'MJ'
pinpis_MJ <- piN_piS_table_filtered_ftest %>% filter(MJ_MN_COL == "MJ") %>% pull(piN_piS)
piN_piS_MJ_median <- median(pinpis_MJ, na.rm = T)
ttest_MJ_COL <- t.test(pinpis_MJ, pinpis_COL, alternative = "less")$p.value
Wilcx_MJ_COL <- wilcox.test(pinpis_MJ, pinpis_COL, alternative = 'less')$p.value
tmp3 <- tibble(strict = strict, cutoff = cutoff, mdn_COL = piN_piS_COL_median, group = group, 
               mdn_group=piN_piS_MJ_median, p_val_ttest= ttest_MJ_COL, p_val_M_W= Wilcx_MJ_COL)


stats <- rbind(stats, tmp1, tmp2, tmp3)

write.table(stats, 'path/to/piN_piS_stats.txt', row.names = TRUE)

#### ADDITIONAL SCRIPT: exprolative statistics 

####PREPARING DISTRIBUTIONS TO PLOT + balanced/imbalanced ####
balanced <- c("Inv10","Inv17", "Inv18","Inv23.1","Inv26","Inv7.1","Inv14.3","Inv14.5","Inv14.1","Inv22.1","Inv22.4","Inv22.2","Inv22.3","Inv16")



piN_piS_table_filtered_only_col_plot <- piN_piS_table_filtered_only_col_filtered %>% 
  select(names,contig, start, end, inversion, MJ_MN_COL, piN_piS)  %>% mutate(distribution = 'collinear') %>% mutate(bal_imb ="collinear")


haplotypes_on_col_dist_merge <- piN_piS_table_filtered_for_col_dist %>%
   select(names,contig, start, end, inversion, MJ_MN_COL, piN_piS) %>%
  mutate(distribution = 'invertion') %>% mutate(bal_imb = ifelse(names %in% balanced, 'balanced', 'imbalanced')) 

piN_piS_table_inv_plot <- piN_piS_table_filtered_ftest_arranged %>% filter(inversion != 'collinear') %>%
  select(names,contig, start, end, inversion, MJ_MN_COL, piN_piS) %>% mutate(distribution = 'invertion') %>%
  mutate(bal_imb = ifelse(names %in% balanced, 'balanced', 'imbalanced')) 

piNpiS_dist_to_plot <- rbind(piN_piS_table_filtered_only_col_plot, piN_piS_table_inv_plot, haplotypes_on_col_dist_merge)
###OR!!!!!
piNpiS_dist_to_plot <- rbind(piN_piS_table_filtered_only_col_plot, piN_piS_table_inv_plot)


#### stats balanced/imbalanced  ####
stats <- tibble(strict = factor(), cutoff = numeric(), mdn_COL = numeric(), group = factor(), 
                mdn_group=numeric(), p_val_ttest= numeric(), p_val_M_W= numeric())


pinpis_COL <- piN_piS_table_filtered_only_col_filtered %>% pull(piN_piS)
piN_piS_COL_median <- median(pinpis_COL, na.rm = T)

group <- 'BAL_INV'
pinpis_bal <- piNpiS_dist_to_plot %>% filter(bal_imb == "balanced") %>% pull(piN_piS)   #length(pinpis_bal)
#pinpis_bal <- pinpis_bal[pinpis_bal < piN_piS_cutoff]
piN_piS_bal_median <- median(pinpis_bal, na.rm = T)
ttest_BAL_COL <- t.test(pinpis_bal, pinpis_COL, alternative = "less")$p.value
Wilcx_BAL_COL <- wilcox.test(pinpis_bal, pinpis_COL, alternative = 'less')$p.value
tmp1 <- tibble(strict = strict, cutoff = cutoff, mdn_COL = piN_piS_COL_median, group = group, 
               mdn_group=piN_piS_bal_median, p_val_ttest= ttest_BAL_COL, p_val_M_W= Wilcx_BAL_COL)

group <- 'BAL_MJ'
piN_piS_bal_MJ <- piNpiS_dist_to_plot %>% filter(bal_imb == "balanced") %>% filter(MJ_MN_COL == "MJ") %>% pull(piN_piS)
#piN_piS_bal_MJ <- piN_piS_bal_MJ[piN_piS_bal_MJ < piN_piS_cutoff] 
piN_piS_bal_MJ_MEDIAN <- median(piN_piS_bal_MJ, na.rm = T)
ttest_BALmj_COL <- t.test(piN_piS_bal_MJ, pinpis_COL, alternative = "less")$p.value
Wilcx_BALmj_COL <- wilcox.test(piN_piS_bal_MJ, pinpis_COL, alternative = 'less')$p.value
tmp2 <- tibble(strict = strict, cutoff = cutoff, mdn_COL = piN_piS_COL_median, group = group, 
               mdn_group=piN_piS_bal_MJ_MEDIAN, p_val_ttest= ttest_BALmj_COL, p_val_M_W= Wilcx_BALmj_COL)


group <- 'BAL_MN'
piN_piS_bal_MN <- piNpiS_dist_to_plot %>% filter(bal_imb == "balanced") %>% filter(MJ_MN_COL == "MN") %>% pull(piN_piS)
#piN_piS_bal_MN <- piN_piS_bal_MN[piN_piS_bal_MN < piN_piS_cutoff] 
piN_piS_bal_MN_MEDIAN <- median(piN_piS_bal_MN, na.rm = T)
ttest_BALmn_COL <- t.test(piN_piS_bal_MN, pinpis_COL, alternative = "less")$p.value
Wilcx_INVmn_COL <- wilcox.test(piN_piS_bal_MN, pinpis_COL, alternative = 'less')$p.value
tmp3 <- tibble(strict = strict, cutoff = cutoff, mdn_COL = piN_piS_COL_median, group = group, 
               mdn_group=piN_piS_bal_MN_MEDIAN, p_val_ttest= ttest_BALmn_COL, p_val_M_W= Wilcx_INVmn_COL)


group <- 'imBAL_INV'
pinpis_imbal <- piNpiS_dist_to_plot %>% filter(bal_imb == "imbalanced") %>% pull(piN_piS)   #length(pinpis_bal)
#pinpis_bal <- pinpis_bal[pinpis_bal < piN_piS_cutoff]
piN_piS_imbal_median <- median(pinpis_imbal, na.rm = T)
ttest_imBAL_COL <- t.test(pinpis_imbal, pinpis_COL, alternative = "less")$p.value
Wilcx_imBAL_COL <- wilcox.test(pinpis_imbal, pinpis_COL, alternative = 'less')$p.value
tmp4 <- tibble(strict = strict, cutoff = cutoff, mdn_COL = piN_piS_COL_median, group = group, 
               mdn_group=piN_piS_imbal_median, p_val_ttest= ttest_imBAL_COL, p_val_M_W= Wilcx_imBAL_COL)

group <- 'imBAL_MJ'
piN_piS_imbal_MJ <- piNpiS_dist_to_plot %>% filter(bal_imb == "imbalanced") %>% filter(MJ_MN_COL == "MJ") %>% pull(piN_piS)
#piN_piS_bal_MJ <- piN_piS_bal_MJ[piN_piS_bal_MJ < piN_piS_cutoff] 
piN_piS_imbal_MJ_MEDIAN <- median(piN_piS_imbal_MJ, na.rm = T)
ttest_imBALmj_COL <- t.test(piN_piS_imbal_MJ, pinpis_COL, alternative = "less")$p.value
Wilcx_imBALmj_COL <- wilcox.test(piN_piS_imbal_MJ, pinpis_COL, alternative = 'less')$p.value
tmp5 <- tibble(strict = strict, cutoff = cutoff, mdn_COL = piN_piS_COL_median, group = group, 
               mdn_group=piN_piS_imbal_MJ_MEDIAN, p_val_ttest= ttest_imBALmj_COL, p_val_M_W= Wilcx_imBALmj_COL)


group <- 'imBAL_MN'
piN_piS_imbal_MN <- piNpiS_dist_to_plot %>% filter(bal_imb == "imbalanced") %>% filter(MJ_MN_COL == "MN") %>% pull(piN_piS)
#piN_piS_bal_MN <- piN_piS_bal_MN[piN_piS_bal_MN < piN_piS_cutoff] 
piN_piS_imbal_MN_MEDIAN <- median(piN_piS_imbal_MN, na.rm = T)
ttest_imBALmn_COL <- t.test(piN_piS_imbal_MN, pinpis_COL, alternative = "less")$p.value
Wilcx_imINVmn_COL <- wilcox.test(piN_piS_imbal_MN, pinpis_COL, alternative = 'less')$p.value
tmp6 <- tibble(strict = strict, cutoff = cutoff, mdn_COL = piN_piS_COL_median, group = group, 
               mdn_group=piN_piS_imbal_MN_MEDIAN, p_val_ttest= ttest_imBALmn_COL, p_val_M_W= Wilcx_imINVmn_COL)


stats <- rbind(stats, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6)

write.table(stats, '/home/anstmykh/Desktop/mut_load/piN_piS_stats_bal_imb.txt', row.names = TRUE)

######permutations#####


n_permutations <- 10000
p_val_perm <- c()

random_v <-c(rep(1, nrow(piN_piS_table_filtered_only_col_plot)), rep(2, nrow(piN_piS_table_filtered_ftest_plot)))


for (i in 1:n_permutations){
  test <- piNpiS_dist_to_plot %>% mutate(random = sample(random_v))
  vec_1 <- test %>% filter(random == '1') %>% pull(piN_piS)
  vec_2 <- test %>% filter(random == '2') %>% pull(piN_piS)
  ttest <- t.test(vec_2, vec_1, alternative = "less")$p.value
  p_val_perm <- c(p_val_perm, ttest)
}
  
  
perc_sign <- length(p_val_perm[p_val_perm < 0.05]) *100/10000



####  PLOTTING DISTRIBUTIONS  ######
piNpiS_dist_to_plot[piNpiS_dist_to_plot == "INT"] <- "MN"  # AS AN OPTION! 
piNpiS_dist_to_plot_cutoff <- piNpiS_dist_to_plot %>% filter(piN_piS < 0.5)

#counts
COL_VS_INV_counts <- ggplot(piNpiS_dist_to_plot) +  
  geom_histogram(aes(x=piN_piS, fill = distribution , alpha = distribution), binwidth = 0.01, position = position_identity()) +
  scale_alpha_manual(values = c( 1, 0.6))

#FREQ!!!!
COL_VS_INV_fr <- ggplot(piNpiS_dist_to_plot_cutoff) +  
  geom_histogram(aes(x=piN_piS, y= after_stat(density*width),fill = distribution , alpha = distribution), binwidth = 0.01, position = position_identity()) +
  scale_alpha_manual(values = c( 1, 0.6)) +    ## IMPORTANT PIECE OF CODE!! 
  scale_fill_discrete(name = '') +
  guides(alpha = 'none') +  
  ylab("frequency") + xlab("piN/piS") +
  ggtitle(paste("piN/piS cutoff -", 0.5 )) +
  theme_bw() +
  theme(axis.text =element_text(size=25), axis.title = element_text(size=25), legend.text = element_text(size=30), legend.position = c(0.8, 0.85)) 


jpeg("/home/anstmykh/Desktop/mut_load/mut_load_col_inv_0.5.jpeg", width = 1000, height = 600)
print(COL_VS_INV_fr)
dev.off()


# freq, MJ, Min
COL_VS_INV_fr <- ggplot(piNpiS_dist_to_plot_cutoff) +  
  geom_histogram(aes(x=piN_piS, y= after_stat(density*width),fill = MJ_MN_COL , alpha = MJ_MN_COL), binwidth = 0.01, position = position_identity()) +
  scale_alpha_manual(values = c( 0.4, 0.2, 0.2)) 

###MJ
COL_VS_INV_fr <- ggplot() +  
  geom_histogram(aes(x=pinpis_COL, y= after_stat(density*width), fill = 'col', alpha = 'col'), binwidth = 0.01) +
  geom_histogram(aes(x=pinpis_MJ, y= after_stat(density*width), fill = 'mj', alpha = 'mj'), binwidth = 0.01) +
  scale_alpha_manual(values = c('col'= 1, 'mj' = 0.6)) 

####MN
COL_VS_INV_fr <- ggplot() +  
  geom_histogram(aes(x=pinpis_COL, y= after_stat(density*width), fill = 'col', alpha = 'col'), binwidth = 0.01) +
  geom_histogram(aes(x=pinpis_MN, y= after_stat(density*width), fill = 'mn', alpha = 'mn'), binwidth = 0.01) +
  scale_alpha_manual(values = c('col'= 1, 'mn' = 0.6)) 


# freq, balanced\imbalanced
COL_VS_INV_fr <- ggplot(piNpiS_dist_to_plot_cutoff) +  
  geom_histogram(aes(x=piN_piS, y= after_stat(density*width),fill = bal_imb , alpha = bal_imb), binwidth = 0.01, position = position_identity()) +
  scale_alpha_manual(values = c( 0.4, 0.2, 0.2)) 


##### balanced ######
pinpis_bal <- piNpiS_dist_to_plot %>% filter(bal_imb == "balanced") %>% pull(piN_piS)   #length(pinpis_bal)
#pinpis_bal <- pinpis_bal[pinpis_bal < piN_piS_cutoff]
piN_piS_bal_median <- median(pinpis_bal, na.rm = T)

piN_piS_bal_MJ <- piNpiS_dist_to_plot %>% filter(bal_imb == "balanced") %>% filter(MJ_MN_COL == "MJ") %>% pull(piN_piS)
#piN_piS_bal_MJ <- piN_piS_bal_MJ[piN_piS_bal_MJ < piN_piS_cutoff] 
piN_piS_bal_MJ_MEDIAN <- median(piN_piS_bal_MJ, na.rm = T)

piN_piS_bal_MN <- piNpiS_dist_to_plot %>% filter(bal_imb == "balanced") %>% filter(MJ_MN_COL == "MN") %>% pull(piN_piS)
#piN_piS_bal_MN <- piN_piS_bal_MN[piN_piS_bal_MN < piN_piS_cutoff] 
piN_piS_bal_MJ_MEDIAN <- median(piN_piS_bal_MN, na.rm = T)


COL_VS_INV_fr <- ggplot() +  
  geom_histogram(aes(x=pinpis_COL, y= after_stat(density*width), fill = 'col', alpha = 'col'), binwidth = 0.01) +
  geom_histogram(aes(x=pinpis_bal, y= after_stat(density*width), fill = 'balanced', alpha = 'balanced'), binwidth = 0.01) +
  scale_alpha_manual(values = c('col'= 1, 'balanced' = 0.6)) 

ttest_BAL_COL <- t.test(pinpis_bal, pinpis_COL, alternative = "less")$p.value  ## significant! 0.0009058513 when strict / when not is 0.0004945551
ttest_BALmj_COL <- t.test(piN_piS_bal_MJ, pinpis_COL, alternative = "less")$p.value ## 0.001290495 when strict / 0.0009994226 when not 
ttest_BALmn_COL <- t.test(piN_piS_bal_MN, pinpis_COL, alternative = "less")$p.value ## 0.02426695 sing when strict / 0.02200701 when not 

######  imbalanced #####
pinpis_imbal <- piNpiS_dist_to_plot %>% filter(bal_imb == "imbalanced") %>% pull(piN_piS)    # length(pinpis_imbal)
#pinpis_imbal <- pinpis_imbal[pinpis_imbal < piN_piS_cutoff]
piN_piS_imbal <- median(pinpis_imbal, na.rm = T)

piN_piS_imbal_MJ <- piNpiS_dist_to_plot %>% filter(bal_imb == "imbalanced") %>% filter(MJ_MN_COL == "MJ") %>% pull(piN_piS)
#piN_piS_imbal_MJ <- piN_piS_imbal_MJ[piN_piS_imbal_MJ < piN_piS_cutoff] 
piN_piS_imbal_MJ_MEDIAN <- median(piN_piS_imbal_MJ, na.rm = T)

piN_piS_imbal_MN <- piNpiS_dist_to_plot %>% filter(bal_imb == "imbalanced") %>% filter(MJ_MN_COL == "MN") %>% pull(piN_piS)
#piN_piS_imbal_MN <- piN_piS_imbal_MN[piN_piS_imbal_MN < piN_piS_cutoff] 
piN_piS_imbal_MN_MEDIAN <- median(piN_piS_imbal_MN, na.rm = T)


COL_VS_INV_fr <- ggplot() +  
  geom_histogram(aes(x=pinpis_COL, y= after_stat(density*width), fill = 'col', alpha = 'col'), binwidth = 0.01) +
  geom_histogram(aes(x=pinpis_imbal, y= after_stat(density*width), fill = 'imbalanced', alpha = 'imbalanced'), binwidth = 0.01) +
  scale_alpha_manual(values = c('col'= 1, 'imbalanced' = 0.6)) 



ttest_IMBAL_COL <- t.test(pinpis_imbal, pinpis_COL, alternative = "less")$p.value  #  (0.03071521 when strict / 0.0004181822 when not) : 5,4,4
                                                                                   #  (0.009358752 when strict / 0.0001522373 when not) : 10,10,4

ttest_IMBALmj_COL <- t.test(piN_piS_imbal_MJ, pinpis_COL, alternative = "less")$p.value #  (5.100997e-08 when strict / 0.0001248394 when not) : 5,4,4
                                                                                        #  (0.000132558 when strict / 5.100997e-08 when not) : 10,10,4

ttest_IMBALmn_COL <- t.test(piN_piS_imbal_MN, pinpis_COL, alternative = "less")$p.value #  (0.131002 when strict / 0.03128847 when not) : 5,4,4
                                                                                        #  (0.5360315 when strict / 0.4132952 when not) : 10,10,4

jpeg("mut_load_col_IMBAL_DIST.jpeg", width = 1000, height = 600)
print(COL_VS_INV_fr)
dev.off()


####final boxplots with statistics ####
piN_piS_table_no_north <- piN_piS_table_no_north %>% mutate(names = replace(names, names == "collinear", "COL"))   # filter(piN_piS_table_no_north, piN_piS < 5)
piN_piS_table_no_north_no_6_17 <-piN_piS_table_no_north %>% filter(names != 'Inv6' 
                                                                   & names != 'Inv17' )
###
piN_piS_table_filtered_ftest_arranged <- piN_piS_table_filtered_ftest  %>% arrange(win_n)

balanced <- c("Inv10","Inv17", "Inv18","Inv23.1","Inv26","Inv7.1","Inv14.3","Inv14.5","Inv14.1","Inv22.1","Inv22.4","Inv22.2","Inv22.3","Inv16")

piN_piS_table_filtered_only_col_plot <- piN_piS_table_filtered_only_col_filtered %>% 
  select(names,contig, start, end, inversion, MJ_MN_COL, piN_piS)  %>% mutate(distribution = 'collinear') %>% mutate(bal_imb ="collinear")

piN_piS_table_inv_plot <- piN_piS_table_filtered_ftest_arranged %>% filter(inversion != 'collinear') %>%
  select(names,contig, start, end, inversion, MJ_MN_COL, piN_piS) %>% mutate(distribution = 'invertion') %>%
  mutate(bal_imb = ifelse(names %in% balanced, 'balanced', 'imbalanced')) 

piNpiS_dist_to_plot <- rbind(piN_piS_table_filtered_only_col_plot, piN_piS_table_inv_plot)
####


piN_piS_table_filtered_ftest_COL <- piNpiS_dist_to_plot %>%
  mutate(names = replace(names, names == "collinear", "COL")) %>% 
  mutate(names = factor(names, levels = unique(names)))


boxplot <- ggplot(filter(piN_piS_table_filtered_ftest_COL, piN_piS < 2 ) , aes(x=MJ_MN_COL, y=piN_piS, fill = MJ_MN_COL)) +  #x=MJ_MN_COL tp order 
  geom_boxplot(na.rm = T) + 
  theme(text = element_text(size = 20)) + 
  facet_grid(~names, scale="free", space = 'free') +
  scale_fill_discrete(name = '') +
  #geom_text(aes(label = gen_fr), check_overlap = T, na.rm = T) +
  #ggtitle(paste("cutoffs: piN/piS -", piN_piS_cutoff, "   ","gene density -", gene_density_cutoff, 
   #             "   ","min_hapl_count -",min_hapl_count, "   ","min_win_n -", min_win_n ,sep = ' ')) +
  ggtitle(paste("cutoffs: piN/piS -", 2, '; gene.dens/min.hapl/min.win-', cutoff )) +
  theme_bw()

jpeg("/home/anstmykh/Desktop/mut_load/mut_load_boxplot_all_genotypes_filtered_codon_arranged_strict_filtering.jpeg", width = 1300, height = 600)
print(boxplot)
dev.off()


#####boxplots MAJ/MIN


boxplot <- ggplot(filter(piN_piS_table_filtered_ftest_COL) , aes(x=MJ_MN_COL, y=piN_piS, fill = MJ_MN_COL)) +  #x=MJ_MN_COL tp order 
  geom_boxplot(na.rm = T) + 
  theme(text = element_text(size = 20)) + 
  #facet_grid(~names, scale="free", space = 'free') +
  scale_fill_discrete(name = '') +
  ggtitle(paste("cutoffs: piN/piS -", 'no', '; gene.dens/min.hapl/min.win-', cutoff )) +
  theme_bw()

jpeg("/home/anstmykh/Desktop/mut_load/mut_load_boxplot_col_vs_MJ_vs_.jpeg", width = 1300, height = 600)
print(boxplot)
dev.off()






######comparing haplotypes ####


names <- read.table("/home/anstmykh/Desktop/mut_load/names.txt", header=F, skip=0, sep="\t")

MN_MJ_R_2s_all_filtered <- tibble(names=character(), contig=character(), start=numeric(),
                             end=numeric(), inversion=character(),piN_piS_MJ_median=numeric(),
                             piN_piS_MN_median=numeric(), piN_piS_R_median=numeric(),p_value_ttest=numeric(), p_value_MW=numeric())



i <- 11
##this is the whole loop
for (i in 1:nrow(names)) {
  
  contig_ref <- names[i,2]
  start_ref <-names[i,3]
  end_ref <- names[i,4]
  inversion_ref <- names[i,5]
  name <- names[i,1]
  
  pinpis_tmp <- piN_piS_table_filtered_ftest %>% filter(contig == contig_ref & start == start_ref & end == end_ref & 
                                                          inversion == inversion_ref)
  
  
  ###cheching for the presense of both alleles
  check_for_allele_n <- unique(pinpis_tmp$MJ_MN_COL) %>% length()
  
  if(check_for_allele_n < 2){
    next;
    ##########
  }
  
  
  pinpis_MJ <- pinpis_tmp %>% filter(MJ_MN_COL == "MJ")  %>% pull(piN_piS)
  pinpis_MN <- pinpis_tmp %>% filter(MJ_MN_COL == "MN")  %>% pull(piN_piS)
  pinpis_R <- pinpis_tmp %>% filter(MJ_MN_COL == "R") %>% pull(piN_piS)
  
  piN_piS_median_MJ <- median(pinpis_MJ, na.rm = T)
  piN_piS_median_MN <- median(pinpis_MN, na.rm = T)
  piN_piS_median_R <- median(pinpis_R, na.rm = T)
  
  
  ### 3 haplotypes ####
  if(length(pinpis_R) != 0) {
    
    
    pinpis_tmp_MJ_MN <- piN_piS_table_filtered_ftest %>% filter(contig == contig_ref & start == start_ref & end == end_ref & 
                                                                  inversion == inversion_ref) %>% filter(MJ_MN_COL == "MJ" | MJ_MN_COL == "MN")
    
    pinpis_tmp_R_MN <- piN_piS_table_filtered_ftest %>% filter(contig == contig_ref & start == start_ref & end == end_ref & 
                                                                 inversion == inversion_ref) %>% filter(MJ_MN_COL == "R" | MJ_MN_COL == "MN")
    
    pinpis_tmp_MJ_R <- piN_piS_table_filtered_ftest %>% filter(contig == contig_ref & start == start_ref & end == end_ref & 
                                                                 inversion == inversion_ref) %>% filter(MJ_MN_COL == "MJ" | MJ_MN_COL == "R")
    
    
    ## f-test between MJ and MN alleles
    if(length(unique(pinpis_tmp_MJ_MN$MJ_MN_COL)) > 1){
      
      
      ttest1 <- t.test(pinpis_MJ,pinpis_MN, alternative = "two.sided")$p.value
      non_param_pval_1 <- wilcox.test(piN_piS ~ MJ_MN_COL, data=pinpis_tmp_MJ_MN, alternative = "two.sided")$p.value
      
      tmp1 <- tibble(names=name, contig=contig_ref, start=start_ref, end=end_ref, inversion=inversion_ref,
                     piN_piS_MJ_median=piN_piS_median_MJ, piN_piS_MN_median=piN_piS_median_MN, piN_piS_R_median=NA, 
                     p_value_ttest=ttest1,p_value_MW=non_param_pval_1)
      
      MN_MJ_R_2s_all_filtered <- rbind(MN_MJ_R_2s_all_filtered, tmp1)
      
    }
    ## f-test between MJ and R alleles
    if(length(unique(pinpis_tmp_MJ_R$MJ_MN_COL)) > 1){
      
      
      ttest2 <- t.test(pinpis_MJ,pinpis_R, alternative = "two.sided")$p.value
      non_param_pval_2 <- wilcox.test(piN_piS ~ MJ_MN_COL, data=pinpis_tmp_MJ_R, alternative = "two.sided")$p.value
      
      tmp2 <- tibble(names=name, contig=contig_ref, start=start_ref, end=end_ref, inversion=inversion_ref, 
                     piN_piS_MJ_median=piN_piS_median_MJ, piN_piS_MN_median=NA, piN_piS_R_median=piN_piS_median_R,
                     p_value_ttest=ttest2,p_value_MW=non_param_pval_2)
      
      MN_MJ_R_2s_all_filtered <- rbind(MN_MJ_R_2s_all_filtered, tmp2)
      
    }
    ## f-test between R and MN alleles
    if(length(unique(pinpis_tmp_R_MN$MJ_MN_COL)) > 1){
      
      ttest3 <- t.test(pinpis_R,pinpis_MN, alternative = "two.sided")$p.value
      non_param_pval_3 <- wilcox.test(piN_piS ~ MJ_MN_COL, data=pinpis_tmp_R_MN, alternative = "two.sided")$p.value
      
      tmp3 <- tibble(names=name, contig=contig_ref, start=start_ref, end=end_ref, inversion=inversion_ref, 
                     piN_piS_MJ_median=NA, piN_piS_MN_median=piN_piS_median_MN, piN_piS_R_median=piN_piS_median_R,
                     p_value_ttest=ttest3,p_value_MW=non_param_pval_3)
      
      MN_MJ_R_2s_all_filtered <- rbind(MN_MJ_R_2s_all_filtered, tmp3)
      
      
    }
    ### 2 haplotypes  ####
  } else {
    
    
    ttest <- t.test(pinpis_MJ,pinpis_MN, alternative = "two.sided")$p.value
    
    non_param_pval <- wilcox.test(piN_piS ~ MJ_MN_COL, data=pinpis_tmp, alternative = "two.sided")$p.value
    
    tmp <- tibble(names=name, contig=contig_ref, start=start_ref, end=end_ref, inversion=inversion_ref, 
                  piN_piS_MJ_median=piN_piS_median_MJ, piN_piS_MN_median=piN_piS_median_MN, piN_piS_R_median=NA, 
                  p_value_ttest=ttest, p_value_MW=non_param_pval )
    
    MN_MJ_R_2s_all_filtered <- rbind(MN_MJ_R_2s_all_filtered, tmp)
    
    
    
    
  }
  
}

write.table(MN_MJ_R_2s_all_filtered, '/home/anstmykh/Desktop/mut_load/hapl_comp_all_stats.txt', row.names = TRUE)

##testing the loop

unique(piN_piS_table_filtered_ftest$names)
unique(MN_MJ_R_2s_all_filtered$names)


setdiff(unique(piN_piS_table_filtered_ftest$names), c(unique(MN_MJ_INT_mann_whit_wilc$names), unique(MN_MJ_INT_2s_ttest$names)))



####### 2s t-test and MW on everything
MN_MJ_R_2s_ttest_all_filtered <- tibble(names=character(), contig=character(), start=numeric(),
                                          end=numeric(), inversion=character(),piN_piS_MJ_median=numeric(),
                                          piN_piS_MN_median=numeric(), piN_piS_R_median=numeric(),
                                        p_value_ttest=numeric(), p_value_MW=numeric())
i <- 29
for (i in 1:nrow(names)) {
  
  contig_ref <- names[i,2]
  start_ref <-names[i,3]
  end_ref <- names[i,4]
  inversion_ref <- names[i,5]
  name <- names[i,1]
  
  pinpis_tmp <- piN_piS_table_filtered_ftest %>% filter(contig == contig_ref & start == start_ref & end == end_ref & 
                                                          inversion == inversion_ref)
  
  
  ###cheching for the presense of both alleles
  check_for_allele_n <- unique(pinpis_tmp$MJ_MN_COL) %>% length()
  
  if(check_for_allele_n < 2){
    next;
    ##########
  }
  
  
  pinpis_MJ <- pinpis_tmp %>% filter(MJ_MN_COL == "MJ")  %>% pull(piN_piS)
  pinpis_MN <- pinpis_tmp %>% filter(MJ_MN_COL == "MN")  %>% pull(piN_piS)
  pinpis_INT <- pinpis_tmp %>% filter(MJ_MN_COL == "INT") %>% pull(piN_piS)
  
  piN_piS_median_MJ <- median(pinpis_MJ, na.rm = T)
  piN_piS_median_MN <- median(pinpis_MN, na.rm = T)
  piN_piS_median_INT <- median(pinpis_INT, na.rm = T)
  
  
  ### 3 haplotypes ####
  if(length(pinpis_INT) != 0) {
    
    
    pinpis_tmp_MJ_MN <- piN_piS_table_filtered_ftest %>% filter(contig == contig_ref & start == start_ref & end == end_ref & 
                                                                  inversion == inversion_ref) %>% filter(MJ_MN_COL == "MJ" | MJ_MN_COL == "MN")
    
    pinpis_tmp_INT_MN <- piN_piS_table_filtered_ftest %>% filter(contig == contig_ref & start == start_ref & end == end_ref & 
                                                                   inversion == inversion_ref) %>% filter(MJ_MN_COL == "INT" | MJ_MN_COL == "MN")
    
    pinpis_tmp_MJ_INT <- piN_piS_table_filtered_ftest %>% filter(contig == contig_ref & start == start_ref & end == end_ref & 
                                                                   inversion == inversion_ref) %>% filter(MJ_MN_COL == "MJ" | MJ_MN_COL == "INT")
    
    
    ## f-test between MJ and MN alleles
    if(length(unique(pinpis_tmp_MJ_MN$MJ_MN_COL)) > 1){
      ## now performing t-test  
      ttest1 <- t.test(pinpis_MJ,pinpis_MN, alternative = "two.sided")$p.value
      
      tmp1.1 <- tibble(names=name, contig=contig_ref, start=start_ref, end=end_ref, inversion=inversion_ref, 
                       piN_piS_MJ_median=piN_piS_median_MJ, piN_piS_MN_median=piN_piS_median_MN, piN_piS_INT_median=NA, p_value=ttest1 )
      
      MN_MJ_INT_2s_ttest_all_filtered <- rbind(MN_MJ_INT_2s_ttest_all_filtered, tmp1.1)
      
    }
    ## f-test between MJ and INT alleles
    if(length(unique(pinpis_tmp_MJ_INT$MJ_MN_COL)) > 1){
      ## now performing t-test 
      ttest2 <- t.test(pinpis_MJ,pinpis_INT, alternative = "two.sided")$p.value
      
      tmp2.1 <- tibble(names=name, contig=contig_ref, start=start_ref, end=end_ref, inversion=inversion_ref, 
                       piN_piS_MJ_median=piN_piS_median_MJ, piN_piS_MN_median=NA, piN_piS_INT_median=piN_piS_median_INT, p_value=ttest2 )
      
      MN_MJ_INT_2s_ttest_all_filtered <- rbind(MN_MJ_INT_2s_ttest_all_filtered, tmp2.1)
    }
    ## f-test between INT and MN alleles
    if(length(unique(pinpis_tmp_INT_MN$MJ_MN_COL)) > 1){
      ## now performing t-test 
      ttest3 <- t.test(pinpis_INT,pinpis_MN, alternative = "two.sided")$p.value
      
      tmp3.1 <- tibble(names=name, contig=contig_ref, start=start_ref, end=end_ref, inversion=inversion_ref, 
                       piN_piS_MJ_median=NA, piN_piS_MN_median=piN_piS_median_MN, piN_piS_INT_median=piN_piS_median_INT, p_value=ttest3 )
      
      MN_MJ_INT_2s_ttest_all_filtered <- rbind(MN_MJ_INT_2s_ttest_all_filtered, tmp3.1)
    }
    ### 2 haplotypes  ####
  } else {
    
    ttest1 <- t.test(pinpis_MJ,pinpis_MN, alternative = "two.sided")$p.value
    
    tmp1.1 <- tibble(names=name, contig=contig_ref, start=start_ref, end=end_ref, inversion=inversion_ref, 
                     piN_piS_MJ_median=piN_piS_median_MJ, piN_piS_MN_median=piN_piS_median_MN, piN_piS_INT_median=NA, p_value=ttest1 )
    
    MN_MJ_INT_2s_ttest_all_filtered <- rbind(MN_MJ_INT_2s_ttest_all_filtered, tmp1.1)
    
  }
  
}


######comparing haplotypes with collinear sample ####
names <- read.table("/home/anstmykh/Desktop/mut_load/names.txt", header=F, skip=0, sep="\t")
MN_MJ_R_COL_1s_ttest <- tibble(names=character(), contig=character(), start=numeric(),
                                 end=numeric(), inversion=character(), MJ_MN_COL=character(), piN_piS_hapl_median=numeric(),
                                 piN_piS_COL_median=numeric() , 
                               p_value_ttest_gr=numeric(), p_value_ttest_ls=numeric(), 
                               p_value_MW_gr=numeric(), p_value_MW_ls=numeric())



i <- 11


pinpis_COL <- piN_piS_table_filtered_only_col_filtered %>% filter(MJ_MN_COL == "collinear") %>% pull(piN_piS)
piN_piS_COL_median <- median(pinpis_COL, na.rm = T)


for (i in 1:nrow(names)) {
  
  contig_ref <- names[i,2]
  start_ref <-names[i,3]
  end_ref <- names[i,4]
  inversion_ref <- names[i,5]
  name <- names[i,1]
  
  pinpis_tmp <- piN_piS_table_filtered_ftest %>% filter(contig == contig_ref & start == start_ref & end == end_ref & 
                                                          inversion == inversion_ref)
  
  joint_tibble <- rbind(piN_piS_table_filtered_only_col_filtered, pinpis_tmp)
  
  
  
  pinpis_MJ <- pinpis_tmp %>% filter(MJ_MN_COL == "MJ")  %>% pull(piN_piS)
  pinpis_MN <- pinpis_tmp %>% filter(MJ_MN_COL == "MN")  %>% pull(piN_piS)
  pinpis_R <- pinpis_tmp %>% filter(MJ_MN_COL == "R") %>% pull(piN_piS)
  
  piN_piS_median_MJ <- median(pinpis_MJ, na.rm = T)
  piN_piS_median_MN <- median(pinpis_MN, na.rm = T)
  piN_piS_median_R <- median(pinpis_R, na.rm = T)
  
  
  ### 3 haplotypes ####
  if(length(pinpis_R) != 0) {
    
    
    pinpis_tmp_MJ_COL <- joint_tibble %>% filter(MJ_MN_COL == "MJ" | MJ_MN_COL == "collinear")
    
    pinpis_tmp_R_COL <- joint_tibble %>% filter(MJ_MN_COL == "R" | MJ_MN_COL == "collinear")
    
    pinpis_tmp_MN_COL <- joint_tibble %>% filter(MJ_MN_COL == "MN" | MJ_MN_COL == "collinear")
    
    
    ## f-test between MJ and collinear
    if(length(unique(pinpis_tmp_MJ_COL$MJ_MN_COL)) == 2){
     

        ttest1g <- t.test(pinpis_MJ, pinpis_COL, alternative = "greater")$p.value
        ttest1l <- t.test(pinpis_MJ, pinpis_COL, alternative = "less")$p.value
        non_param_pval_1g <- wilcox.test(pinpis_MJ, pinpis_COL, data=pinpis_tmp_MJ_COL,  alternative = "greater")$p.value
        non_param_pval_1l <- wilcox.test(pinpis_MJ, pinpis_COL, data=pinpis_tmp_MJ_COL,  alternative = "less")$p.value
        
        tmp1 <- tibble(names=name, contig=contig_ref, start=start_ref, end=end_ref, inversion=inversion_ref, 
                         MJ_MN_COL='MJ', piN_piS_hapl_median=piN_piS_median_MJ,
                         piN_piS_COL_median= piN_piS_COL_median,
                       p_value_ttest_gr=ttest1g, p_value_ttest_ls=ttest1l, 
                       p_value_MW_gr=non_param_pval_1g, p_value_MW_ls=non_param_pval_1l)
        
        MN_MJ_R_COL_1s_ttest <- rbind(MN_MJ_R_COL_1s_ttest, tmp1)

        
    }
    ## f-test between R and collinear
    if(length(unique(pinpis_tmp_R_COL$MJ_MN_COL)) == 2){
      
      ttest2g <- t.test(pinpis_R, pinpis_COL, alternative = "greater")$p.value
      ttest2l <- t.test(pinpis_R, pinpis_COL, alternative = "less")$p.value
      non_param_pval_2g <- wilcox.test(pinpis_R, pinpis_COL, data=pinpis_tmp_R_COL,  alternative = "greater")$p.value
      non_param_pval_2l <- wilcox.test(pinpis_R, pinpis_COL, data=pinpis_tmp_R_COL,  alternative = "less")$p.value
      
      tmp2 <- tibble(names=name, contig=contig_ref, start=start_ref, end=end_ref, inversion=inversion_ref, 
                     MJ_MN_COL='R', piN_piS_hapl_median=piN_piS_median_R,
                     piN_piS_COL_median= piN_piS_COL_median,
                     p_value_ttest_gr=ttest2g, p_value_ttest_ls=ttest2l, 
                     p_value_MW_gr=non_param_pval_2g, p_value_MW_ls=non_param_pval_2l)
      
      MN_MJ_R_COL_1s_ttest <- rbind(MN_MJ_R_COL_1s_ttest, tmp2)
  
    }
    ## f-test between MN and collinear 
    if(length(unique(pinpis_tmp_MN_COL$MJ_MN_COL)) == 2){
      
      ttest3g <- t.test(pinpis_MN, pinpis_COL, alternative = "greater")$p.value
      ttest3l <- t.test(pinpis_MN, pinpis_COL, alternative = "less")$p.value
      non_param_pval_3g <- wilcox.test(pinpis_MN, pinpis_COL, data=pinpis_tmp_MN_COL,  alternative = "greater")$p.value
      non_param_pval_3l <- wilcox.test(pinpis_MN, pinpis_COL, data=pinpis_tmp_MN_COL,  alternative = "less")$p.value
      
      tmp3 <- tibble(names=name, contig=contig_ref, start=start_ref, end=end_ref, inversion=inversion_ref, 
                     MJ_MN_COL='MN', piN_piS_hapl_median=piN_piS_median_MN,
                     piN_piS_COL_median= piN_piS_COL_median,
                     p_value_ttest_gr=ttest3g, p_value_ttest_ls=ttest3l, 
                     p_value_MW_gr=non_param_pval_3g, p_value_MW_ls=non_param_pval_3l)
      
      MN_MJ_R_COL_1s_ttest <- rbind(MN_MJ_R_COL_1s_ttest, tmp3)
    }
    ### 2 haplotypes  ####
  } else {
    
    ## 2 haplotypes #######
    
    pinpis_MJ <- pinpis_tmp %>% filter(MJ_MN_COL == "MJ")  %>% pull(piN_piS)
    pinpis_MN <- pinpis_tmp %>% filter(MJ_MN_COL == "MN")  %>% pull(piN_piS)
    
    pinpis_tmp_MJ_COL <- joint_tibble %>% filter(MJ_MN_COL == "MJ" | MJ_MN_COL == "collinear")
    pinpis_tmp_MN_COL <- joint_tibble %>% filter(MJ_MN_COL == "MN" | MJ_MN_COL == "collinear")
    
    if(length(pinpis_MJ) != 0) {
        
        ttest_MJ_COLg <- t.test(pinpis_MJ,pinpis_COL, alternative = "greater")$p.value
        ttest_MJ_COLl <- t.test(pinpis_MJ,pinpis_COL, alternative = "less")$p.value
        non_param_pval_MJ_COLg <- wilcox.test(pinpis_MJ,pinpis_COL, data=pinpis_tmp_MJ_COL,  alternative = "greater")$p.value
        non_param_pval_MJ_COLl <- wilcox.test(pinpis_MJ,pinpis_COL, data=pinpis_tmp_MJ_COL,  alternative = "less")$p.value
        
        tmp1 <- tibble(names=name, contig=contig_ref, start=start_ref, end=end_ref, inversion=inversion_ref, 
                       MJ_MN_COL='MJ', piN_piS_hapl_median=piN_piS_median_MJ,
                       piN_piS_COL_median= piN_piS_COL_median,
                       p_value_ttest_gr=ttest_MJ_COLg, p_value_ttest_ls=ttest_MJ_COLl, 
                       p_value_MW_gr=non_param_pval_MJ_COLg, p_value_MW_ls=non_param_pval_MJ_COLl)
        
        MN_MJ_R_COL_1s_ttest <- rbind(MN_MJ_R_COL_1s_ttest, tmp1)
        

    }# else {
      
      if(length(pinpis_MN) != 0) {
        
        ttest_MN_COLg <- t.test(pinpis_MN,pinpis_COL, alternative = "greater")$p.value
        ttest_MN_COLl <- t.test(pinpis_MN,pinpis_COL, alternative = "less")$p.value
        non_param_pval_MN_COLg <- wilcox.test(pinpis_MN,pinpis_COL, data=pinpis_tmp_MN_COL,  alternative = "greater")$p.value
        non_param_pval_MN_COLl <- wilcox.test(pinpis_MN,pinpis_COL, data=pinpis_tmp_MN_COL,  alternative = "less")$p.value
        
        tmp2 <- tibble(names=name, contig=contig_ref, start=start_ref, end=end_ref, inversion=inversion_ref, 
                       MJ_MN_COL='MN', piN_piS_hapl_median=piN_piS_median_MN,
                       piN_piS_COL_median= piN_piS_COL_median,
                       p_value_ttest_gr=ttest_MN_COLg, p_value_ttest_ls=ttest_MN_COLl, 
                       p_value_MW_gr=non_param_pval_MN_COLg, p_value_MW_ls=non_param_pval_MN_COLl)
        
        MN_MJ_R_COL_1s_ttest <- rbind(MN_MJ_R_COL_1s_ttest, tmp2)
      } else {
        next;
      }
  }
  
}








MN_MJ_R_COL_1s_ttest %>% filter(p_value_ttest_gr < 0.05 | p_value_ttest_ls < 0.05 |
                                  p_value_MW_gr < 0.05 | p_value_MW_ls < 0.05) -> significant_MN_MJ_R_COL
write.table(MN_MJ_R_COL_1s_ttest, '/home/anstmykh/Desktop/mut_load/hapl_to_col_all_stats.txt', row.names = TRUE)

####### only  1s t-test on everything

MN_MJ_INT_COL_1s_ttest_all_filtered <- tibble(names=character(), contig=character(), start=numeric(),
                                 end=numeric(), inversion=character(), MJ_MN_COL=character(), piN_piS_hapl_median=numeric(),
                                 piN_piS_COL_median=numeric() ,p_value=numeric())


pinpis_COL <- piN_piS_table %>% filter(MJ_MN_COL == "collinear") %>% pull(piN_piS)
piN_piS_COL_median <- median(pinpis_COL, na.rm = T)
for (i in 1:nrow(names)) {
  
  contig_ref <- names[i,2]
  start_ref <-names[i,3]
  end_ref <- names[i,4]
  inversion_ref <- names[i,5]
  name <- names[i,1]
  
  pinpis_tmp <- piN_piS_table_filtered_ftest %>% filter(contig == contig_ref & start == start_ref & end == end_ref & 
                                                          inversion == inversion_ref)
  
  joint_tibble <- rbind(piN_piS_table_filtered_only_col, pinpis_tmp)
  
  
  
  pinpis_MJ <- pinpis_tmp %>% filter(MJ_MN_COL == "MJ")  %>% pull(piN_piS)
  pinpis_MN <- pinpis_tmp %>% filter(MJ_MN_COL == "MN")  %>% pull(piN_piS)
  pinpis_INT <- pinpis_tmp %>% filter(MJ_MN_COL == "INT") %>% pull(piN_piS)
  
  piN_piS_median_MJ <- median(pinpis_MJ, na.rm = T)
  piN_piS_median_MN <- median(pinpis_MN, na.rm = T)
  piN_piS_median_INT <- median(pinpis_INT, na.rm = T)
  
  
  ### 3 haplotypes ####
  if(length(pinpis_INT) != 0) {
    
    
    pinpis_tmp_MJ_COL <- joint_tibble %>% filter(MJ_MN_COL == "MJ" | MJ_MN_COL == "collinear")
    
    pinpis_tmp_INT_COL <- joint_tibble %>% filter(MJ_MN_COL == "INT" | MJ_MN_COL == "collinear")
    
    pinpis_tmp_MN_COL <- joint_tibble %>% filter(MJ_MN_COL == "MN" | MJ_MN_COL == "collinear")
    
    
    ## f-test between MJ and collinear
    if(length(unique(pinpis_tmp_MJ_COL$MJ_MN_COL)) == 2){
      ## now performing t-test  
      ttest1 <- t.test(pinpis_MJ, pinpis_COL, alternative = "greater")$p.value
      
      tmp1.1 <- tibble(names=name, contig=contig_ref, start=start_ref, end=end_ref, inversion=inversion_ref, 
                       MJ_MN_COL='MJ', piN_piS_hapl_median=piN_piS_median_MJ,
                       piN_piS_COL_median= piN_piS_COL_median, p_value=ttest1 )
      
      MN_MJ_INT_COL_1s_ttest_all_filtered <- rbind(MN_MJ_INT_COL_1s_ttest_all_filtered, tmp1.1)
      
    }
    ## f-test between INT and collinear
    if(length(unique(pinpis_tmp_INT_COL$MJ_MN_COL)) == 2){
      ## now performing t-test 
      ttest2 <- t.test(pinpis_INT, pinpis_COL,alternative = "greater")$p.value
      
      tmp2.1 <- tibble(names=name, contig=contig_ref, start=start_ref, end=end_ref, inversion=inversion_ref, 
                       MJ_MN_COL='INT', piN_piS_hapl_median=piN_piS_median_INT,
                       piN_piS_COL_median= piN_piS_COL_median, p_value=ttest2 )
      
      MN_MJ_INT_COL_1s_ttest_all_filtered <- rbind(MN_MJ_INT_COL_1s_ttest_all_filtered, tmp2.1)
    }
    ## f-test between MN and collinear 
    if(length(unique(pinpis_tmp_MN_COL$MJ_MN_COL)) == 2){
      ## now performing t-test 
      ttest3 <- t.test(pinpis_MN, pinpis_COL,alternative = "greater")$p.value
      
      tmp3.1 <- tibble(names=name, contig=contig_ref, start=start_ref, end=end_ref, inversion=inversion_ref, 
                       MJ_MN_COL='MN', piN_piS_hapl_median=piN_piS_median_MN,
                       piN_piS_COL_median= piN_piS_COL_median, p_value=ttest3 )
      
      MN_MJ_INT_COL_1s_ttest_all_filtered <- rbind(MN_MJ_INT_COL_1s_ttest_all_filtered, tmp3.1)
    }
    ### 2 haplotypes  ####
  } else {
    
    ## 2 haplotypes #######
    
    pinpis_MJ <- pinpis_tmp %>% filter(MJ_MN_COL == "MJ")  %>% pull(piN_piS)
    pinpis_MN <- pinpis_tmp %>% filter(MJ_MN_COL == "MN")  %>% pull(piN_piS)
    
    pinpis_tmp_MJ_COL <- joint_tibble %>% filter(MJ_MN_COL == "MJ" | MJ_MN_COL == "collinear")
    pinpis_tmp_MN_COL <- joint_tibble %>% filter(MJ_MN_COL == "MN" | MJ_MN_COL == "collinear")
    
    if(length(pinpis_MJ) != 0) {
      
      ttest_MJ_COL <- t.test(pinpis_MJ,pinpis_COL, alternative = "greater")$p.value
      
      tmp1.1 <- tibble(names=name, contig=contig_ref, start=start_ref, end=end_ref, inversion=inversion_ref, 
                       MJ_MN_COL='MJ', piN_piS_hapl_median=piN_piS_median_MJ,
                       piN_piS_COL_median= piN_piS_COL_median, p_value=ttest_MJ_COL )
      
      MN_MJ_INT_COL_1s_ttest_all_filtered <- rbind(MN_MJ_INT_COL_1s_ttest_all_filtered, tmp1.1)
    }# else {
    
    if(length(pinpis_MN) != 0) {
      
      ttest_MN_COL <- t.test(pinpis_MN,pinpis_COL, alternative = "greater")$p.value
      
      tmp1.1 <- tibble(names=name, contig=contig_ref, start=start_ref, end=end_ref, inversion=inversion_ref, 
                       MJ_MN_COL='MN', piN_piS_hapl_median=piN_piS_median_MN,
                       piN_piS_COL_median= piN_piS_COL_median, p_value=ttest_MN_COL )
      
      MN_MJ_INT_COL_1s_ttest_all_filtered <- rbind(MN_MJ_INT_COL_1s_ttest_all_filtered, tmp1.1)
      
    } else {
      next;
    }
    #}
  }
  
}

######## placing the sample in a collinear distribution 



#Preparing data for plotting: significant pinpis outliers
sign_outl_sample_for_col_dist <- rbind(significant_MN_MJ_INT_COL_mnn_whtn_wlc %>% select(piN_piS_hapl_median, names, MJ_MN_COL), 
                                                        significant_MN_MJ_INT_COL_1s_ttest  %>% select(piN_piS_hapl_median, names, MJ_MN_COL))
#Preparing data for plotting: good sample, non-significant
MN_MJ_INT_COL_mnn_whtn_wlc_nonsign_plot <- MN_MJ_INT_COL_mnn_whtn_wlc %>% filter(p_value > 0.05) %>% select(piN_piS_hapl_median, names, MJ_MN_COL)
MN_MJ_INT_COL_mnn_whtn_wlc_plot <- MN_MJ_INT_COL_mnn_whtn_wlc %>% select(piN_piS_hapl_median, names, MJ_MN_COL)
#Preparing data for plotting: piN/piS medians for all genotypes
piNpiS_all_genotypes_plot <- rbind(MN_MJ_INT_COL_mnn_whtn_wlc_plot, haplotypes_on_col_dist_merge[,1:3])





# palette_12<-brewer.pal(21,"Set3")
# names(palette_12) <- levels(significant_MN_MJ_INT_COL_mnn_whtn_wlc$piN_piS)
# custom_colors <- scale_colour_manual(name = "pop", values = colors)

boxplot_col2 <- ggplot(piN_piS_table_filtered_only_col_plot) +  
  geom_histogram(aes(x=piN_piS), binwidth = 0.01) +
  geom_vline(data = haplotypes_on_col_dist, aes(xintercept = piN_piS, color = "bad_sample"), linetype = "dashed") + 
  geom_vline(data = significant_MN_MJ_INT_COL_mnn_whtn_wlc, aes(xintercept = piN_piS_hapl_median, color = "sign_outliers"), linetype = "solid") +
  geom_vline(data = MN_MJ_INT_COL_mnn_whtn_wlc_nonsign_plot, aes(xintercept = piN_piS_hapl_median, color = "nonsign"), linetype = "solid") +
  scale_color_manual(name = "genotypes", values = c(bad_sample = "blue", sign_outliers = "red", nonsign = "green"))

jpeg("mut_load_col_dist_with_sing_and_nonsign_and_lim_gen.jpeg", width = 1000, height = 600)
print(boxplot_col2)
dev.off()






          
