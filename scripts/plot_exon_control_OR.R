
# READ REPEATMASKER

################################################################################

library(dplyr)
library(tidyr)
library(purrr)
library(writexl)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

################################################################################

# INPUT


rm_cryptics <- "data/repeatmasker/junction_windows.fa_rm.bed"
rm_ctrl <- "data/repeatmasker/control_all_controls_as_exons.fa_rm.bed"

dir_meta <- "data/repeatmasker/metadata_mapping_exon_controls.txt"
sr_table <- "data/repeatmasker/simplerepeats_table.txt"



source("~/Documents/GitHub/online/cryptic_biology/scripts/repeatmasker_analyses.R")
################################################################################

# CODE

metadata <- fread(dir_meta, header = T)
cand_ctrl <- fread(dir_meta, header = T)

cryp_rep <- fread(rm_cryptics)
colnames(cryp_rep) <- c("id", "rep_start", "rep_end", "repeats", "rep_length", "rep_strand", "rep_class", "rep_family", "V9", "rep_ID")
cryp_rep$name <- sapply(strsplit(as.character(cryp_rep$id), "::"), `[`, 1)
cryp_rep = cryp_rep %>% 
    left_join(cand_ctrl %>% distinct(cryptic_name,junc_cat,new_cate),by = c("name" = 'cryptic_name'))

ctrl_rep <- fread(rm_ctrl)
colnames(ctrl_rep) <- c("name", "rep_start", "rep_end", "repeats", "rep_length", "rep_strand", "rep_class", "rep_family", "V9", "rep_ID")
ctrl_rep$id = ctrl_rep$name
ctrl_rep = ctrl_rep %>% 
    left_join(cand_ctrl %>% distinct(name,junc_cat,new_cate))

acc_early <- cand_ctrl %>% filter(new_cate == "Early" & junc_cat == "acceptor")
cryp_acc_early <- unique(acc_early$cryptic_name)
ctrl_acc_early <- unique(acc_early$name)

don_early <- cand_ctrl %>% filter(new_cate == "Early" & junc_cat == "donor")
cryp_don_early <- unique(don_early$cryptic)
ctrl_don_early <- unique(don_early$ctrl)

acc_late <- cand_ctrl %>% filter(new_cate == "Late" & junc_cat == "acceptor")
cryp_acc_late <- unique(acc_late$cryptic_name)
ctrl_acc_late <- unique(acc_late$name)

don_late <- cand_ctrl %>% filter(new_cate == "Late" & junc_cat == "donor")
cryp_don_late <- unique(don_late$cryptic_name)
ctrl_don_late <- unique(don_late$name)




rep_classes <- unique(c(ctrl_rep$rep_class, cryp_rep$rep_class))



for(i in 1:length(rep_classes)){
    
    repeat_c <- rep_classes[i]
    print(repeat_c)
    
    # Early Acceptor
    cryp_acc_early_df <- cryp_rep %>% filter(name %in% cryp_acc_early & rep_class == repeat_c)
    ctrl_acc_early_df <- ctrl_rep %>% filter(name %in% ctrl_acc_early & rep_class == repeat_c)
    
    cryp_acc_early_df_rep <- rep_ranges(cryp_acc_early_df, "cryptic")
    ctrl_acc_early_df_rep <- rep_ranges(ctrl_acc_early_df, "controls")
    
    acc_early_df <- merge(cryp_acc_early_df_rep, ctrl_acc_early_df_rep, by="position")
    acc_early_df$n_cry <- length(cryp_acc_early)
    acc_early_df$n_ctrl <- length(ctrl_acc_early)
    
    acc_early_df <- rowwise_fisher_test(
        acc_early_df, 
        cryptic_col = "cryptic", 
        set1_col = "controls", 
        n_cry_col = "n_cry", 
        n_ctrl_col = "n_ctrl",
        pval_out = "pv",
        or_out = "OR"
    )
    
    acc_early_df$new_cate <- "Early"
    acc_early_df$junc_cat <- "acceptor"
    
    
    
    # Early Donor
    cryp_don_early_df <- cryp_rep %>% filter(name %in% cryp_don_early & rep_class == repeat_c)
    ctrl_don_early_df <- ctrl_rep %>% filter(name %in% ctrl_don_early & rep_class == repeat_c)
    
    cryp_don_early_df_rep <- rep_ranges(cryp_don_early_df, "cryptic")
    ctrl_don_early_df_rep <- rep_ranges(ctrl_don_early_df, "controls")
    
    don_early_df <- merge(cryp_don_early_df_rep, ctrl_don_early_df_rep, by="position")
    don_early_df$n_cry <- length(cryp_don_early)
    don_early_df$n_ctrl <- length(ctrl_don_early)
    
    don_early_df <- rowwise_fisher_test(
        don_early_df, 
        cryptic_col = "cryptic", 
        set1_col = "controls", 
        n_cry_col = "n_cry", 
        n_ctrl_col = "n_ctrl",
        pval_out = "pv",
        or_out = "OR"
    )
    
    don_early_df$new_cate <- "Early"
    don_early_df$junc_cat <- "donor"
    
    
    
    # Late Acceptor
    cryp_acc_late_df <- cryp_rep %>% filter(name %in% cryp_acc_late & rep_class == repeat_c)
    ctrl_acc_late_df <- ctrl_rep %>% filter(name %in% ctrl_acc_late & rep_class == repeat_c)
    
    cryp_acc_late_df_rep <- rep_ranges(cryp_acc_late_df, "cryptic")
    ctrl_acc_late_df_rep <- rep_ranges(ctrl_acc_late_df, "controls")
    
    acc_late_df <- merge(cryp_acc_late_df_rep, ctrl_acc_late_df_rep, by="position")
    acc_late_df$n_cry <- length(cryp_acc_late)
    acc_late_df$n_ctrl <- length(ctrl_acc_late)
    
    acc_late_df <- rowwise_fisher_test(
        acc_late_df, 
        cryptic_col = "cryptic", 
        set1_col = "controls", 
        n_cry_col = "n_cry", 
        n_ctrl_col = "n_ctrl",
        pval_out = "pv",
        or_out = "OR"
    )
    
    acc_late_df$new_cate <- "Late"
    acc_late_df$junc_cat <- "acceptor"
    
    
    
    # Late Donor
    cryp_don_late_df <- cryp_rep %>% filter(name %in% cryp_don_late & rep_class == repeat_c)
    ctrl_don_late_df <- ctrl_rep %>% filter(name %in% ctrl_don_late & rep_class == repeat_c)
    
    cryp_don_late_df_rep <- rep_ranges(cryp_don_late_df, "cryptic")
    ctrl_don_late_df_rep <- rep_ranges(ctrl_don_late_df, "controls")
    
    don_late_df <- merge(cryp_don_late_df_rep, ctrl_don_late_df_rep, by="position")
    don_late_df$n_cry <- length(cryp_don_late)
    don_late_df$n_ctrl <- length(ctrl_don_late)
    
    don_late_df <- rowwise_fisher_test(
        don_late_df, 
        cryptic_col = "cryptic", 
        set1_col = "controls", 
        n_cry_col = "n_cry", 
        n_ctrl_col = "n_ctrl",
        pval_out = "pv",
        or_out = "OR"
    )
    
    don_late_df$new_cate <- "Late"
    don_late_df$junc_cat <- "donor"
    
    
    head(acc_early_df)
    head(don_early_df)
    head(acc_late_df)
    head(don_late_df)
    
    final_df <- rbind(acc_early_df, don_early_df)
    final_df <- rbind(final_df, acc_late_df)
    final_df <- rbind(final_df, don_late_df)
    final_df$position <- final_df$position -250
    
    
    
    final_df$repeat_cat <- repeat_c
    
    if(i == 1){
        repeats_df <- final_df
    } else {
        repeats_df <- rbind(repeats_df, final_df) 
    }
    
    
    
}



# write.table(repeats_df, paste0(dir_plot, "repeats_class_table.txt"), col.names = T, row.names = F, sep = "\t", quote = F)



repeats_df2 <- repeats_df
repeats_df2$logOR <- log((repeats_df2$OR)+0.01)
repeats_df2[repeats_df2$pv>0.05,]$logOR <- 0
repeats_df2$logOR[is.infinite(repeats_df2$logOR)] <- 0.5
repeats_df2$rows <- paste0(repeats_df2$new_cate, "_", repeats_df2$repeat_cat)
repeats_df2 <- repeats_df2[,c(1,11,12,9,10)]

acceptor <- repeats_df2 %>% filter(junc_cat == "acceptor")
donor <- repeats_df2 %>% filter(junc_cat == "donor")

acceptor <- acceptor %>%
    group_by(repeat_cat) %>%
    filter(sum(logOR, na.rm = TRUE) != 0) %>%
    ungroup()

donor <- donor %>%
    group_by(repeat_cat) %>%
    filter(sum(logOR, na.rm = TRUE) != 0) %>%
    ungroup()





col_fun = colorRamp2(c(-4,0, 4), c("blue", "gray95", "red"))


p1 <- repeats_heatmap(acceptor, title = "Acceptor", col_fun=col_fun)
p2 <- repeats_heatmap(donor, title = "Donor", col_fun=col_fun)

cowplot::plot_grid(plotlist = list(as.grob(p1),as.grob(p2)))
