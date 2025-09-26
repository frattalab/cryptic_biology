source("~/Documents/GitHub/online/cryptic_biology/scripts/repeatmasker_analyses.R")
library(tidyverse)
library(data.table)
rm_cryptics <- "data/repeatmasker/junction_windows.fa_rm.bed"
cand_ctrl_dir <- "data/repeatmasker/cryp_set1.txt"
dir_meta <- "data/repeatmasker/metadata.txt"
sr_table <- "data/repeatmasker/simplerepeats_table.txt"

metadata <- read.table(dir_meta, header = T)
cand_ctrl <- read.table(cand_ctrl_dir, header = T)
cand_ctrl <- merge(cand_ctrl,metadata,  by.x="cryptic", by.y="name")
cand_ctrl$junc_type <- sapply(strsplit(as.character(cand_ctrl$cryptic), "_"), `[`, 2)

cryp_rep <- read.table(rm_cryptics)
colnames(cryp_rep) <- c("id", "rep_start", "rep_end", "repeats", "rep_length", "rep_strand", "rep_class", "rep_family", "V9", "rep_ID")
cryp_rep$name <- sapply(strsplit(as.character(cryp_rep$id), "::"), `[`, 1)
cryp_rep <- merge(cryp_rep, distinct(cand_ctrl[,c("cryptic", "junc_type", "new_cate")]), by.x="name", by.y = "cryptic")


# Folder path
control_folder <- "data/output_rm_different_exons/"

# Get all files that match your pattern
files <- list.files(control_folder, pattern = "\\.fa_rm\\.bed$", full.names = TRUE)

# Read and bind all files
dt_all <- rbindlist(lapply(files, fread))


possible_exon_inclusion_rates = unique(gsub("/donor","",gsub("/acceptor","",gsub(control_folder,"",files))))
combinations_df <- expand.grid(
    timing = c("Early", "Late"),
    junc_type = c("acceptor", "donor"),
    exon_inclusion_type = possible_exon_inclusion_rates,
    stringsAsFactors = FALSE
)

results_list <- data.table()

for(i in 1:nrow(combinations_df)) {

    timing <- combinations_df$timing[i]
    type <- combinations_df$junc_type[i]
    exon_inclusion_type <- combinations_df$exon_inclusion_type[i]
    
    
    cat("Processing:", timing, type, exon_inclusion_type,"\n")
    
    # Filter cryptic set for this combination
    this_cryptic_set <- cryp_rep %>% 
        filter(new_cate == timing, 
               junc_type == type)
    # Read in the control
    control_file_path <- file.path(control_folder, paste0(type, exon_inclusion_type))
    control_set <- fread(control_file_path)
    
    
    colnames(control_set) <- c("id", "rep_start", "rep_end", "repeats", "rep_length", 
                               "rep_strand", "rep_class", "rep_family", "V9", "rep_ID")
    # Get all the possible repeat classes
    repeat_clases = unique(c(this_cryptic_set$rep_class))
    cat("This event type has these classes:", timing, type, repeat_clases,"\n")

    for(rep in repeat_clases){

        cryptic_set_this_repeat = this_cryptic_set %>% 
            filter(rep_class == rep)
        # Generate the position dataframe for cryptics
        cryptic_postion_df <- rep_ranges(cryptic_set_this_repeat, "cryptic")
        control_set_this_rep <- control_set %>% filter(rep_class == rep)
        control_postion_df <- rep_ranges(control_set_this_rep, "controls")
        # Merge and create final combination
        this_combination <- cryptic_postion_df %>% 
            left_join(control_postion_df) %>% 
            mutate(n_cry = length(unique(this_cryptic_set$name)),
                   n_ctrl = length(unique(control_set$id)),
                   timing = timing,
                   rep_class = rep,
                   junc_type = type,
                   exon_inclusion_type = exon_inclusion_type)
        
        # Store result (assuming you want to collect all results)
        results_list = rbind(results_list,this_combination)
    }

}

final_results <- bind_rows(results_list)

test_out <- rowwise_fisher_test(
    results_list,
    cryptic_col = "cryptic",
    set1_col = "controls",
    n_cry_col = "n_cry",
    n_ctrl_col = "n_ctrl",
    pval_out = "pv",
    or_out = "OR"
)
test_out$adj.p = p.adjust(test_out$pv)

repeats_df2 <- test_out
repeats_df2$logOR <- log((repeats_df2$OR)+0.01)
repeats_df2 = repeats_df2 %>% 
    mutate(logOR = ifelse(adj.p > 0.05,0,logOR)) %>% 
    mutate(logOR = ifelse(is.infinite(logOR),0.5,logOR))

repeats_df2$rows <- paste0(repeats_df2$timing, "_", repeats_df2$rep_class)


filtered_output = repeats_df2 %>% 
    group_by(junc_type,timing,rep_class) %>%
    filter(sum(logOR, na.rm = TRUE) != 0) %>%
    ungroup() %>% 
    dplyr::rename(repeat_cat = rep_class)

filtered_output %>% 
    mutate(position = position - 250) %>% 
    filter(grepl("In",exon_inclusion_type)) %>%  
    mutate(exon_inclusion_type = case_match(exon_inclusion_type,
                                            "_99___In_frame.fa_rm.bed" ~ "99%+",
                                            "1_25___In_frame.fa_rm.bed" ~ "1-25%",
                                            "25_85___In_frame.fa_rm.bed" ~ "25-85%",
                                            "85_99___In_frame.fa_rm.bed" ~ "85-99%")) %>% 
    ggplot(aes(x = position, y = logOR, color = exon_inclusion_type,linetype = timing)) + 
    facet_grid(rows = vars(repeat_cat),col = vars(junc_type)) +
    geom_line(size = 1.3) + 
    theme_classic() +
    ggtitle("In-frame") + 
    scale_color_manual(values = c("#BE6FD2", "#B0E17C", "#DC977F", "#A8CCD5"))



filtered_output %>% 
    mutate(position = position - 250) %>% 
    filter(grepl("Out",exon_inclusion_type)) %>%  
    mutate(exon_inclusion_type = case_match(exon_inclusion_type,
                                            "_99___Out_frame.fa_rm.bed" ~ "99%+",
                                            "1_25___Out_frame.fa_rm.bed" ~ "1-25%",
                                            "25_85___Out_frame.fa_rm.bed" ~ "25-85%",
                                            "85_99___Out_frame.fa_rm.bed" ~ "85-99%")) %>% 
    ggplot(aes(x = position, y = logOR, color = exon_inclusion_type,linetype = timing)) + 
    facet_grid(rows = vars(repeat_cat),col = vars(junc_type)) +
    geom_line(size = 1.3) + 
    theme_classic() +
    ggtitle("Out-frame") + 
    scale_color_manual(values = c("#BE6FD2", "#B0E17C", "#DC977F", "#A8CCD5"))

# Enrichment of simple repeat types ---------------------------------------
sr <- read.table(sr_table, header = T)


possible_simples <- expand.grid(
    timing = c("Early", "Late"),
    junc_type = c("acceptor", "donor"),
    exon_inclusion_type = possible_exon_inclusion_rates,
    rep_type = unique(sr$rep_type),
    stringsAsFactors = FALSE
)



simple_results_list <- data.table()

for(i in 1:nrow(possible_simples)) {
    
    timing <- possible_simples$timing[i]
    type <- possible_simples$junc_type[i]
    exon_inclusion_type <- possible_simples$exon_inclusion_type[i]
    this_rep_type <- possible_simples$exon_inclusion_type[i]
    
    cat("Processing:", timing, type, this_rep_type,"\n")
    
    # Filter cryptic set for this combination
    this_cryptic_set <- cryp_rep %>% 
        filter(new_cate == timing, 
               junc_type == type)
    # Read in the control
    control_file_path <- file.path(control_folder, paste0(type, exon_inclusion_type))
    control_set <- fread(control_file_path)
    
    
    colnames(control_set) <- c("id", "rep_start", "rep_end", "repeats", "rep_length", 
                               "rep_strand", "rep_class", "rep_family", "V9", "rep_ID")
    # Get all the possible repeat classes
    simple_cryptics = this_cryptic_set %>% 
        left_join(sr) %>% 
        filter(!is.na(rep_type))
    repeat_clases = unique(c(this_cryptic_set$rep_class))
    cat("This event type has these classes:", timing, type, repeat_clases,"\n")
    
    for(rep in repeat_clases){
        
        cryptic_set_this_repeat = this_cryptic_set %>% 
            filter(rep_class == rep)
        # Generate the position dataframe for cryptics
        cryptic_postion_df <- rep_ranges(cryptic_set_this_repeat, "cryptic")
        control_set_this_rep <- control_set %>% filter(rep_class == rep)
        control_postion_df <- rep_ranges(control_set_this_rep, "controls")
        # Merge and create final combination
        this_combination <- cryptic_postion_df %>% 
            left_join(control_postion_df) %>% 
            mutate(n_cry = length(unique(this_cryptic_set$name)),
                   n_ctrl = length(unique(control_set$id)),
                   timing = timing,
                   rep_class = rep,
                   junc_type = type,
                   exon_inclusion_type = exon_inclusion_type)
        
        # Store result (assuming you want to collect all results)
        results_list = rbind(results_list,this_combination)
    }
    
}