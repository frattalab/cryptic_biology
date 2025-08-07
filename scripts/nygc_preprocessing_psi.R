library(tidyverse)
library(arrow)
library(pROC)
library(ggrepel)
library(data.table)
# plot_function -----------------------------------------------------------

plot_junction = function(junc,plotin_table = complete_table,measure = 'spliced_reads'){
    
    vals = c(`ALS-TDP` = "#E1BE6A", Control = "#40B0A6", `FTLD-TDP` = "#E1BE6A", 
             `ALS-non-TDP` = "#408A3E", `FTLD-non-TDP` = "#408A3E")
    

    plotin_table = plotin_table %>% 
        filter(junctions_coords == junc)
    plotin_table = plotin_table %>% 
        mutate(tdp_proteinopathy = ifelse(tdp_proteinopathy %in% c("FTLD-TAU","FTLD-FUS"),"FTLD-non-TDP",tdp_proteinopathy))
    
    plot_gene_name = plotin_table %>% pull(gene_name) %>% unique()
    plot_gene_name = plot_gene_name[1]
    
    plotin_table = plotin_table %>% filter(gene_name == plot_gene_name)
    
    # junc_type = plotin_table[junctions_coords == junc,unique(type)]
    plot_title = glue::glue("{plot_gene_name} - {junc}")
    
    if(measure == 'spliced_reads'){
        plt = plotin_table |>
            # filter(!tissue_clean %in% c("Other","Choroid","Hippocampus","Liver","Sensory_Cortex","Occipital_Cortex")) |> 
            # mutate(disease =  gsub("-n","\nn",disease)) |> 
            mutate(tdp_proteinopathy = fct_relevel(tdp_proteinopathy,"Control", "ALS-non-TDP","ALS-TDP", "FTLD-non-TDP", "FTLD-TDP")) |>
            ggplot(aes(x = tdp_proteinopathy, y = junction_count, fill = tdp_proteinopathy)) +
            geom_boxplot(outlier.colour = NA) +
            geom_jitter(height = 0,alpha = 0.7, pch = 21) +
            facet_wrap(~rna_tissue_source_simplified, scales = "free_y") +
            scale_fill_manual(values = vals)  +
            ylab("N spliced reads") +
            xlab("") +
            ggtitle(plot_title) +
            ggpubr::theme_pubr() +
            theme(legend.position = 'none') +
            theme(text = element_text(size = 10)) +
            scale_x_discrete(guide = guide_axis(n.dodge=2)) 

    }else{
        plt = plotin_table |>
            # filter(!tissue_clean %in% c("Other","Choroid","Hippocampus","Liver","Sensory_Cortex","Occipital_Cortex")) |> 
            # mutate(disease =  gsub("-n","\nn",disease)) |> 
            mutate(tdp_proteinopathy = fct_relevel(tdp_proteinopathy,"Control", "ALS-non-TDP","ALS-TDP", "FTLD-non-TDP", "FTLD-TDP")) |>
            ggplot(aes(x = tdp_proteinopathy, y = psi, fill = tdp_proteinopathy)) +
            geom_boxplot(outlier.colour = NA) +
            geom_jitter(height = 0,alpha = 0.7, pch = 21) +
            facet_wrap(~rna_tissue_source_simplified, scales = "free_y") +
            scale_fill_manual(values = vals)  +
            ylab("PSI") +
            xlab("") +
            ggtitle(plot_title) +
            ggpubr::theme_pubr() +
            theme(legend.position = 'none') +
            theme(text = element_text(size = 10)) +
            scale_x_discrete(guide = guide_axis(n.dodge=2)) +
            scale_y_continuous(labels = scales::percent_format())
        
    }

    
    
    
    return(plt)
    
}
# roc_function ------------------------------------------------------------

process_psi_data <- function(data, not_included = c("sample", "numeric_prediction")) {
    # Scale the data and filter out columns with NaN values
    wide_psi_scale <- data %>% 
        mutate(across(starts_with("chr"), ~ scale(.) %>% as.vector())) %>% 
        select(where(~ !any(is.nan(.))))
    
    # Prepare the ROC vector, excluding certain columns
    ROCvector <- colnames(wide_psi_scale)
    ROCvector <- ROCvector[!(ROCvector %in% not_included )]
    
    my_df <- data.table()
    
    # Loop through each predictor to calculate ROC metrics
    for (predictor in ROCvector) {
        pROC_obj <- roc_(data = wide_psi_scale, "numeric_prediction", predictor, quiet = TRUE)
        
        # Calculate the FPR and TPR values
        # This would be the median PSI of the event (e.g.anything greater than this is an 'above normal' level of expression across the dataset)
        
        predictor_value <- 0
        fpr_value <- coords(pROC_obj, x = predictor_value, input = "threshold", ret = "specificity", 
                            best.method = "closest.topleft", transpose = FALSE)
        fpr_value <- 1 - fpr_value
        tpr_value <- coords(pROC_obj, x = predictor_value, input = "threshold", ret = "sensitivity", 
                            best.method = "closest.topleft", transpose = FALSE)
        
        # Compile results into a list
        lil_list <- list(paste_into_igv_junction = predictor,
                         auc = as.numeric(pROC_obj$auc),
                         avg_fpr = 1 - mean(pROC_obj$specificities[-c(1, length(pROC_obj$sensitivities))]),
                         avg_tpr = mean(pROC_obj$sensitivities[-c(1, length(pROC_obj$sensitivities))]),
                         fpr_value = as.numeric(fpr_value),
                         tpr_value = as.numeric(tpr_value))
        
        # Bind results to the my_df data.table
        my_df <- rbind(my_df, lil_list)
    }
 return(my_df)   
}
summary_post <- fread("data/cryptic_detection_summary_long_version.csv")

wide_chx = fread('~/Desktop/wide_chx.csv') %>% 
    mutate(group = 'CHX -SHY5Y')

wide_upf1 = fread('~/Desktop/wide_upf1.csv') %>% 
    mutate(group = 'UPF1 - i3Cortic')

wide_smg1 = fread('~/Desktop/wide_smg1.csv') %>% 
    mutate(group = 'SMG1i -SHY5Y')

these_genes = c("AARS1", "ACTL6B", "ARHGAP22", "CELSR3",
                "DNM1", "EPB41L4A", "HDGFL2", "IGLON5", "MYO18A", 
                "NECAB2", "PHF2", "PTPRZ1", "PXDN", 
                "SLC24A3", "SLC24A3", "STMN2", 
                "SYNJ2", "XPO4", "ZNF423", "ACSF2")

nygc_metadata = fread("/Users/annaleigh/Desktop/nygc_metadata_done.csv") %>% 
    dplyr::rename(sample = V4) %>% 
    select(-bam_path,-bam_name,-id)


nygc_events = arrow::open_dataset("/Users/annaleigh/Desktop/parquets/")


master_table <- fread(here::here("data/master_table.csv"))
master_table <- master_table %>%
    mutate(new_cate = case_when(tdp43_sensitivity_be2 == "" ~ tdp43_sensitivity_sh,
                                tdp43_sensitivity_sh == "" ~ tdp43_sensitivity_be2,
                                tdp43_sensitivity_sh == tdp43_sensitivity_be2 ~ tdp43_sensitivity_sh,
                                tdp43_sensitivity_be2 != 'Late' & tdp43_sensitivity_sh == 'Early' ~ 'Early',
                                tdp43_sensitivity_be2 != 'Early' & tdp43_sensitivity_sh == 'Late' ~ 'Late',
                                tdp43_sensitivity_sh == 'Intermediate'~  'Intermediate')) %>% 
    mutate(new_cate = ifelse(is.na(new_cate),"Ambiguous",new_cate)) %>% 
    mutate(comparison = 'curve_experiments')

compendium_events = fread('data/compendium_cryptics.csv')
bind_table = master_table %>% 
    distinct(paste_into_igv_junction,gene,strand) %>% 
    dplyr::rename(gene_name = gene) %>% 
    rbind((compendium_events %>% distinct(gene_name,paste_into_igv_junction,strand))) %>% 
    dplyr::rename(junctions_coords = paste_into_igv_junction) %>% 
    unique()  %>% 
    left_join(annotables::grch38 %>% select(symbol,strand),by = c("gene_name" = "symbol")) %>% 
    mutate(strand.y = ifelse(strand.y == 1, "+","-")) %>% 
    mutate(strand_final = ifelse(is.na(strand.x),strand.y,strand.x)) %>% 
    select(-strand.x,-strand.y) %>% 
    dplyr::rename(strand = strand_final) %>% 
    unique() %>% 
    separate(junctions_coords,remove = FALSE, convert = TRUE, into = c("chr","start","end")) %>% 
    mutate(paste_into_igv_junction = junctions_coords)




bind_table = bind_table %>% 
    rbind(wide_chx %>% 
              distinct(gene_name,strand,paste_into_igv_junction) %>% 
              separate(paste_into_igv_junction,remove = FALSE, convert = TRUE, into = c("chr","start","end")) %>% 
              mutate(junctions_coords = paste_into_igv_junction))
bind_table = bind_table %>% 
    rbind(wide_smg1 %>% 
              distinct(gene_name,strand,paste_into_igv_junction) %>% 
              separate(paste_into_igv_junction,remove = FALSE, convert = TRUE, into = c("chr","start","end")) %>% 
              mutate(junctions_coords = paste_into_igv_junction))

bind_table = bind_table %>% 
    rbind(wide_upf1 %>% 
              distinct(gene_name,strand,paste_into_igv_junction) %>% 
              separate(paste_into_igv_junction,remove = FALSE, convert = TRUE, into = c("chr","start","end")) %>% 
              mutate(junctions_coords = paste_into_igv_junction))
nygc_cryptic = nygc_events %>%
    filter(chrom %in% bind_table$chr, start %in% bind_table$start, end %in% bind_table$end) %>%
    collect()
# 
# 
nygc_cryptic_filtered = nygc_cryptic %>%
    mutate(junctions_coords = as.character(glue::glue("{chrom}:{start}-{end}"))) %>%
    filter(junctions_coords %in% bind_table$junctions_coords) %>%
    left_join(bind_table,by = c("start","end",'strand','junctions_coords','chrom' = 'chr'))

nygc_cryptic_filtered %>% 
    arrow::write_parquet('~/Desktop/nygc_cryptic_filtered_with_nmd.parquet')

nygc_cryptic_filtered = read_parquet('~/Desktop/nygc_cryptic_filtered_with_nmd.parquet')

# which cryptic events are specific to TDP path? --------------------------
nygc_cryptic_filtered = read_parquet('data/nygc_cryptic_filtered.parquet')
neural_meta = nygc_metadata %>% 
    filter(rna_tissue_source_simplified %in% c("Spinal_Cord_Thoracic", "Spinal_Cord_Cervical", "Spinal_Cord_Lumbar", 
                                               "Cortex_Frontal", "Cortex_Motor", "Cortex_Temporal", 
                                               "Medulla", "Hippocampus", "Choroid", 
                                               "Cortex_Unspecified", 
                                               "Spinal_Cord_Unspecified","Cerebellum")) %>% 
    filter(tdp_proteinopathy != "Unknown")
    
nygc_cryptic_filtered = nygc_cryptic_filtered %>% 
    filter(sample %in% neural_meta$sample) %>% 
    left_join(neural_meta %>% select(sample,rna_tissue_source_simplified,tdp_proteinopathy)) %>% 
    unique()
    
nygc_cryptic_filtered = nygc_cryptic_filtered %>% 
    mutate(group = case_when(grepl("Cord",rna_tissue_source_simplified) ~ 'Cord',
                             grepl("Cortex_Motor",rna_tissue_source_simplified) & 
                                 tdp_proteinopathy %in% c("Control","ALS-TDP","ALS-non-TDP") ~ "Motor Cortex",
                             rna_tissue_source_simplified %in% c("Cortex_Frontal","Cortex_Temporal") &
                                 tdp_proteinopathy %in% c("Control","FTLD-TAU","FTLD-FUS","FTLD-TDP") ~ "FTD frontal/temporal")) %>% 
    filter(rna_tissue_source_simplified %in% c("Cortex_Frontal","Cortex_Temporal","Cortex_Motor") |
               grepl("Cord",rna_tissue_source_simplified)) 

nygc_cryptic_filtered %>% arrow::write_parquet('~/Desktop/nygc_cryptic_filtered_metadata.parquet')
# what's the actual observation rate? -------------------------------------


neural_meta = neural_meta %>% 
    mutate(group = case_when(grepl("Cord",rna_tissue_source_simplified) ~ 'Cord',
                             grepl("Cortex_Motor",rna_tissue_source_simplified) & 
                                 tdp_proteinopathy %in% c("Control","ALS-TDP","ALS-non-TDP") ~ "Motor Cortex",
                             rna_tissue_source_simplified %in% c("Cortex_Frontal","Cortex_Temporal") &
                                 tdp_proteinopathy %in% c("Control","FTLD-TAU","FTLD-FUS","FTLD-TDP") ~ "FTD frontal/temporal")) %>% 
    filter(rna_tissue_source_simplified %in% c("Cortex_Frontal","Cortex_Temporal","Cortex_Motor") |
               grepl("Cord",rna_tissue_source_simplified)) %>% 
    mutate(variable = ifelse(tdp_proteinopathy %in% c("ALS-TDP","FTLD-TDP"), 'tdp_path', 'not_path'))

denom_tbl_nygc = neural_meta %>% 
    filter(!is.na(group)) %>% 
    group_by(variable,group) %>% 
    summarize(n_samples = n_distinct(sample))


n_sample_junc_nygc = nygc_cryptic_filtered %>% 
    filter(group %in% denom_tbl_nygc$group) %>% 
    # filter(paste_into_igv_junction %in% compendium_events$paste_into_igv_junction | 
    #            paste_into_igv_junction %in% master_table$paste_into_igv_junction) %>% 
    filter(junction_count >=2) %>%
    mutate(numeric_prediction = ifelse(tdp_proteinopathy %in% c("ALS-TDP","FTLD-TDP"), 'tdp_path', 'not_path')) %>% 
    group_by(paste_into_igv_junction,numeric_prediction,group) %>% 
    summarize(n_obs = n_distinct(sample))

n_sample_junc_nygc = n_sample_junc_nygc %>% 
    ungroup() %>% 
    pivot_wider(names_from = 'numeric_prediction', values_from = 'n_obs',values_fill = 0) %>% 
    reshape2::melt(id.vars = c("paste_into_igv_junction","group")) %>% 
    left_join(denom_tbl_nygc) %>% 
    mutate(frac_ob = value / n_samples)


# convert to wide psi -----------------------------------------------------


wide_psi = nygc_cryptic_filtered %>% 
    mutate(numeric_prediction = ifelse(tdp_proteinopathy %in% c("ALS-TDP","FTLD-TDP"), 1, 0)) %>% 
    distinct(sample,junctions_coords,numeric_prediction,psi,rna_tissue_source_simplified,tdp_proteinopathy) %>% 
    pivot_wider(values_from = 'psi',
                names_from = 'junctions_coords',values_fill = 0) %>% 
    as.data.table() 

wide_counts = nygc_cryptic_filtered %>%
    mutate(numeric_prediction = ifelse(tdp_proteinopathy %in% c("ALS-TDP","FTLD-TDP"), 1, 0)) %>%
    distinct(sample,junctions_coords,numeric_prediction,rna_tissue_source_simplified,tdp_proteinopathy,junction_count) %>%
    pivot_wider(values_from = 'junction_count',
                names_from = 'junctions_coords',values_fill = 0) %>%
    as.data.table()

complete_table = arrow::read_parquet("/Users/annaleigh/Documents/GitHub/cryptic_biology/data/complete_nygc_counts.parqet")
# Prediction with scaling psi - cortex only ---------------------------------------------
cortex_df = wide_psi %>% 
    filter(grepl("Cortex",rna_tissue_source_simplified)) %>% 
    process_psi_data(not_included = c("sample",'numeric_prediction',"rna_tissue_source_simplified","tdp_proteinopathy"))

ftd_cortex_df = wide_psi %>% 
    filter(rna_tissue_source_simplified %in% c("Cortex_Frontal","Cortex_Temporal")) %>% 
    filter(tdp_proteinopathy %in% c("Control","FTLD-TAU","FTLD-FUS","FTLD-TDP")) %>%     
    process_psi_data(not_included = c("sample",'numeric_prediction',"rna_tissue_source_simplified","tdp_proteinopathy"))

als_cortex_df = wide_psi %>% 
    filter(grepl("Cortex_Motor",rna_tissue_source_simplified)) %>% 
    filter(tdp_proteinopathy %in% c("Control","ALS-TDP","ALS-non-TDP")) %>% 
    process_psi_data(not_included = c("sample",'numeric_prediction',"rna_tissue_source_simplified","tdp_proteinopathy"))


cerebellum_psi = complete_table %>% 
    filter(rna_tissue_source_simplified == 'Cerebellum') %>% 
    filter(tdp_proteinopathy %in% c("Control","FTLD-TAU","FTLD-FUS","FTLD-TDP")) %>% 
    mutate(numeric_prediction = ifelse(tdp_proteinopathy %in% c("ALS-TDP","FTLD-TDP"), 1, 0)) %>% 
    distinct(sample,junctions_coords,numeric_prediction,psi,rna_tissue_source_simplified,tdp_proteinopathy) %>% 
    pivot_wider(values_from = 'psi',
                names_from = 'junctions_coords',values_fill = 0) %>% 
    as.data.table() 



cerebellum_df = process_psi_data(cerebellum_psi, not_included = c("sample",'numeric_prediction',"rna_tissue_source_simplified","tdp_proteinopathy"))

cerebellum_df %>% 
    left_join(complete_table %>% select(strand,gene_name,junctions_coords), by = c("paste_into_igv_junction" = "junctions_coords")) %>% 
    unique() %>% 
    arrange(fpr_value) %>% 
    fwrite('~/Desktop/cerbellum_predictions.csv')

# comparing predictive ability in cortex ----------------------------------
bind_table = complete_table %>% distinct(strand,gene_name,junctions_coords) %>% 
    dplyr::rename(paste_into_igv_junction = junctions_coords)
ftd_cortex_df %>% 
    left_join(als_cortex_df,by = c("paste_into_igv_junction"),
              suffix = c("_ftd_cortex","_als_motor")) %>% 
    filter(fpr_value_als_motor < 0.1 & fpr_value_ftd_cortex < 0.1) %>% 
    filter(tpr_value_als_motor > 0.1 | tpr_value_ftd_cortex > 0.1) %>% 
    left_join(bind_table) %>%  
    ggplot(aes(x = tpr_value_ftd_cortex,y = tpr_value_als_motor)) + 
    geom_point() + 
    geom_abline() + 
    geom_text_repel(aes(label = gene_name))


# comparison frontal and cerebellum ---------------------------------------
ftd_cortex_df %>% 
    left_join(cerebellum_df,by = c("paste_into_igv_junction"),
              suffix = c("_ftd_cortex","_ftd_cerebellum")) %>% 
    filter(fpr_value_ftd_cerebellum < 0.1 & fpr_value_ftd_cortex < 0.1) %>% 
    filter(tpr_value_ftd_cerebellum > 0.1 | tpr_value_ftd_cortex > 0.1) %>% 
    left_join(bind_table) %>%  
    ggplot(aes(x = tpr_value_ftd_cortex,y = tpr_value_ftd_cerebellum)) + 
    geom_point() + 
    geom_abline() + 
    geom_text_repel(aes(label = gene_name))

ftd_cortex_df %>% 
    left_join(cerebellum_df,by = c("paste_into_igv_junction"),
              suffix = c("_ftd_cortex","_ftd_cerebellum")) %>% 
    filter(fpr_value_ftd_cerebellum < 0.1 & fpr_value_ftd_cortex < 0.1) %>% 
    filter(tpr_value_ftd_cerebellum > 0.1 | tpr_value_ftd_cortex > 0.1) %>% 
    left_join(bind_table) %>%  
    ggplot(aes(x = auc_ftd_cortex,y = auc_ftd_cerebellum)) + 
    geom_point() + 
    geom_abline() + 
    geom_text_repel(aes(label = gene_name)) +
    ylab("AUC FTD Cerebellum") +
    xlab("AUC FTD Frontal Cortex") +
    theme_bw() +
    geom_hline(yintercept = 0.6,linetype = 'dotted') +
    geom_vline(xintercept = 0.6,linetype = 'dotted')


ftd_cortex_df %>% 
    full_join(cerebellum_df,by = c("paste_into_igv_junction"),
              suffix = c("_ftd_cortex","_ftd_cerebellum")) %>% 
    filter(fpr_value_ftd_cerebellum < 0.1) %>% 
    filter(tpr_value_ftd_cerebellum > 0.1 | tpr_value_ftd_cortex > 0.1) %>% 
    left_join(bind_table) %>%  
    ggplot(aes(x = auc_ftd_cortex,y = auc_ftd_cerebellum)) + 
    geom_point() + 
    geom_abline() + 
    geom_text_repel(aes(label = gene_name)) +
    ylab("AUC FTD Cerebellum") +
    xlab("AUC FTD Frontal Cortex") +
    theme_bw() +
    geom_hline(yintercept = 0.6,linetype = 'dotted') +
    geom_vline(xintercept = 0.6,linetype = 'dotted')

pdf('~/Desktop/cerebellum_auc_g60_fpr_l10.pdf')
dem = cerebellum_df %>% filter(auc > 0.6 & fpr_value < 0.1)
for(d in dem$paste_into_igv_junction){
    p = plot_junction(junc = d,measure = 'psi')
    p2 = plot_junction(junc = d)
    print(p)
    print(p2)
}
dev.off()
# Prediction with scaling psi - cord only ---------------------------------------------
cord_df = wide_psi %>% 
    filter(grepl("Cord",rna_tissue_source_simplified)) %>% 
    process_psi_data(not_included = c("sample",'numeric_prediction',"rna_tissue_source_simplified","tdp_proteinopathy"))
# Prediction with scaling psi - whole tissue ---------------------------------------------
total_df = wide_psi %>% 
    process_psi_data(not_included = c("sample",'numeric_prediction',"rna_tissue_source_simplified","tdp_proteinopathy"))
# Prediction with scaling psi - hippo only ---------------------------------------------
hippo_df = wide_psi %>% 
    filter(grepl("Hippo",rna_tissue_source_simplified)) %>% 
    process_psi_data(not_included = c("sample",'numeric_prediction',"rna_tissue_source_simplified","tdp_proteinopathy"))
cerebellum_df = wide_psi %>% 
    filter(grepl("Cerebellum",rna_tissue_source_simplified)) %>% 
    process_psi_data(not_included = c("sample",'numeric_prediction',"rna_tissue_source_simplified","tdp_proteinopathy"))
# wilcoxin test upregulated cortex ----------------------------------------
wilcox_test_res_cortex = complete_table %>% 
    filter(junctions_coords %in% bind_table$junctions_coords) %>% 
    filter(grepl("Cortex",rna_tissue_source_simplified)) %>% 
    mutate(condition = ifelse(tdp_proteinopathy %in% c("ALS-TDP","FTLD-TDP"),"TDP-43 proteinopathy","non-TDP-43-proteinopathy")) %>% 
    group_by(junctions_coords) %>% 
    nest() |> 
    mutate(wicox_test_result = map(data, ~ broom::tidy(wilcox.test(psi ~ condition, data = .x))),
           mean_path = map_dbl(data, ~ mean(.x$psi[.x$condition == "TDP-43 proteinopathy"], na.rm = TRUE)),
           mean_nonpath = map_dbl(data, ~ mean(.x$psi[.x$condition == "non-TDP-43-proteinopathy"], na.rm = TRUE)),
           ratio = log2(mean_path / mean_nonpath)) |> 
    unnest() %>% 
    distinct(junctions_coords,p.value,ratio,mean_path,mean_nonpath,ratio) %>% 
    left_join(bind_table,by = c("junctions_coords" = "junctions_coords")) %>% 
    arrange(-ratio) %>% 
    mutate(adj_pvalue = p.adjust(p.value))


upregulated_late_events = wilcox_test_res_cortex %>%
    filter(adj_pvalue < 0.01 & ratio >0) %>% 
    left_join(master_table) %>% 
    filter(!is.na(new_cate)) %>% 
    filter(new_cate == 'Late'&cryptic == "yes") %>% 
    arrange(-ratio)

# Cord - do cryptic events cluster together? Do ALS patients clust --------
pred_cord = cord_df %>% 
    filter(fpr_value < 0.1 & tpr_value > 0.1)
spinal_cluster = wide_psi %>% 
    filter(rna_tissue_source_simplified == 'Spinal_Cord_Cervical') %>% 
    filter(tdp_proteinopathy == 'ALS-TDP') %>% 
    select(sample,pred_cord$paste_into_igv_junction) %>% 
    column_to_rownames('sample') 
library(caret)
nzv <- nearZeroVar(spinal_cluster)
spinal_cluster_filtered <- spinal_cluster[, -nzv]

pred_cor = cortex_df %>% 
    filter(fpr_value < 0.1 & tpr_value > 0.1)
cortex_cluster = wide_psi %>% 
    # filter(rna_tissue_source_simplified == 'Cortex_Motor') %>% 
    # filter(tdp_proteinopathy == 'ALS-TDP') %>% 
    filter(rna_tissue_source_simplified == 'Cortex_Frontal') %>%
    filter(tdp_proteinopathy == 'FTLD-TDP') %>%
    select(sample,pred_cord$paste_into_igv_junction) %>% 
    column_to_rownames('sample') 
library(caret)
nzv <- nearZeroVar(cortex_cluster)
cortex_cluster_filtered <- cortex_cluster[, -nzv]
cortex_cluster_filtered = cortex_cluster_filtered %>% 
    rownames_to_column('sample') %>% 
    melt() %>% 
    left_join(master_table %>% select(new_cate,paste_into_igv_junction,ID_exon,gene),
              by = c("variable" = "paste_into_igv_junction")) %>% 
    filter(!is.na(new_cate)) %>% 
    group_by(sample,ID_exon) %>% 
    summarize(value = mean(value),new_cate,gene) %>% 
    unique() %>% 
    mutate(ident = paste(gene,new_cate,ID_exon))

cortex.cor <- cortex_cluster_filtered %>% 
    ungroup %>% 
    select(sample,ident,value) %>% 
    pivot_wider(names_from = ident) %>%    # t() for matrix transpose
    correlate() %>%    # correlate() is equivalent to cor() but                                                put NA as its diagonal entry and different class
    shave(upper = TRUE) %>%            # Shave the data frame to lower triangular matrix
    stretch(na.rm = TRUE) %>%           
    filter(abs(r) > 0.2) 

swiss.graph <- as_tbl_graph(cortex.cor, directed = FALSE)


set.seed(100)
ggraph(swiss.graph) + 
    geom_edge_link(aes(width = r), alpha = 0.2) + 
    scale_edge_width(range = c(0.2, 3)) +
    geom_node_point(, size = 3) +
    geom_node_text(aes(label = name), size = 3, repel = TRUE) +
    theme_graph()
# Do samples with upregulated late events have more cryptic splicig? --------
late_up_detection_cortex = complete_table %>% 
    filter(junctions_coords %in% upregulated_late_events$junctions_coords) %>% 
    filter(grepl("Cortex",rna_tissue_source_simplified)) %>% 
    mutate(detected = junction_count >=2) %>% 
    group_by(sample) %>% 
    summarize(n_late_up = sum(detected))
    
complete_table %>% 
    filter(junctions_coords %in% (cortex_df %>% filter(fpr_value < 0.1 & tpr_value > 0.1) %>% pull(paste_into_igv_junction))) %>% 
    filter(grepl("Cortex",rna_tissue_source_simplified)) %>% 
    left_join(late_up_detection_cortex) %>% 
    mutate(detected = junction_count >=2) %>% 
    group_by(sample) %>% 
    mutate(n_detected_selective = sum(detected)) %>% 
    ungroup() %>% 
    group_by(tdp_proteinopathy,n_late_up) %>% 
    summarize(mean_psi = mean(psi),mean_detected = mean(n_detected_selective)) %>% 
    ggplot(aes(x = as.factor(n_late_up),fill = tdp_proteinopathy,y = mean_detected)) + 
    geom_col(position = 'dodge')

# wilcoxin upregulated cord -----------------------------------------------
wilcox_test_res_cord = complete_table %>% 
    filter(junctions_coords %in% bind_table$junctions_coords) %>% 
    filter(grepl("Cord",rna_tissue_source_simplified)) %>% 
    mutate(condition = ifelse(tdp_proteinopathy %in% c("ALS-TDP","FTLD-TDP"),"TDP-43 proteinopathy","non-TDP-43-proteinopathy")) %>% 
    group_by(junctions_coords) %>% 
    nest() |> 
    mutate(wicox_test_result = map(data, ~ broom::tidy(wilcox.test(psi ~ condition, data = .x))),
           mean_path = map_dbl(data, ~ mean(.x$psi[.x$condition == "TDP-43 proteinopathy"], na.rm = TRUE)),
           mean_nonpath = map_dbl(data, ~ mean(.x$psi[.x$condition == "non-TDP-43-proteinopathy"], na.rm = TRUE)),
           ratio = log2(mean_path / mean_nonpath)) |> 
    unnest() %>% 
    distinct(junctions_coords,p.value,ratio,mean_path,mean_nonpath,ratio) %>% 
    left_join(bind_table,by = c("junctions_coords" = "junctions_coords")) %>% 
    arrange(-ratio) %>% 
    mutate(adj_pvalue = p.adjust(p.value))


# Boolean table construction ---------------------------------------------
library(eulerr)
detected_in_cortex_als_motor = complete_table %>% 
    filter(grepl("Cortex_Motor",rna_tissue_source_simplified)) %>% 
    filter(tdp_proteinopathy %in% c("Control","ALS-TDP","ALS-non-TDP")) %>%   
    filter(junction_count >=2) %>% 
    pull(junctions_coords) %>% unique()
detected_in_cortex_ftd_ft = complete_table %>% 
    filter(rna_tissue_source_simplified %in% c("Cortex_Frontal","Cortex_Temporal")) %>% 
    filter(tdp_proteinopathy %in% c("Control","FTLD-TAU","FTLD-FUS","FTLD-TDP")) %>%     
    filter(junction_count >=2) %>% 
    pull(junctions_coords) %>% unique()

bool_cortex = bind_table %>% 
    mutate(detected_in_cortex_als_motor = junctions_coords %in% detected_in_cortex_als_motor) %>% 
    mutate(detected_in_cortex_ftd_ft = junctions_coords %in% detected_in_cortex_ftd_ft) %>% 
    left_join(als_cortex_df,by = c("junctions_coords" = "paste_into_igv_junction")) %>% 
    mutate(specific_by_fpr_als_cortex = fpr_value < 0.1 & tpr_value > 0.10) %>% 
    select(junctions_coords,detected_in_cortex_als_motor,detected_in_cortex_ftd_ft,specific_by_fpr_als_cortex) %>% 
    left_join(ftd_cortex_df,by = c("junctions_coords" = "paste_into_igv_junction")) %>% 
    mutate(specific_by_fpr_ftd_ft = fpr_value < 0.1 & tpr_value > 0.10) %>% 
    replace(is.na(.), FALSE) %>% 
    select(junctions_coords,detected_in_cortex_als_motor,detected_in_cortex_ftd_ft,specific_by_fpr_als_cortex,specific_by_fpr_ftd_ft) %>% 
    unique()
    

bool_cord = bind_table %>% 
    left_join(cord_df,by = c("junctions_coords" = "paste_into_igv_junction")) %>% 
    mutate(specific_by_fpr_cord = fpr_value < 0.1 & tpr_value > 0.10) %>% 
    select(junctions_coords,specific_by_fpr_cord) %>% 
    replace(is.na(.), FALSE) %>% 
    unique()

bool_cortex %>% 
    left_join(bool_cord) %>% 
    select(-contains("detect")) %>% 
    left_join(rimod,by = c("junctions_coords" = "paste_into_igv_junction")) %>% 
    mutate(specific_by_fpr_rimod = fpr_value < 0.1 & tpr_value > 0.10) %>% 
    replace(is.na(.), FALSE) %>% 
    select(junctions_coords,contains('pecific')) %>% 
    remove_rownames() %>% 
    column_to_rownames('junctions_coords') %>% 
    euler(shape = "ellipse") %>% 
    plot(quantities = TRUE)


# Comparing cord and cortex ---------------------------------------------
tmp_2 = bool_cortex %>% 
    left_join(bool_cord) %>% 
    select(-contains("detect")) %>% 
    left_join(rimod,by = c("junctions_coords" = "paste_into_igv_junction")) %>% 
    mutate(specific_by_fpr_rimod = fpr_value < 0.1 & tpr_value > 0.10) %>% 
    filter(specific_by_fpr_ftd_ft == TRUE) %>% 
    arrange(-auc) %>% 
    left_join(master_table,by = c("junctions_coords" = "paste_into_igv_junction")) %>% 
    select(junctions_coords,gene,new_cate,max_psi_be2,max_psi_sh) %>% 
    filter(!is.na(new_cate)) %>% 
    mutate(new_max = case_when(max_psi_be2 == max_psi_sh ~ max_psi_sh,
                               max_psi_be2 == "" ~ max_psi_sh,
                               max_psi_sh == "" ~ max_psi_be2,
                               max_psi_be2 == 'Medium' ~ max_psi_sh,
                               max_psi_sh == 'Medium' ~ max_psi_be2)) %>% 
    filter(!is.na(new_max))

complete_table %>% 
    filter(rna_tissue_source_simplified == 'Cortex_Frontal' & tdp_proteinopathy == 'Control') %>% 
    filter(junctions_coords %in% tmp_2$junctions_coords) %>% 
    left_join(tmp_2) %>% 
    left_join(nygc_metadata %>% select(sample,age_at_death,cause_of_death,subject_group)) %>% 
    ungroup() %>% 
    mutate(age_at_death = ifelse(age_at_death == '90 or Older',100, age_at_death)) %>% 
    mutate(age_at_death = as.numeric(age_at_death)) %>% 
    filter(!is.na(age_at_death)) %>% 
    mutate(age_at_death = cut_number(age_at_death,3)) %>% 
    group_by(age_at_death) %>% 
    mutate(n_samp_total = n_distinct(sample)) %>% 
    ungroup() %>% 
    mutate(observed = junction_count >=2) %>% 
    group_by(junctions_coords,age_at_death,observed) %>% 
    summarize(n_samp = n_distinct(sample),gene,new_cate,new_max,age_at_death,n_samp_total) %>% 
    ungroup() %>% 
    filter(observed == TRUE) %>% 
    mutate(prop = n_samp/n_samp_total) %>% 
    arrange(-n_samp) %>% 
    unique() %>% 
    filter(new_cate != "Intermediate") %>% 
    mutate(plot_gene = ifelse(prop > 0.14,gene,NA_character_)) %>% 
    ggplot(aes(x = age_at_death,y = prop,fill = new_max)) + geom_boxplot() +
    ylab("N control samples junction detected") + 
    theme_minimal() +
    # scale_fill_manual(values = plot_sensitivity_colors) +
    # scale_color_manual(values = plot_sensitivity_colors) +
    geom_text_repel(aes(label = plot_gene,color = new_cate))

complete_table %>% 
    filter(rna_tissue_source_simplified == 'Cortex_Frontal' & tdp_proteinopathy == 'Control') %>% 
    filter(junctions_coords %in% tmp_2$junctions_coords) %>% 
    left_join(tmp_2) %>% 
    left_join(nygc_metadata %>% select(sample,age_at_death,cause_of_death,subject_group)) %>% 
    filter(gene_name == 'PFKP') %>% 
    mutate(age_at_death = ifelse(age_at_death == '90 or Older',100, age_at_death)) %>% 
    mutate(age_at_death = as.numeric(age_at_death)) %>% 
    mutate(pfkp_cut = junction_count > 2) %>% 
    ggplot(aes(x = pfkp_cut,y = age_at_death)) + 
    geom_boxplot()
# Comparing cord and cortex ---------------------------------------------
cerebellum_df %>% 
    left_join(cortex_df,by = 'paste_into_igv_junction',suffix = c("_cerebellum","_cortex")) %>% 
    left_join(bind_table) %>% 
    filter(fpr_value_cerebellum < 0.1) %>% arrange(-auc_cerebellum)
    ggplot(aes(x = fpr_value_cortex,
               y = fpr_value_cord)) + 
    geom_point() +
    geom_abline() + 
    geom_text_repel(aes(label = gene_name)) +
    theme_classic()

cortex_df %>% 
    left_join(cord_df,by = 'paste_into_igv_junction',suffix = c("_cortex","_cord")) %>% 
    left_join(bind_table) %>% 
    filter(tpr_value_cortex > 0.1 & tpr_value_cord > 0.1) %>% 
    filter(fpr_value_cortex < 0.1 | fpr_value_cord < 0.1) %>% 
    ggplot(aes(x = fpr_value_cortex,
               y = fpr_value_cord)) + 
    geom_point() +
    geom_abline() + 
    geom_text_repel(aes(label = gene_name)) +
    theme_classic()

cortex_df %>% 
    left_join(cord_df,by = 'paste_into_igv_junction',suffix = c("_cortex","_cord")) %>% 
    left_join(bind_table) %>% 
    filter(fpr_value_cortex < 0.1 & fpr_value_cord < 0.1) %>% 
    filter(auc_cortex > 0.6)

plot_junction("chr1:42456499-42457247",complete_table,"psi")
plot_junction("chr19:50506572-50511117",complete_table,measure = 'psi')

wide_counts %>% 
    filter(rna_tissue_source_simplified == '')
    select(sample,`chr19:50506572-50511117`) %>% 
    janitor::clean_names() %>% 
    filter(chr19_50506572_50511117 > 1)
cerebellum_df %>% 
    left_join(bind_table) %>% 
    left_join(wid)

# How which is more predictive, early or late? ------------------
cortex_predictive = cortex_df %>% 
    left_join(master_table %>% select(new_cate,paste_into_igv_junction,cryptic,gene)) %>% 
    filter(!is.na(new_cate)) %>% 
    filter(cryptic == 'yes') %>% 
    mutate(predictive = fpr_value < 0.1 & tpr_value > 0.1) %>% 
    filter(predictive == TRUE) %>% 
    pull(paste_into_igv_junction)

simple_cortex = complete_table %>% 
    filter(grepl("Cortex",rna_tissue_source_simplified)) %>% 
    filter(junctions_coords %in% cortex_predictive) %>% 
    group_by(sample) %>% 
    reframe(simple_score = sum(psi),tdp_proteinopathy) %>% 
    unique()
simple_cortex %>% 
    ggplot(aes(y = simple_score,x = tdp_proteinopathy)) + 
    geom_violin()
library(skimr)

simple_cortex %>%
    group_by(tdp_proteinopathy) %>%
    summarize(
        mean = mean(simple_score, na.rm = TRUE),
        sd = sd(simple_score, na.rm = TRUE),
        median = median(simple_score, na.rm = TRUE),
        min = min(simple_score, na.rm = TRUE),
        max = max(simple_score, na.rm = TRUE),
        n = n()
    )

simple_cortex = simple_cortex %>% 
    mutate(numeric_prediction = ifelse(tdp_proteinopathy %in% c("ALS-TDP","FTLD-TDP"), 1, 0))
    
pROC_obj <- roc_(data = simple_cortex, "numeric_prediction", 'simple_score', quiet = TRUE)

low_expressors = simple_cortex %>% filter(tdp_proteinopathy == 'FTLD-TDP') %>% slice_min(simple_score,n = 32)
high_expressors = simple_cortex %>% filter(tdp_proteinopathy == 'FTLD-TDP') %>% slice_max(simple_score,n = 32)



# How often are early/late events dectected across NYGC? ------------------

these_counts = tmp %>% 
    filter(gene_name %in% these_genes) %>% 
    mutate(sample_id_rna = gsub(".SJ.out","",V4)) %>% 
    mutate(junction_coords = as.character(glue::glue("{V1}:{V2}-{V3}"))) %>% 
    left_join(nygc_metadata,by ='sample_id_rna') %>% 
    filter(!is.na(tdp_proteinopathy)) %>% 
    distinct(junction_coords,V5,sample_id_rna)

these_counts = these_counts %>%
    complete(
        junction_coords, sample_id_rna,
        fill = list(V5 = 0)
    ) %>% 
    left_join(subset_sneak) %>% 
    left_join((nygc_metadata %>% select(sample_id_rna,tdp_proteinopathy,rna_tissue_source_simplified)))

plot_junction = function(junc,plotin_table = these_counts){

    tmp_tbl = plotin_table[junction_coords == junc]
    if(nrow(tmp_tbl) >0){
        tmp_tbl = tmp_tbl %>% 
            mutate(disease = tdp_proteinopathy) %>% 
            mutate(disease = case_when(disease %in% c("FTLD-TAU", "FTLD-FUS") ~ "FTLD-non-TDP",
                                       TRUE ~ disease))
        vals = c(`ALS-TDP` = "#E1BE6A", Control = "#40B0A6", `FTLD-TDP` = "#E1BE6A", 
                 `ALS\nnon-TDP` = "#408A3E", `FTLD\nnon-TDP` = "#408A3E")
        
        # gene_name = plotin_table[junction_coords == junc,unique(gene_name)]
        # junc_type = plotin_table[junction_coords == junc,unique(type)]
        # plot_title = glue::glue("{gene_name} - {junc} - {junc_type}")
        
        plt = tmp_tbl |>
            filter(disease != "Unknown") %>% 
            filter(!rna_tissue_source_simplified %in% c("Other","Choroid","Hippocampus","Liver","Cortex_Sensory","Cortex_Occipital","Medulla")) |> 
            mutate(disease =  gsub("-n","\nn",disease)) |> 
            mutate(disease = fct_relevel(disease,"Control", "ALS\nnon-TDP","ALS-TDP", "FTLD\nnon-TDP", "FTLD-TDP")) |>
            mutate(rna_tissue_source_simplified = fct_relevel(rna_tissue_source_simplified, "Cerebellum", after = Inf)) %>% 
            ggplot(aes(x = disease, y = psi, fill = disease)) +
            geom_boxplot(outlier.colour = NA) +
            geom_jitter(height = 0,alpha = 0.7, pch = 21) +
            facet_wrap(~rna_tissue_source_simplified) +
            scale_fill_manual(values = vals)  +
            ylab("N spliced reads") +
            xlab("") +
            # ggtitle(plot_title) +
            ggpubr::theme_pubr() +
            theme(legend.position = 'none') 
        
        
        return(plt) 
    }else{
        print("no found")
        print(junc)
    }

    
}

# How often are dark cryptics detected in NYGC? ---------------------------------------------
table_setup = function(df){
    wide_psi = df %>% 
        mutate(numeric_prediction = ifelse(tdp_proteinopathy %in% c("ALS-TDP","FTLD-TDP"), 1, 0)) %>% 
        distinct(sample,junction_coords,numeric_prediction,psi,rna_tissue_source_simplified,tdp_proteinopathy) %>% 
        pivot_wider(values_from = 'psi',
                    names_from = 'junction_coords',values_fill = 0) %>% 
        as.data.table() 
    
    wide_counts = df %>% 
        mutate(numeric_prediction = ifelse(tdp_proteinopathy %in% c("ALS-TDP","FTLD-TDP"), 1, 0)) %>% 
        distinct(sample,junction_coords,numeric_prediction,rna_tissue_source_simplified,tdp_proteinopathy,value) %>% 
        pivot_wider(values_from = 'value',
                    names_from = 'junction_coords',values_fill = 0) %>% 
        as.data.table()
    
    complete_psi = wide_psi %>% select(-numeric_prediction) %>% melt() %>% 
        dplyr::rename(junction_coords = variable,psi = value)
    
    complete_counts = wide_counts %>% select(-numeric_prediction) %>% melt() %>% 
        dplyr::rename(junction_coords = variable,junction_count = value)
    
    complete_table = complete_psi %>% 
        left_join(complete_counts) 
    
    return(complete_table)
}

big_data_chx <- read.csv(here::here('data','nmd_or_not.csv'))
big_data_chx = big_data_chx %>% 
    separate(paste_into_igv_junction,remove = FALSE,into = c("chr",'start','end'),
             convert = TRUE)

chx_dark = wide_chx |> 
    filter(normal_cryptic_chx == FALSE & chx_cryptic == TRUE) |> 
    pull(paste_into_igv_junction) |> unique()
wide_chx = big_data_chx |> 
    select(paste_into_igv_junction,.id,mean_psi_per_lsv_junction,color_gene_name) |> 
    pivot_wider(names_from = '.id',values_from = 'mean_psi_per_lsv_junction',
                values_fn = max) %>% 
    unique() %>% 
    mutate(normal_cryptic_chx = Control_Control < 0.05 &  Control_TDP43KD > 0.1) |> 
    mutate(chx_cryptic = Control_Control < 0.05 &  Cycloheximide_TDP43KD > 0.1) |> 
    mutate(chx_effect = case_when((Cycloheximide_TDP43KD - Control_TDP43KD) > 0.05 ~ "NMD rescued",
                                  Control_TDP43KD > 0.94 ~ "PSI > 95% at baseline",
                                  TRUE ~ "non-NMD rescued")) %>% 
   left_join(big_data_chx %>% distinct(paste_into_igv_junction,gene_name))


nygc_chx = nygc_reg %>% 
    filter(chrom %in% big_data_chx$chr, start %in% big_data_chx$start, end %in% big_data_chx$end) %>% 
    collect()

nygc_chx = nygc_chx %>% 
    filter(sample %in% neural_meta$sample) %>% 
    left_join(neural_meta %>% select(sample,rna_tissue_source_simplified,tdp_proteinopathy)) %>% 
    unique()
nygc_chx = nygc_chx %>% 
    mutate(junction_coords = as.character(glue::glue("{chrom}:{start}-{end}"))) 

nygc_chx_complete = table_setup(nygc_chx)
nygc_chx_wide = nygc_chx_complete %>% 
    mutate(numeric_prediction = ifelse(tdp_proteinopathy %in% c("ALS-TDP","FTLD-TDP"), 1, 0)) %>% 
    distinct(sample,junction_coords,numeric_prediction,psi,rna_tissue_source_simplified,tdp_proteinopathy) %>% 
    pivot_wider(values_from = 'psi',
                names_from = 'junction_coords',values_fill = 0) %>% 
    as.data.table() 

cortex_dark = nygc_chx_wide %>% 
    filter(grepl("Cortex",rna_tissue_source_simplified)) %>% 
    process_psi_data(not_included = c("sample",'numeric_prediction',"rna_tissue_source_simplified","tdp_proteinopathy"))


nygc_chx_complete %>% 
    mutate(tdp_path = ifelse(tdp_proteinopathy %in% c("ALS-TDP","FTLD-TDP"), 'path', 'not_path')) %>% 
    group_by(tdp_path) %>% 
    
    mutate(n_sample = n_distinct(sample)) %>% 
    mutate(detected = junction_count > 2) %>% 
    group_by(tdp_path,junction_coords) |> 
    summarise(n_obs = sum(detected)) |> 
    ungroup() %>% 
    unique() %>% 
    pivot_wider(values_from = 'n_obs',
                names_from = 'tdp_path')
# Compare NMD vs Not cryptic detection NYGC ---------------------------------------------


# Ranking the cortex samples by their expression of the selective  --------
selective_cortex = cortex_df %>% 
    filter(fpr_value < 0.1 & tpr_value > 0.1) %>% 
    pull(paste_into_igv_junction)

complete_table %>% 
    filter(junctions_coords %in% selective_cortex) %>% 
    group_by(sample) %>% 
    mutate(sum_psi = sum(psi)) %>% 
    ungroup() %>% 
    left_join(nygc_metadata[,.(sample,id,simplified_mutations)]) %>% 
    filter(id != "") %>% 
    arrange(-sum_psi) %>% 
    distinct(sample,rna_tissue_source_simplified,tdp_proteinopathy,sum_psi,id,simplified_mutations) %>% 
    slice_max(prop = 0.15,sum_psi)

complete_table %>% 
    filter(junctions_coords %in% selective_cortex) %>% 
    filter(rna_tissue_source_simplified == 'Cerebellum') %>% 
    arrange(-psi) %>% 
    filter(tdp_proteinopathy == 'FTLD-TDP')


# generate postmortem summary ---------------------------------------------
n_sample_junc_nygc = n_sample_junc_nygc %>% 
    pivot_wider(names_from = 'variable',
                values_from = c("value","frac_ob","n_samples")) %>% 
    filter(!is.na(value_tdp_path))

rbind(cord_df %>% 
    mutate(group = "Cord"),
    ftd_cortex_df %>% mutate(group = 'FTD frontal/temporal'),
    als_cortex_df %>% mutate(group = 'Motor Cortex')) %>% 
    left_join(n_sample_junc_nygc) %>% 
    unique() %>% 
    fwrite('~/Desktop/nygc_postmortem_summary.csv')
