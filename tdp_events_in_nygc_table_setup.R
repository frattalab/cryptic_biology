library(data.table)
library(tidyverse)
library(tidyr)
library(ggpubr)
library(ggplot2)
#lil function that will read in the counts or tpm table, select the samples given
#as a list and then remove it out
gene_count_reader = function(samples, type = "tpm", genes = ""){
    if(type == "counts"){
        full_count_file = fread(file.path(here::here(),"data/rsem_counts_nygc.csv"))
    }else if(type == "tpm"){
        full_count_file = fread(file.path("/Users/annaleigh/Documents/GitHub/tdp_43_psi_rankings/","data/rsem_tpm_nygc.csv"))
        full_count_file = full_count_file |> mutate(gene = gsub("\\..*", "", gene))
    }
    return_me = full_count_file |> 
        dplyr::select(gene,samples) 
    if(genes != ""){
        return_me = return_me |> 
            filter(gene %in% genes)
    }
    return(return_me)
}

plot_junction = function(junc,plotin_table = spliced_counts){

    vals = c(`ALS-TDP` = "#E1BE6A", Control = "#40B0A6", `FTD-TDP` = "#E1BE6A", 
      `ALS\nnon-TDP` = "#408A3E", `FTD\nnon-TDP` = "#408A3E")
    
    gene_name = plotin_table[paste_into_igv_junction == junc,unique(gene)]
    junc_type = plotin_table[paste_into_igv_junction == junc,unique(type)]
    plot_title = glue::glue("{gene_name} - {junc} - {junc_type}")
    
    plt = plotin_table[paste_into_igv_junction == junc] |>
        filter(!tissue_clean %in% c("Other","Choroid","Hippocampus","Liver","Sensory_Cortex","Occipital_Cortex")) |> 
        mutate(disease =  gsub("-n","\nn",disease)) |> 
        mutate(disease = fct_relevel(disease,"Control", "ALS\nnon-TDP","ALS-TDP", "FTD\nnon-TDP", "FTD-TDP")) |>
        ggplot(aes(x = disease, y = spliced_reads, fill = disease)) +
        geom_boxplot(outlier.colour = NA) +
        geom_jitter(height = 0,alpha = 0.7, pch = 21) +
        facet_wrap(~tissue_clean, scales = "free_y") +
        scale_fill_manual(values = vals)  +
        ylab("N spliced reads") +
        xlab("") +
        ggtitle(plot_title) +
        ggpubr::theme_pubr() +
        theme(legend.position = 'none') +
        theme(text = element_text(size = 16)) 
        
        
    
    return(plt)
    
}
create_formated_metadata = function(dt){
    dt[,disease := ifelse(disease == "FTD" & pathology %in% c("FTD-TAU","FTD-FUS"),"FTD-non-TDP",disease)]
    dt[,disease := ifelse(disease == "FTD","FTD-TDP",disease)]
    
    dt[,disease := ifelse(grepl("ALS",disease), "ALS-TDP",disease)]
    dt[,disease := ifelse(grepl("ALS",disease) & mutations %in% c("FUS","SOD1"),"ALS-non-TDP",disease)]
    

    dt = dt |>  mutate(disease_tissue = case_when(grepl("FTD",disease_full) & grepl("Front|Temp",tissue_clean)  ~ T,
                                                   (grepl("ALS",disease) & grepl("Cord|Motor|Front|Temp",tissue_clean))  ~ T,
                                                   (grepl("Occipital|Sensory",tissue_clean)) ~ F,
                                                   (grepl("Control",disease) & grepl("Cord|Cortex",tissue_clean)) ~ T,
                                                   TRUE ~ F)) 
    dt = dt |> 
        mutate(tdp_path = case_when((disease %in% c("ALS-TDP","FTD-TDP") & disease_tissue == T) ~ 'path',
                                    T ~ "not_path"))
    dt = unique(dt[disease != "Other"])
    return(dt)
}
meta_data_full = fread(here::here('NYGC_all_RNA_samples_support.tsv'))
# meta_data_full = fread("/Users/annaleigh/Documents/GitHub/tdp_43_psi_rankings/tdp_43_psi_rankings/NYGC_all_RNA_samples_support.tsv")
meta_data = meta_data_full[,.SD,.SDcols = c("sample", "individual", "region", "tissue", "tissue_clean",
                                       "disease", "disease_full", "age",
                                       "onset", "mutations", "pathology")]
raw_counts = unique(fread(here::here("data/everything_except_liu_ferguson_nygc_counts.aggregated.clean.annotated.bed")))
# raw_counts = unique(fread("/Users/annaleigh/Documents/GitHub/cryptic_biology/data/everything_except_liu_ferguson_nygc_counts.aggregated.clean.annotated.bed"))
# raw_counts = unique(fread("/Users/annaleigh/Documents/GitHub/tdp_43_psi_rankings/tdp_43_psi_rankings/everything_and_the_kitchen_sink_nygc_counts.aggregated.clean.annotated.bed"))
#raw_counts = unique(fread("/Users/annaleigh/Documents/GitHub/tdp_43_psi_rankings/only_the_neuronal_linesaggregated.clean.annotated.bed"))
# raw_counts = unique(fread("/Users/annaleigh/Documents/GitHub/tdp_43_psi_rankings/nygc_all_recent_cryptics.bedaggregated.clean.annotated.bed"))
# raw_counts = unique(fread("/Users/annaleigh/Documents/GitHub/tdp_43_psi_rankings/raw_counts_cryptic_in_any_kdaggregated.clean.annotated.bed"))
raw_counts[,sample := gsub(".SJ.out","",V4)]
raw_counts[, gene := tstrsplit(V7,"\\|")[[1]]]
raw_counts[, type := tstrsplit(V7,"\\|")[[2]]]
raw_counts[, n_invitro_observed := tstrsplit(V7,"\\|")[[3]]]

raw_counts[,paste_into_igv_junction := paste0(V1,":",V2,"-",V3)]

junction_info = unique(raw_counts[,.(paste_into_igv_junction,gene,type,V6,n_invitro_observed)])
setnames(junction_info, "V6","strand")
temp = complete(raw_counts[,.(paste_into_igv_junction,sample,V5)], 
         paste_into_igv_junction, sample,fill = list(V5 = 0)) |> as.data.table()

setnames(temp,"V5","spliced_reads")

temp = temp |> left_join(junction_info)

temp = temp |> left_join(meta_data)

spliced_counts = create_formated_metadata(temp)

# Observed >= 2 spliced reads
# Highly expressed >= 10 TDP-43 path tissues
# Selective <= 5 non-TDP path
total_path = spliced_counts |> group_by(tdp_path) |> 
  summarize(n_samp = n_distinct(sample)) |> filter(tdp_path == "path") |> pull()

total_notpath = spliced_counts |> group_by(tdp_path) |> 
  summarize(n_samp = n_distinct(sample)) |> filter(tdp_path == "not_path") |> pull()

expression_by_pathology = spliced_counts |> 
  dplyr::select(disease_tissue,spliced_reads,tdp_path,sample,paste_into_igv_junction,n_invitro_observed) |> 
  unique() |> 
  filter(disease_tissue == TRUE) |> 
  mutate(observed = spliced_reads >= 2) |> 
  group_by(tdp_path,paste_into_igv_junction) |> 
  summarise(n_obs = sum(observed)) |> 
  ungroup() |> 
  pivot_wider(values_from = 'n_obs',
              names_from = 'tdp_path') |> 
  left_join(junction_info)  |> 
  as.data.table() |> 
  mutate(fraction_not_path = not_path / total_notpath) |> 
  mutate(fraction_path = path / total_path)


# Find how many of the 11K events were detected ---------------------------


# fwrite(expression_by_pathology,here::here('expression_by_pathology_counts_just_neuronal.csv'))

expression_by_pathology_tissue_sep = spliced_counts |> 
    select(disease_tissue,spliced_reads,tdp_path,sample,paste_into_igv_junction,n_invitro_observed) |> 
    unique() |> 
    filter(disease_tissue == TRUE) |> 
    mutate(observed = spliced_reads >= 2) |> 
    group_by(tdp_path,paste_into_igv_junction) |> 
    summarise(n_obs = sum(observed)) |> 
    ungroup() |> 
    pivot_wider(values_from = 'n_obs',
                names_from = 'tdp_path') |> 
    left_join(junction_info)  |> 
    as.data.table() |> 
    mutate(fraction_not_path = not_path / total_notpath) |> 
    mutate(fraction_path = path / total_path)

# 
# # remove cerebellum splicing ----------------------------------------------
# cerebellum_splicing_fraction = spliced_counts |> 
#     filter(tissue == "Cerebellum") |> 
#     filter(spliced_reads >= 2) |> 
#     mutate(total_sample = n_distinct(sample)) |> 
#     group_by(paste_into_igv_junction) |> 
#     summarize(n_sample_obs = n_distinct(sample),total_sample) |> 
#     ungroup() |> unique() |> 
#     mutate(fraction_cerebellum_observed = n_sample_obs / total_sample) |> 
#     arrange(-fraction_cerebellum_observed)
# 
# 
potential_new_selective = expression_by_pathology |>
    filter(path >= 10) |>
    filter(not_path <= 5)
    # filter(fraction_path >= 0.01) |>
    # filter(fraction_not_path <= 0.002)

selective_junctions  = potential_new_selective$paste_into_igv_junction

potential_new_selective |> 
  mutate(obs_tot = path + not_path) |> 
  mutate(gene_name = glue::glue("{gene}-{type}")) |> 
  # select(gene_name,not_path,path,obs_tot) |> 
  select(gene_name,fraction_not_path,fraction_path,obs_tot) |> 
  slice_max(obs_tot,n = 50) |> 
  melt(id.vars = c('obs_tot',"gene_name")) |> 
  mutate(gene_name = fct_reorder(gene_name,obs_tot)) |> 
  ggplot(aes(x = gene_name,y = value,fill = variable)) +
  ggpubr::theme_pubr() + 
  geom_col(position = 'dodge') + 
  coord_flip() + 
  scale_fill_manual(values = c("#408A3E","#E1BE6A"),
                    labels = c("non-TDP Path", "TDP-43 Path"))  +
  # ylab("N tissue samples observed (>= 3 spliced reads)") + 
  ylab("Fraction tissue samples observed (>= 2 spliced reads)") + 
  xlab(NULL) +
  labs(fill = "Tissue pathology type observed in") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) 


  
not_nmd = unique(c("AARS1", "ACTL6B", "CDO1", "CELSR3", "CERS4", 
                   "DENND2B", "EPB41L4A", "HDGFL2", "ICA1", "KCNIP2", 
                   "MNAT1", "MYO18A", "NECAB2", "PAOX", "PRELID3A", 
                   "PTPRD", "PTPRD", "PTPRZ1", "PTPRZ1", "PXDN", "SEMA6D", 
                   "SLC24A3", "SLC24A3", "UNC13B",
                   "ZNF423", "PXDN", "DNM1", "XPO4"))  

selective_genes = unique(potential_new_selective$gene)



controls_with_pfkp = spliced_counts[paste_into_igv_junction == 'chr10:3099819-3101365' & 
                   spliced_reads > 2 & disease == "Control",individual] |> unique()

meta_data |> 
    filter(disease == "Control") |> 
    select(individual,age) |> 
    unique() |> 
    mutate(has_pfkp = ifelse(individual %in% controls_with_pfkp,"PFKP cryptic found", "PFKP cryptic not found")) |> 
    ggplot(aes(x = has_pfkp, y = age,fill = has_pfkp)) + 
    geom_boxplot() + 
    stat_cor() + 
    ggpubr::stat_compare_means() + 
    theme_pubr() + 
    ggtitle("Control samples with PFKP cryptic exon are older") + 
    xlab(NULL) +
    ylab("Individual age") + 
    theme(legend.title = element_blank()) +
    scale_fill_brewer(palette = 2)




library(clusterProfiler)
library('org.Hs.eg.db')
ens_selective = annotables::grch38 |> 
    filter(symbol %in% c("AARS",selective_genes)) |> 
    arrange(symbol) |> 
    pull(ensgene)

selective_go = clusterProfiler::enrichGO(ens_selective,
                                              keyType = 'ENSEMBL',
                                              OrgDb = org.Hs.eg.db,
                                              ont = 'ALL',
                                              readable = TRUE)

dotplot(selective_go)
cnetplot(selective_go)


# sj = potential_new_selective |> arrange(-path) |> pull(paste_into_igv_junction)
# pdf(width = 18, height =  8.27, file = 'selective_junctions_nygc_expression.pdf')
# for (j in sj){
#     p = plot_junction(j)
#     print(p)
# }
# dev.off()
# unique(cryptics_sig[n_significant_data == 1,.(paste_into_igv_junction,comparison)]) |> 
#     left_join(potential_new_selective) |> 
#     filter(!is.na(fraction_path)) |> 
#     arrange(-fraction_path)

  
# # reading in the liu events -----------------------------------------------
# liu_counts = unique(fread("/Users/annaleigh/Documents/GitHub/tdp_43_psi_rankings/patients_liu_parseaggregated.clean.annotated.bed"))
# liu_counts[,sample := gsub(".SJ.out","",V4)]
# liu_counts[,junction := paste0(V7,"|",V1,":",V2,"-",V3)]
# 
# temp = dcast(unique(liu_counts), junction ~ sample, value.var = "V5",fill = 0)
# liu_counts = unique(melt(temp, id.vars = 'junction'))
# setnames(liu_counts,"variable","sample")
# liu_counts = unique(liu_counts[meta_data, on="sample"])
# liu_counts = create_formated_metadata(liu_counts)
# 
# full_counts = unique(rbind(liu_counts, raw_counts))
# full_counts[, junction_coords := tstrsplit(junction,"\\|")[[3]]]
# full_counts[, gene := tstrsplit(junction,"\\|")[[1]]]
# 
# full_counts[, type := tstrsplit(gene,"\\:")[[2]]]
# full_counts[is.na(type),type := tstrsplit(junction,"\\|")[[2]]]
# 
# 
# full_counts$junction = NULL
# full_counts[,gene := gsub("\\:.*","\\1",gene)]
# full_counts = unique(full_counts)
# 
# full_counts[,observed := value >= 3]
# full_counts[,n_observed := sum(observed), by = .(junction_coords)]
# double_juncs = unique(full_counts[,.(junction_coords,gene,type)])[,.N, by = junction_coords][N > 1]$junction_coord |> unique()
# full_counts[junction_coords %in% double_juncs, n_observed := n_observed / 2]
# # widen to find selective events ------------------------------------------
# 
# wide_counts_pathology = dcast(full_counts[tissue_clean != "Cerebellum" &
#                                               disease_tissue == T,
#                                           mean(value,na.rm = T), by = .(tdp_path,junction_coords)], junction_coords ~ tdp_path, value.var = "V1", fill = 0)
# wide_counts_pathology[,log2readsTDPpath :=log2((path + 1) / (not_path + 1))]
# wide_counts_pathology = unique(wide_counts_pathology[unique(full_counts[,.(gene,type,junction_coords,n_observed)]), on = "junction_coords"])
# 
# wide_counts_pathology[n_observed >= 10 &log2readsTDPpath > 0  & not_path <0.5][order(-log2readsTDPpath)]
# 
# wide_selective = fread("/Users/annaleigh/Documents/GitHub/tdp_43_psi_rankings/specific_hits.csv")
# selective_junctions  = wide_selective$junction_coords
# 
# 
# # wide table by number of tissue samples observed at least 3 junct --------
# 
# selective_observation  = spliced_counts[!tissue_clean %in% c("Hippocampus","Cerebellum","Occipital_Cortex")][,sum(observed),by = .(junction_coords,tdp_path,gene)]
# 
# selective_observation = dcast(selective_observation[,mean(V1,na.rm = T), by = .(tdp_path,junction_coords)], junction_coords ~ tdp_path, value.var = "V1", fill = 0)
# 
# 
# selective_observation = selective_observation |>
#     left_join(unique(full_counts[,.(junction_coords,gene)])) |>
#     arrange(-path) |> unique()


# count how many selective junctions appears per sample -------------------
double_juncs = unique(spliced_counts[,.(paste_into_igv_junction,gene,type)])[,.N, by = paste_into_igv_junction][N > 1]$paste_into_igv_junction |> unique()

selective_counts = spliced_counts[paste_into_igv_junction %in% potential_new_selective$paste_into_igv_junction]
selective_counts[,observed := spliced_reads >= 2]
cryptics_observed = selective_counts |> 
    filter(spliced_reads >=2) |> 
    group_by(sample) |> 
    mutate(n_cryptics_observed = n_distinct(paste_into_igv_junction)) |> 
    as.data.table() |> 
    dplyr::select(sample,tissue_clean,disease,disease_tissue,age,n_cryptics_observed) |> 
    unique() |> 
    arrange(-n_cryptics_observed) 

selective_counts[,observed_in_patients := sum(observed),by = paste_into_igv_junction]
selective_counts[paste_into_igv_junction %in% double_juncs, observed_in_patients := observed_in_patients / 2]

# What correlates with the number of cryptic events observed? -------------

cryptics_observed[disease_tissue == TRUE] |> 
    ggplot(aes(x = age, y = n_cryptics_observed)) + 
    geom_point() + 
    geom_rug() + 
    facet_wrap(~tissue_clean,scales = 'free') + 
    stat_cor(method = 'spearman',
             cor.coef.name = c("rho"),
    ) + 
    ggpubr::theme_pubr() + 
    geom_smooth(method = 'lm')


cryptics_observed[disease_tissue == TRUE] |> 
    left_join(meta_data_full) |> 
    filter(!is.na(rin)) |> 
    ggplot(aes(x = rin, y = n_cryptics_observed)) + 
    geom_point() + 
    geom_rug() + 
    # facet_wrap(~tissue_clean,scales = 'free') +
    stat_cor(method = 'spearman',
             cor.coef.name = c("rho"),
    ) + 
    ggpubr::theme_pubr() + 
    geom_smooth(method = 'lm')


cryptics_observed[disease_tissue == TRUE] |> 
    left_join(meta_data_full) |> 
    filter(!is.na(rin)) |> 
    ggplot(aes(x = pmi, y = n_cryptics_observed)) + 
    geom_point() + 
    geom_rug() + 
    facet_wrap(~tissue_clean,scales = 'free') +
    stat_cor(method = 'spearman',
             cor.coef.name = c("rho"),
    ) + 
    ggpubr::theme_pubr() + 
    geom_smooth(method = 'lm')


cryptics_observed[disease_tissue == TRUE] |> 
    left_join(meta_data_full) |> 
    filter(!is.na(platform)) |> 
    ggplot(aes(x = platform, y = n_cryptics_observed)) + 
    geom_boxplot() + 
    geom_rug() + 
    facet_wrap(~tissue_clean) +
    ggpubr::stat_compare_means() + 
    ggpubr::theme_pubr()

cell_type = fread("https://raw.githubusercontent.com/frattalab/unc13a_cryptic_splicing/main/data/All_1917_samples_Darmanis_dtangle_deconv.tsv")
library_sizes = fread('https://raw.githubusercontent.com/frattalab/cryptic_biology/f7f01522ef8d3285ac0b952d5aee6affaabbba75/data/nygc_library_sizes.csv?token=GHSAT0AAAAAACAKYNAWLBRELBCM3F76VVD2ZBCY2LQ')


# Library depth is the most correlated with N cryptic observed  --------

cryptics_observed[disease_tissue == TRUE] |> 
    left_join(library_sizes,by = c("sample" = "sample_id")) |> 
    ggplot(aes(x = library_size,
               y = n_cryptics_observed)) + 
    facet_wrap(~tissue_clean,scales = 'free') +
    geom_point() + 
    geom_rug() + 
    stat_cor() + 
    ggpubr::theme_pubr() + 
    geom_smooth(method = 'lm')
    
cryptics_observed[disease_tissue == TRUE] |> 
    left_join(meta_data_full,by = 'sample') |> 
    select(tissue,individual,disease.x,n_cryptics_observed) |> unique() |> 
    pivot_wider(names_from = "tissue",
                values_from = "n_cryptics_observed",
                id_cols = c("individual","disease.x"),values_fn = mean)



