library(readxl)
library(tidyverse)
library(arrow)
library_depth = read_parquet('data/nygc_library_sizes.parquet')
rsem_tpm = fread(file.path(here::here(),"data/rsem_tpm_nygc.csv"))
df_cat_shsy5y_curve <- read_excel("data/df_cat_shsy5y_curve.xlsx")
selective_counts = read_parquet("data/selective_counts.parquet")
df_cat_shsy5y_curve |> 
    mutate(paste_into_igv_junction = glue::glue("{chr}:{start}-{end}")) |> 
    filter(is_cryptic == 'Cryptic') |> 
    select(paste_into_igv_junction,Type) |> 
    left_join(selective_counts) |> 
    filter(!is.na(spliced_reads)) |>
    filter(disease_tissue == TRUE) |> 
    ggplot(aes(y = spliced_reads,
               x = Type)) + 
    geom_boxplot() + 
    facet_wrap(~tissue_clean)



selective_curve = df_cat_shsy5y_curve |> 
    mutate(paste_into_igv_junction = glue::glue("{chr}:{start}-{end}")) |> 
    filter(is_cryptic == 'Cryptic') |> 
    select(paste_into_igv_junction,Type) |> 
    left_join(selective_counts) |> 
    filter(!is.na(spliced_reads)) |>
    filter(disease_tissue == TRUE) |> 
    mutate(observed = spliced_reads >= 2)
    


selective_curve |>  group_by(sample,Type) |> 
    summarize(n_ob = sum(observed)) |> 
    ungroup() |> 
    pivot_wider(names_from = 'Type',
                values_from = 'n_ob') |> 
    mutate(has_late = ifelse(Late >=1,'late detected','no late')) |> 
    mutate(has_early = ifelse(Early >=1,'early detected','no early')) |> 
    filter(Early + Late != 0) |> 
    select(has_early,has_late) |> 
    table() |> as.data.frame() |> 
    ggplot(aes(y = Freq,
               fill = has_early,
               x = has_late)) + 
    geom_col()
    filter(Late >=1 & Early <1)

have_early_detected = selective_curve |>  group_by(sample,Type) |> 
    summarize(n_ob = sum(observed)) |> 
    ungroup() |> 
    pivot_wider(names_from = 'Type',
                values_from = 'n_ob') |> 
    filter(Early >=1) |> pull(sample)  

have_late_detected = selective_curve |>  group_by(sample,Type) |> 
    summarize(n_ob = sum(observed)) |> 
    ungroup() |> 
    pivot_wider(names_from = 'Type',
                values_from = 'n_ob') |> 
    filter(Late >=1) |> pull(sample)  
    
have_late_no_early = selective_curve |>  group_by(sample,Type) |> 
    summarize(n_ob = sum(observed)) |> 
    ungroup() |> 
    pivot_wider(names_from = 'Type',
                values_from = 'n_ob') |> 
    filter(Late >=1 & Early <1) |> pull(sample)
    
library_depth |> 
    mutate(hlne = sample_id %in% have_late_no_early) |> 
    ggplot(aes(y = library_size,
               x = hlne)) + 
    geom_violin()

meta_data |> 
    filter(sample %in% have_late_detected) |> 
    mutate(hlne = sample %in% have_late_no_early) |> 
    count(tissue_clean, hlne) |> 
    group_by(tissue_clean) %>%
    mutate(freq = n / sum(n)) |> 
    ungroup() |> 
        filter(tissue_clean %in% c("Frontal_Cortex", "Lumbar_Spinal_Cord", "Cervical_Spinal_Cord", 
                                   "Motor_Cortex", "Thoracic_Spinal_Cord", "Temporal_Cortex")) |> 
    ggplot(aes(x = tissue_clean,
               fill = hlne,
               y = freq)) +
    geom_col() + 
    coord_flip() +
    scale_fill_manual('Have a late cryptic but no early',values=c('blue', 'red')) +
    ylab("Tissue") + 
    xlab("Frequency")
    
meta_data |> 
    filter(sample %in% have_late_detected) |> 
    mutate(hlne = sample %in% have_late_no_early) |> 
    count(tissue_clean, hlne) |> 
    group_by(tissue_clean) %>%
    mutate(freq = n / sum(n)) |> 
    ungroup() |> 
    filter(tissue_clean %in% c("Frontal_Cortex", "Lumbar_Spinal_Cord", "Cervical_Spinal_Cord", 
                               "Motor_Cortex", "Thoracic_Spinal_Cord", "Temporal_Cortex")) |> 
    ggplot(aes(x = tissue_clean,
               fill = hlne,
               y = n)) +
    geom_col() + 
    coord_flip() +
    scale_fill_manual('Have a late cryptic but no early',values=c('blue', 'red')) +
    ylab("Tissue") + 
    xlab("N")


meta_data |> 
    filter(sample %in% have_early_detected) |> 
    mutate(hlne = sample %in% have_late_detected) |> 
    count(tissue_clean, hlne) |> 
    group_by(tissue_clean) %>%
    mutate(freq = n / sum(n)) |> 
    ungroup() |> 
    filter(tissue_clean %in% c("Frontal_Cortex", "Lumbar_Spinal_Cord", "Cervical_Spinal_Cord", 
                               "Motor_Cortex", "Thoracic_Spinal_Cord", "Temporal_Cortex")) |> 
    ggplot(aes(x = tissue_clean,
               fill = hlne,
               y = n)) +
    geom_col() + 
    coord_flip() +
    scale_fill_manual('Have a late cryptic',values=c('blue', 'red')) +
    ylab("Tissue") + 
    xlab("N") + 
    ggtitle("Samples with an early cryptic")

meta_data |> 
    filter(sample %in% have_early_detected) |> 
    mutate(hlne = sample %in% have_late_detected) |> 
    count(tissue_clean, hlne) |> 
    group_by(tissue_clean) %>%
    mutate(freq = n / sum(n)) |> 
    ungroup() |> 
    filter(tissue_clean %in% c("Frontal_Cortex", "Lumbar_Spinal_Cord", "Cervical_Spinal_Cord", 
                               "Motor_Cortex", "Thoracic_Spinal_Cord", "Temporal_Cortex")) |> 
    ggplot(aes(x = tissue_clean,
               fill = hlne,
               y = freq)) +
    geom_col() + 
    coord_flip() +
    scale_fill_manual('Have a late cryptic',values=c('blue', 'red')) +
    ylab("Tissue") + 
    xlab("Frequency") +
    ggtitle("Samples with an early cryptic")

early_normed = s_cryptic_psi |> 
    filter(disease_tissue == TRUE) |> 
    left_join(selective_curve) |> 
    filter(Type == "Early") |> 
    pivot_wider(names_from = 'paste_into_igv_junction',
                values_from = 'psi',
                id_cols = 'sample') |> 
    tibble::column_to_rownames('sample') |> 
    scale() %>%
    as.data.frame()  %>%
    mutate(early_cryptic = rowSums(.)) |> 
    select(early_cryptic) |> 
    tibble::rownames_to_column('sample') |> 
    left_join(meta_data) 

late_normed = s_cryptic_psi |> 
    filter(disease_tissue == TRUE) |> 
    left_join(selective_curve) |> 
    filter(Type == "Late") |> 
    pivot_wider(names_from = 'paste_into_igv_junction',
                values_from = 'psi',
                id_cols = 'sample') |> 
    tibble::column_to_rownames('sample') |> 
    scale() %>%
    as.data.frame()  %>%
    mutate(late_cryptic = rowSums(.)) |> 
    select(late_cryptic) |> 
    tibble::rownames_to_column('sample') 

early_normed |> 
    left_join(late_normed) |>
    select(sample,early_cryptic,late_cryptic) |> 
    unique() |> 
    ggplot(aes(x = early_cryptic,
               y = late_cryptic)) +
    geom_point()


early_normed |> 
    left_join(late_normed) |> 
    filter(early_cryptic < 0 & late_cryptic > 5)
    ggplot(aes(x = early_cryptic,
               y = late_cryptic)) +
    geom_point()