 # doi: 10.1038/s41467-022-35494-w
library(stringr)
conflicted::conflict_prefer('select','dplyr')
conflicted::conflicts_prefer(dplyr::filter)
hi = fread(here::here("data/HiSeq_RobustSubtypeAssignment_10-11-21.csv"))
nova = fread(here::here("data/NovaSeq_RobustSubtypeAssignment_10-11-21_Tiebreaker.csv"))
all_data = fread(here::here("data/AllSubjects_coldata_SUBTYPES.csv"))
rsem_tpm = fread(file.path(here::here(),"data/rsem_tpm_nygc.csv"))
# What correlates with the number of cryptic events observed? -------------
cryptic_psi = fread("data/everything_except_liu_ferguson_nygc_psi.aggregated.clean.annotated.bed")

# Define the regex

cryptic_psi[,sample := str_extract(V4, "^[^_]+")]
cryptic_psi[,psi := as.numeric(str_extract(V5, "^[^_]+"))]
cryptic_psi[,paste_into_igv_junction := paste0(V1,":",V2,"-",V3)]

s_cryptic_psi = cryptic_psi|> filter(paste_into_igv_junction %in% selective_junctions)
s_cryptic_psi = meta_data |> 
    left_join(s_cryptic_psi |> 
                  dplyr::select(paste_into_igv_junction,psi,sample)) |> 
    unique()  |> 
    dplyr::select(paste_into_igv_junction,psi,sample) |> 
    unique() |> 
    complete(paste_into_igv_junction, sample,fill = list(psi = 0)) |> unique()

s_cryptic_psi = s_cryptic_psi |> 
    left_join(meta_data) |> 
    as.data.table() |> 
    create_formated_metadata()

long_score = s_cryptic_psi |> 
    filter(disease_tissue == TRUE) |> 
    pivot_wider(names_from = 'paste_into_igv_junction',
                values_from = 'psi',
                id_cols = 'sample') |> 
    tibble::column_to_rownames('sample') |> 
    scale() %>%
    as.data.frame()  %>%
    mutate(total_cryptic = rowSums(.)) |> 
    select(total_cryptic) |> 
    tibble::rownames_to_column('sample') |> 
    left_join(meta_data)  |> 
    as.data.table() |> 
    create_formated_metadata()
    

# How does the meta-score vary across tissue and disease? ------------------
long_score |> 
    ggplot(aes(fill = disease,
               x = total_cryptic),alpha = 0.3) + 
    geom_density() + 
    facet_wrap(~tissue_clean)

long_score |> 
    ggplot(aes(color = disease,
               x = total_cryptic)) + 
    geom_hline(yintercept = 0.5) + 
    geom_vline(xintercept = 0) +
    stat_ecdf(linewidth = 1.2) + 
    facet_wrap(~tissue_clean) + 
    ggpubr::theme_pubr()  +
    xlab("Meta-score selective cryptic")

    



hi |> 
    filter(V1 %in% c("", "Majority")) |> 
    mutate(V1 = ifelse(V1 == "","sample",V1)) |> 
    tibble::column_to_rownames('V1') |> 
    t() |> 
    as.data.table() |> 
    mutate(sample = gsub("\\.","-",sample)) |> 
    left_join(cryptics_observed[disease_tissue == TRUE] |> 
                  left_join(meta_data_full,by = 'sample')) |> 
    as.data.table() |> 
    ggplot(aes(x = Majority,
               y = n_cryptics_observed)) + 
    geom_boxplot() + 
    ggsignif::geom_signif(comparisons = list(c("GLIA","OX")))


all_data |> 
    mutate(sample = gsub("\\.","-",Subject)) |> 
    left_join(cryptics_observed[disease_tissue == TRUE]) |> 
    as.data.table() |> 
    filter(!is.na(n_cryptics_observed) & Subtype != "Control")
    group_by(Subtype) |> 
    mutate(mean_cryptic = mean(n_cryptics_observed)) |> 
    ungroup() |> 
    mutate(Subtype = fct_reorder(Subtype,-mean_cryptic)) |> 
    ggplot(aes(x = Subtype,
               y = n_cryptics_observed)) + 
    geom_boxplot(outlier.shape = NULL) + 
    ggsignif::geom_signif(comparisons = list(c("GLIA","OX"),
                                             c("GLIA","TE"),
                                             c("OX","TE")),
                          map_signif_level = TRUE,
                          tip_length  = 0,
                          step_increase = 0.1) +
    ggpubr::theme_pubr() +
    ylab("N selective cryptic observed")


all_data = all_data |> 
    mutate(sample = gsub("\\.","-",Subject)) 

all_data |> 
    left_join(long_score,by = 'sample') |> 
    as.data.table() |> 
    filter(!is.na(total_cryptic) & Subtype != "Control") |> 
    group_by(Subtype) |> 
    mutate(mean_cryptic = mean(total_cryptic)) |> 
    ungroup() |> 
    mutate(Subtype = fct_reorder(Subtype,-mean_cryptic)) |> 
    ggplot(aes(x = Subtype,
               y = total_cryptic)) + 
    geom_boxplot(outlier.shape = NULL) + 
    ggsignif::geom_signif(comparisons = list(c("GLIA","OX"),
                                             c("GLIA","TE"),
                                             c("OX","TE")),
                          map_signif_level = TRUE,
                          tip_length  = 0,
                          step_increase = 0.1) +
    ggpubr::theme_pubr() +
    ylab("Meta-score selective cryptic")

g = annotables::grch38 |> filter(symbol == 'TRIM23') |> pull(ensgene)
# g = 'ENSG00000109654'
s = annotables::grch38 |> filter(grepl(g,ensgene)) |> pull(symbol) |> unique()

rsem_tpm |> 
    filter(grepl(g,gene)) |> 
    melt() |> 
    left_join(cryptics_observed[disease_tissue == TRUE],by = c('variable' = 'sample')) |> 
    filter(grepl("Cortex",tissue_clean)) |>
    mutate(tpm = as.numeric(value)) |> 
    mutate(scaled = scale(tpm)) |> 
    mutate(scaled_c = scale(n_cryptics_observed)) |> 
    ggplot(aes(y = n_cryptics_observed, x = tpm)) + 
    geom_point() + 
    stat_cor() +
    geom_rug() + 
    geom_smooth(method = 'lm') + 
    ggtitle(s) + 
    ggpubr::theme_pubr() + 
    ylab("N selective cryptic observed") + 
    xlab("gene TPM")

rsem_tpm |> 
    filter(grepl(g,gene)) |> 
    melt() |> 
    left_join(long_score[disease_tissue == TRUE],by = c('variable' = 'sample')) |> 
    filter(grepl("Cortex",tissue_clean)) |> 
    mutate(tpm = as.numeric(value)) |> 
    mutate(scaled = scale(tpm)) |> 
    ggplot(aes(y = total_cryptic, x = tpm,color = disease)) + 
    geom_point() + 
    stat_cor() +
    geom_rug() + 
    geom_smooth(method = 'lm') + 
    ggtitle(s) + 
    ggpubr::theme_pubr() + 
    ylab("Scaled selective psi meta-score") + 
    xlab("gene TPM")

only_one = cryptics_observed[disease_tissue == TRUE & n_cryptics_observed < 4,sample]
selective_counts[sample %in% only_one] |> 
    filter(spliced_reads >= 2) |> 
    as.data.table() |> 
    dplyr::count(gene) |> 
    arrange(-n) |> View()


rsem_flip_scale = rsem_tpm %>%  
    mutate(xmean = rowMeans(select(., starts_with("C")))) |> 
    filter(xmean > 5) |> 
    select(-xmean) |> tibble::column_to_rownames('gene') |> t() |> scale() |> 
    as.data.frame() |> 
    tibble::rownames_to_column('sample') |> 
    as.data.table()

rsem_flip_scale_in_tdp = rsem_flip_scale |> filter(sample %in% (selective_counts[tdp_path == "path" & disease_tissue == TRUE,sample] |> unique()))
# high_brain = rsem_flip_scale |> filter(gene == 'CGND-HRA-00969') |> melt() |> 
#     as.data.table() |> arrange(-value) |> 
#     mutate(ensgene = gsub("\\..*","",variable)) |> 
#     left_join(annotables::grch38) 
#     
#     
    
melt_flip = rsem_flip_scale_in_tdp[] |> 
    melt(id.cols = 'sample') |> 
    left_join(cryptics_observed,by = c('sample'))

melt_flip_nest = melt_flip[!is.na(n_cryptics_observed) & disease_tissue == TRUE] |> 
    select(sample,variable,value,n_cryptics_observed) |> unique()

tested = melt_flip_nest[,stats::cor.test(value,n_cryptics_observed),by = 'variable']
sig = tested[p.value < 0.01,.(estimate,variable,p.value)] |> unique() |> mutate(ensgene = gsub("\\..*","",variable)) |> left_join(annotables::grch38)
