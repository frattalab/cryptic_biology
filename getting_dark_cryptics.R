library(ggplot2)
splicing_full = fread("data/splicing_full_delta_psi_tables.csv")

merged_nmd_ce = fread('https://raw.githubusercontent.com/frattalab/cryptic_biology/f7f01522ef8d3285ac0b952d5aee6affaabbba75/data/merged_nmd_CE.csv?token=GHSAT0AAAAAACAKYNAXRYLQSFZY7NGDVHAEZBGVUAQ')
merged_nmd_ce = merged_nmd_ce |> 
    mutate(gene_id.x =  gsub("\\..*","",gene_id.x),
           gene_id.y =  gsub("\\..*","",gene_id.y))
# 
# potential_new_selective |> 
#     mutate(obs_tot = path + not_path) |> 
#     mutate(gene_name = glue::glue("{gene}-{type}")) |> 
#     left_join(merged_nmd_ce,by = c("paste_into_igv_junction")) |> 
#     filter(!is.na(color_gene_name.x) | !is.na(color_gene_name.y)) |> 
#     select(paste_into_igv_junction,not_path,path,gene_name,color_gene_name.x,color_gene_name.y) |> 
#     dplyr::rename(nmd_call_chx = color_gene_name.x,
#                   nmd_call_upf1 = color_gene_name.y) |> 
#     unique() |> 
#     ggplot(aes(x = nmd_call_upf1,
#                y = path)) + 
#     geom_boxplot() + 
#     scale_y_continuous(trans = scales::pseudo_log_trans())



affected_by_chx = merged_nmd_ce |> 
    select(gene_name.x,paste_into_igv_junction,Control_Control:Cycloheximide_TDP43KD,color_gene_name.x) |> 
    filter(Control_Control < 0.05 & 
               color_gene_name.x == "Delta PSI > 0.05") |> 
    pull(paste_into_igv_junction)

dark_chx = merged_nmd_ce |> 
    select(gene_name.x,paste_into_igv_junction,Control_Control:Cycloheximide_TDP43KD,color_gene_name.x) |> 
    filter(Control_Control < 0.05 & 
               Control_TDP43KD < 0.1 & 
               Cycloheximide_TDP43KD > 0.1 & 
           color_gene_name.x == "Delta PSI > 0.05") |> pull(paste_into_igv_junction) |> unique()

dark_upf1 =  merged_nmd_ce |> 
    select(gene_name.y,paste_into_igv_junction,ctrl_ctrl:TDP43_UPF1,color_gene_name.y) |> 
    filter(ctrl_ctrl < 0.05 & 
               TDP43_ctrl < 0.1 & 
               TDP43_UPF1 > 0.1 & 
               color_gene_name.y == "Delta PSI > 0.05") |> pull(paste_into_igv_junction) |> unique()

everything_but_chx = splicing_full |> 
    filter(!(comparison %in% c("controlfergusonhela-tdp43kdfergusonhela",
                               "controlliufacsneurons-tdp43kdliufacsneurons",
                               "cycloheximidecontrol-cycloheximidetdp43kd", 
                               "controltdp43kd-cycloheximidetdp43kd"))) |> 
    select(comparison,baseline_PSI,contrast_PSI,junc_cat,paste_into_igv_junction) |> 
    unique() |> 
    group_by(paste_into_igv_junction) |>  
    mutate(quantified_n_datasets = n_distinct(comparison)) |> 
    ungroup()

everything_but_upf1 = splicing_full |> 
    filter(!(comparison %in% c("controlfergusonhela-tdp43kdfergusonhela",
                               "controlliufacsneurons-tdp43kdliufacsneurons",
                               "ctrlupf1-tdp43upf1" ,
                               "tdp43ctrl-tdp43upf1"))) |> 
    select(comparison,baseline_PSI,contrast_PSI,junc_cat,paste_into_igv_junction) |> 
    unique() |> 
    group_by(paste_into_igv_junction) |>  
    mutate(quantified_n_datasets = n_distinct(comparison)) |> 
    ungroup()

c2_chx = everything_but_chx |> 
    filter(baseline_PSI < 0.05 & contrast_PSI > 0.1) |> 
    group_by(paste_into_igv_junction) |>  
    mutate(cryptic_n_datasets = n_distinct(comparison)) |> 
    ungroup() |> 
    select(paste_into_igv_junction,cryptic_n_datasets) |> 
    unique()


c2_upf1 = everything_but_upf1 |> 
    filter(baseline_PSI < 0.05 & contrast_PSI > 0.1) |> 
    group_by(paste_into_igv_junction) |>  
    mutate(cryptic_n_datasets = n_distinct(comparison)) |> 
    ungroup() |> 
    select(paste_into_igv_junction,cryptic_n_datasets) |> 
    unique()

everything_but_chx = everything_but_chx |> 
    left_join(c2_chx) |> 
    mutate(cryptic_n_datasets = ifelse(is.na(cryptic_n_datasets),0,cryptic_n_datasets)) |> 
    unique()

merged_nmd_ce |> 
    select(gene_name.x,paste_into_igv_junction,Control_Control:Cycloheximide_TDP43KD,color_gene_name.x) |> 
    filter(Control_Control < 0.05 & 
               Control_TDP43KD < 0.1 & 
               Cycloheximide_TDP43KD > 0.1 & 
               color_gene_name.x == "Delta PSI > 0.05") |> 
    left_join(splicing_full2) |> 
    filter(baseline_PSI < 0.05) |> unique() |> View()


everything_but_chx |>
    filter(paste_into_igv_junction %in% dark_chx) |> 
    filter(cryptic_n_datasets == 0) |> 
    filter(baseline_PSI < 0.05) |> 
    left_join(expression_by_pathology) |> 
    filter(gene == "CYFIP2")
    filter(path >= 10 & not_path <=5) |> 
    View()

chx_rna_full = fread("https://docs.google.com/spreadsheets/d/e/2PACX-1vRKe8lDDucaHQWNwbdHBeY3P6FeHa8Q14mZ3o24dem37vZLG0DfubVnTkT5l5jAjy8BvHUuCYOWsDY8/pub?output=csv")    




# Are dark cryptic genes downregulated in normal KD? ----------------------
dark_chx_genes = merged_nmd_ce |> 
    filter(paste_into_igv_junction %in% affected_by_chx) |> 
    pull(gene_id.x) |> unique() 

dark_chx_genes = gsub("\\..*","",dark_chx_genes)
dark_chx_genes = gsub("\\..*","",dark_chx_genes)

chx_rna_full |> 
    filter(experiment %in% c("Control_TDP43KD|Control_Control",
                             "Cycloheximide_TDP43KD|Control_TDP43KD")) |> 
    mutate(dark_cryptic_containing = ensgene %in% dark_chx_genes) |> 
    filter(dark_cryptic_containing) |> 
    ggplot(aes(x = experiment,
               y = log2fold_change,
               group =ensgene)) + 
    geom_point() +
    geom_line(alpha = 0.3) + 
    geom_hline(yintercept = 0)


chx_dark_rna = chx_rna_full |> 
    filter(padj < 0.1) |> 
    filter(experiment %in% c("Control_TDP43KD|Control_Control",
                             "Cycloheximide_TDP43KD|Control_TDP43KD",
                             "Cycloheximide_Control|Control_Control")) |> 
    mutate(dark_cryptic_containing = ensgene %in% dark_chx_genes) |> 
    filter(dark_cryptic_containing) |> 
    pivot_wider(names_from = 'experiment',
                values_from = c("log2fold_change"),
                id_cols = c('gene_name',"ensgene")) |> 
    drop_na() |> 
    janitor::clean_names() |> 
    mutate(gene_movement = case_when(control_tdp43kd_control_control < 0 & cycloheximide_tdp43kd_control_tdp43kd > 0 ~ "down_rescued_level",
                                     cycloheximide_tdp43kd_control_tdp43kd < 0 & control_tdp43kd_control_control > 0 ~ "direction_swap",
                                     TRUE ~ "other")) |> 
    melt(id.vars = c("gene_name","ensgene","gene_movement"))
    

chx_dark_rna |> 
    mutate(variable = fct_relevel(variable, "cycloheximide_control_control_control")) |> 
    ggplot(aes(x = variable,
               y = value,
               group =ensgene)) + 
    geom_point() +
    geom_line(alpha = 0.3) + 
    geom_hline(yintercept = 0) + 
    facet_wrap(~gene_movement) +
    scale_x_discrete(labels=c("cycloheximide_tdp43kd_control_tdp43kd" = "CHX-KD vs KD", 
                              "control_tdp43kd_control_control" = "KD vs Control",
                              'cycloheximide_control_control_control' = "CHX ctrl")) + 
    ylab("Log2Fold Change") + 
    ggtitle("Genes that contain CHX rescued splicing (Delta PSI > 0.05)") + 
    ggpubr::theme_pubr() +
    xlab(element_blank())
    
    
chx_dark_rna |> 
    filter(gene_movement == "swapped_but_why")
