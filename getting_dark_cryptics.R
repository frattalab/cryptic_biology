library(ggplot2)
merged_nmd_ce = fread('/Users/annaleigh/Downloads/merged_nmd_CE (1).csv')

    

merged_nmd_ce = fread('/Users/annaleigh/Downloads/merged_nmd_CE (1).csv')
merged_nmd_ce |> 
    mutate(paste_into_igv_junction)

potential_new_selective |> 
    mutate(obs_tot = path + not_path) |> 
    mutate(gene_name = glue::glue("{gene}-{type}")) |> 
    left_join(merged_nmd_ce,by = c("paste_into_igv_junction")) |> 
    filter(!is.na(color_gene_name.x) | !is.na(color_gene_name.y)) |> 
    select(paste_into_igv_junction,not_path,path,gene_name,color_gene_name.x,color_gene_name.y) |> 
    dplyr::rename(nmd_call_chx = color_gene_name.x,
                  nmd_call_upf1 = color_gene_name.y) |> 
    unique() |> 
    ggplot(aes(x = nmd_call_upf1,
               y = path)) + 
    geom_boxplot() + 
    scale_y_continuous(trans = scales::pseudo_log_trans())




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

c2 = splicing_full2 |> 
    filter(baseline_PSI < 0.05 & contrast_PSI > 0.1) |> 
    group_by(paste_into_igv_junction) |>  
    mutate(cryptic_n_datasets = n_distinct(comparison)) |> 
    ungroup() |> 
    select(paste_into_igv_junction,cryptic_n_datasets) |> 
    unique()

splicing_full2 = splicing_full2 |> 
    left_join(c2) |> 
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


splicing_full2 |> 
    filter(paste_into_igv_junction %in% dark_chx) |> 
    mutate(baseline_PSI_any =  baseline_PSI > 0.05) |> unique() |> 
    group_by(paste_into_igv_junction) |> 
    mutate(baseline_PSI_any =  sum(baseline_PSI_any)) |> 
    ungroup() |> 
    filter(baseline_PSI_any == 0) |> 
    View()
