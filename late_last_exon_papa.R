my_clean_reader <- function(file){
    df = data.table::as.data.table(janitor::clean_names(data.table::fread(file)))
    return(df)
}

fix_column_names = function(df){
    ori_colnames = df |>  colnames()
    ori_colnames[19] = "mean_ppau_baseline"
    ori_colnames[20] = "mean_ppau_contrast"
    ori_colnames[21] = "deltaPPAU"
    colnames(df) = ori_colnames
    
    return(df)
    
}
all_curves = list.files('/Users/annaleigh/cluster/sbs_projects/data/PAPA/runs/gencode_v40/quantify/i3_cortical_zanovello_novel',pattern = 'zanovello_shsy5y_curve_',full.names = TRUE)
samp_ids = gsub("/Users/annaleigh/cluster/sbs_projects/data/PAPA/runs/gencode_v40/quantify/i3_cortical_zanovello_novel/","",all_curves)

all_curves = paste0(all_curves,'/differential_apa/saturn_apa.results.processed.tsv')

last_curve_full = purrr::map(all_curves,my_clean_reader)
##add on the name as an additional column
last_curve_full = purrr::map2(last_curve_full, samp_ids, ~cbind(.x, experiment = .y))
last_curve_full = purrr::map(last_curve_full,fix_column_names)

last_curve_full = data.table::rbindlist(last_curve_full,fill=TRUE) |> 
    unique()

last_curve_full |>
    select(le_id,gene_name,event_type,annot_status,mean_ppau_baseline,mean_ppau_contrast,experiment) |> 
    unique() |> 
    # filter(event_type == "last_exon_extension") |> 
    unique() |> 
    group_by(le_id) |> 
    mutate(mean_ppau_baseline = mean(mean_ppau_baseline)) |> 
    ungroup() |> 
    pivot_wider(names_from = 'experiment',
                values_from = 'mean_ppau_contrast',
                values_fill = 0,
                id_cols = c("le_id",
                            "gene_name",
                            "event_type",
                            "annot_status",
                            "mean_ppau_baseline")) |> 
    unique() |>
    filter(mean_ppau_baseline < 0.15) |>
    filter(zanovello_shsy5y_curve_00125 < 0.1) |>
    filter(zanovello_shsy5y_curve_00187 < 0.1) |>
    filter(zanovello_shsy5y_curve_0075 - mean_ppau_baseline > 0.1) |>
    filter(zanovello_shsy5y_curve_0025 - mean_ppau_baseline > 0.1) |>
    filter(zanovello_shsy5y_curve_00125 < zanovello_shsy5y_curve_0075) |>
    # pull(le_id) |> unique() |> clipr::write_clip()
    melt(id.vars = c("le_id",
                     "gene_name",
                     "event_type",
                     "annot_status",
                     "mean_ppau_baseline") ) |> 
    ggplot(aes(x = variable,
               y = value,
               group = le_id)) + 
    geom_point() + 
    geom_line() + 
    geom_text_repel(aes(label = gene_name))


possible_late_last_exons = c("ENSG00000132849.22_1",
  "ENSG00000150760.13_2",
  "ENSG00000166046.11_3",
  "ENSG00000174989.13_2",
  "ENSG00000184923.12_2",
  "ENSG00000185158.12_3",
  "ENSG00000244405.8_3")


last_curve_full |> 
    filter(le_id %in% possible_late_last_exons) |> 
    select(le_id,chromosome,start,end,strand) |> 
    unique() |> 
    separate_rows(start,end) |> 
    unique() |> 
    mutate(score = 0) |> 
    makeGRangesFromDataFrame(,keep.extra.columns = TRUE) |> 
    rtracklayer::export("possible_late_le.bed")


last_curve_full |> 
    filter(le_id %in% possible_late_last_exons) |> 
    select(le_id,chromosome,start,end,strand,event_type,gene_name)
