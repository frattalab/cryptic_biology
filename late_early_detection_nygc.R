library(readxl)
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
