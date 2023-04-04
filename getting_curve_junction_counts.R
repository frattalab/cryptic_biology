library(tidyverse)
library(data.table)
library(tidytable)
df_cat_sushi = readxl::read_excel(here::here('data/df_cat_shsy5y_curve.xlsx')) |> as.data.table() |> rowid_to_column()
#used the parsing tool to get junction counts - let's read that in 
junction_counts = fread(here::here('data/validation_check_one.aggregated.bed'))
#let's just check the correct samples
junction_counts = junction_counts |> 
    filter(grepl("dox-conc",V4)) |> 
    mutate(sample = gsub(".SJ.out","",V4))
#let's check that all the df_cat_sushi junctions were picked up by the parsing tool
match_me_junctions = df_cat_sushi |> 
    mutate(m = glue::glue("{chr}:{start}-{end}:{strand}")) |> 
    pull(m) |> 
    unique()

found_junctions = junction_counts |> 
    mutate(m = glue::glue("{V1}:{V2}-{V3}:{V6}")) |> 
    pull(m) |> 
    unique()

table(match_me_junctions %in% found_junctions)

#some 33 junctions weren't found in the STAR output files
not_found_by_parser = df_cat_sushi |> 
    mutate(m = glue::glue("{chr}:{start}-{end}:{strand}")) |> 
    filter(!(m %in% found_junctions)) |> 
    tibble::column_to_rownames('rowid')

#these appear to be one offs
#chr15 74180243 74180782 - parsed 
#chr15 74180244 74180781 - in table
#chr10  3120044  3125097 - parsed
#chr10  3120045  3125096 - in table

not_found_by_parser = not_found_by_parser |> 
    mutate(start = start - 1) |> 
    mutate(end = end + 1)

lost_id_map = not_found_by_parser |> 
    select(m) |> 
    tibble::rownames_to_column('rowid')

pattern <- "DOX_[0-9.]+|NT_[0-9]+"

not_found_counts = not_found_by_parser |> 
    left_join(junction_counts,
              by = c("chr" = "V1",
                     "start" = "V2",
                     "end" = "V3")) |> 
    select(gene_id,m,sample,V5) |> 
    complete(m, sample, fill = list(V5 = 0)) |> 
    group_by(m) %>%
    fill(gene_id, .direction = "downup") %>%
    ungroup() |> 
    unique() |> 
    mutate(condition =  str_extract(sample, pattern)) |> 
    filter(!is.na(condition)) |> 
    group_by(condition,m) |> 
    summarise(mean_count = mean(V5),
              median_count = median(V5),
              total_count = sum(V5),
              gene_id) |> 
    ungroup() |> 
    unique() |> 
    pivot_wider(names_from = 'condition',
                values_from = c('mean_count', 'median_count', 'total_count')) |> 
    janitor::clean_names() |> 
    mutate(found_in_star_outputs = FALSE) |> 
    left_join(lost_id_map)

#do the same analysis for those found from star's parser
found_id_map = df_cat_sushi |> 
    mutate(m = glue::glue("{chr}:{start}-{end}:{strand}")) |> 
    select(rowid,m)

found_counts = df_cat_sushi |> 
    mutate(m = glue::glue("{chr}:{start}-{end}:{strand}")) |> 
    filter((m %in% found_junctions)) |> 
    left_join(junction_counts,
              by = c("chr" = "V1",
                     "start" = "V2",
                     "end" = "V3")) |> 
    select(gene_id,m,sample,V5) |> 
    complete(m, sample, fill = list(V5 = 0)) |> 
    group_by(m) %>%
    fill(gene_id, .direction = "downup") %>%
    ungroup() |> 
    unique() |> 
    mutate(condition =  str_extract(sample, pattern)) |> 
    filter(!is.na(condition)) |> 
    group_by(condition,m) |> 
    summarise(mean_count = mean(V5),
              median_count = median(V5),
              total_count = sum(V5),
              gene_id) |> 
    ungroup() |> 
    unique() |> 
    pivot_wider(names_from = 'condition',
                values_from = c('mean_count', 'median_count', 'total_count')) |> 
    janitor::clean_names() |> 
    mutate(found_in_star_outputs = TRUE) |> 
    left_join(found_id_map)



junction_raw_counts = rbind(found_counts,not_found_counts) |> mutate(rowid = as.numeric(rowid))
fwrite(junction_raw_counts,'data/junction_raw_counts_sushi.csv')
