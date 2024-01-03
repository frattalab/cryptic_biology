library(tidyverse)
library(data.table)
# Clustering
library(cluster) 
library(factoextra)
set.seed(420)

# read in the data cleaned by Matteo ----------------------------------------

manual_validation_sy5y_distinct = fread('data/manual_validation_sy5y_distinct.csv')
manual_validation_dz_distinct = fread('data/manual_validation_dz_distinct.csv')
# Cluster SH-SY5Y ----------------------------------------

sushi_curve_cryptic_matrix = manual_validation_sy5y_distinct |> 
    janitor::clean_names() |> 
    select(gene_id,paste_into_igv_junction,dox0,dox0125:dox075) |> 
    mutate(name = glue::glue("{gene_id}_{paste_into_igv_junction}")) |> 
    tibble::column_to_rownames('name') |> 
    select(-c(gene_id,paste_into_igv_junction)) |> 
    as.matrix()

set.seed(420)
p = pheatmap::pheatmap(sushi_curve_cryptic_matrix, cluster_cols = FALSE,kmeans_k = 3)

sushi_curve_cryptic_matrix |> 
    as.data.frame() |> 
    tibble::rownames_to_column('name') |> 
    mutate(cluster_assignments = p$kmeans$cluster) |> 
    melt(id.vars = c('name','cluster_assignments')) |> 
    separate(name, into = c("gene",'paste_into_igv_junction'),sep = '_',remove = FALSE) |> 
    ggplot(aes(x = variable,
               y = value,
               group = name)) + 
    geom_line() + 
    facet_wrap(~cluster_assignments)


dz_curve_cryptic_matrix = manual_validation_dz_distinct |> 
    janitor::clean_names() |> 
    select(gene_id,paste_into_igv_junction,dox0,dz002_vs_nt:dz1_vs_nt) |> 
    mutate(name = glue::glue("{gene_id}_{paste_into_igv_junction}")) |> 
    tibble::column_to_rownames('name') |> 
    select(-c(gene_id,paste_into_igv_junction)) |> 
    as.matrix()

set.seed(420)
p2 = pheatmap::pheatmap(dz_curve_cryptic_matrix, cluster_cols = FALSE,kmeans_k = 3)




dz_curve_cryptic_matrix |> 
    as.data.frame() |> 
    tibble::rownames_to_column('name') |> 
    mutate(cluster_assignments = p2$kmeans$cluster) |> 
    left_join()
    melt(id.vars = c('name','cluster_assignments')) |> 
    separate(name, into = c("gene",'paste_into_igv_junction'),sep = '_',remove = FALSE) |> 
    ggplot(aes(x = variable,
               y = value,
               group = name)) + 
    geom_line() + 
    facet_wrap(~cluster_assignments)



# Comparing cluster SK & SH -----------------------------------------------


dz_compare = dz_curve_cryptic_matrix |> 
    as.data.frame() |> 
    tibble::rownames_to_column('name') |> 
    mutate(cluster_assignments = p2$kmeans$cluster) |> 
    select(name,cluster_assignments) |> 
    dplyr::rename(cluster_dz = cluster_assignments)
    
sushi_compare = sushi_curve_cryptic_matrix |> 
    as.data.frame() |> 
    tibble::rownames_to_column('name') |> 
    mutate(cluster_assignments = p$kmeans$cluster) |> 
    select(name,cluster_assignments) |> 
    dplyr::rename(cluster_sushi = cluster_assignments)


df = full_join(dz_compare,sushi_compare) |> 
    mutate(cluster_dz = ifelse(is.na(cluster_dz),"not found",cluster_dz)) |>
    mutate(cluster_sushi = ifelse(is.na(cluster_sushi),"not found",cluster_sushi)) |>
    filter(!is.na(cluster_dz)) |> 
    filter(!is.na(cluster_sushi)) |> 
    mutate(cluster_dz = case_when(cluster_dz == 1 ~ 'early',
                                  cluster_dz == 2 ~ 'intermediate',
                                  cluster_dz == 3 ~ 'late',
                                  TRUE ~ cluster_dz)) |> 
    mutate(cluster_sushi = case_when(cluster_sushi == 1 ~ 'early',
                                     cluster_sushi == 2 ~ 'intermediate',
                                     cluster_sushi == 3 ~ 'late',
                                     TRUE ~ cluster_sushi)) |> 
    dplyr::rename(`SH-SY-5Y` = cluster_sushi) |> 
    dplyr::rename(`SK-N-BE(2)` = cluster_dz) |> 
    ggsankey::make_long(`SK-N-BE(2)`,`SH-SY-5Y`)



# Chart 1
pl <- ggplot(df, aes(x = x
                     , next_x = next_x
                     , node = node
                     , next_node = next_node
                     , fill = factor(node)
                     , label = node)
)
pl <- pl +geom_sankey(
                      , node.color = "black"
                      ,show.legend = FALSE)
pl <- pl +geom_sankey_label(size = 3, color = "black", fill= "white", hjust = -0.5)
pl <- pl +  theme_bw()
pl <- pl + theme(legend.position = "none")
pl <- pl +  theme(axis.title = element_blank()
                  , axis.text.y = element_blank()
                  , axis.ticks = element_blank()  
                  , panel.grid = element_blank())
pl <- pl + scale_fill_viridis_d(option = "inferno")

pl <- pl + labs(fill = 'Nodes')
pl


# Sankey with manual ------------------------------------------------------


sank_dz = manual_validation_dz_distinct |> 
    select(gene_id,paste_into_igv_junction,Type) |> 
    janitor::clean_names() |> 
    mutate(name = glue::glue("{gene_id}_{paste_into_igv_junction}")) |> 
    select(-c(gene_id,paste_into_igv_junction)) |> 
    dplyr::rename(`SK-N-BE(2)` = type)  

sank_sushi = manual_validation_sy5y_distinct |> 
    select(gene_id,paste_into_igv_junction,Type) |> 
    janitor::clean_names() |> 
    mutate(name = glue::glue("{gene_id}_{paste_into_igv_junction}")) |> 
    select(-c(gene_id,paste_into_igv_junction)) |> 
    dplyr::rename(`SH-SY-5Y` = type) 
    

df3 = full_join(sank_dz,sank_sushi) |> 
    mutate(`SK-N-BE(2)` = ifelse(is.na(`SK-N-BE(2)`),"not found",`SK-N-BE(2)`)) |>
    mutate(`SH-SY-5Y` = ifelse(is.na(`SH-SY-5Y`),"not found",`SH-SY-5Y`)) |>
    ggsankey::make_long(`SK-N-BE(2)`,`SH-SY-5Y`)



# Chart 1
pl2 <- ggplot(df3, aes(x = x
                     , next_x = next_x
                     , node = node
                     , next_node = next_node
                     , fill = factor(node)
                     , label = node)
)
pl2 <- pl2 +geom_sankey(
    , node.color = "black"
    ,show.legend = FALSE)
pl2 <- pl2 +geom_sankey_label(size = 3, color = "black", fill= "white", hjust = -0.5)
pl2 <- pl2 +  theme_bw()
pl2 <- pl2 + theme(legend.position = "none")
pl2 <- pl2 +  theme(axis.title = element_blank()
                  , axis.text.y = element_blank()
                  , axis.ticks = element_blank()  
                  , panel.grid = element_blank())

pl2 <- pl2 + labs(fill = 'Nodes')
pl2
