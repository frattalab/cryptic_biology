###removed from supplementary


print("S1D")
master_table <- master_table %>%
    mutate(tdp43_sensitivity_be2 = ifelse(tdp43_sensitivity_be2 == "","Not found",tdp43_sensitivity_be2)) %>% 
    mutate(tdp43_sensitivity_sh = ifelse(tdp43_sensitivity_sh == "","Not found",tdp43_sensitivity_sh)) %>% 
    mutate(new_cate = case_when(tdp43_sensitivity_be2 == "Not found" ~ tdp43_sensitivity_sh,
                                tdp43_sensitivity_sh == "Not found" ~ tdp43_sensitivity_be2,
                                tdp43_sensitivity_sh == tdp43_sensitivity_be2 ~ tdp43_sensitivity_sh,
                                tdp43_sensitivity_be2 != 'Late' & tdp43_sensitivity_sh == 'Early' ~ 'Early',
                                tdp43_sensitivity_be2 != 'Early' & tdp43_sensitivity_sh == 'Late' ~ 'Late',
                                tdp43_sensitivity_sh == 'Intermediate'~  'Intermediate')) %>% 
    mutate(new_cate = ifelse(is.na(new_cate),"Ambiguous",new_cate)) 

print("S1C - UMAP clustering of events in SH and BE2 together")
master_table %>% 
    filter(tdp43_sensitivity_sh != "") %>% 
    ggplot(aes(x = UMAP_1_together, y = UMAP_2_together, color = as.character(new_cate))) +
    geom_point(show.legend = F) +
    xlab("UMAP 1 together") +
    ylab("UMAP 2 together") +
    labs(color = "cluster") +
    scale_color_manual(values = plot_sensitivity_colors) +
    theme_bw()


print("S1D - UMAP clustering of events in SH and in BE2 separated")

g <- master_table %>% 
    filter(tdp43_sensitivity_be2 != "") %>% 
    ggplot(aes(x = UMAP_1_sh, y = UMAP_2_sh, color = as.character(tdp43_sensitivity_sh))) +
    geom_point(show.legend = F) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    labs(color = "TDP43 sensitivity", title = "SH-SY5Y") +
    scale_color_manual(values = plot_sensitivity_colors) +
    theme_bw()
h <- master_table %>% 
    filter(tdp43_sensitivity_be2 != "") %>% 
    ggplot(aes(x = UMAP_1_be2, y = UMAP_2_be2, color = as.character(tdp43_sensitivity_be2))) +
    geom_point() +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    labs(color = "TDP43 sensitivity", title = "SK-N-BE(2)") +
    scale_color_manual(values = plot_sensitivity_colors) +
    theme_bw()
g+h


print("S1E - differences in differential splicing are not due to changes in gene expression")

print("S1F")
#apa_sh
print("S1G")
#apa_sk



gene_list_slope_upf1 <- c("STMN2", "UNC13A", "AARS1", "HDGFL2", "CYFIP2", "SYNE1")
input_list <- c("ctrl_ctrl",
                "ctrl_UPF1", 
                "TDP43_ctrl", 
                "TDP43_UPF1")
big_data_upf1 <- read.csv(here::here('data','nmd_or_not_upf1.csv'))
a <- slope_plot_nmd(big_data_upf1, 
                    gene_list_slope_upf1)

print("S3B")

