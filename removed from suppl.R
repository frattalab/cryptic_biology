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




### targeted panel
#```{r targeted_panel, echo=FALSE, warning=FALSE, message=FALSE}

## Figure 4 - compendium / in vivo
print("4A")
##load targeted seq table and compute CE ratio / rename samples in agreement to UNC13A paper
maxc <- read.csv(here::here('data',"count_table_of_aligned_remove_200bp_plus_product_reads.csv")) %>%
    mutate(fraction_cryptic = n_cryptic / n_counts_for_gene) %>%
    mutate(sample_name_ari = case_when(sample_name == "P64_11" ~ "HC1",
                                       sample_name == "P17_07" ~ "HC2",
                                       sample_name == "P47_11" ~ "HC3",
                                       sample_name == "P35_07" ~ "HC4",
                                       sample_name == "P56_13" ~ "FTD1",
                                       sample_name == "P45_15" ~ "FTD2",
                                       sample_name == "P63_05" ~ "FTD3",
                                       sample_name == "P28_07" ~ "FTD4",
                                       sample_name == "P40_04" ~ "FTD5",
                                       sample_name == "P07_15" ~ "FTD6",
                                       sample_name == "P16_09" ~ "FTD7",
                                       sample_name == "P86_08" ~ "FTD8",
                                       sample_name == "P13_13" ~ "FTD9",
                                       sample_name == "P11_07" ~ "FTE10",
                                       sample_name == "ddH2O" ~ "HO2O")) #%>%
#pivot_wider(names_from = sample_name_ari, values_from = ratio_cryptic)
table(maxc$sample_name_ari)

##use only TDP43path-specific genes (use list from AL) - this is exluding more than 2/3 of genes!!
al <- read.csv(here::here('data',"expression_by_pathology_updated.csv")) %>%
    dplyr::filter(fraction_path > 0.01 & fraction_not_path < 0.005) %>% ####need to carefully check this with AL
    distinct(gene, .keep_all = T) ##check when binding with annotated validation table - should be ok since left_join
names(al)[4] <- "gene_name"
maxc_al <- maxc %>%
    left_join(al, by = "gene_name") %>%
    dplyr::filter(!is.na(paste_into_igv_junction))
table(maxc_al$sample_name_ari)
table(maxc_al$gene_name)

##check igv files for these genes
maxc_type <- maxc_al %>%
    select(-type) %>%
    mutate(gene = gene_name) %>%
    filter(gene_name != "SETD5") %>%
    left_join(master_table, by = "paste_into_igv_junction") %>%
    filter(!is.na(tdp43_sensitivity_sh)) %>% 
    mutate(group = ifelse(grepl("FT", sample_name_ari), "FTD",
                          ifelse(grepl("HC", sample_name_ari), "HC",
                                 "H2O"))) %>% 
    group_by(gene_name) %>%
    filter(sum(n_cryptic) > 0) %>%
    ungroup() %>% 
    filter(group != "H2O")

maxc_type %>%
    filter(tdp43_sensitivity_sh != "") %>%
    filter(gene_name != "SETD5") %>%
    #filter(n_cryptic>0) %>% 
    #filter(detected == "Yes") %>%
    filter(group != "H2O") %>%
    ggplot(aes(x = tdp43_sensitivity_sh, y = n_cryptic)) +
    geom_col(aes(fill = tdp43_sensitivity_sh)) +
    facet_wrap(facets = vars(group), ncol = 2) +
    scale_fill_brewer(palette = "Set1") +
    xlab("") +
    ylab("Count of cryptic exons detected by targeted panel") +
    labs(fill = "Category") +
    theme_bw()
#ggsave("~/Desktop/grouped_panel.png", width = 8, height = 4)

maxc_earli <- maxc_type %>%
    filter(tdp43_sensitivity_sh == "Earli") %>%
    filter(gene_name != "SETD5") %>% 
    ggplot(aes(x = sample_name_ari, y = fraction_cryptic)) +
    geom_bar(aes(fill = gene_name), stat = "identity", position = "stack", color="black") + #position = position_dodge(width = 0.8)) +
    #facet_grid(facets = vars(gene_name)) +
    scale_fill_brewer(palette = "Reds") +
    scale_y_continuous(limits = c(0,0.3)) +
    labs(x = "", y = "", fill = "Earli") +
    #scale_fill_manual(values = randomcoloR::distinctColorPalette(k = 30)) +
    theme_bw()
maxc_early <- maxc_type %>%
    filter(tdp43_sensitivity_sh == "Early") %>%
    filter(gene_name != "SETD5") %>% 
    ggplot(aes(x = sample_name_ari, y = fraction_cryptic)) +
    geom_bar(aes(fill = gene_name), stat = "identity", position = "stack", color="black") +#position = position_dodge(width = 0.8)) +
    #facet_grid(facets = vars(gene_name)) +
    scale_fill_brewer(palette = "Reds") +
    scale_y_continuous(limits = c(0,0.3)) +
    labs(x = "", y = "", fill = "Early") +
    #scale_fill_manual(values = randomcoloR::distinctColorPalette(k = 30)) +
    theme_bw()
maxc_intermediate <- maxc_type %>%
    filter(tdp43_sensitivity_sh == "Intermediate") %>%
    filter(gene_name != "SETD5") %>%
    ggplot(aes(x = sample_name_ari, y = fraction_cryptic)) +
    geom_bar(aes(fill = gene_name), stat = "identity", position = "stack", color="black") +#position = position_dodge(width = 0.8)) +
    #facet_grid(facets = vars(gene_name)) +
    scale_fill_brewer(palette = "Greens") +
    scale_y_continuous(limits = c(0,0.3)) +
    labs(x = "", y = "", fill = "Intermediate") +
    #scale_fill_manual(values = randomcoloR::distinctColorPalette(k = 30)) +
    theme_bw()
maxc_late <- maxc_type %>%
    filter(tdp43_sensitivity_sh == "Late") %>%
    filter(gene_name != "SETD5") %>%
    ggplot(aes(x = sample_name_ari, y = fraction_cryptic)) +
    geom_bar(aes(fill = gene_name), stat = "identity", position = "stack", color="black") +#position = position_dodge(width = 0.8)) +
    #facet_grid(facets = vars(gene_name)) +
    scale_fill_brewer(palette = "Blues") +
    scale_y_continuous(limits = c(0,0.3)) +
    labs(x = "", y = "", fill = "Late") +
    #scale_fill_manual(values = randomcoloR::distinctColorPalette(k = 30)) +
    theme_bw()
maxc_early / maxc_intermediate / maxc_late #/ maxc_na


### Supplementary figure 4 - in vivo
#```{r, echo=FALSE, warning=FALSE, message=FALSE}
print("S4A")

print("S4B")

print("S4C")

maxc <- read.csv(here::here('data',"count_table_of_aligned_remove_200bp_plus_product_reads.csv")) %>%
    mutate(fraction_cryptic = n_cryptic / n_counts_for_gene) %>%
    mutate(sample_name_ari = case_when(sample_name == "P64_11" ~ "HC1",
                                       sample_name == "P17_07" ~ "HC2",
                                       sample_name == "P47_11" ~ "HC3",
                                       sample_name == "P35_07" ~ "HC4",
                                       sample_name == "P56_13" ~ "FTD1",
                                       sample_name == "P45_15" ~ "FTD2",
                                       sample_name == "P63_05" ~ "FTD3",
                                       sample_name == "P28_07" ~ "FTD4",
                                       sample_name == "P40_04" ~ "FTD5",
                                       sample_name == "P07_15" ~ "FTD6",
                                       sample_name == "P16_09" ~ "FTD7",
                                       sample_name == "P86_08" ~ "FTD8",
                                       sample_name == "P13_13" ~ "FTD9",
                                       sample_name == "P11_07" ~ "FTE10",
                                       sample_name == "ddH2O" ~ "HO2O")) #%>%
#pivot_wider(names_from = sample_name_ari, values_from = ratio_cryptic)
table(maxc$sample_name_ari)

##use only TDP43path-specific genes (use list from AL) - this is exluding more than 2/3 of genes!!
al <- read.csv(here::here('data',"expression_by_pathology_updated.csv")) %>%
    dplyr::filter(fraction_path > 0.01 & fraction_not_path < 0.005) %>% ####need to carefully check this with AL
    distinct(gene, .keep_all = T) ##check when binding with annotated validation table - should be ok since left_join
names(al)[4] <- "gene_name"
maxc_al <- maxc %>%
    left_join(al, by = "gene_name") %>%
    dplyr::filter(!is.na(paste_into_igv_junction))
table(maxc_al$sample_name_ari)
table(maxc_al$gene_name)

##check igv files for these genes
maxc_type <- maxc_al %>%
    select(-type) %>%
    mutate(gene = gene_name) %>%
    filter(gene_name != "SETD5") %>%
    left_join(master_table, by = "paste_into_igv_junction") %>%
    filter(!is.na(tdp43_sensitivity_sh)) %>% 
    mutate(group = ifelse(grepl("FT", sample_name_ari), "FTD",
                          ifelse(grepl("HC", sample_name_ari), "HC",
                                 "H2O"))) %>% 
    group_by(gene_name) %>%
    filter(sum(n_cryptic) > 0) %>%
    ungroup() %>% 
    filter(group != "H2O")


maxc_type %>%
    filter(group != "H2O") %>%
    filter(tdp43_sensitivity_sh != "") %>% 
    ggplot(aes(x = sample_name_ari, y = gene_name, fill = fraction_cryptic)) +
    geom_tile() +
    #scale_fill_viridis_d() +
    scale_fill_gradient(low = "white", high = "black") +
    facet_grid(rows = vars(tdp43_sensitivity_sh), cols = vars(group), 
               scales = "free", space = "free") +
    ylab("") + xlab("Sample name") + labs(fill = "CE fraction") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90))










chx_g_e = chx_slope_plot_table$filtered_table %>% 
    ungroup() %>%  
    distinct(color_gene_name,gene_name) %>% 
    filter(gene_name %in% master_table$gene) %>% group_by(gene_name) %>%
    # Create a new column indicating if any row in the group has "NMD rescued"
    mutate(nmd_status = if_else(any(color_gene_name == "NMD rescued"), 
                                "NMD rescued", 
                                "non-NMD rescued")) %>%
    # Remove the grouping
    ungroup() %>% 
    select(-color_gene_name) %>% 
    unique()
smgi_g_e = pseudo_bulk_plot_slope$filtered_table %>% 
    ungroup() %>%  
    distinct(color_gene_name,gene_name) %>% 
    filter(gene_name %in% master_table$gene) %>% group_by(gene_name) %>%
    # Create a new column indicating if any row in the group has "NMD rescued"
    mutate(nmd_status = if_else(any(color_gene_name == "NMD rescued"), 
                                "NMD rescued", 
                                "non-NMD rescued")) %>%
    # Remove the grouping
    ungroup() %>% 
    select(-color_gene_name) %>% 
    unique()
chx_g_e %>% 
    left_join(smgi_g_e,by = c('gene_name')) %>% 
    filter(nmd_status.x != nmd_status.y)
####read in data!!
big_delta_de <- read.csv(here::here('data','big_delta_de.csv'))

for_plot <- big_delta_de[,c(4,8,10,11)] %>%
    # dplyr::filter(symbol == "STMN2" | symbol == "AARS1" | symbol == "HDGFL2" | symbol == "UNC13A" |
    #                 symbol == "SYNE1" | symbol == "CYFIP2")
    filter(symbol %in% master_table$gene)
for_plot_0 <- for_plot %>% 
    distinct(symbol, .keep_all = T) %>% 
    mutate(log2FoldChange = 0, padj = 0, source = 0)
for_plot <- rbind(for_plot_0,for_plot) %>% 
    mutate(symbol = factor(symbol, levels = c("AARS1", "HDGFL2", "SYNE1", "CYFIP2", "UNC13A", "STMN2")))

big_delta_de_tdp <- big_delta_de[,c(4,8,10,11)] %>% 
    dplyr::filter(symbol == "TARDBP")
for_plot_0_tdp <- big_delta_de_tdp[1,] %>% 
    mutate(log2FoldChange = 0, padj = 0, source = 0)
big_delta_de_tdp <- rbind(for_plot_0_tdp,big_delta_de_tdp)

for_plot$source <- as.factor(for_plot$source)
for_plot %>%   
    left_join(chx_g_e,by = c("symbol" = 'gene_name')) %>% 
    left_join(master_table %>% distinct(gene,new_cate,cryptic),by = c("symbol" = 'gene')) %>% 
    filter(new_cate %in% c("Early","Late")) %>% 
    filter(!is.na(nmd_status)) %>% 
    filter(cryptic == 'yes' & nmd_status == 'NMD rescued' & log2FoldChange > 1)
ggplot(aes(x = source, y= log2FoldChange,fill = new_cate)) + 
    facet_grid(vars(nmd_status),vars(cryptic)) + 
    geom_boxplot()
ggplot() +
    # geom_col(aes(x = as.character(source), y = log2FoldChange, fill = symbol), 
    #          data = big_delta_de_tdp, alpha = 0.3, show.legend = F) +
    geom_point(aes(x = as.character(source), y = log2FoldChange, color = symbol, group = symbol), ) +
    geom_line(aes(x = as.character(source), y = log2FoldChange, color = symbol, group = symbol), show.legend = F) +
    scale_fill_manual(values = "black") +
    facet_wrap(~nmd_status)
# scale_color_manual(values = c("#A6CEE3", "#1F78B4", "#FF7F00", "#FDBF6F", "#FB9A99", "#E31A1C"))+ # "#B2DF8A" "#33A02C" "#CAB2D6" "#6A3D9A")) +
scale_y_continuous(limits = c(-6.2,0.6)) +
    labs(x = "Doxycycline concentration", y = "Log2 Fold Change compared to control", color = "", fill = "Gene") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90))
a
#ggsave("~/Desktop/curves_nmd.png", width = 6, height = 4)



big_delta_de_dz <- read.csv(here::here('data','big_delta_de_dz.csv'))


big_delta_de_dz_tdp <- big_delta_de_dz[,c(4,8,10,11)] %>% 
    dplyr::filter(symbol == "TARDBP")
for_plot_dz_0_tdp <- big_delta_de_dz_tdp[1,] %>% 
    mutate(log2FoldChange = 0, padj = 0, source = 0)
big_delta_de_dz_tdp <- rbind(for_plot_dz_0_tdp,big_delta_de_dz_tdp)

for_plot_dz <- big_delta_de_dz[,c(4,8,10,11)] %>%
    dplyr::filter(symbol == "STMN2" | symbol == "AARS1" | symbol == "HDGFL2" | symbol == "UNC13A" |
                      symbol == "SYNE1" | symbol == "CYFIP2")
for_plot_dz_0 <- for_plot_dz %>% 
    distinct(symbol, .keep_all = T) %>% 
    mutate(log2FoldChange = 0, padj = 0, source = 0)
for_plot_dz <- rbind(for_plot_dz_0,for_plot_dz) %>% 
    mutate(symbol = factor(symbol, levels = c("AARS1", "HDGFL2", "SYNE1", "CYFIP2", "UNC13A", "STMN2")))

for_plot_dz$source <- as.factor(for_plot_dz$source)
b <- for_plot_dz %>%   
    ggplot() +
    geom_col(aes(x = as.character(source), y = log2FoldChange, fill = symbol), data = big_delta_de_dz_tdp, alpha = 0.3) +
    geom_point(aes(x = source, y = log2FoldChange, color = symbol, group = symbol)) +
    geom_line(aes(x = source, y = log2FoldChange, color = symbol, group = symbol)) +
    scale_fill_manual(values = "black") +
    # scale_color_manual(values = c("#A6CEE3", "#1F78B4", "#FF7F00", "#FDBF6F", "#FB9A99", "#E31A1C"))+ # "#B2DF8A" "#33A02C" "#CAB2D6" "#6A3D9A")) +
    scale_y_continuous(limits = c(-5.2,0.5)) +
    labs(x = "Doxycycline concentration", y = "Log2 Fold Change compared to control", color = "", fill = "Gene") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90))
b
#ggsave("~/Desktop/curves_nmd_dz.png", width = 6, height = 4)

#library(patchwork)

a+b

#ggsave("~/Desktop/caio.png", caio, width = 8, height = 4)




print("3G")




