slope_plot_nmd <- function(big_data, gene_list){
big_data_filtereds <- big_data %>%
  group_by(lsv_junc) %>%
  dplyr::filter(n() == 4) %>% #???
  dplyr::filter(de_novo_junctions == 1) %>%
  dplyr::filter(mean_psi_per_lsv_junction[.id == input_list[1]] < 0.05) %>%
  dplyr::filter(mean_psi_per_lsv_junction[.id == input_list[2]] < 0.05) %>%
  #dplyr::filter(mean_psi_per_lsv_junction[.id == input_list[3]] < 0.05) %>% 
  dplyr::filter(mean_psi_per_lsv_junction[.id == input_list[4]] > 0.1) %>%
  dplyr::filter((mean_psi_per_lsv_junction[.id == input_list[3]] - mean_psi_per_lsv_junction[.id == input_list[1]] > 0.05 |
           mean_psi_per_lsv_junction[.id == input_list[4]] - mean_psi_per_lsv_junction[.id == input_list[1]] > 0.05)) %>%
  #dplyr::filter(mean_psi_per_lsv_junction[.id == input_list[4]] - mean_psi_per_lsv_junction[.id == input_list[3]] > 0) %>%# &
  mutate(color_gene_name = as.character(ifelse(mean_psi_per_lsv_junction[.id == input_list[4]] - mean_psi_per_lsv_junction[.id == input_list[3]] > 0.1, 
                                               "NMD rescued", "non-NMD rescued"))) %>%
  mutate(alpha_gene_name = as.character(ifelse(mean_psi_per_lsv_junction[.id == input_list[4]] > 0.25 &
                                                 gene_name %in% gene_list &
                                                 paste_into_igv_junction %in% c("chr19:17641556-17642414",
                                                                                "chr6:152247944-152249161",
                                                                                "chr8:79611214-79616822",
                                                                                "chr19:4492152-4493703",
                                                                                "chr5:157361336-157361468",
                                                                                "chr16:70271972-70272796"), 1, 0))) %>%
  
  
  #"chr19:17641556-17642414","chr8:79611214-79616822"), #UNC13A, STMN2
                                                #                                "chr19:7169831-7170537", "chr19:7168102-7169727", "chr16:70271972-70272796") #AARS1
                                                                                #"chr5:157361336-157361468", #CYFIP2
                                                                                #"chr6:152247944-152249161" #, "chr6:152244656-152247823" #SYNE1
                                                                              
                                                                                #"chr20:63439708-63444659","chr10:3099819-3101365",
                                                                                #"chr2:241668985-241670726","chr19:7169831-7170537"
                                                #                                1, 0))) %>% 
  mutate(label_junction = case_when(.id == input_list[4] & mean_psi_per_lsv_junction[.id == input_list[4]] > 0.2 &
                                      alpha_gene_name == 1 ~ gene_name, T ~ ""))
                                           #gene_name %in% gene_list ~ gene_name, T ~ ""))
table(big_data_filtereds$color_gene_name)

#big_data_filtereds_list <- big_data_filtereds %>% filter(color_gene_name == "Delta PSI > 0.05") %>% pull(gene_name) %>% unique()

#rite.table(big_data_filtereds, "~/Desktop/nmd_or_not.csv", quote = F, row.names = F, sep = ",")

ploss <- big_data_filtereds %>%
  dplyr::filter(.id %in% c(input_list[3],input_list[4])) %>%
  ggplot(mapping = aes(x = .id, y = mean_psi_per_lsv_junction)) +
  facet_wrap(facets = vars(color_gene_name)) +
  geom_point(aes(color = color_gene_name, alpha = alpha_gene_name, group = lsv_junc), show.legend = F) + 
  geom_line(aes(color = color_gene_name, alpha = alpha_gene_name, group = lsv_junc), show.legend = F) +
  geom_text_repel(aes(x = .id, y = mean_psi_per_lsv_junction, label = label_junction), #point.padding = 0.3,
                  #nudge_y = 0.2, min.segment.length = 0.5, box.padding  = 2, 
                  max.overlaps = Inf
                  #, size=4, show.legend = F
                  ) +
  scale_color_manual(values = c(#"#ABDDA4", 
                                "#3288BD",
                                "#D53E4F")) +
  scale_alpha_manual(values = c(0.02,1)) +
  xlab("") +
  ylab("PSI") +
  scale_x_discrete(labels = c("TDP43-KD", "TDP43-KD + \n NMD inhibition")) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1), labels = c("0%", "25%", "50%", "75%", "100%")) +
  theme_classic() +
  theme(text = element_text(size = 18))
plot(ploss)
#ggsave(filename = "~/Desktop/upf1_slope.png", ploss)
return(ploss)
}

#avg_read_counts <- featureCounts %>%
#  mutate(Control_Control = rowMeans(dplyr::select(featureCounts, contains("CTRL_ctrl")), na.rm = TRUE)) %>%
#  mutate(Cycloheximide_Control = rowMeans(dplyr::select(featureCounts, contains("CTRL_chx")), na.rm = TRUE)) %>%
#  mutate(Control_TDP43KD = rowMeans(dplyr::select(featureCounts, contains("DOX_ctrl")), na.rm = TRUE)) %>%
#  mutate(Cycloheximide_TDP43KD = rowMeans(dplyr::select(featureCounts, contains("DOX_chx")), na.rm = TRUE)) %>%
#  dplyr::select(c(18:22)) %>%
#  pivot_longer(cols = c("Control_Control", "Cycloheximide_Control", "Control_TDP43KD", "Cycloheximide_TDP43KD"), 
#               names_to = ".id", values_to = "avg_count") %>%
#  group_by(gene_name, .id) %>%
#  summarise(avg_count = max(avg_count))

#big_data_avg_counts <- left_join(big_data, avg_read_counts)
#big_data_avg_counts_annot <- left_join(big_data_avg_counts, big_delta_filter)

#big_data_avg_counts_annot %>%
#  filter(.id %in% c("Control_TDP43KD","Cycloheximide_TDP43KD")) %>%
#  ggplot(mapping = aes(x = .id, y = avg_count)) +
  #facet_wrap(facets = vars(junc_cat)) +
#  geom_point(aes(color = color_gene_name, alpha = alpha_gene_name, group = lsv_junc), show.legend = T) + 
  #geom_violin(aes(x = .id, y = mean_psi_per_lsv_junction_normal_fake, group = color_gene_name, color = color_gene_name), show.legend = F) +
#  geom_line(aes(color = color_gene_name, alpha = alpha_gene_name, group = lsv_junc), show.legend = F) +
#  geom_text_repel(aes(label = label_junction), point.padding = 0.3,
#                  nudge_y = 0.2, min.segment.length = 0.5, box.padding  = 2, max.overlaps = Inf, size=4, show.legend = F) +
#  scale_color_manual(values = c(#"#ABDDA4", 
#                                "#3288BD",
#                                "#D53E4F")) +
#  scale_alpha_manual(values = c(0.1,1)) +
  #scale_x_log10() +
#  scale_y_log10() +
#  theme_classic()