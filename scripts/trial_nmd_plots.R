


slope_plot_nmd <- function(big_data, gene_list, input_list = c("Control_Control",
                                                                             "Cycloheximide_Control", 
                                                                             "Control_TDP43KD", 
                                                                             "Cycloheximide_TDP43KD"),
                       nmd_label = "CHX treatment",
                       true_dark_junctions = c("chr7:1898390-1957629")){

big_data_filtereds <- big_data %>%
    filter(.id %in% input_list) %>% 
  group_by(lsv_junc) %>%
  dplyr::filter(n() == 4) %>% #Only junctions that are 
  dplyr::filter(de_novo_junctions == 1) %>%
  dplyr::filter(mean_psi_per_lsv_junction[.id == input_list[1]] < 0.05) %>% #control has a PSI less than 5%
  dplyr::filter(mean_psi_per_lsv_junction[.id == input_list[2]] < 0.05) %>% #control with NMDi has a PSI less than 5%
  #dplyr::filter(mean_psi_per_lsv_junction[.id == input_list[3]] < 0.05) %>% 
  dplyr::filter(mean_psi_per_lsv_junction[.id == input_list[4]] > 0.1) %>% #KD with NMDi has a PSI greater than 10%
  dplyr::filter((mean_psi_per_lsv_junction[.id == input_list[3]] - mean_psi_per_lsv_junction[.id == input_list[1]] > 0.05 |
           mean_psi_per_lsv_junction[.id == input_list[4]] - mean_psi_per_lsv_junction[.id == input_list[1]] > 0.05)) %>% #KD is more than 5% higher than control OR KD with NMDi is more than 5% higher
  #dplyr::filter(mean_psi_per_lsv_junction[.id == input_list[4]] - mean_psi_per_lsv_junction[.id == input_list[3]] > 0) %>%# &
    mutate(color_gene_name = as.character(ifelse(
        mean_psi_per_lsv_junction[.id == input_list[4]] - mean_psi_per_lsv_junction[.id == input_list[3]] > 0.1 |
            (mean_psi_per_lsv_junction[.id == input_list[3]] < 0.09 & mean_psi_per_lsv_junction[.id == input_list[4]] > 0.1), 
        "NMD rescued", "non-NMD rescued")))
  # mutate(alpha_gene_name = as.character(ifelse(mean_psi_per_lsv_junction[.id == input_list[4]] > 0.25 &
  #                                                gene_name %in% gene_list,
  #                                                # paste_into_igv_junction %in% c("chr19:17641556-17642414",
  #                                                #                                "chr6:152247944-152249161",
  #                                                #                                "chr8:79611214-79616822",
  #                                                #                                "chr19:4492152-4493703",
  #                                                #                                "chr5:157361336-157361468",
  #                                                #                                "chr16:70271972-70272796"), 
  #                                              1, 0))) %>%
  # 
  # mutate(label_junction = case_when(.id == input_list[4] & mean_psi_per_lsv_junction[.id == input_list[4]] > 0.2 &
  #                                     alpha_gene_name == 1 ~ gene_name, T ~ NA_character_))
                                           #gene_name %in% gene_list ~ gene_name, T ~ ""))
highlighting_table = big_data_filtereds %>% 
    filter(gene_name %in% gene_list) %>% 
    filter(.id == c(input_list[4])) %>% 
    group_by(gene_name) %>% 
    slice_max(mean_psi_per_lsv_junction) %>% 
    mutate(alpha_gene_name = 1) %>% 
    mutate(label_junction = gene_name) %>% 
    distinct(paste_into_igv_junction,gene_name,alpha_gene_name,label_junction)


higlight_dark = big_data_filtereds %>% 
    filter(gene_name %in% gene_list) %>% 
    filter(.id == c(input_list[4])) %>% 
    group_by(gene_name) %>% 
    slice_max(mean_psi_per_lsv_junction) %>% 
    mutate(alpha_gene_name = 1) %>% 
    mutate(label_junction = gene_name) %>% 
    distinct(paste_into_igv_junction,gene_name,alpha_gene_name,label_junction) %>% 
    filter(paste_into_igv_junction %in% true_dark_junctions)
#big_data_filtereds_list <- big_data_filtereds %>% filter(color_gene_name == "Delta PSI > 0.05") %>% pull(gene_name) %>% unique()

#rite.table(big_data_filtereds, "~/Desktop/nmd_or_not.csv", quote = F, row.names = F, sep = ",")
true_dark_table = big_data_filtereds %>%
    dplyr::filter(.id %in% c(input_list[3],input_list[4])) %>%
    left_join(highlighting_table) %>% 
    mutate(label_junction = ifelse(.id == input_list[3], NA_character_,label_junction)) %>% 
    mutate(alpha_gene_name = ifelse(is.na(alpha_gene_name),0.4,alpha_gene_name)) %>% 
    group_by(paste_into_igv_junction,.id) %>% 
    mutate(label_junction = ifelse(mean_psi_per_lsv_junction == max(mean_psi_per_lsv_junction),label_junction, NA_character_)) %>% 
    mutate(alpha_gene_name = ifelse(mean_psi_per_lsv_junction == max(mean_psi_per_lsv_junction),alpha_gene_name,0.4)) %>% 
    ungroup() %>% 
    mutate(with_dark_color = color_gene_name) %>% 
    filter(paste_into_igv_junction %in% true_dark_junctions)

ploss = big_data_filtereds %>%
  dplyr::filter(.id %in% c(input_list[3],input_list[4])) %>%
     left_join(highlighting_table) %>% 
     mutate(label_junction = ifelse(.id == input_list[3], NA_character_,label_junction)) %>% 
     mutate(alpha_gene_name = ifelse(is.na(alpha_gene_name),0.4,alpha_gene_name)) %>% 
    group_by(paste_into_igv_junction,.id) %>% 
    mutate(label_junction = ifelse(mean_psi_per_lsv_junction == max(mean_psi_per_lsv_junction),label_junction, NA_character_)) %>% 
    mutate(alpha_gene_name = ifelse(mean_psi_per_lsv_junction == max(mean_psi_per_lsv_junction),alpha_gene_name,0.4)) %>% 
    ungroup() %>% 
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
  # scale_alpha_manual(values = c(0.2,1)) +
  xlab("") +
  ylab("PSI") +
  scale_x_discrete(labels = c("TDP43-KD", paste0("TDP43-KD\n+\n",nmd_label))) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1), labels = c("0%", "25%", "50%", "75%", "100%")) +
  theme_classic() +
  theme(text = element_text(size = 18))
ploss = ploss + 
    geom_point(data = true_dark_table,color = '#4A4747') + 
    geom_line(data = true_dark_table, aes(group = lsv_junc), color = '#4A4747',alpha = 0.3)
plot(ploss)

return_list = list(ploss,big_data_filtereds)
names(return_list) = c('plot','filtered_table')
#ggsave(filename = "~/Desktop/upf1_slope.png", ploss)
return(return_list)
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



amigoingmad = function (package = "dplyr", fix = TRUE, iteration = 0) 
{
    if (iteration > 1) {
        stop("Can't fix.")
    }
    conf = unique(conflicts())
    want_package = paste0("package:", package)
    conflicts_desired_package = conf[conf %in% ls(want_package)]
    conflict_envs = sapply(conflicts_desired_package, FUN = function(x) {
        environmentName(pryr::where(x, globalenv()))
    })
    is_good = conflict_envs == want_package
    potentially_bad_confs = conflicts_desired_package[!is_good]
    potentially_bad_envs = conflict_envs[!is_good]
    have_to_fix = rep(FALSE, length(potentially_bad_confs))
    for (i in seq_along(potentially_bad_confs)) {
        if (!identical(body(get(potentially_bad_confs[i], pos = want_package)), 
                       body(get(potentially_bad_confs[i])))) {
            have_to_fix[i] = TRUE
        }
    }
    if (any(have_to_fix)) {
        message("The following functions don't have the environment you want.")
        print(data.frame(function. = potentially_bad_confs[have_to_fix], 
                         environment = potentially_bad_envs[have_to_fix]), 
              row.names = F)
        if (fix) {
            base::detach(name = want_package, character.only = TRUE)
            base::library(package, character.only = TRUE)
            message("Tried to fix this, calling myself again to make sure...")
            amigoingmad(package, fix, iteration + 1)
            message("Sanity restored!")
        }
    }
    else if (iteration == 0) {
        message("Everything looks normal. Maybe it's you.")
    }
}
