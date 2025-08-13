library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(rlang)
generate_summary_table <- function(
        summary_post, 
        junction_table, 
        category_column = "color_gene_name", 
        detected_threshold_for_observed = 2,
        detected_threshold_for_bio = 5,
        selective_fpr_threshold = 0.1,
        selective_tpr_threshold = 0.1
) {
    # Convert column names to symbols for tidy evaluation
    category_col <- sym(category_column)
    
    # Step 1: Get detected and selective junctions from postmortem summary
    detected_junctions <- summary_post %>% 
        filter(n_samples_tdp_path >= detected_threshold_for_observed) %>% 
        select(paste_into_igv_junction, group) %>% 
        mutate(detected = TRUE)
    
    selective_junctions <- summary_post %>% 
        filter(n_samples_tdp_path >= detected_threshold_for_bio) %>% 
        filter(fpr_value <= selective_fpr_threshold & tpr_value >= selective_tpr_threshold) %>% 
        select(paste_into_igv_junction, group) %>% 
        mutate(selective = TRUE)
    
    # Step 2: Get unique junctions with their category classification
    unique_junctions <- junction_table %>% 
        ungroup() %>% 
        distinct(paste_into_igv_junction, !!category_col) %>%
        group_by(!!category_col) %>%
        mutate(n_type = n_distinct(paste_into_igv_junction)) %>%
        ungroup()
    
    # Step 3: Get all unique groups (excluding NA)
    all_groups <- detected_junctions %>%
        filter(!is.na(group)) %>%
        pull(group) %>%
        unique()

    # Step 4: Create all combinations of junctions and groups
    all_junction_group_combinations <- unique_junctions %>%
        crossing(group = all_groups)
    
    # Step 5: Join with detection data (set FALSE for missing combinations)
    complete_junction_dataset <- all_junction_group_combinations %>%
        left_join(
            detected_junctions,
            by = c("paste_into_igv_junction", "group")
        ) %>%
        mutate(detected = ifelse(is.na(detected),FALSE,detected)) %>% 
        left_join(selective_junctions) %>% 
        mutate(selective = ifelse(is.na(selective),FALSE,selective)) %>% 
        group_by(!!category_col,group) %>% 
        summarize(n_type,
                  n_detected = sum(detected), 
                  n_selective = sum(selective)) %>% 
        ungroup() %>% 
        unique() %>% 
        mutate(n_undetected = n_type - n_detected) %>% 
        mutate(n_detected_not_selective = n_detected - n_selective) %>% 
        mutate(detection_rate = n_detected / n_type) %>% 
        mutate(selectivity_rate = n_selective / n_detected)
        
        
    return(complete_junction_dataset)
    
}



run_fisher_test = function(summary_table, category_column = "color_gene_name"){

    category_col <- sym(category_column)

    detection_nested = summary_table %>%
        dplyr::rename(comparison = !!category_col) %>% 
        distinct(comparison, group, n_detected, n_undetected) %>%
        nest(data = -group) %>% 
        mutate(
            crosstab = map(data, ~{
                .x %>%
                    select(comparison, n_detected, n_undetected) %>%
                    column_to_rownames("comparison") %>%
                    as.matrix()
            })
        )
    
    detection_nested$fisher = purrr::map(detection_nested$crosstab, ~ broom::tidy(fisher.test(.x)))
    detection_nested = detection_nested %>% distinct(group,fisher) %>% unnest(fisher) %>% mutate(rate = 'detection')
    
    
    
    selective_nested = summary_table %>%
        dplyr::rename(comparison = !!category_col) %>% 
        distinct(comparison, group, n_selective, n_detected_not_selective) %>%
        nest(data = -group) %>% 
        mutate(
            crosstab = map(data, ~{
                .x %>%
                    select(comparison, n_selective, n_detected_not_selective) %>%
                    column_to_rownames("comparison") %>%
                    as.matrix()
            })
        )
    
    selective_nested$fisher = purrr::map(selective_nested$crosstab, ~ broom::tidy(fisher.test(.x)))
    selective_nested = selective_nested %>% distinct(group,fisher) %>% unnest(fisher) %>% mutate(rate = 'selectivity')
    
    return(list(detection_fisher = detection_nested,selective_fisher = selective_nested))
    


}


make_final_plot = function(summary_table, 
                           fisher_result, 
                           category_column = "color_gene_name",
                           color_choices = c("#3288BD","#D53E4F")){

    category_col <- sym(category_column)
    
    det_fish = fisher_result$detection_fisher %>% 
        mutate(p_label = case_when(
            p.value < 0.001 ~ "***",
            p.value < 0.01 ~ "**",
            p.value < 0.05 ~ "*",
            TRUE ~ "ns"
        )) 
    

    det_plot_data <- summary_table %>% 
        mutate(fill_col = paste0(!!category_col, ' (n=', n_type, ')')) %>% 

        group_by(group) %>%
        mutate(max_detection_rate = max(detection_rate),
               p_y_pos = max_detection_rate + 0.02) %>%
        ungroup()
    
    det_fish = det_fish %>% 
        left_join(det_plot_data %>% distinct(group,p_y_pos))

    
    detection_plot = det_plot_data %>% 
        ggplot( aes(x = group, y = detection_rate, fill = fill_col)) +
        geom_bar(stat = "identity", position = "dodge") +
        geom_text(aes(label = n_detected), 
                  position = position_dodge(width = 0.9), 
                  vjust = -0.5, 
                  size = 3.5) +
        geom_text(data = det_fish, aes(y = p_y_pos, label = p_label,x = group),inherit.aes = FALSE, 
                  position = position_dodge(width = 0.9),
                  vjust = -0.5,
                  size = 10,
                  fontface = "bold") +
        labs(title = "Detection Rate by Group",
             x = "Group",
             y = "Detection Rate",
             fill = "") +
        scale_y_continuous(limits = c(0, max(det_plot_data$detection_rate) * 1.1)) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom",
              plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = color_choices) +
    scale_y_continuous(labels = scales::percent)
    
    
    
    
    
    sel_fish = fisher_result$selective_fisher %>% 
        mutate(p_label = case_when(
            p.value < 0.001 ~ "***",
            p.value < 0.01 ~ "**",
            p.value < 0.05 ~ "*",
            TRUE ~ "ns"
        )) 
    
    sel_plot_data <- summary_table %>% 
        mutate(fill_col = paste0(!!category_col, ' (n=', n_type, ')')) %>% 
        group_by(group) %>%
        mutate(max_detection_rate = max(selectivity_rate),
               p_y_pos = max_detection_rate + 0.02) %>%
        ungroup()
    
    sel_fish = sel_fish %>% 
        left_join(sel_plot_data %>% distinct(group,p_y_pos))
    
    
    selectivity_plot = sel_plot_data %>% 
        ggplot( aes(x = group, y = selectivity_rate, fill = fill_col)) +
        geom_bar(stat = "identity", position = "dodge") +
        geom_text(aes(label = n_selective), 
                  position = position_dodge(width = 0.9), 
                  vjust = -0.5, 
                  size = 3.5) +
        geom_text(data = sel_fish, aes(y = p_y_pos, label = p_label,x = group),inherit.aes = FALSE, 
                  position = position_dodge(width = 0.9),
                  vjust = -0.5,
                  size = 10,
                  fontface = "bold") +
        labs(title = "Selectivity Rate by Group",
             x = "Group",
             y = "Selectivity Rate",
             fill = "") +
        scale_y_continuous(limits = c(0, max(sel_plot_data$selectivity_rate) * 1.1)) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom",
              plot.title = element_text(hjust = 0.5)) +
        scale_fill_manual(values = color_choices) + 
        scale_y_continuous(labels = scales::percent)
    
    
    return_list = list(detection_plot = detection_plot, selectivity_plot = selectivity_plot)
    return(return_list)
    
    
}



