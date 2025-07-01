library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(rlang)
analyze_junction_detection_and_selectivity <- function(
        summary_post, 
        junction_table, 
        category_column = "color_gene_name", 
        baseline_category = "non-NMD rescued",
        detected_threshold_for_observed = 2,
        detected_threshold_for_bio = 5,
        selective_fpr_threshold = 0.1,
        selective_tpr_threshold = 0.1,
        color_choices = c("#1F77B4", "#FF7F0E")
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
        mutate(detected = coalesce(detected, FALSE))
    
    # Step 6: Summarize for detection analysis
    detection_summary <- complete_junction_dataset %>% 
        group_by(!!category_col, group, detected) %>% 
        summarise(
            n_type = first(n_type),
            n_junc = n_distinct(paste_into_igv_junction),
            .groups = "drop"
        ) %>% 
        mutate(detection_ratio = n_junc / n_type)
    
    # Step 7: Create detection rate summary for plotting
    detection_rate_summary <- detection_summary %>%
        group_by(!!category_col, group) %>%
        summarize(
            detected_count = sum(n_junc[detected == TRUE]),
            not_detected_count = sum(n_junc[detected == FALSE]),
            total_count = sum(n_junc),
            detection_rate = detected_count / total_count,
            .groups = "drop"
        )
    
    # Step 8: Perform Fisher's exact test for detection
    detection_fisher_results <- data.frame(
        group = character(),
        p_value = numeric(),
        odds_ratio = numeric(),
        ci_lower = numeric(),
        ci_upper = numeric(),
        stringsAsFactors = FALSE
    )
    
    for (g in all_groups) {
        # Filter data for current group
        group_data <- detection_summary %>% filter(group == g)
        
        # Get all categories
        categories <- unique(group_data[[category_column]])
        
        # Only proceed if we have exactly 2 categories (for 2x2 contingency table)
        if (length(categories) == 2) {
            # Identify the comparison category (not baseline)
            comparison_category <- setdiff(categories, baseline_category)
            
            # Create contingency table for detected vs not detected
            cont_table <- matrix(
                c(
                    sum(group_data$n_junc[group_data[[category_column]] == baseline_category & group_data$detected == FALSE]),
                    sum(group_data$n_junc[group_data[[category_column]] == baseline_category & group_data$detected == TRUE]),
                    sum(group_data$n_junc[group_data[[category_column]] == comparison_category & group_data$detected == FALSE]),
                    sum(group_data$n_junc[group_data[[category_column]] == comparison_category & group_data$detected == TRUE])
                ),
                nrow = 2,
                byrow = TRUE,
                dimnames = list(
                    c(baseline_category, comparison_category),
                    c("Not Detected", "Detected")
                )
            )
            
            # Perform Fisher's exact test
            fisher_result <- fisher.test(cont_table)
            
            # Store results
            detection_fisher_results <- rbind(detection_fisher_results, data.frame(
                group = g,
                p_value = fisher_result$p.value,
                odds_ratio = fisher_result$estimate,
                ci_lower = fisher_result$conf.int[1],
                ci_upper = fisher_result$conf.int[2],
                stringsAsFactors = FALSE
            ))
        }
    }
    
    # Format p-values for display
    detection_fisher_results$p_formatted <- ifelse(
        detection_fisher_results$p_value < 0.001, 
        "p < 0.001",
        paste0("p = ", sprintf("%.3f", detection_fisher_results$p_value))
    )
    
    # Step 9: Selectivity analysis (only on detected junctions)
    selectivity_analysis <- complete_junction_dataset %>%
        # Only analyze junctions that are detected
        filter(detected == TRUE) %>%
        # Add selectivity information
        left_join(selective_junctions, by = c("paste_into_igv_junction", "group")) %>%
        # Handle NAs in selective column
        mutate(selective = coalesce(selective, FALSE))
    
    # Step 10: Summarize for selectivity analysis
    selectivity_summary <- selectivity_analysis %>%
        group_by(!!category_col, group) %>%
        summarize(
            total_detected = n(),
            selective_count = sum(selective),
            non_selective_count = total_detected - selective_count,
            selectivity_rate = selective_count / total_detected,
            .groups = "drop"
        )
    
    # Step 11: Perform Fisher's exact test for selectivity
    selectivity_fisher_results <- data.frame(
        group = character(),
        p_value = numeric(),
        odds_ratio = numeric(),
        ci_lower = numeric(),
        ci_upper = numeric(),
        stringsAsFactors = FALSE
    )
    
    for (g in all_groups) {
        # Filter data for current group
        group_data <- selectivity_summary %>% filter(group == g)
        
        # Get all categories
        categories <- unique(group_data[[category_column]])
        
        # Only proceed if we have exactly 2 categories (for 2x2 contingency table)
        if (length(categories) == 2) {
            # Identify the comparison category (not baseline)
            comparison_category <- setdiff(categories, baseline_category)
            
            # Create contingency table for selective vs non-selective among detected junctions
            cont_table <- matrix(
                c(
                    group_data$selective_count[group_data[[category_column]] == baseline_category],
                    group_data$non_selective_count[group_data[[category_column]] == baseline_category],
                    group_data$selective_count[group_data[[category_column]] == comparison_category],
                    group_data$non_selective_count[group_data[[category_column]] == comparison_category]
                ),
                nrow = 2,
                byrow = TRUE,
                dimnames = list(
                    c(baseline_category, comparison_category),
                    c("Selective", "Non-selective")
                )
            )
            
            # Perform Fisher's exact test
            fisher_result <- fisher.test(cont_table)
            
            # Store results
            selectivity_fisher_results <- rbind(selectivity_fisher_results, data.frame(
                group = g,
                p_value = fisher_result$p.value,
                odds_ratio = fisher_result$estimate,
                ci_lower = fisher_result$conf.int[1],
                ci_upper = fisher_result$conf.int[2],
                stringsAsFactors = FALSE
            ))
        }
    }
    
    # Format p-values for display
    selectivity_fisher_results$p_formatted <- ifelse(
        selectivity_fisher_results$p_value < 0.001, 
        "p < 0.001",
        paste0("p = ", sprintf("%.3f", selectivity_fisher_results$p_value))
    )

    # Step 12: Create detection visualization
    detection_plot <- create_detection_visualization(
        detection_rate_summary, 
        detection_fisher_results,
        category_column = category_column,
        baseline_category = baseline_category,
        color_choices
    )
    
    # Step 13: Create selectivity visualization
    selectivity_plot <- create_selectivity_visualization(
        selectivity_summary, 
        selectivity_fisher_results,
        category_column = category_column,
        baseline_category = baseline_category,
        color_choices
    )
    
    # Step 14: Create combined visualizations
    combined_rate_plot <- create_combined_visualization(
        detection_rate_summary,
        selectivity_summary,
        detection_fisher_results,
        selectivity_fisher_results,
        category_column = category_column,
        baseline_category = baseline_category,
        color_choices = color_choices
    )
    
    combined_odds_plot <- create_combined_odds_ratio_plot(
        detection_fisher_results,
        selectivity_fisher_results,
        category_column = category_column,
        baseline_category = baseline_category
    )
    
    # Step 15: Create summary tables
    detection_summary_table <- detection_fisher_results %>%
        mutate(
            odds_ratio_ci = paste0(sprintf("%.2f", odds_ratio),
                                   " (", sprintf("%.2f", ci_lower), "-",
                                   sprintf("%.2f", ci_upper), ")"),
            interpretation = ifelse(
                p_value < 0.05,
                ifelse(
                    odds_ratio > 1,
                    paste0(comparison_category, " more likely to be detected"),
                    paste0(baseline_category, " more likely to be detected")
                ),
                "No significant difference"
            )
        ) %>%
        select(group, p_value, odds_ratio_ci, interpretation)
    
    selectivity_summary_table <- selectivity_fisher_results %>%
        mutate(
            odds_ratio_ci = paste0(sprintf("%.2f", odds_ratio), 
                                   " (", sprintf("%.2f", ci_lower), "-", 
                                   sprintf("%.2f", ci_upper), ")"),
            interpretation = ifelse(
                p_value < 0.05,
                ifelse(
                    odds_ratio > 1,
                    paste0(baseline_category, " more likely to be selective (among detected)"),
                    paste0(comparison_category, " more likely to be selective (among detected)")
                ),
                "No significant difference in selectivity"
            )
        ) %>%
        select(group, p_value, odds_ratio_ci, interpretation)
    
    # Return all results in a list
    return(list(
        # Data summaries
        detection_summary = detection_summary,
        detection_rate_summary = detection_rate_summary,
        selectivity_analysis = selectivity_analysis,
        selectivity_summary = selectivity_summary,
        
        # Statistical results
        detection_fisher_results = detection_fisher_results,
        selectivity_fisher_results = selectivity_fisher_results,
        
        # Visualization plots
        detection_plot = detection_plot,
        selectivity_plot = selectivity_plot,
        combined_rate_plot = combined_rate_plot,
        combined_odds_plot = combined_odds_plot,
        
        # Summary tables
        detection_summary_table = detection_summary_table,
        selectivity_summary_table = selectivity_summary_table,
        
        # Parameters used
        parameters = list(
            category_column = category_column,
            baseline_category = baseline_category,
            detected_threshold_for_bio = detected_threshold_for_bio,
            detected_threshold_for_observed = detected_threshold_for_observed,
            selective_fpr_threshold = selective_fpr_threshold,
            selective_tpr_threshold = selective_tpr_threshold
        )
    ))
}

# Updated helper functions with category_column and baseline_category parameters

# Updated helper function to create selectivity visualization
create_selectivity_visualization <- function(selectivity_summary, selectivity_fisher_results, 
                                             category_column = "color_gene_name", 
                                             baseline_category = "non-NMD rescued",
                                             color_choices) {
    # Reorder the groups for better visualization
    selectivity_summary$group <- factor(selectivity_summary$group, 
                                        levels = c("Cord", "Motor Cortex", "FTD frontal/temporal", "RiMOD"))
    
    # Get comparison category
    categories <- unique(selectivity_summary[[category_column]])
    comparison_category <- setdiff(categories, baseline_category)
    
    # Create the selectivity rate plot
    plot <- ggplot(selectivity_summary, aes_string(x = "group", y = "selectivity_rate", fill = category_column)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
        geom_text(
            aes(
                label = paste0(round(selectivity_rate * 100, 1), "%"),
                y = selectivity_rate + 0.03
            ),
            position = position_dodge(width = 0.9),
            size = 3
        ) +
        scale_fill_manual(values = setNames(color_choices, c(baseline_category, comparison_category))) +
        scale_y_continuous(
            labels = percent_format(),
            limits = c(0, 1),
            expand = expansion(mult = c(0, 0.1))
        ) +
        labs(
            title = "Selectivity Rates Among Detected Junctions",
            subtitle = paste("Comparing", baseline_category, "vs.", comparison_category, "(among detected junctions only)"),
            x = "Group",
            y = "Proportion of Detected Junctions that are Selective",
            fill = "Junction Type"
        ) +
        theme_minimal() +
        theme(
            legend.position = "bottom",
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.minor = element_blank()
        )
    
    # Add Fisher's exact test p-values to the plot
    plot_with_stats <- plot + 
        geom_text(
            data = selectivity_fisher_results,
            aes(x = group, y = 0.9, label = p_formatted),
            size = 3.5,
            inherit.aes = FALSE
        )
    
    return(plot_with_stats)
}

# Updated helper function to create detection visualization
create_detection_visualization <- function(detection_rate_summary, detection_fisher_results, 
                                           category_column = "color_gene_name", 
                                           baseline_category = "non-NMD rescued",
                                           color_choices) {
    # Reorder the groups for better visualization
    detection_rate_summary$group <- factor(detection_rate_summary$group, 
                                           levels = c("Cord", "Motor Cortex", "FTD frontal/temporal", "RiMOD"))
    
    # Get comparison category
    categories <- unique(detection_rate_summary[[category_column]])
    comparison_category <- setdiff(categories, baseline_category)
    
    # Create the detection rate plot
    plot <- ggplot(detection_rate_summary, aes_string(x = "group", y = "detection_rate", fill = category_column)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
        geom_text(
            aes(
                label = paste0(round(detection_rate * 100, 1), "%"),
                y = detection_rate + 0.03
            ),
            position = position_dodge(width = 0.9),
            size = 3
        ) +
        scale_fill_manual(values = setNames(color_choices, c(baseline_category, comparison_category))) +
        scale_y_continuous(
            labels = percent_format(),
            limits = c(0, 1),
            expand = expansion(mult = c(0, 0.1))
        ) +
        labs(
            title = "Detection Rates by Junction Type",
            x = "Group",
            y = "Detection Rate",
            fill = "Junction Type"
        ) +
        theme_minimal() +
        theme(
            legend.position = "bottom",
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.minor = element_blank()
        )
    
    # Add Fisher's exact test p-values to the plot
    plot_with_stats <- plot + 
        geom_text(
            data = detection_fisher_results,
            aes(x = group, y = 0.9, label = p_formatted),
            size = 3.5,
            inherit.aes = FALSE
        )
    
    return(plot_with_stats)
}

# Updated helper function to create combined visualization
create_combined_visualization <- function(detection_data, selectivity_data, 
                                          detection_fisher, selectivity_fisher,
                                          category_column = "color_gene_name", 
                                          baseline_category = "non-NMD rescued",
                                          color_choices) {
    # Ensure groups are ordered consistently
    detection_data$group <- factor(detection_data$group, 
                                   levels = c("Cord", "Motor Cortex", "FTD frontal/temporal", "RiMOD"))
    selectivity_data$group <- factor(selectivity_data$group, 
                                     levels = c("Cord", "Motor Cortex", "FTD frontal/temporal", "RiMOD"))
    
    # Get comparison category
    categories <- unique(detection_data[[category_column]])
    comparison_category <- setdiff(categories, baseline_category)
    
    # Prepare data for combined plot
    detection_for_plot <- detection_data %>%
        mutate(metric = "Detection Rate") %>%
        rename(rate = detection_rate)
    
    selectivity_for_plot <- selectivity_data %>%
        mutate(metric = "Selectivity Rate\n(among detected)") %>%
        rename(rate = selectivity_rate)
    
    # Combine the datasets
    combined_data <- bind_rows(detection_for_plot, selectivity_for_plot)
    
    # Similarly combine fisher results
    detection_fisher$metric <- "Detection Rate"
    selectivity_fisher$metric <- "Selectivity Rate\n(among detected)"
    combined_fisher <- bind_rows(detection_fisher, selectivity_fisher)
    
    # Create the combined plot
    combined_plot <- ggplot(combined_data, aes_string(x = "group", y = "rate", fill = category_column)) +
        facet_wrap(~ metric, scales = "free_y") +
        geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
        geom_text(
            aes(
                label = paste0(round(rate * 100, 1), "%"),
                y = rate + 0.02
            ),
            position = position_dodge(width = 0.9),
            size = 3
        ) +
        scale_fill_manual(values = setNames(color_choices, c(baseline_category, comparison_category))) +
        labs(
            title = paste("Comparison of", baseline_category, "vs", comparison_category, "Rates"),
            subtitle = "Detection rates (all junctions) and Selectivity rates (detected junctions only)",
            x = "Group",
            y = "Rate",
            fill = "Junction Type"
        ) +
        theme_minimal() +
        theme(
            legend.position = "bottom",
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = "lightgray", color = NA),
            strip.text = element_text(face = "bold")
        )
    
    # Add p-values to each facet
    combined_plot <- combined_plot +
        geom_text(
            data = combined_fisher,
            aes(x = group, y = Inf, label = p_formatted),
            vjust = 1.5,
            position = position_dodge(width = 0.9),
            inherit.aes = FALSE
        )
    
    return(combined_plot)
}

# Updated helper function to create combined odds ratio plot
create_combined_odds_ratio_plot <- function(detection_fisher, selectivity_fisher,
                                            category_column = "color_gene_name",
                                            baseline_category = "non-NMD rescued") {
    # Add a type column to distinguish detection from selectivity data
    detection_fisher$metric <- "Detection"
    selectivity_fisher$metric <- "Selectivity\n(among detected)"
    
    # Combine the datasets
    combined_fisher <- bind_rows(detection_fisher, selectivity_fisher)
    
    # For detection, we want to interpret odds ratio for comparison vs. baseline
    # For selectivity, we want to interpret odds ratio for baseline vs. comparison
    # We'll flip the detection odds ratios to maintain consistent interpretation
    combined_fisher <- combined_fisher %>%
        mutate(
            odds_ratio = ifelse(metric == "Detection", 1/odds_ratio, odds_ratio),
            ci_lower_adj = ifelse(metric == "Detection", 1/ci_upper, ci_lower),
            ci_upper_adj = ifelse(metric == "Detection", 1/ci_lower, ci_upper)
        )
    
    # Create the combined odds ratio plot
    combined_or_plot <- ggplot(combined_fisher, aes(x = group, y = odds_ratio, color = metric)) +
        geom_point(position = position_dodge(width = 0.5), size = 3) +
        geom_errorbar(
            aes(ymin = ci_lower_adj, ymax = ci_upper_adj),
            position = position_dodge(width = 0.5),
            width = 0.3
        ) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
        scale_y_continuous(trans = "log10") +
        scale_color_manual(values = c("Detection" = "#E41A1C", "Selectivity\n(among detected)" = "#4DAF4A")) +
        labs(
            title = paste("Odds Ratios for", baseline_category, "vs. Comparison"),
            subtitle = paste("Values > 1 indicate", baseline_category, "junctions are more likely to be detected/selective"),
            x = "Group",
            y = "Odds Ratio (log scale)",
            color = "Comparison"
        ) +
        theme_minimal() +
        theme(
            legend.position = "bottom",
            axis.text.x = element_text(angle = 45, hjust = 1)
        )
    
    return(combined_or_plot)
}
# 
# 
# # Original usage with NMD/non-NMD
# upf1_result_d <- analyze_junction_detection_and_selectivity(
#     summary_post, 
#     upf1_slope_plot_table$filtered_table,
#     category_column = "color_gene_name",
#     baseline_category = "non-NMD rescued",
#     color_choices = c("#D53E4F","#3288BD")
# )
# 

# 
# 
# combined_detection_table = rbind((smg1_d$detection_fisher_results %>% 
#                                       mutate(experiment = 'SMG1 inhibition')),
#                                  (chx_result_d$detection_fisher_results %>% 
#                                       mutate(experiment = 'CHX')),
#                                  (upf1_result_d$detection_fisher_results %>% 
#                                       mutate(experiment = 'UPF1 KD'))) %>% 
#     mutate(metric = "Detection")
# 
# 
# combined_selectivity_table = rbind((smg1_d$selectivity_fisher_results %>% 
#                                         mutate(experiment = 'SMG1 inhibition')),
#                                    (chx_result_d$selectivity_fisher_results %>% 
#                                         mutate(experiment = 'CHX')),
#                                    (upf1_result_d$selectivity_fisher_results %>% 
#                                         mutate(experiment = 'UPF1 KD'))) %>% 
#     mutate(metric = "Selectivity\n(among detected)")
# 
# 
# rbind(combined_detection_table, combined_selectivity_table) %>% 
#     mutate(
#         odds_ratio = ifelse(metric == "Detection", 1/odds_ratio, odds_ratio),
#         ci_lower_adj = ifelse(metric == "Detection", 1/ci_upper, ci_lower),
#         ci_upper_adj = ifelse(metric == "Detection", 1/ci_lower, ci_upper),
#         # Create a variable to determine fill based on p-value
#         significant = p_value < 0.005
#     ) %>% 
#     ggplot(aes(x = group, y = odds_ratio, color = metric)) + 
#     geom_point(
#         aes(fill = ifelse(significant, metric, NA)),  # Fill only if significant
#         position = position_dodge(width = 0.5), 
#         size = 3,
#         shape = 21,  # Use shape 21 which allows both color (outline) and fill
#         stroke = 1.5  # Make the outline a bit thicker
#     ) + 
#     geom_errorbar(
#         aes(ymin = ci_lower_adj, ymax = ci_upper_adj),
#         position = position_dodge(width = 0.5),
#         width = 0.3
#     ) + 
#     geom_hline(yintercept = 1, linetype = "dashed", color = "red") + 
#     scale_y_continuous(trans = "log10") + 
#     scale_color_manual(values = c("Detection" = "#E41A1C", "Selectivity\n(among detected)" = "#4DAF4A")) +
#     scale_fill_manual(
#         values = c("Detection" = "#E41A1C", "Selectivity\n(among detected)" = "#4DAF4A"),
#         na.value = "white",  # Unfilled (white) for non-significant points
#         guide = "none"  # Don't show fill legend since it's redundant with color
#     ) +
#     facet_wrap(~experiment) +
#     labs(
#         subtitle = paste("Values > 1 indicate", 'NMD-insensitive', "junctions are more likely to be detected/selective"),
#         x = "Group",
#         y = "Odds Ratio",
#         color = "Comparison"
#     ) +
#     theme_bw() + 
#     theme(
#         legend.position = "bottom",
#         axis.text.x = element_text(angle = 45, hjust = 1)) 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# another_table = master_table %>% filter(cryptic == 'yes' & new_cate %in% c("Early","Late")) %>% select(paste_into_igv_junction,new_cate)
# # New usage with early/late categories
# early_late_results <- analyze_junction_detection_and_selectivity(
#     summary_post, 
#     another_table,
#     category_column = "new_cate",
#     baseline_category = "Early",
#     color_choices = c("#de53b0", "#173dd3")
# )
# 
# chx_result_d$combined_rate_plot + ggtitle('CHX junctions')
# smg1_d$combined_rate_plot + ggtitle('SMG1 junctions')
# upf1_result_d$combined_rate_plot + ggtitle('UPF1 juncions')
# 
# early_late_results$combined_rate_plot + theme(axis.text.x = element_text(size = 15)) + ggtitle('',subtitle = '')
# 
# 
# early_late_results$selectivity_analysis %>% 
#     filter(detected) %>% 
#     left_join(master_table %>% select(paste_into_igv_junction,junc_cat,strand,gene)) %>% 
#     distinct(paste_into_igv_junction,junc_cat,new_cate,selective,strand,gene) %>% 
#     group_by(paste_into_igv_junction) %>% 
#     add_count(paste_into_igv_junction) %>% filter(n == 2)
# separate(paste_into_igv_junction,remove = FALSE,into = c("chr",'start','end'),
#          convert = TRUE) %>% 
#     makeGRangesFromDataFrame(,keep.extra.columns = TRUE)
# 
# sel_er_gr$longest_ug_region <- Biostrings::vcountPattern('TGTGTG', seq_tmp, max.mismatch = 1)
# 
# sel_er_gr %>% 
#     as.data.table() %>% 
#     mutate(ug_den = longest_ug_region / width) %>% 
#     ggplot(aes(x = selective, y = ug_den, fill = selective)) + 
#     geom_boxplot() +
#     facet_wrap(~junc_cat)
# 
# 
# chx_result_d$detection_summary %>% 
#     mutate(experiment = 'chx') %>% 
#     
#     