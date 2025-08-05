library(ggplot2)
library(ggsignif)
library(dplyr)
library(tibble)
library(dunn.test)

plot_significant_comparisons <- function(data, column_name, alpha = 0.05,ylab_chosen = 'Age at death') {
    # Ensure the data has the required columns
    if (!("dataset" %in% colnames(data))) {
        stop("Data must contain a 'dataset' column")
    }
    
    if (!(column_name %in% colnames(data))) {
        stop(paste("Column", column_name, "not found in data"))
    }
    
    # Remove rows with missing values for the specified column
    data_clean <- data %>%
        filter(!is.na(.data[[column_name]]))
    
    # Get unique datasets
    datasets <- unique(data_clean$dataset)
    
    if (length(datasets) < 2) {
        stop("Need at least 2 datasets for comparisons")
    }
    
    # First perform Kruskal-Wallis test
    kw_test <- kruskal.test(data_clean[[column_name]], data_clean$dataset)
    
    cat("Kruskal-Wallis Test Results:\n")
    cat("Chi-squared =", round(kw_test$statistic, 4), "\n")
    cat("df =", kw_test$parameter, "\n")
    cat("p-value =", format.pval(kw_test$p.value), "\n\n")
    
    comparisons <- list()
    comparison_results <- NULL
    
    # Only proceed with post-hoc tests if Kruskal-Wallis is significant
    if (kw_test$p.value < alpha) {
        cat("Kruskal-Wallis test is significant. Proceeding with Dunn's test for post-hoc comparisons.\n\n")
        
        # Perform Dunn's test for post-hoc pairwise comparisons
        dunn_result <- dunn.test(data_clean[[column_name]], 
                                 data_clean$dataset, 
                                 method = "BH",  # Benjamini-Hochberg correction
                                 alpha = alpha)
        
        # Extract comparison results
        comparison_results <- data.frame(
            comparison = dunn_result$comparisons,
            Z = dunn_result$Z,
            p.value = dunn_result$P,
            p.adjusted = dunn_result$P.adjusted,
            significant = dunn_result$P.adjusted < alpha,
            stringsAsFactors = FALSE
        )
        
        # Create comparisons list for plotting (only significant ones)
        for (i in 1:nrow(comparison_results)) {
            if (comparison_results$significant[i]) {
                groups <- strsplit(comparison_results$comparison[i], " - ")[[1]]
                comparisons[[length(comparisons) + 1]] <- c(groups[1], groups[2])
            }
        }
        
    } else {
        cat("Kruskal-Wallis test is not significant (p =", format.pval(kw_test$p.value), ").\n")
        cat("No post-hoc comparisons will be performed.\n\n")
    }
    
    # Create the base plot
    p <- data_clean %>%
        rowid_to_column() %>%
        ggplot(aes(x = dataset, y = .data[[column_name]])) +
        geom_boxplot() +
        ylab(ylab_chosen) +
        xlab("") +
        theme_minimal()
    
    # Add significance annotations using geom_signif if there are significant comparisons
    if (length(comparisons) > 0) {
        # Prepare significance labels from Dunn's test results
        sig_labels <- character()
        sig_comparisons <- list()
        
        for (i in 1:nrow(comparison_results)) {
            if (comparison_results$significant[i]) {
                groups <- strsplit(comparison_results$comparison[i], " - ")[[1]]
                sig_comparisons[[length(sig_comparisons) + 1]] <- groups
                
                # Create significance label based on adjusted p-value
                p_adj <- comparison_results$p.adjusted[i]
                if (p_adj < 0.001) {
                    sig_labels <- c(sig_labels, "***")
                } else if (p_adj < 0.01) {
                    sig_labels <- c(sig_labels, "**")
                } else if (p_adj < 0.05) {
                    sig_labels <- c(sig_labels, "*")
                } else {
                    sig_labels <- c(sig_labels, "ns")
                }
            }
        }
        
        # Add geom_signif with the results
        p <- p + geom_signif(
            comparisons = sig_comparisons,
            annotations = sig_labels,
            step_increase = 0.1,
            tip_length = 0.01,
            vjust = 0.5
        )
    }
    
    # Add overall Kruskal-Wallis test result
    kw_label <- paste0("Kruskal-Wallis, p = ", format.pval(kw_test$p.value))
    p <- p + labs(subtitle = kw_label)
    
    # Print summary of post-hoc comparisons if performed
    if (!is.null(comparison_results)) {
        cat("Dunn's Test Results (Benjamini-Hochberg corrected):\n")
        print(comparison_results)
        
        cat("\nNumber of significant pairwise comparisons:", sum(comparison_results$significant), "\n")
    }
    
    return(p)
}
