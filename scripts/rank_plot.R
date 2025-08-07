# Load required libraries
library(ggplot2)
library(ggrepel)  # for better label positioning


# Create the rank plot
create_rank_plot <- function(data, expression_col, gene_name_col, cord_only_col) {
    
    # Sort data by expression value (descending order for typical rank plots)
    data_sorted <- data[order(data[[expression_col]], decreasing = TRUE), ]
    
    # Add rank column
    data_sorted$rank <- 1:nrow(data_sorted)
    
    # Create a logical vector for labeling (TRUE where cord-only has a value)
    # Adjust this condition based on what "value in column 'cord-only'" means
    # This assumes you want to label genes where cord-only is not NA or not empty
    data_sorted$label_gene <- !is.na(data_sorted[[cord_only_col]]) & 
        data_sorted[[cord_only_col]] != "" &
        data_sorted[[cord_only_col]] != 0  # adjust condition as needed
    
    # Create the plot
    p <- ggplot(data_sorted, aes(x = rank, y = .data[[expression_col]])) +
        geom_point(alpha = 0.6, size = 1.5) +
        geom_text_repel(
            data = data_sorted[data_sorted$label_gene, ],
            aes(label = .data[[gene_name_col]]),
            size = 3,
            max.overlaps = 20,  # adjust as needed
            box.padding = 0.3,
            point.padding = 0.3
        ) +
        labs(
            x = "Rank",
            y = "Log10 Gene Expression"
        ) +
        theme_minimal() +
        theme(
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5)
        ) +
        theme(text = element_text(size = 18))
    
    return(p)
}

