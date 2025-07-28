# FUNCTION - author Flaminia Pellegrini

rowwise_fisher_test <- function(
        
    df, 
    cryptic_col, 
    set1_col, 
    n_cry_col, 
    n_ctrl_col,
    pval_out = "pvalue_fisher",
    or_out = "oddsratio_fisher"
) {
    
    
    pvalues <- numeric(nrow(df))
    odds_ratio <- numeric(nrow(df))
    
    for (i in 1:nrow(df)) {
        a <- df[[cryptic_col]][i]
        b <- df[[set1_col]][i]
        c <- df[[n_cry_col]][i] - a
        d <- df[[n_ctrl_col]][i] - b
        
        
        mat <- matrix(c(a, b, c, d), nrow=2, byrow=FALSE)
        ft <- fisher.test(mat)
        pvalues[i] <- ft$p.value
        odds_ratio[i] <- ft$estimate
    }
    
    df[[pval_out]] <- pvalues
    df[[or_out]] <- odds_ratio
    
    return(df)
    
}



rep_ranges <- function(df, colname) {
    
    output_df <- data.frame(position = 0:499)
    
    if(nrow(df) == 0){
        
        output_df[[colname]] <- 0
        
    } else {
        
        unique_ids <- unique(df$id)
        
        for (id in unique_ids) {
            output_df[[id]] <- 0
        }
        
        
        for (i in 1:nrow(df)) {
            id <- df$id[i]
            start <- as.numeric(df$rep_start[i])
            end <- as.numeric(df$rep_end[i])
            output_df[[id]][output_df$position >= start & output_df$position <= end] <- 1
        }
        
        if(ncol(output_df)>2){
            output_df[[colname]] <- rowSums(output_df[ , !(names(output_df) %in% "position") ])
        }else{
            output_df[[colname]] <- output_df[ ,2]
        }
        
        output_df <- output_df[,c("position", colname)]
        
        
    }
    
    return(output_df)
    
}


repeats_heatmap <- function(
        df,
        title = "heatmap",
        col_fun = circlize::colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
        column_label_interval = 50,
        border = TRUE,
        cluster_rows = FALSE,
        cluster_columns = FALSE
) {
    # Pivot data to wide format
    acceptor_wide <- df %>%
        select(rows, position, logOR, repeat_cat) %>%
        pivot_wider(
            id_cols = c(rows, repeat_cat),
            names_from = position,
            values_from = logOR
        )
    
    # Prepare matrix for heatmap
    rownames_mat <- acceptor_wide$rows
    row_split <- acceptor_wide$repeat_cat
    heatmap_mat <- as.matrix(acceptor_wide[ , !(names(acceptor_wide) %in% c("rows", "repeat_cat")) ])
    rownames(heatmap_mat) <- rownames_mat
    
    # Custom column labels: label every interval, blanks for others
    all_colnames <- colnames(heatmap_mat)
    column_labels <- ifelse((seq_along(all_colnames) - 1) %% column_label_interval == 0, all_colnames, "")
    
    # Draw heatmap
    p <- Heatmap(
        heatmap_mat,
        name = "logOR",
        column_title = title,
        show_row_names = TRUE,
        row_names_gp = grid::gpar(fontsize = 10),
        show_column_names = TRUE,
        column_labels = column_labels,
        row_split = row_split,
        row_title = NULL,
        border = border,
        cluster_rows = cluster_rows,
        cluster_columns = cluster_columns,
        col = col_fun
    )
    
    return(p)
}

