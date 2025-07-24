create_coverage_vectors <- function(range_a, rbp_grange) {
    
    # Use IRanges for maximum efficiency
    coverage_list <- lapply(seq_along(range_a), function(i) {
        current_range <- range_a[i]
        
        #findOverlaps
        hits <- findOverlaps(current_range, rbp_grange)
        
        if (length(hits) == 0) {
            return(rep(0L, width(current_range)))
        }
        
        # Get overlapping ranges and restrict to current range
        overlapping <- rbp_grange[subjectHits(hits)]
        intersected <- pintersect(rep(current_range, length(overlapping)), overlapping,strict.strand = TRUE)
        
        # Convert to relative coordinates (1-based)
        rel_ranges <- IRanges(
            start = start(intersected) - start(current_range) + 1L,
            end = end(intersected) - start(current_range) + 1L
        )
        
        # Use coverage for ultra-fast binary vector creation
        cov <- coverage(rel_ranges, width = width(current_range))
        return(as.integer(cov > 0))
    })
    
    return(coverage_list)
}

# coverage_numeric_tdp_flank <- create_coverage_vectors(flank_region, tdp_3)
# names(coverage_numeric_tdp_flank) = flank_region$name
# coverage_numeric_not_flank <- create_coverage_vectors(flank_region, not_tdp)
# names(coverage_numeric_not_flank) = flank_region$name
# 
# coverage_matrix_tdp<- do.call(rbind, coverage_numeric_tdp_flank) %>% 
#     as.tibble() %>% 
#     mutate(name = flank_region$name,
#            ID_exon = flank_region$ID_exon) %>% 
#     dplyr::relocate(name,ID_exon)
# coverage_matrix_not_tdp<- do.call(rbind, coverage_numeric_tdp_flank) %>% 
#     as.tibble() %>% 
#     mutate(name = flank_region$name,
#            ID_exon = flank_region$ID_exon) %>% 
#     dplyr::relocate(name,ID_exon)
# coverage_matrix_tdp %>% 
#     as.data.table() %>% 
#     rowwise() %>% 
#     mutate(sumVar = sum(c_across(V1:V501))) %>% 
#     filter(sumVar > 0) %>% 
#     select(-sumVar) %>% 
#     filter(grepl("donor",name)) %>% 
#     column_to_rownames('name') %>% 
#     pheatmap::pheatmap(,cluster_cols = FALSE)
# 
# coverage_matrix_tdp %>% 
#     as.data.table() %>% 
#     rowwise() %>% 
#     mutate(sumVar = sum(c_across(V1:V501))) %>% 
#     filter(sumVar > 0) %>% 
#     select(-sumVar) %>% 
#     filter(grepl("acceptor",name)) %>% 
#     column_to_rownames('name') %>% 
#     pheatmap::pheatmap(,cluster_cols = FALSE,cutree_rows = 6)
# # coverage plot in region -------------------------------------------------
# 
# library(zoo)
# coverage_matrix_tdp %>% 
#     mutate(cryptic = ifelse(grepl("_yes",name),"cryptic","alternative")) %>% 
#     mutate(type = ifelse(grepl("acceptor",name),"acceptor","donor")) %>% 
#     mutate(timing = gsub(".*_","",name)) %>% 
#     reshape2::melt(id.vars = c("name","cryptic","timing","type")) %>% 
#     mutate(pos = parse_number(as.character(variable))) %>% 
#     left_join(data.frame(name = flank_region$name, strand = strand(flank_region))) %>% 
#     mutate(pos = ifelse(strand == "-", max(pos) - pos + 1, pos)) %>% 
#     group_by(pos,timing,type,cryptic) %>% 
#     summarise(cov_val = mean(value)) %>% 
#     mutate(pos = pos - (flank_size)) %>% 
#     ungroup() %>% 
#     group_by(timing,type,cryptic) %>% 
#     filter(timing %in% c("early","late")) %>% 
#     filter(cryptic == 'cryptic') %>% 
#     mutate(
#         smoothed_weight = rollmean(
#             cov_val,
#             k = 5,
#             fill = NA,
#             align = "center"
#         )) %>% 
#     ungroup() %>% 
#     # mutate(type = fct_relevel(type,"donor")) %>% 
#     ggplot(aes(x = pos, y = smoothed_weight,color = timing)) + 
#     geom_line(size = 1.6) + 
#     facet_wrap(~type) + 
#     theme_bw() +
#     scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
#     ylab("Fraction covered TDP-43 Postar3") + 
#     geom_vline(xintercept = 0,linetype = 'dotted') +
#     xlab("Distance from cryptic splice site (bp)") +
#     scale_color_manual(values = c("#de53b0" ,   "#173dd3")) + 
#     theme(text = element_text(size = 18))