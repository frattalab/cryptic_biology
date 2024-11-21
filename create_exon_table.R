source("scripts/splicejame_closestExonToJunctions.R")
build_all_potential_exon_table = function(grange_object,exon_annotation){
    
    ##check that exon_coords exists as a metadata column
    
    
    
    
    if(!("exon_coords" %in% (colnames(mcols(exon_annotation))))){
        exon_annotation$exon_coords = paste0(as.character(seqnames(exon_annotation)),
                                             ":",
                                             as.character(start(exon_annotation)),
                                             "-",
                                             as.character(end(exon_annotation)),
                                             ":",
                                             as.character(strand(exon_annotation)))
    }
    ##check that the exons have names, splice jam needs it
    if(is.null(names(exon_annotation))){
        names(exon_annotation) = exon_annotation$exon_coords
    }
    ####build a table of all potential exons for all junctions####
    splicejamobject = splicejam_closestExonToJunctions(
        grange_object,
        exon_annotation,
        flipNegativeStrand = TRUE,
        sampleColname = "sample_id",
    )
    
    ####build a table of all potential exons for all junctions####
    spliceStartExonEndD = splicejamobject$spliceStartExonEndD
    spliceEndExonStartD = splicejamobject$spliceEndExonStartD
    
    #retrieve the junctions that we queried, the splice start
    inputSpliceStarts = as.data.table(grange_object[spliceStartExonEndD$queryHits,]) %>%
        dplyr::select(seqnames,
                      start,
                      end,
                      strand)
    
    # get all the exons taht end on those splice starts
    inputExonsEnds = as.data.table(exon_annotation[spliceStartExonEndD$subjectHits,]) %>%
        dplyr::select(exon_coords)
    
    allPossibleSourceExons = cbind(inputSpliceStarts,inputExonsEnds,spliceStartExonEndD$strandedDistance) %>%
        filter(V3 == 0) %>%
        dplyr::select(-V3) %>%
        unique()
    
    ####now for the TargetExons
    ####all possible target exons ####
    inputSpliceEnds = as.data.table(grange_object[spliceEndExonStartD$queryHits,]) %>%
        dplyr::select(seqnames,
                      start,
                      end,
                      strand)
    
    inputExonsStarts = as.data.table(exon_annotation[spliceEndExonStartD$subjectHits,]) %>%
        dplyr::select(exon_coords)
    
    allPossibleTargetExons = cbind(inputSpliceEnds,inputExonsStarts,spliceEndExonStartD$strandedDistance) %>%
        filter(V3 == 0) %>% #only take distance zero
        dplyr::select(-V3) %>%
        unique()
    
    ###merge together the source and target
    allPossibleJunctionExon = allPossibleSourceExons %>%
        rbind(allPossibleTargetExons) %>%
        unique() %>% 
        as.data.table()
    
    allPossibleJunctionExon[,paste_into_igv_junction := paste0(seqnames, ":",start, "-",end)]
    
    return(allPossibleJunctionExon)
}
parse_coords <- function(coord_str) {
    parts <- str_match(coord_str, "(chr\\w+):(\\d+)-(\\d+)")
    list(chrom = parts[, 2], start = as.numeric(parts[, 3]), end = as.numeric(parts[, 4]))
}

merge_junctions <- function(df, distance_threshold = 2000) {
    df <- df %>% arrange(paste_into_igv_junction)
    
    result <- list()
    i <- 1
    while (i <= nrow(df)) {
        current <- parse_coords(df$paste_into_igv_junction[i])
        
        if (i < nrow(df)) {
            next_junction <- parse_coords(df$paste_into_igv_junction[i + 1])
            
            # Check if chromosomes are the same and distance is less than threshold
            if (current$chrom == next_junction$chrom &&
                (next_junction$start - current$end) < distance_threshold &&
                (next_junction$start - current$end) > 0) {
                
                # Merge the junctions
                merged_coord <- sprintf("%s:%d-%d", current$chrom, current$end, next_junction$start)
                result[[length(result) + 1]] <- c(
                    list(exon_coord = merged_coord),
                    new_cate = df$new_cate[i],
                    cryptic = df$cryptic[i],
                    gene = df$gene[i],
                    paste_into_igv_junction = df$paste_into_igv_junction[i]
                )
                
                i <- i + 2  # Skip the next row as it's been merged
                next
            }
        }
        
        # If not merged, don't add to the result
        i <- i + 1
    }
    
    # Convert result to data frame
    result_df <- do.call(rbind, lapply(result, data.frame))
    
    # Remove rows with same start or end coordinates
    result_df <- result_df %>%
        mutate(
            start = sapply(paste_into_igv_junction, function(x) parse_coords(x)$start),
            end = sapply(paste_into_igv_junction, function(x) parse_coords(x)$end)
        ) %>%
        group_by(start, end) %>%
        filter(n() == 1) %>%
        ungroup()
}

exons_combine <- merge_junctions(master_table) %>%
    mutate(
        symmetric = (end-start+1) %% 3==0
    )

exons1 <- rtracklayer::import("data/noDox-dox0075.cryptic_exons.bed") %>% as.data.table()
exons1 <- exons1 %>%
    tidyr::separate(name, into = c("gene", "paste_into_igv_junction"), sep = "\\|", remove = FALSE) %>% 
    mutate(exon_coord = paste0(seqnames,":",start,"-",end))


exons1 <- merge(exons1,
                master_table %>%
                    select(paste_into_igv_junction, new_cate, cryptic),
                by = "paste_into_igv_junction", all.x=T) %>%
    distinct(paste_into_igv_junction, .keep_all = T) %>%
    mutate(
        symmetric = (end-start) %% 3==0
    )

exons2 <- rtracklayer::import("data/all_the_ce_length_manual_filtered.bed")

master_table_gr <- master_table %>%
    separate(paste_into_igv_junction,remove = FALSE,into = c("seqnames",'start','end')) %>% 
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)

exons2_gr = exons2 %>% 
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)

exons2_gr$exon_coords = paste0(as.character(seqnames(exons2_gr)),
                               ":",
                               as.character(start(exons2_gr)),
                               "-",
                               as.character(end(exons2_gr)),
                               ":",
                               as.character(strand(exons2_gr)))

exons2_junced = build_all_potential_exon_table(master_table_gr,exons2_gr)
exons2_junced = exons2_junced %>% 
    separate(exon_coords,into = c("chr","start",'end','strand'),convert = TRUE,remove = FALSE) %>% 
    mutate(symmetric = (end-start) %% 3==0) %>% 
    dplyr::rename(exon_coord = exon_coords)


exons2_junced = exons2_junced %>% 
    left_join(master_table %>%
                  select(gene, paste_into_igv_junction, new_cate, cryptic,last_exon)) %>% 
    select(gene, exon_coord,paste_into_igv_junction, new_cate, cryptic, symmetric) %>% 
    arrange(gene)
# At this point we manually verfied exons bounds ---------------------
# In rare circumstances, e.g. lowly expressed last exons like ZNF681 the method is imperfect ---------------------
# Therefore the exon bounds were determined manually by inspected coverage in RNA-seq with IGV

exons_combine = fread("data/check.csv")
result_mrg <- rbind(
    exons_combine %>% select(gene, exon_coord, paste_into_igv_junction, new_cate, cryptic),
    exons1 %>% select(gene,exon_coord, paste_into_igv_junction, new_cate, cryptic),
    exons2_junced %>% select(gene, exon_coord,paste_into_igv_junction, new_cate, cryptic)) %>%
    na.omit() %>%
    mutate(exon_coord = str_remove(exon_coord, ":[+-]")) %>% 
    left_join(master_table %>% select(paste_into_igv_junction,last_exon,strand)) %>% 
    distinct() 

fwrite(result_mrg,'data/result_exon_table_merged.csv')

