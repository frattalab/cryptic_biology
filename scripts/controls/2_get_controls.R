
# GET INTRONIC CONTROLS

################################################################################

# LIBRARIES

library(dplyr)

################################################################################

# INPUT

dir_in <- "D:/Documenti/FlaminiaD/LABORATORIO FRATTA/PROJECTS/cryptics/input/bed_files/"

intronic_window <- 250
exonic_window <- 250
step <- 50

n_controls <- 10
n_control_sets <- 5

################################################################################

# FUNCTIONS

Dirs <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}



get_nonoverlapping <- function(df, windows_per_junc) {
  all_selected <- list()
  gene_list <- unique(df$gene)
  
  for (gene in gene_list) {
    subdf_gene <- df[df$gene == gene, ]
    # Extract all unique junctions for this gene
    all_juncs <- unique(unlist(strsplit(paste(subdf_gene$id_junc, collapse = ";"), ";")))
    
    for (junc in all_juncs) {
      # Pool: all windows for this gene containing this junction
      pool <- subdf_gene[grepl(junc, subdf_gene$id_junc), ]
      selected <- data.frame()
      
      while (nrow(pool) > 0 && nrow(selected) < windows_per_junc) {
        idx <- sample(seq_len(nrow(pool)), 1)
        picked <- pool[idx, , drop = FALSE]
        selected <- rbind(selected, picked)
        
        # Remove all overlapping windows from pool
        overlap_idx <- which(
          (pool$start <= picked$end & pool$end >= picked$start)
        )
        pool <- pool[-overlap_idx, ]
      }
      
      if (nrow(selected) < windows_per_junc) {
        warning(sprintf("Gene %s, Junction %s: Only %d non-overlapping windows found.", gene, junc, nrow(selected)))
      }
      
      if (nrow(selected) > 0) {
        selected$junction <- junc
        all_selected[[paste0(gene, "_", junc)]] <- selected
        # Remove selected windows from subdf_gene (so they are not reused for next junction)
        subdf_gene <- subdf_gene[!subdf_gene$name %in% selected$name, ]
      }
    }
  }
  do.call(rbind, all_selected)
}



################################################################################

# CODE

############################
# 0. MAKE OUTPUT DIRECTORY #
############################

parameters <- paste0("intw_", intronic_window, "_exw_", exonic_window, "_step_", step, "_nctrl_", n_controls)
dir_out <- paste0(dir_in, parameters, "/")
Dirs(dir_out)


###################
# 1. IMPORT FILES #
###################

cryptic_bed <- read.table(paste0(dir_in, "cryptic_junctions.bed"), sep = "\t", header = F)
colnames(cryptic_bed) <- c("chr", "start", "end", "name", "score", "strand")

introns_bed <- read.table(paste0(dir_in, "matched_introns.bed"), sep = "\t", header = F)
colnames(introns_bed) <- c("chr", "start", "end", "name", "score", "strand")



####################################
# 2. GET CRYPTIC JUNCTIONS WINDOWS #
####################################

cryptic_bed$junc_type <- sapply(strsplit(as.character(cryptic_bed$name), "_"), `[`, 2)

cryptic_bed$splice_site <- with(cryptic_bed,
                                ifelse(junc_type == "acceptor" & strand == "+", end,
                                       ifelse(junc_type == "donor" & strand == "+", start,
                                              ifelse(junc_type == "acceptor" & strand == "-", start,
                                                     ifelse(junc_type == "donor" & strand == "-", end, NA)))))


cryptic_bed <- cryptic_bed %>% filter(junc_type != "exon")

cryptic_bed <- cryptic_bed %>%
  mutate(
    junc_start = case_when(
      junc_type == "acceptor" & strand == "+" ~ splice_site - intronic_window,
      junc_type == "donor"    & strand == "+" ~ splice_site - exonic_window,
      junc_type == "acceptor" & strand == "-" ~ splice_site - exonic_window,
      junc_type == "donor"    & strand == "-" ~ splice_site - intronic_window,
      TRUE ~ NA_real_
    ),
    junc_end = case_when(
      junc_type == "acceptor" & strand == "+" ~ splice_site + exonic_window,
      junc_type == "donor"    & strand == "+" ~ splice_site + intronic_window,
      junc_type == "acceptor" & strand == "-" ~ splice_site + intronic_window,
      junc_type == "donor"    & strand == "-" ~ splice_site + exonic_window,
      TRUE ~ NA_real_
    )
  )


cryptic_bed2 <- cryptic_bed[,c("chr", "junc_start", "junc_end", "name", "score", "strand")]

write.table(
  cryptic_bed2,
  file = paste0(dir_out, "junction_windows.bed"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)



#################################
# 3. SPLIT INTRONS INTO WINDOWS #
#################################

window_size <- intronic_window + exonic_window

windows_list <- lapply(1:nrow(introns_bed), function(i) {
  row <- introns_bed[i, ]
  max_start <- row$end - window_size + 1
  if (max_start >= row$start) {
    starts <- seq(row$start, max_start, by = step)
    ends <- starts + window_size - 1
    n <- length(starts)
    name_with_num <- paste0(row$name, "_", seq_len(n))
    data.frame(
      chr = row$chr,
      start = starts,
      end = ends,
      name = name_with_num,
      score = row$score,
      strand = row$strand
    )
  } else {
    NULL
  }
})

sliding_windows_bed <- do.call(rbind, windows_list)

write.table(
  sliding_windows_bed,
  file = paste0(dir_out, "introns_sliding_windows.bed"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


#################################
# 4. RANDOMLY PICK CONTROL SETS #
#################################

cryptic_bed2$gene <- sapply(strsplit(as.character(cryptic_bed$name), "_"), `[`, 1)
n_cryp <- cryptic_bed2[,c("gene", "name")]
n_cryp <- n_cryp %>%
  group_by(gene) %>%
  summarize(
    n_junc = n(),
    id_junc = paste(name, collapse = ";"))



sliding_windows_bed$gene <- sub("_.*", "", sliding_windows_bed$name)
sliding_windows_bed <- merge(sliding_windows_bed, n_cryp, by="gene")


for(i in 1:n_control_sets){
  
  print(paste0("set ", i))
  
  set <- get_nonoverlapping(sliding_windows_bed, n_controls)
  warnings()
  rownames(set) <- NULL
  
  set_info <- set %>% group_by(gene, junction) %>% summarise(n_ctrl = n())
  set_info <- merge(set[,c("junction", "name")], set_info, by="junction")
  set_info$set_n <- i
  
  set$gene <- NULL
  set <- set[,c("chr", "start", "end", "name", "score", "strand")]
  set$start <- as.integer(set$start)
  set$end <- as.integer(set$end)
  
  write.table(
    set,
    file = paste0(dir_out, "set", i,".bed"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  if(i==1){
    set_info_total <- set_info
  } else {
    set_info_total <- rbind(set_info_total, set_info)
  }
  
}








