
# CREATE BED INTRONS

################################################################################

# LIBRARIES

library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(dplyr)
library(tidyr)

################################################################################

# INPUT

dir_cryptics <- "D:/Documenti/FlaminiaD/LABORATORIO FRATTA/PROJECTS/cryptics/input/master_table.csv"
dir_gtf <- "D:/Documenti/FlaminiaD/LABORATORIO FRATTA/Roba Utile/GRCh38.113.gtf"

dir_out <- "D:/Documenti/FlaminiaD/LABORATORIO FRATTA/PROJECTS/cryptics/input/bed_files/"

################################################################################

# CODE


###################
# 1. IMPORT FILES #
###################

# 1.1 IMPORT GTF AND GET INTRONS COORDINATES

gtf <- import(dir_gtf)

txdb <- makeTxDbFromGFF(dir_gtf, format = "gtf")
introns_by_tx <- intronsByTranscript(txdb, use.names = TRUE)
introns <- unlist(introns_by_tx, use.names = FALSE)

rm(txdb, introns_by_tx)



# 1.2 IMPORT CRYPTIC TABLE
cryptic <- read.csv(dir_cryptics)
cryptic <- cryptic %>% filter(cryptic == "yes" & junc_cat != "exon_skip")

list_genes <- unique(cryptic$gene)

cryptic <- cryptic %>%
  group_by(gene, junc_cat) %>%
  mutate(
    name = paste0(gene, "_", junc_cat, "_", row_number())
  ) %>%
  ungroup()

write.table(cryptic[,c("name", "strand", "new_cate")], "D:/Documenti/FlaminiaD/LABORATORIO FRATTA/PROJECTS/cryptics/input/metadata.txt", col.names = T, row.names = F, sep = "\t", quote = F)


################################################
# 2. CREATE BED FILES FROM CRYPTIC ANNOTATIONS #
################################################

# 2.1 CRYPTIC JUNCTIONS BED

bed_cryptic <- cryptic[,c(2,4,46)]

bed_cryptic <- bed_cryptic %>%
  separate(paste_into_igv_junction, into = c("seqnames", "positions"), sep = ":") %>%
  mutate(seqnames = sub("^chr", "", seqnames)) %>%
  separate(positions, into = c("start", "end"), sep = "-") %>%
  mutate(
    start = as.integer(start),
    end = as.integer(end),
    score = "."
  ) %>%
  dplyr::select(seqnames, start, end, name, score, strand)



# 2.2 CRYPTIC INTRONS BED

bed_junctions <- distinct(cryptic[,c(1,3,4)])

bed_junctions <- bed_junctions %>%
  group_by(gene) %>%
  mutate(
    n = if(n() > 1) row_number() else NA_integer_,
    name = ifelse(!is.na(n), paste0(gene, "_", n), gene)
  ) %>%
  ungroup() %>%
  dplyr::select(-n, -gene)

bed_junctions <- bed_junctions %>%
  separate(paste_into_igv_junction_annotated, into = c("seqnames", "positions"), sep = ":") %>%
  mutate(seqnames = sub("^chr", "", seqnames)) %>%
  separate(positions, into = c("start", "end"), sep = "-") %>%
  mutate(
    start = as.integer(start),
    end = as.integer(end),
    score = "."
  ) %>%
  dplyr::select(seqnames, start, end, name, score, strand)

bed_junctions_gr <- GRanges(
  seqnames = bed_junctions$seqnames,
  ranges = IRanges(start = bed_junctions$start, end = bed_junctions$end),
  strand = bed_junctions$strand
)




######################
# 3. FILTER GTF FILE #
######################

# 3.1 KEEP ONLY THE INTRONS NON OVERLAPPING WITH EXONS

exons   <- gtf[mcols(gtf)$type == "exon"]

hits <- findOverlaps(introns, exons, ignore.strand = TRUE)

introns_with_overlap <- unique(queryHits(hits))
introns_idx <- seq_along(introns)
gtf2 <- introns[!(introns_idx %in% introns_with_overlap)]



# 3.2 KEEP GENES OF INTEREST ONLY

gtf_genes <- gtf[mcols(gtf)$gene_name %in% list_genes]
gtf_genes   <- gtf_genes[mcols(gtf_genes)$type == "gene"]

hits <- findOverlaps(gtf2, gtf_genes, type = "within", ignore.strand = FALSE)
gtf2_genes <- gtf2[unique(queryHits(hits))]

gtf3 <- gtf2[mcols(gtf2)$gene_name %in% list_genes]



# 3.3 REMOVE INTRONS INVOLVED IN CRYPTIC SPLICING

hits <- findOverlaps(gtf2_genes, bed_junctions_gr, ignore.strand = FALSE)
introns_gr <- gtf2_genes[-unique(queryHits(hits))]



# 3.4 ASSIGN GENE IDENTITY TO INTRONS

hits <- findOverlaps(introns_gr, gtf_genes, ignore.strand = FALSE)

gene_names <- rep(NA_character_, length(introns_gr))

hit_df <- data.frame(
  intron_idx = queryHits(hits),
  gene_idx = subjectHits(hits),
  gene_name = mcols(gtf_genes)$gene_name[subjectHits(hits)]
)

gene_map <- hit_df %>% 
  group_by(intron_idx) %>%
  summarise(gene_name = paste(unique(gene_name), collapse = ","))

gene_names[gene_map$intron_idx] <- gene_map$gene_name

mcols(introns_gr)$gene_name <- gene_names

introns_bed <- data.frame(
  chrom = as.character(seqnames(introns_gr)),
  chromStart = start(introns_gr) - 1,
  chromEnd = end(introns_gr),
  name = if ("gene_name" %in% colnames(mcols(introns_gr))) mcols(introns_gr)$gene_name else ".",
  score = ".",
  strand = as.character(strand(introns_gr))
)

introns_bed <- distinct(introns_bed)

introns_bed$name <- with(introns_bed, paste0(name, "_", ave(seq_along(name), name, FUN=seq_along)))




#################
# 4. SAVE FILES #
#################

write.table(
  bed_cryptic,
  file = paste0(dir_out, "cryptic_junctions.bed"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


write.table(
  introns_bed,
  file = paste0(dir_out, "matched_introns.bed"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)








