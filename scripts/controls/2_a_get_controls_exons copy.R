
# GET Exon CONTROLS

################################################################################

# LIBRARIES

library(tidyverse)
library(plyranges)
library(txdbmaker)
library(GenomicFeatures)
library(annotatr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(data.table)

################################################################################

# INPUT
source("scripts/splicejame_closestExonToJunctions.R")
master_table <- read.csv(here::here("data/master_table.csv"))
master_table <- master_table %>%
    mutate(tdp43_sensitivity_be2 = ifelse(tdp43_sensitivity_be2 == "","Not found",tdp43_sensitivity_be2)) %>% 
    mutate(tdp43_sensitivity_sh = ifelse(tdp43_sensitivity_sh == "","Not found",tdp43_sensitivity_sh)) %>% 
    mutate(new_cate = case_when(tdp43_sensitivity_be2 == "Not found" ~ tdp43_sensitivity_sh,
                                tdp43_sensitivity_sh == "Not found" ~ tdp43_sensitivity_be2,
                                tdp43_sensitivity_sh == tdp43_sensitivity_be2 ~ tdp43_sensitivity_sh,
                                tdp43_sensitivity_be2 != 'Late' & tdp43_sensitivity_sh == 'Early' ~ 'Early',
                                tdp43_sensitivity_be2 != 'Early' & tdp43_sensitivity_sh == 'Late' ~ 'Late',
                                tdp43_sensitivity_sh == 'Intermediate'~  'Intermediate')) %>% 
    mutate(new_cate = ifelse(is.na(new_cate),"Ambiguous",new_cate)) 

gtf = txdbmaker::makeTxDbFromGFF('/Users/annaleigh/Downloads/gencode.v42.primary_assembly.annotation.gtf')
gr_master = master_table %>% 
    filter(cryptic == 'yes') %>% 
    filter(junc_cat %in% c("donor",'acceptor')) %>% 
    separate(paste_into_igv_junction,convert = TRUE, into = c("chrom",'start','end'),remove = FALSE) %>% 
    group_by(gene, junc_cat) %>%
    mutate(
        name = paste0(gene, "_", junc_cat, "_", row_number())
    ) %>%
    makeGRangesFromDataFrame(,keep.extra.columns = TRUE)

## get the genes that overlap 
all_txp = unlist(transcriptsBy(gtf,'gene'))
overlapping_txp = findOverlaps(all_txp,gr_master,ignore.strand = FALSE)
overlapping_txp = all_txp[queryHits(overlapping_txp)]
overlapping_txp$gene_id = names(overlapping_txp)
names(overlapping_txp) = NULL
## get the exons that overlap 

all_exons = exonsBy(gtf,'tx',use.names = TRUE)
all_exons = unlist(all_exons)
all_exons$tx = names(all_exons)
names(all_exons) = NULL

## get the exons in transcripts which overlap
possible_control_exons = all_exons %>%
    filter(tx %in% overlapping_txp$tx_name) 

## assign genes and transcripts to cryptic junctions

annotated_genes_transcripts = annotate_regions(gr_master,annotations =overlapping_txp ) %>% 
    as.data.table() %>% 
    distinct(gene,paste_into_igv_junction,annot.tx_name,annot.gene_id)


# remove flanking exons of the cryptic itself -----------------------------

find_flank_gr = all_exons %>% 
    as.data.frame() %>% 
    distinct(seqnames,start,end,strand) %>% 
    mutate(name = glue::glue("{seqnames}:{start}-{end}:{strand}")) %>% 
    makeGRangesFromDataFrame(,keep.extra.columns = TRUE)
names(find_flank_gr) = find_flank_gr$name

assigned_closest_exon = splicejam_closestExonToJunctions(gr_master,exonsGR = find_flank_gr)$spliceGRgene 
assigned_closest_exon = assigned_closest_exon %>% 
    filter(distFrom == 0 | distTo == 0)

assigned_closest_exon = assigned_closest_exon %>% 
    as.data.table() %>% 
    distinct(gene,strand,junc_cat,nameTo,nameFrom,paste_into_igv_junction,distFrom,distTo) %>% 
    mutate(flank_exon = ifelse(distFrom == 0,nameFrom,nameTo)) %>% 
    separate(flank_exon,into = c("chrom",'start','end'),convert = TRUE) %>% 
    makeGRangesFromDataFrame(,keep.extra.columns = TRUE) %>% 
    unique()



## for each cryptic event, generate 5 control exons
all_control_exons = GRanges()
set.seed(69)
for(g in 1:length(gr_master)){
    print(g)
    this_event = gr_master[g]

    this_event_strand = as.data.frame(this_event) %>% pull(strand) %>% unique()
    this_events_junc_cate =  as.data.frame(this_event) %>% pull(junc_cat) %>% unique()
    this_events_transcripts = annotated_genes_transcripts %>% filter(paste_into_igv_junction == this_event$paste_into_igv_junction)
    
    these_possible_exons = all_exons %>% 
        filter(tx %in% this_events_transcripts$annot.tx_name)
    these_possible_exons$tx = NULL
    these_possible_exons$exon_id = NULL
    these_possible_exons$exon_name = NULL
    these_possible_exons$exon_rank = NULL
    these_possible_exons = these_possible_exons %>% unique()
    
    if(length(these_possible_exons) > 2){
        these_possible_exons = subsetByOverlaps(these_possible_exons,assigned_closest_exon,invert = TRUE)
    }else{
        print(this_event)

    }

    
    if(this_events_junc_cate == 'acceptor'){
            
            chosen_exons = resize(these_possible_exons,fix = 'start',width = 1)

            chosen_exons = resize(chosen_exons,fix = 'center',width = 500)

  
    }else{
        chosen_exons = resize(these_possible_exons,fix = 'end',width = 1)
        chosen_exons = resize(chosen_exons,fix = 'center',width = 500)
        

    }

    exon_ids = these_possible_exons %>% as.data.table() %>% mutate(exon_id = glue::glue("{seqnames}:{start}-{end}")) %>% pull(exon_id)
    
    chosen_exons$gene_name = this_event$gene
    chosen_exons$cryptic_name = this_event$name
    chosen_exons$junc_cat = this_events_junc_cate
    chosen_exons$exon_ids = exon_ids
    chosen_exons = chosen_exons %>% unique()
    
    chosen_exons = chosen_exons %>% 
        group_by(gene_name, junc_cat) %>%
        mutate(
            name = paste0(gene_name, "_", junc_cat, "_", row_number())
        ) %>%
        ungroup()
    #sample 5 random exons
    # if(length(chosen_exons) >=5){
    # 
    #     rows_to_sample = sample(1:length(chosen_exons),size = 5,replace = FALSE)
    # }else{
    # 
    #     rows_to_sample = sample(1:length(chosen_exons),size = length(chosen_exons),replace = FALSE)
    # }
    
    # chosen_exons = chosen_exons[rows_to_sample]
    
    all_control_exons = c(all_control_exons,chosen_exons)

}

sampled_controls = all_control_exons %>% 
    as.data.table() %>% 
    group_by(gene_name,junc_cat) %>% 
    slice_sample(n = 5,replace = FALSE) %>% 
    ungroup() %>% 
    add_count(exon_ids) %>% 
    makeGRangesFromDataFrame(,keep.extra.columns = TRUE)


all_control_exons_seq = getSeq(BSgenome.Hsapiens.UCSC.hg38,sampled_controls)
names(all_control_exons_seq) = sampled_controls$name
writeXStringSet(all_control_exons_seq, file = "~/cluster/vyplab/sbs_projects/SORRY_AL_MADE_THIS_WILL_DELETE/repeatmasker/control_all_controls_as_exons_420seed2.fa")
rtracklayer::export(sampled_controls,"~/cluster/vyplab/sbs_projects/SORRY_AL_MADE_THIS_WILL_DELETE/repeatmasker/control_all_controls_as_exons_420seed2.bed" )

meta_table = sampled_controls %>% 
    as.data.table() %>% 
    distinct(cryptic_name,junc_cat,name)

meta_table %>% 
    left_join(gr_master %>% as.data.frame() %>% distinct(name,paste_into_igv_junction,new_cate),by = c("cryptic_name" = 'name')) %>% 
    fwrite('data/repeatmasker/metadata_mapping_exon_controls.txt')
