
# BEDTOOLS GETFASTA

################################################################################

# INPUT

bedtools_bin <- "/SAN/vyplab/alb_projects/tools/bedtools"

dir_in <- "/SAN/vyplab/sbs_projects/flaminia/projects/cryptics/input/bed_files/"
dir_out <- "/SAN/vyplab/sbs_projects/flaminia/projects/cryptics/input/fasta_files/"

input_fasta <- "/SAN/vyplab/sbs_projects/flaminia/data_file/index/STAR/GRCh38.fa"

################################################################################

# CODE

files <- list.files(dir_in)

for(i in 1:length(files)){
  
  line <- gsub(".bed", "", files[i])
  print(line)
  
  bed_in <- paste0(dir_in, files[i])
  fasta_out <- paste0(dir_out, line, ".fa")
  
  system(paste0(bedtools_bin, " getfasta -fo ", fasta_out, " -s -name -fi ", input_fasta, " -bed ", bed_in))
  
}
