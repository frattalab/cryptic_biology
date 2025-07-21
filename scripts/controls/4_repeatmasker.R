
# REPEATMASKER
# To run this script you need to activate the repeatmasker environment

################################################################################

# INPUT

species <- "human"

dir_fasta <- "/SAN/vyplab/sbs_projects/flaminia/projects/cryptics/input/fasta_files/"
dir_out <- "/SAN/vyplab/sbs_projects/flaminia/projects/cryptics/output/RepeatMasker/"

python3_bin <- "/SAN/neuroscience/cryptic_circ/tools/miniconda3/envs/analysis/bin/python3.8"
repeatmasker_bin <- "/SAN/vyplab/sbs_projects/flaminia/tools/RepeatMasker/RepeatMasker"
RM2Bed_bin <- "/SAN/vyplab/sbs_projects/flaminia/tools/RepeatMasker/util/RM2Bed.py"

################################################################################

# CODE

files <- list.files(dir_fasta)

for(i in 1:length(files)){
  
  fasta_file <- files[i]
  print(paste0(fasta_file))
  
  print("Running RepeatMasker...")
  command1 <- paste0(repeatmasker_bin, " -species ", species," -dir ",  dir_out, " ", dir_fasta, fasta_file)
  print(command1)
  system(command1)
  
  print("Done!")
  
  print("Converting RM to BED...")
  basename <- fasta_file
  rm_file <- paste0(dir_out, basename, ".out")
  command2 <- paste0(python3_bin, " ", RM2Bed_bin, " -d ", dir_out, " ", rm_file)
  print(command2)
  system(command2)
  
  print("Done!")
  
}
