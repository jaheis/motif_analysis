# Install and load rentrez
if (!requireNamespace("rentrez", quietly = TRUE)) {
  install.packages("rentrez")
}
library(rentrez)
library(Biostrings)
library(rtracklayer)

setwd("C:/Users/jahei/OneDrive/Työpöytä/Power_bi_project/R")

# Specify the organism and download genome
genome_id <- "BA000019.2"  # Example: *Nostoc* PCC 7210
fasta_file <- entrez_fetch(db = "nuccore", id = genome_id, rettype = "fasta")

# Write to a file
writeLines(fasta_file, "nostoc_7120_genome.fasta")

# Load the genome
genome <- readDNAStringSet("nostoc_7120_genome.fasta")

# Define a function to extract upstream sequences (800bp)
extract_upstream <- function(gene, sequence, upstream_length = 800) {
  start_pos <- start(gene)
  end_pos <- end(gene)
  strand_info <- as.character(strand(gene))
  
  # Adjust the start position based on the strand
  if (strand_info == "+") {
    upstream_start <- max(1, start_pos - upstream_length)  # Avoid negative start
    upstream_end <- start_pos - 1
  } else {
    upstream_start <- end_pos + 1
    upstream_end <- end_pos + upstream_length
  }
  
  # Extract the upstream sequence
  upstream_sequence <- subseq(sequence, upstream_start, upstream_end)
  return(upstream_sequence)
}

# Search for binding sites using a consensus sequence
fura <- "AAATAAATTCTCAATAAAT"
ntca <- "TGTNNNNNNNNNACA"
fura_binding_sites <- matchPattern(fura, genome[[1]])
ntca_binding_sites <- matchPattern(ntca, genome[[1]])

# Load a GFF file
gtf_file <- "nostoc_7120.gtf"
gtf_data <- import(gtf_file)
#gtf <- read.table(gtf_file, sep = "\t", header = FALSE, comment.char = "#")

# Example: Map binding sites to features
binding_positions <- start(fura_binding_sites)
annotated_sites <- gtf[apply(gtf, 1, function(row) {
  any(binding_positions >= as.numeric(row[4]) & binding_positions <= as.numeric(row[5]))
}), ]
