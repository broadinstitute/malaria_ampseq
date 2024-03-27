#!/bin/R env

########################################
# LOAD LIBRARIES AND PREPARE VARIABLES #
########################################

if (!require("argparse")) {
  install.packages("argparse", repos="http://cran.rstudio.com/")
  library("argparse")
}

if (!require("stringdist")) {
  install.packages("stringdist", repos="http://cran.rstudio.com/")
  library("stringdist")
}

# Custom filtering, denoising parameters (if not default) can be provided as a separate config file?
parser <- ArgumentParser()
parser$add_argument("-p", "--path_to_fastq", help="Path to merged fastq file (required)")
parser$add_argument("-d", "--dir", help="Working directory path for writing all bbmerge contamination output files")
parser$add_argument("-b", "--barcodes", help="The path to a csv file with the sample_id,Forward,Reverse, where Forward and Reverse are columns with the barcodes for the sample")
parser$add_argument("--terra", action='store_true', help="State whether the scripts are running in Terra")
args <- parser$parse_args()

# Universal parameters
fastq_file <- args$path_to_fastq
work_dir <- args$dir
path_to_flist <- args$barcodes
terra = args$terra

print(paste0("Processing file: ", fastq_file))
write(paste0("Processing file: ", fastq_file), stderr())
########################################
#             BARCODE REPORT           #
########################################

if (terra) { 
	source(file.path("/Code", "matching_functions.R"))
} else {
	source(file.path("Code", "matching_functions.R"))
}
#source(paste0(file.path(dirname(dirname(work_dir)), "Code", "matching_functions.R")))
barcodes = read.csv(path_to_flist, sep = ",", header = TRUE)
dist = 2


# Open the Fastq file and count the number of reads
con <- file(fastq_file, "r")
num_lines <- length(readLines(con))
close(con)
num_reads <- num_lines / 4
print(paste("Number of reads in the FASTQ file:", num_reads))

seq_names <- vector("character", length = num_reads)
forward_barcodes <- vector("character", length = num_reads)
reverse_barcodes <- vector("character", length = num_reads)
forward_distances <- vector("numeric", length = num_reads)
reverse_distances <- vector("numeric", length = num_reads)
five_prime_ends <- vector("numeric", length = num_reads)
three_prime_ends <- vector("numeric", length = num_reads)
match_categories <- vector("numeric", length = num_reads)
insert_sizes <- vector("numeric", length = num_reads)

# Reopen the Fastq file
con <- file(fastq_file, "r")
sample <- sub(".*/([^/]+)_merged\\.fastq", "\\1", fastq_file)
print(sample)

i <- 0
while (length(line <- readLines(con, n = 4)) > 0) {
  i <- i + 1
  seq_name <- line[1]
  sequence <- line[2]

  match_tmp <- match_fun(sample, sequence, barcodes, dist)

  # Assign values directly to preallocated vectors
  seq_names[i] <- seq_name
  forward_barcodes[i] <- match_tmp[[1]]
  reverse_barcodes[i] <- match_tmp[[2]]
  forward_distances[i] <- match_tmp[[3]]
  reverse_distances[i] <- match_tmp[[4]]
  five_prime_ends[i] <- match_tmp[[5]]
  three_prime_ends[i] <- match_tmp[[6]]
  match_categories[i] <- match_tmp[[7]]
  insert_sizes[i] <- match_tmp[[8]]
}

close(con)

# Create data frame
df <- data.frame(
  sequence = seq_names,
  forward_barcodes = forward_barcodes,
  reverse_barcodes = reverse_barcodes,
  forward_distances = forward_distances,
  reverse_distances = reverse_distances,
  five_prime_ends = five_prime_ends,
  three_prime_ends = three_prime_ends,
  match_status = match_categories,
  insert_sizes = insert_sizes
)
              
output_filename_o = file.path(dirname(dirname(work_dir)), "Report", "Merge", paste0(sample, "_final.tsv"))
write.table(df, file=output_filename_o, quote = FALSE, sep = "\t", row.names = FALSE)
print(paste0("BBmerge contamination report written to ", output_filename_o))
