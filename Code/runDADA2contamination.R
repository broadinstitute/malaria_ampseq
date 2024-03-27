#!/bin/R env

########################################
# LOAD LIBRARIES AND PREPARE VARIABLES #
########################################

if (!require("dada2")) {
  install.packages("dada2", repos="http://cran.rstudio.com/")
  library("dada2")
}
if (!require("limma")) {
  install.packages("limma", repos="http://cran.rstudio.com/")
  library("limma")
}
if (!require("data.table")) {
  install.packages("data.table", repos="http://cran.rstudio.com/")
  library("data.table")
}
if (!require("argparse")) {
  install.packages("argparse", repos="http://cran.rstudio.com/")
  library("argparse")
}
if (!require("stringdist")) {
  install.packages("stringdist", repos="http://cran.rstudio.com/")
  library("stringdist")
}
if (!require("viridisLite")) {
  install.packages("viridis", repos="http://cran.rstudio.com/")
  library("viridis")
}

# Custom filtering, denoising parameters (if not default) can be provided as a separate config file?
parser <- ArgumentParser()
parser$add_argument("-p", "--path_to_meta", help="Path to input meta file listing fastqs (required)")
parser$add_argument("-r", "--path_to_data_repo", help="Path to original fastq files (required)")
parser$add_argument("-b", "--path_to_flist", help="The path to a one column file with the sample_id of the samples in the experiment (required)")
parser$add_argument("-c", "--class", help="Class specifying 'parasite' or 'vector' (required if '--default' is specified)")
parser$add_argument("-d", "--dir", help="Working directory path for writing all dada2 output files")
parser$add_argument("-o", "--output_filename", help="output tab-separated filename (required)")
parser$add_argument("-s", "--save_run", help="save Run as R workspace image")
parser$add_argument("-ee", "--maxEE", help="Maximum expected errors for filtering forward and reverse read")
parser$add_argument("-tR", "--trimRight", help="Length for trimming from right for both forward and reverse reads")
parser$add_argument("-mL", "--minLen", type="integer", help="Minimum length required for reads on both end. Shorter reads are discarded")
parser$add_argument("-tQ", "--truncQ", help="Truncate reads to first occurence of truncQ. All filtered reads have quality >= truncQ")
parser$add_argument("-id", "--matchIDs", type="integer", help="Match ids on fastqs to make sure reads on forward and reverse end are in same order")
parser$add_argument("-mC", "--max_consist", type="integer", help="Maximum cycles for error model until consistency. If no convergence, error values at max_consist cycle are used")
parser$add_argument("-wA", "--omega_a", type="double", help="P-value threshold in sample inference for forming a new partition")
parser$add_argument("-jC", "--justConcatenate", type="integer", help="Specify whether ASVs need to be concatenated with Ns instead of merging")
parser$add_argument("-mM", "--maxMismatch", type="integer", help="Specify the maximum number of mismatches allowed during merging")
parser$add_argument("--bimera", action='store_true', help="Optionally output list of sequences identified as bimeras")
parser$add_argument("--terra", action='store_true', help="State whether the scripts are running in Terra")
args <- parser$parse_args()

# work_dir = '/Users/jorgeamaya/Desktop/Terra_Development/amplicon_decontamination_pipeline/Results/DADA2_OP_Contamination'
# path_to_meta = '/Users/jorgeamaya/Desktop/Terra_Development/amplicon_decontamination_pipeline/Results/PrimerRem_OP/mixed_op_prim_meta.tsv'
# output_filename = '/Users/jorgeamaya/Desktop/Terra_Development/amplicon_decontamination_pipeline/Results/DADA2_OP_Contamination/seqtab.tsv'
# path_to_data_repo = '/Users/jorgeamaya/Desktop/Broad_Test/Benchmarking/iSeq_ci_Benchmarking'
# path_to_flist = '/Users/jorgeamaya/Desktop/Terra_Development/amplicon_decontamination_pipeline/Data/barcodes_matches_iseq_benchmarking.csv'
# 
# matchIDs <- '0'
# maxEE <- '5,5'
# trimRight <- '0,0'
# minLen <- '30'
# truncQ <- '5,5'
# max_consist <- '10'
# omega_a <- '1e-120'
# class = 'parasite'
# justConcatenate <- '0'
# maxMismatch = '0'
# bimera = TRUE
# save_run = "dada2.RData"

# Universal parameters

path_to_meta <- args$path_to_meta
path_to_data_repo = args$path_to_data_repo
work_dir <- args$dir
path_to_flist <- args$path_to_flist
output_filename = args$output_filename
save_run = args$save_run
terra = args$terra

# obtain/initialize Parameters
# (Universal) Parameters
randomize=TRUE
selfConsist=TRUE
filter = TRUE
matchIDs <- args$matchIDs
class = args$class

# Parameters for merging
justConcatenate <- args$justConcatenate
bimera = args$bimera
maxMismatch = args$maxMismatch

# DADA2 and Filtering parameters
maxEE <- args$maxEE
trimRight <- args$trimRight
minLen <- args$minLen
truncQ <- args$truncQ
max_consist <- args$max_consist
omega_a <- args$omega_a

if (file.exists(path_to_meta)) {
  metafile <- fread(path_to_meta, sep = "\t", header=FALSE)
  sample.names <- metafile$V1
  fnFs <- metafile$V2
  fnRs <- metafile$V3
} else {
  stop(paste("metafile",path_to_meta,"not found!"))
}

if (is.null(matchIDs)||matchIDs == '') {
  matchIDs = TRUE
} else {
  matchIDs = as.logical(as.numeric(matchIDs))
}
if (is.null(trimRight)||trimRight == '') {
  trimRight = c(0,0)
} else {
  trimRight <- as.numeric(strsplit(trimRight,',')[[1]])
}
if (is.null(truncQ)||truncQ == '') {
  truncQ = c(5,5)
} else {
  truncQ <- as.numeric(strsplit(truncQ,',')[[1]])
}

if (class == "parasite") {
  # Parameters for filtering
  if (is.null(maxEE)||maxEE == '') {
    maxEE = c(5,5)
  } else {
    maxEE <- as.numeric(strsplit(maxEE,',')[[1]])
  }
  if (is.null(minLen)||minLen == '') {
    minLen=30
  } else {
    minLen = as.numeric(minLen)
  }
  
  # Parameters for Denoising
  if (is.null(max_consist)||max_consist == '') {
    max_consist=10
  } else {
    max_consist = as.numeric(max_consist)
  }
  if (is.null(omega_a)||omega_a == '') {
    omega_a=1e-120
  } else {
    omega_a = as.numeric(omega_a)
  }
  
  #Parameters for merging
  if (is.null(justConcatenate)||justConcatenate == '') {
    justConcatenate = FALSE
  } else {
    justConcatenate = as.logical(as.numeric(justConcatenate))
  }
} else if (class == "vector") {
  # Parameters for filtering
  if (is.null(maxEE)||maxEE == '') {
    maxEE = c(2,2)
  } else {
    maxEE <- as.numeric(strsplit(maxEE,',')[[1]])
  }
  if (is.null(minLen)||minLen == '') {
    minLen=75
  } else {
    minLen = as.numeric(minLen)
  }
  
  # Parameters for Denoising
  if (is.null(max_consist)||max_consist == '') {
    max_consist=20
  } else {
    max_consist = as.numeric(max_consist)
  }
  if (is.null(omega_a)||omega_a == '') {
    omega_a=1e-40
  } else {
    omega_a = as.numeric(omega_a)
  }
  
  #Parameters for merging
  if (is.null(justConcatenate)||justConcatenate == '') {
    justConcatenate = TRUE
  } else {
    justConcatenate = as.logical(as.numeric(justConcatenate))
  }
  if (is.null(maxMismatch)||maxMismatch == '') {
    maxMismatch=30
  } else {
    maxMismatch = as.numeric(maxMismatch)
  }
  
} else {
  stop("Please provide valid option for the '--class' argument")
}

#Output parameters
if (dirname(output_filename) != ".") {
  output_filename <- output_filename
} else {
  output_filename <- paste0(work_dir, "/", output_filename)
}

#Datatable to summarize parmeters
parameter_df <- data.frame(maxEE=maxEE,
                           trimRight=trimRight,
                           minLen=minLen,
                           truncQ=truncQ,
                           matchIDs=matchIDs,
                           max_consist=max_consist,
                           randomize=randomize,
                           selfConsist=selfConsist,
                           OMEGA_A=omega_a,
                           justConcatenate=justConcatenate)

print(parameter_df)

# List files and sample names
if (length(fnFs) == 0 || length(fnFs) != length(fnRs)) {
  stop("fastq files incomplete or not found")
}

########################################
#                DADA2                 #
########################################

# Create paths for filtered fastq
filtFs <- file.path(work_dir, "filtered", paste0(sample.names, "_filt_R1.fastq.gz"))
filtRs <- file.path(work_dir, "filtered", paste0(sample.names, "_filt_R2.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter read
if (filter == TRUE) {
  print("filtering samples...")
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                       maxN=0, maxEE=maxEE, trimRight=trimRight, truncQ=truncQ, minLen=minLen,
                       rm.phix=TRUE, compress=TRUE, multithread=TRUE, verbose=TRUE,
                       matchIDs=matchIDs)
  print("filtering done!")
} else {
  print("skipping filter except mandatory removal of N's... ")
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncQ=c(0,0), maxN=0, rm.phix=TRUE,
                       compress=TRUE, multithread=TRUE, verbose=TRUE, matchIDs=matchIDs)
}

# Report and Correct for samples with zero reads after filter
zeros <- row.names(out)[out[,2] == 0]
filtFs <- filtFs[out[,2] != 0]
filtRs <- filtRs[out[,2] != 0]
sample.names <- sample.names[out[,2] != 0]

# Update Out table
out <- out[(out[,2] != 0),]

#Compute the error model
print("starting error model learning for forward reads...")
errF <- learnErrors(filtFs, multithread=TRUE, verbose=2, randomize=randomize, MAX_CONSIST=max_consist)
print("starting error model learning for reverse reads...")
errR <- learnErrors(filtRs, multithread=TRUE, verbose=2, randomize=randomize, MAX_CONSIST=max_consist)

#DeReplicate the reads
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Run core DADA2 algorithm
print("starting dada2 for forward reads...")
dadaFs <- dada(derepFs, err=errF, selfConsist=selfConsist, multithread=TRUE, verbose=TRUE, OMEGA_A=omega_a)
print("starting dada2 for reverse reads...")
dadaRs <- dada(derepRs, err=errR, selfConsist=selfConsist, multithread=TRUE, verbose=TRUE, OMEGA_A=omega_a)

# Merge reads
print("merging paird ends...")
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, justConcatenate=justConcatenate, trimOverhang = TRUE, maxMismatch = maxMismatch)

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

for(sample in sample.names){
  forward_barcodes = c()
  reverse_barcodes = c()
  forward_distances = c()
  reverse_distances = c()
  five_prime_ends = c()
  three_prime_ends = c()
  match_categories = c()
  insert_sizes = c()
  print(sample)
  for (sequence in mergers[[sample]]$sequence){
    match_tmp = match_fun(sample, sequence, barcodes, dist)
      
    forward_barcodes = c(forward_barcodes, match_tmp[[1]])
    reverse_barcodes = c(reverse_barcodes, match_tmp[[2]])
    forward_distances = c(forward_distances, match_tmp[[3]])
    reverse_distances = c(reverse_distances, match_tmp[[4]])
    five_prime_ends = c(five_prime_ends, match_tmp[[5]])
    three_prime_ends = c(three_prime_ends, match_tmp[[6]])
    match_categories = c(match_categories, match_tmp[[7]])
    insert_sizes = c(insert_sizes, match_tmp[[8]])
  }
  
  if(justConcatenate) {
    #DADA2 produces a dataframe with the sequences column as a list when it runs in justConcatenate mode.
    #(There is no clear explanation as to why). write.table returns an error if it tries to save a list.
    #A future version of dada2 may correct this particularity. For the moment, this if statement is 
    #necessary to correct said issue. 
    mergers[[sample]]$sequence = as.character(mergers[[sample]]$sequence)
  } 
  
  mergers[[sample]]$forward_barcodes = forward_barcodes 
  mergers[[sample]]$reverse_barcodes = reverse_barcodes 
  mergers[[sample]]$forward_distances = forward_distances 
  mergers[[sample]]$reverse_distances = reverse_distances 
  mergers[[sample]]$five_prime_end = five_prime_ends 
  mergers[[sample]]$three_prime_end = three_prime_ends 
  mergers[[sample]]$match_status = match_categories
  mergers[[sample]]$insert_size = insert_sizes 

  output_filename_o = file.path(dirname(dirname(work_dir)), "Report", "DADA2_Contamination", paste0(sample, "_final.tsv"))
  write.table(mergers[[sample]], file=output_filename_o, quote = FALSE, sep = "\t", row.names = FALSE)
  print(paste0("Contamination report written to ", output_filename_o))
}

#Generate sequence table
print("generating sequence table...")
seqtab <- makeSequenceTable(mergers)
print("Number of sequences in table")
print(dim(seqtab))
# Inspect distribution of sequence lengths
print(table(nchar(getSequences(seqtab))))

#Remove Chimeras
if(bimera) {
  print("identifying bimeric sequences...")
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  print("Number of non-bimeric sequences:")
  print(dim(seqtab.nochim)[2])
  print("Percentage of reads which are non-bimeric:")
  print(sum(seqtab.nochim)/sum(seqtab))
  bimeras <- !(colnames(seqtab) %in% colnames(seqtab.nochim))
  write.table(data.frame(sequence = colnames(seqtab), bimera = bimeras), file=paste0(work_dir,"/ASVBimeras.txt"),
              quote=FALSE, sep="\t", row.names=FALSE)
} else {
  print("skipping Bimera identification..")
}

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

#Get reads from original files
reads_report = data.frame()
count_reads <- function(file_path) {
  # Open the gzipped fastq file
  gz_conn <- gzfile(file_path, "r")
  
  # Initialize a variable to count the reads
  read_count <- 0
  
  # Iterate over the lines in the file
  while (length(line <- readLines(gz_conn, n = 4)) > 0) {
    if (length(line) == 4 && substr(line[1], 1, 1) == "@") {
      read_count <- read_count + 1
    }
  }
  
  # Close the connection to the file
  close(gz_conn)
  
  # Return the total number of reads
  return(read_count)
}

for(sample.name in names(track[,1])) {
  orig_f <- file.path(path_to_data_repo, paste0(sample.name, "_L001_R1_001.fastq.gz"))
  orig_r <- file.path(path_to_data_repo, paste0(sample.name, "_L001_R2_001.fastq.gz"))
  
  orig_f_reads = count_reads(orig_f)
  orig_r_reads = count_reads(orig_r)
  
  adap_f <- file.path(dirname(work_dir), "AdaptorRem", paste0(sample.name, "_val_1.fq.gz"))
  adap_r <- file.path(dirname(work_dir), "AdaptorRem", paste0(sample.name, "_val_2.fq.gz"))
  
  adap_f_reads = count_reads(adap_f)
  adap_r_reads = count_reads(adap_r)
  
  if (grepl('adaptorrem_meta.tsv', path_to_meta)) {
    prim = "PrimerRem"
    ext_1 = '_prim_1.fq.gz'
    ext_2 = '_prim_2.fq.gz'
  } else if (grepl('PrimerRem_NOP', path_to_meta)) {
    prim = "PrimerRem_NOP"
    ext_1 = '_mixed_nop_1.fq.gz'
    ext_2 = '_mixed_nop_2.fq.gz'
  } else if (grepl('PrimerRem_OP', path_to_meta)) {
    prim = "PrimerRem_OP"
    ext_1 = '_mixed_op_1.fq.gz'
    ext_2 = '_mixed_op_2.fq.gz'
  }
  
  prim_f <- file.path(dirname(work_dir), prim, paste0(sample.name, ext_1))
  prim_r <- file.path(dirname(work_dir), prim, paste0(sample.name, ext_2))
  
  prim_f_reads = count_reads(prim_f)
  prim_r_reads = count_reads(prim_r)
  
  reads_report = rbind(reads_report,
                       data.frame(sample.name,
                                  orig_f_reads,
                                  orig_r_reads,
                                  adap_f_reads,
                                  adap_r_reads,
                                  prim_f_reads,
                                  prim_r_reads))
}

stopifnot(names(track[,1]) == reads_report$sample.name)
rownames(reads_report) = reads_report$sample.name
reads_report = data.matrix(reads_report)

track = cbind(track, data.matrix(reads_report))
original_discarded = track[,"orig_f_reads"] - track[,"adap_f_reads"] 
adaptor_discarded = track[,"adap_f_reads"] - track[,"filtered"]
filtered_discarded = track[,"filtered"] - track[,"denoisedF"]
merged_discarded = track[,"denoisedF"] - track[,"merged"]
track = cbind(track, original_discarded, adaptor_discarded, filtered_discarded, merged_discarded)

# sink summary from stdout to a file
sink(paste0(work_dir,"/reads_summary.txt"))
print(track)
sink()

track_plot = as.data.frame(track[, c("merged", "merged_discarded", "filtered_discarded", "adaptor_discarded", "original_discarded")])

#Subset the table to the desired experiments
samples_order = read.csv(file.path(path_to_flist), sep = ",", header = FALSE)$V1
track_plot = track_plot[row.names(track_plot) %in% samples_order,] 

track_plot <- track_plot[order(-track_plot[,1], 
                               -track_plot[,2], 
                               -track_plot[,3],
                               -track_plot[,4],
                               -track_plot[,5]),]

track_plot_per = sweep(track_plot, 1, rowSums(track_plot), "/")

track_plot_per <- track_plot_per[order(-track_plot_per[,1], 
                                       -track_plot_per[,2], 
                                       -track_plot_per[,3],
                                       -track_plot_per[,4],
                                       -track_plot_per[,5]),]

color_vector = viridis(nrow(t(as.matrix(track_plot))), option = "D")

pdf(paste0(work_dir,"/stacked_barplot.pdf"), width = 12)
par(mar = c(14, 4, 6, 2), xpd = TRUE) # increase the bottom margin
barplot(
  t(as.matrix(track_plot)),
  col = c(color_vector, "red"),
  border = NA,
  ylab = "Reads Count",
  xlab = "",
  main = "DADA2 performance - Absolute",
  las = 2,
  cex.names = 0.5,
  cex.axis = 0.7,
  cex.main = 2
)

mtext("Experiment ID", side = 1, line = 12) # lower x-axis label

legend(
  "top",
  inset=c(0,-0.1),
  horiz = TRUE,
  legend = c("Merged", "Merged discarded", "Filtered discarded", "Adaptor discarded", "Original discarded"),
  fill = c(color_vector, "red"),
  bty = "n",
  cex = 1
)

dev.off()

pdf(paste0(work_dir,"/stacked_barplot_per.pdf"), width = 12)
par(mar = c(14, 4, 6, 2), xpd = TRUE) # increase the bottom margin
barplot(
  t(as.matrix(track_plot_per)),
  col = c(color_vector, "red"),
  border = NA,
  ylab = "Reads Percentage",
  xlab = "",
  main = "DADA2 performance - Percentage",
  las = 2,
  cex.names = 0.5,
  cex.axis = 0.7,
  cex.main = 2
)

mtext("Experiment ID", side = 1, line = 12) # lower x-axis label

legend(
  "top",
  inset=c(0,-0.1),
  horiz = TRUE,
  legend = c("Merged", "Merged discarded", "Filtered discarded", "Adaptor discarded", "Original discarded"),
  fill = c(color_vector, "red"),
  bty = "n",
  cex = 1
)

dev.off()

#Show the barplot of length distribution
pdf(paste0(work_dir,"/sequences_barplot.pdf"))
print(barplot(table(nchar(getSequences(seqtab)))))
dev.off()

#Add Zero Read Samples
zeros.df = data.frame(matrix(ncol = ncol(seqtab), nrow = length(zeros)))
colnames(zeros.df) = colnames(seqtab)

# Define a function to extract everything before the second underscore
extract_before_second_underscore <- function(input_string) {
  parts <- strsplit(input_string, "_")[[1]]
  result <- paste(parts[1], parts[2], sep = "_")
  return(result)
}
zeros_corrected <- sapply(zeros, extract_before_second_underscore)
rownames(zeros.df) = unname(zeros_corrected)
seqtab = rbind(seqtab, zeros.df)

#Generate output: sequence table to a tsv
write.table(seqtab, file=output_filename, quote = FALSE, sep = "\t")

# Save Run as R workspace image (Optional)
if (is.null(save_run)||save_run == '') {
  print("--save_run not found or empty. skip saving Rdata image to a file")
} else {
  save.image(paste0(work_dir,"/",save_run))
}

