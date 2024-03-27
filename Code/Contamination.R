#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
data_dir = args[1]
out_dir = args[2]
path_to_flist = args[3]
joined_threshold = as.numeric(args[4])
contamination_threshold = as.numeric(args[5])

#data_dir = "/Users/jorgeamaya/Desktop/Terra_Development/amplicon_decontamination_pipeline/Report/DADA2_Contamination/" 
#out_dir = "/Users/jorgeamaya/Desktop/Terra_Development/amplicon_decontamination_pipeline/Report/"
#path_to_flist = "/Users/jorgeamaya/Desktop/Terra_Development/amplicon_decontamination_pipeline/Data/barcodes_matches_iseq_benchmarking.csv"
#joined_threshold = 1000
#contamination_threshold = 0.5

if (!require("ggplot2")) {
  install.packages("ggplot2", repos="http://cran.rstudio.com/")
  library("ggplot2")
}
if (!require("viridis")) {
  install.packages("viridis", repos="http://cran.rstudio.com/")
  library("viridis")
}
if (!require("reshape")) {
  install.packages("reshape", repos="http://cran.rstudio.com/")
  library("reshape")
}
if (!require("gridExtra")) {
  install.packages("gridExtra", repos="http://cran.rstudio.com/")
  library("gridExtra")
}

df = data.frame(sequence = character(),
                abundance = integer(),
                forward = integer(),
                reverse = integer(),
                nmatch = integer(),
                nmismatch = integer(),
                nindel = integer(),
                prefer = integer(),
                accept = character(),
                forward_barcode = character(),
                reverse_barcode = character(),
                forward_distance = numeric(),
                reverse_distance = numeric(),
                five_prime_end = character(),
                three_prime_end = character(),
                match_status = numeric(),
                insert_size = numeric(),
                stringsAsFactors = FALSE)

file_list <- list.files(path = data_dir, pattern = "\\.tsv$")
file_list = file_list[file_list != "bbmergefields.tsv"]

# Loop through each file and append to the data frame row-wise
# This has to be done differently depending on whether the contamination
# report is being written from BBMerge or DADA2

if(basename(data_dir) == 'Merge'){
  for (file in file_list) {
    file_path = paste0(path = data_dir, file) # current directory assumed
    print("file_path")
    print(file_path)
    if (file.info(file_path)$size == 1) {
      next
    }
    data <- read.csv(file_path, sep="\t")     # read file as data frame
    data$sample_id <- sub("_[^_]*$", "", file)
    df <- rbind(df, data)           # append data frame to df
  }
} else {
  #For DADA2_Contamination
  empty_samples = c()
  for (file in file_list) {
    file_path = paste0(path = data_dir, file) # current directory assumed
    print(file_path)
    if (file.info(file_path)$size == 1) {
      next
    }
    data <- read.csv(file_path, sep="\t")     # read file as data frame
    if (nrow(data) == 0) {
      empty_samples = c(empty_samples, sub("_[^_]*$", "", file))
    } else {
      data$sample_id <- sub("_[^_]*$", "", file)
      df <- rbind(df, data)           # append data frame to df
    }
  }
}

barcodes_list = read.csv(path_to_flist, sep = ",", header = TRUE)
samples_order = as.character(barcodes_list$sample_id)

#Handle empty dada2 file
#if(basename(data_dir) == 'DADA2_Contamination'){
#  write.csv(empty_samples, file=file.path(out_dir, "drop_samples_dada2.tsv"), 
#            row.names = FALSE,
#            quote = FALSE)
#  samples_order = samples_order[!samples_order %in% empty_samples]
#}

###################################
###MATCH STATUS PER SAMPLE ID #####
###################################
#Subset the table to the desired experiments
df = df[df$sample_id %in% samples_order,] 

m_sample_status = table(df[,c("sample_id", "match_status")])
m_sample_status_melted = melt(m_sample_status, id.vars=c("sample_id", "match_status"))

m_sample_status_p = round(m_sample_status/rowSums(m_sample_status), 2)
m_sample_status_p_melted = melt(m_sample_status_p, id.vars=c("sample_id", "match_status"))

m_sample_status_melted$sample_id = factor(m_sample_status_melted$sample_id,
                                          levels=rev(samples_order))

m_sample_status_p_melted$sample_id = factor(m_sample_status_p_melted$sample_id,
                                          levels=rev(samples_order))

match_order = c("Match",
                "Missmatch",
                "Forward_Match_Reverse_Missmatch",
                "Forward_Match_Reverse_Missing",
                "Forward_Missmatch_Reverse_Match",
                "Forward_Missing_Reverse_Match",
                "Forward_Missmatch_Reverse_Missing",
                "Reverse_Missmatch_Forward_Missing",
                "No_barcodes")

m_sample_status_melted$match_status = factor(m_sample_status_melted$match_status, levels=match_order)

m_sample_status_p_melted$match_status = factor(m_sample_status_p_melted$match_status, levels=match_order)

#Limit the plots to the match categories that are in the dataset
match_order = match_order[match_order %in% unique(m_sample_status_melted$match_status)]

plotabs = ggplot(m_sample_status_melted, aes(match_status, sample_id)) +    
  geom_tile(aes(fill = value)) +
  guides(fill = guide_colourbar(title = "Merged reads count")) +
  scale_fill_viridis() + theme_bw() +
  xlab("Match Status") + 
  ylab("Experiment/Sample ID") +
  theme(legend.position="top",
        legend.key.width= unit(0.5, 'in'),
        axis.text.x = element_text(angle = 45, hjust=1)) + 
  scale_x_discrete(labels = gsub("_", " ", match_order))

plotper = ggplot(m_sample_status_p_melted, aes(match_status, sample_id)) +   
  geom_tile(aes(fill = value)) +
  guides(fill = guide_colourbar(title = "Proportion of merged reads")) +
  scale_fill_viridis(option = "inferno") + theme_bw() +
  xlab("Match Status") + 
  ylab("Experiment/Sample ID") +
  theme(legend.position="top",
        legend.key.width= unit(0.5, 'in'),
        axis.text.x = element_text(angle = 45, hjust=1)) + 
  scale_x_discrete(labels = gsub("_", " ", match_order))

svg(file.path(out_dir, "Match_report_abs.svg"), width = 11, height = 20)
print(plotabs)
dev.off()

svg(file.path(out_dir, "Match_report_per.svg"), width = 11, height = 20)
print(plotper)
dev.off()

###################################
###MATCH MATRIX BETWEEN BARCODES###
###################################
m_sample_status = table(df[,c("forward_barcodes", "reverse_barcodes")])
m_sample_status_melted = melt(m_sample_status, id.vars=c("forward_barcodes", "reverse_barcodes"))

m_sample_status_p = round(m_sample_status/rowSums(m_sample_status), 2)
m_sample_status_p_melted = melt(m_sample_status_p, id.vars=c("forward_barcodes", "reverse_barcodes"))

plotabs = ggplot(m_sample_status_melted, aes(forward_barcodes, reverse_barcodes)) +    
  geom_tile(aes(fill = value)) +
  scale_fill_viridis() + theme_bw() +
  guides(fill = guide_colourbar(title = "Match counts")) +
  xlab("Forward Barcode") + 
  ylab("Reverse Barcode") +
  theme(legend.position="top",
        legend.key.width= unit(0.5, 'in'),
        axis.text.x = element_text(angle = 45, hjust=1)) 

plotper = ggplot(m_sample_status_p_melted, aes(forward_barcodes, reverse_barcodes)) +   
  geom_tile(aes(fill = value)) +
  scale_fill_viridis(option = "inferno") + theme_bw() +
  guides(fill = guide_colourbar(title = "Match proportion against forward barcode")) +
  xlab("Forward Barcode") + 
  ylab("Reverse Barcode") +
  theme(legend.position="top",
        legend.key.width= unit(0.5, 'in'),
        axis.text.x = element_text(angle = 45, hjust=1)) 

svg(file.path(out_dir, "Barcode_report_abs.svg"), width = 18, height = 18)
print(plotabs)
dev.off()

svg(file.path(out_dir, "Barcode_report_per.svg"), width = 18, height = 18)
print(plotper)
dev.off()

###################################
###BARCODES HAMMING DISTANCES   ###
###################################
m_sample_status = table(df[,c("forward_barcodes", "forward_distances")])
m_sample_status_p = m_sample_status/rowSums(m_sample_status)
m_sample_status_p_melted = melt(m_sample_status_p, id.vars=c("forward_barcodes", "forward_distances"))

m_sample_status_p_melted = m_sample_status_p_melted[-which(m_sample_status_p_melted == 'NULL'),]

m_sample_status_p_melted_reshaped <- reshape(m_sample_status_p_melted, direction = "wide", idvar = "forward_barcodes", timevar = "forward_distances")

if (!'value.0' %in% colnames(m_sample_status_p_melted_reshaped)){
  m_sample_status_p_melted_reshaped$value.0 = rep(0,nrow(m_sample_status_p_melted_reshaped))
}

if (!'value.1' %in% colnames(m_sample_status_p_melted_reshaped)){
  m_sample_status_p_melted_reshaped$value.1 = rep(0,nrow(m_sample_status_p_melted_reshaped))
}

if (!'value.2' %in% colnames(m_sample_status_p_melted_reshaped)){
  m_sample_status_p_melted_reshaped$value.2 = rep(0,nrow(m_sample_status_p_melted_reshaped))
}

m_sample_status_p_melted_reshaped = m_sample_status_p_melted_reshaped[, c("forward_barcodes",
                                                                          "value.0",
                                                                          "value.1",
                                                                          "value.2")]
m_sample_status_p_melted_reshaped = m_sample_status_p_melted_reshaped[which(m_sample_status_p_melted_reshaped$forward_barcode != 'NULL'),]

colnames(m_sample_status_p_melted_reshaped) = c("Forward Barcode", "0", "1", "2")

write.table(m_sample_status_p_melted_reshaped, 
            file=file.path(out_dir, "hamming_forward.tsv"), 
            quote = FALSE, 
            sep = "\t", 
            row.names = FALSE,
            col.names = TRUE)

m_sample_status = table(df[,c("reverse_barcodes", "reverse_distances")])
m_sample_status_p = m_sample_status/rowSums(m_sample_status)
m_sample_status_p_melted = melt(m_sample_status_p, id.vars=c("reverse_barcodes", "reverse_distances"))

m_sample_status_p_melted = m_sample_status_p_melted[-which(m_sample_status_p_melted == 'NULL'),]

m_sample_status_p_melted_reshaped <- reshape(m_sample_status_p_melted, direction = "wide", idvar = "reverse_barcodes", timevar = "reverse_distances")

if (!'value.0' %in% colnames(m_sample_status_p_melted_reshaped)){
  m_sample_status_p_melted_reshaped$value.0 = rep(0,nrow(m_sample_status_p_melted_reshaped))
}

if (!'value.1' %in% colnames(m_sample_status_p_melted_reshaped)){
  m_sample_status_p_melted_reshaped$value.1 = rep(0,nrow(m_sample_status_p_melted_reshaped))
}

if (!'value.2' %in% colnames(m_sample_status_p_melted_reshaped)){
  m_sample_status_p_melted_reshaped$value.2 = rep(0,nrow(m_sample_status_p_melted_reshaped))
}

m_sample_status_p_melted_reshaped = m_sample_status_p_melted_reshaped[, c("reverse_barcodes",
                                                                          "value.0",
                                                                          "value.1",
                                                                          "value.2")]
m_sample_status_p_melted_reshaped = m_sample_status_p_melted_reshaped[which(m_sample_status_p_melted_reshaped$reverse_barcodes != 'NULL'),]

colnames(m_sample_status_p_melted_reshaped) = c("Reverse Barcode", "0", "1", "2")

write.table(m_sample_status_p_melted_reshaped, 
            file=file.path(out_dir, "hamming_reverse.tsv"), 
            quote = FALSE, 
            sep = "\t", 
            row.names = FALSE,
            col.names = TRUE)
        
###################################
###INSERT SIZE DISTRIBUTION     ###
###################################

insert_size = sort(df$insert_size)

png(file.path(out_dir, "Insert_size.png"))
barplot(insert_size, xaxt='n',
        ylab = "Insert size",
        xlab = "Read number",
        main = "Insert size distribution")
axis(side=1, 
     at=seq(0, length(df$insert_size), by=500000), 
     labels=seq(0, length(df$insert_size), by=500000), 
     cex.axis=1,
     ) 
dev.off()

###################################
########   FLAG ALGORITHM       ###
###################################
#DISAGGREGATED REPORT

if(basename(data_dir) == 'Merge'){
  missmatch_df = df
  missmatch_df$barcode_pair = paste0(missmatch_df$forward_barcode, "/", missmatch_df$reverse_barcode)

  m_sample_status = table(missmatch_df[,c("sample_id", "barcode_pair")])
  m_sample_status_p = round(m_sample_status/rowSums(m_sample_status), 3)
  m_sample_status_p_melted = melt(m_sample_status_p, id.vars=c("sample_id", "barcode_pair"))

  m_sample_status_p_melted_reshaped = as.data.frame.matrix(xtabs(value ~ sample_id + barcode_pair, data = m_sample_status_p_melted))

  m_sample_status_p_melted_reshaped$sample_id = rownames(m_sample_status_p_melted_reshaped)
  m_sample_status_p_melted_reshaped = m_sample_status_p_melted_reshaped[match(samples_order, m_sample_status_p_melted_reshaped$sample_id),]

  mergedata = read.csv(file.path(data_dir, "bbmergefields.tsv"), sep="\t", header=TRUE)
  mergedata = mergedata[match(samples_order, mergedata$SampleID),]

  m_sample_status_p_melted_reshaped$Well_Productivity = mergedata$Joined
  m_sample_status_p_melted_reshaped$Productivity_Flag = mergedata$Joined < joined_threshold
  
  #Adjust table to include samples with no reads
  rows = rownames(m_sample_status_p_melted_reshaped) 
  no_mergers = samples_order[!rows %in% samples_order]
  rows[grep('^NA', rows)] = no_mergers
  rownames(m_sample_status_p_melted_reshaped) = rows
  
  #Fill the threshold and wrong barcode combinations flags
  Below_Threshold_Flag = vector()
  Wrong_Barcode_Comb_Flag = vector()
  
  for (sample in rownames(m_sample_status_p_melted_reshaped)) {
    print(sample)
    if (!sample %in% no_mergers) {
      tmp = m_sample_status_p_melted_reshaped[sample, 
                                              !colnames(m_sample_status_p_melted_reshaped) %in% c("sample_id", "Productivity_Flag", "Well_Productivity")]
      abs_values <- sapply(tmp, function(x) abs(x))
      column_with_max_abs <- names(tmp)[which.max(abs_values)]
      max_freq_barcodes = unlist(barcodes_list[barcodes_list$sample_id == sample, c("Forward", "Reverse")])
      names(max_freq_barcodes) = NULL
      
      below_thres = tmp[which.max(abs_values)] < 0.5
      are_equal = !identical(strsplit(column_with_max_abs, "/")[[1]], max_freq_barcodes)
      
      Below_Threshold_Flag = c(Below_Threshold_Flag, below_thres)
      Wrong_Barcode_Comb_Flag = c(Wrong_Barcode_Comb_Flag, are_equal)
    }else{
      Below_Threshold_Flag = c(Below_Threshold_Flag, NA)
      Wrong_Barcode_Comb_Flag = c(Wrong_Barcode_Comb_Flag, NA)
    }
  }
  
  m_sample_status_p_melted_reshaped$Below_Threshold_Flag = Below_Threshold_Flag
  m_sample_status_p_melted_reshaped$Wrong_Barcode_Comb_Flag = Wrong_Barcode_Comb_Flag
  m_sample_status_p_melted_reshaped$Flagged_well = m_sample_status_p_melted_reshaped$Productivity_Flag |
    m_sample_status_p_melted_reshaped$Below_Threshold_Flag |
    m_sample_status_p_melted_reshaped$Wrong_Barcode_Comb_Flag
  
  m_sample_status_p_melted_reshaped = m_sample_status_p_melted_reshaped[,
                                                                        !colnames(m_sample_status_p_melted_reshaped) %in% c("sample_id")]
  
  write.table(m_sample_status_p_melted_reshaped, 
              file=file.path(out_dir, "barcodes_report_bbmerge.tsv"), 
              quote = FALSE, 
              sep = "\t", 
              row.names = TRUE,
              col.names = TRUE)
} else { #Process DADA2 contamination report
  missmatch_df = df
  missmatch_df$barcode_pair = paste0(missmatch_df$forward_barcode, "/", missmatch_df$reverse_barcode)
  
  m_sample_status = table(missmatch_df[,c("sample_id", "barcode_pair")])
  m_sample_status_p = round(m_sample_status/rowSums(m_sample_status), 3)
  m_sample_status_p_melted = melt(m_sample_status_p, id.vars=c("sample_id", "barcode_pair"))
  
  m_sample_status_p_melted_reshaped = as.data.frame.matrix(xtabs(value ~ sample_id + barcode_pair, data = m_sample_status_p_melted))
  
  m_sample_status_p_melted_reshaped$sample_id = rownames(m_sample_status_p_melted_reshaped)
  m_sample_status_p_melted_reshaped = m_sample_status_p_melted_reshaped[match(samples_order, m_sample_status_p_melted_reshaped$sample_id),]
  
  mergedata = read.table(file.path(dirname(dirname(data_dir)), 
                                 "Results", "DADA2_Contamination", "seqtab.tsv"), 
                       sep="\t", header=TRUE, row.names = 2)
  mergedata <- mergedata[match(samples_order, row.names(mergedata)), ]
  mergedata <- mergedata[, !(names(mergedata) == "X")]
  
  m_sample_status_p_melted_reshaped$Well_Productivity = rowSums(mergedata)
  m_sample_status_p_melted_reshaped$Productivity_Flag = m_sample_status_p_melted_reshaped$Well_Productivity < joined_threshold
  
  #Adjust table to include samples with no reads
  rows = rownames(m_sample_status_p_melted_reshaped)
  no_mergers = samples_order[!rows %in% samples_order]
  rows[grep('^NA', rows)] = no_mergers  
  rownames(m_sample_status_p_melted_reshaped) = rows
  
  Below_Threshold_Flag = vector()
  Wrong_Barcode_Comb_Flag = vector() 
  
  for (sample in rownames(m_sample_status_p_melted_reshaped)) {
    print(sample)
    if (!sample %in% no_mergers) {
      tmp = m_sample_status_p_melted_reshaped[sample, 
                                              !colnames(m_sample_status_p_melted_reshaped) %in% c("sample_id", "Productivity_Flag", "Well_Productivity")]
      abs_values <- sapply(tmp, function(x) abs(x))
      column_with_max_abs <- names(tmp)[which.max(abs_values)]
      max_freq_barcodes = unlist(barcodes_list[barcodes_list$sample_id == sample, c("Forward", "Reverse")])
      names(max_freq_barcodes) = NULL
      
      below_thres = tmp[which.max(abs_values)] < 0.5
      are_equal = !identical(strsplit(column_with_max_abs, "/")[[1]], max_freq_barcodes)
      
      Below_Threshold_Flag = c(Below_Threshold_Flag, below_thres)
      Wrong_Barcode_Comb_Flag = c(Wrong_Barcode_Comb_Flag, are_equal)
    }else{
      Below_Threshold_Flag = c(Below_Threshold_Flag, NA)
      Wrong_Barcode_Comb_Flag = c(Wrong_Barcode_Comb_Flag, NA)
    }
  }
  
  m_sample_status_p_melted_reshaped$Below_Threshold_Flag = Below_Threshold_Flag
  m_sample_status_p_melted_reshaped$Wrong_Barcode_Comb_Flag = Wrong_Barcode_Comb_Flag
  m_sample_status_p_melted_reshaped$Flagged_well = m_sample_status_p_melted_reshaped$Productivity_Flag |
    m_sample_status_p_melted_reshaped$Below_Threshold_Flag |
    m_sample_status_p_melted_reshaped$Wrong_Barcode_Comb_Flag
  
  m_sample_status_p_melted_reshaped = m_sample_status_p_melted_reshaped[,
                                                                        !colnames(m_sample_status_p_melted_reshaped) %in% c("sample_id")]
  
  m_sample_status_p_melted_reshaped$samples <- rownames(m_sample_status_p_melted_reshaped)
  rownames(m_sample_status_p_melted_reshaped) <- NULL
  m_sample_status_p_melted_reshaped <- m_sample_status_p_melted_reshaped[, c("samples", setdiff(names(m_sample_status_p_melted_reshaped), "samples"))]
  
  write.table(m_sample_status_p_melted_reshaped, 
              file=file.path(out_dir, "barcodes_report_dada2.tsv"), 
              quote = FALSE, 
              sep = "\t", 
              row.names = FALSE,
              col.names = TRUE)
}

