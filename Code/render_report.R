library(argparse)
library(rmarkdown)

# Create ArgumentParser object
parser <- ArgumentParser()

# Define arguments
parser$add_argument("-d", "--data_dir", help="Path to directory with Mergers Report")
parser$add_argument("-o", "--out_dir", help="Path to output dir")
parser$add_argument("-p", "--path_to_flist", help="Path to file with list of samples and barcode pairs") 
parser$add_argument("-m", "--minreads_threshold", help="Minimum threshold to flag wells for low productivity")
parser$add_argument("-c", "--contamination_threshold", help="Percentage of reads of foreign origin necessary to classify a well as contaminated")
parser$add_argument("-mf", "--missing_files", help="List of missing files")

# Parse the command-line arguments
args <- parser$parse_args()
print(args)

# Assign variables based on command-line arguments
render("/Code/ci_report_layouting.Rmd", params = list(
  data_dir = args$data_dir,
  out_dir = args$out_dir,
  path_to_flist = args$path_to_flist,
  minreads_threshold = args$minreads_threshold,
  contamination_threshold = args$contamination_threshold,
  missing_files = args$missing_files
),
output_dir = "Results/")

print("Leaving render script")
