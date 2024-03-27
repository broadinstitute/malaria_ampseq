#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
data_dir = args[1]
out_dir =args[2]

if (!require("viridis")) {
  install.packages("viridis", repos="http://cran.rstudio.com/")
  library("viridis")
}

###################################
###BB MERGE PERFORMANCE REPORT ####
###################################

mergedata <- read.csv(file.path(data_dir, "bbmergefields.tsv"), sep = "\t", header = TRUE)

mergedata_join = subset(mergedata, select = c(3, 5, 7))
rownames(mergedata_join) = mergedata$SampleID
color_vector = viridis(nrow(t(as.matrix(mergedata_join))), option = "D")
 
mergedata_join <- mergedata_join[, order(colSums(mergedata_join), 
                                         decreasing = TRUE)]
mergedata_join <- mergedata_join[order(-mergedata_join[,1], 
                                       -mergedata_join[,2], 
                                       -mergedata_join[,3]),]
 
mergedata_per = subset(mergedata, select = c(4, 6, 8))
rownames(mergedata_per) = mergedata$SampleID
 
mergedata_per <- mergedata_per[, order(colSums(mergedata_per), 
                                         decreasing = TRUE)]
mergedata_per <- mergedata_per[order(-mergedata_per[, 1], 
                                        -mergedata_per[, 2], 
                                        -mergedata_per[, 3]),]
 
#Add discarded count AFTER the columns have been organized
rownames(mergedata) = mergedata$SampleID
mergedata = mergedata[rownames(mergedata_join), ]
mergedata_join$Discarded = mergedata$Pairs - mergedata$Joined - mergedata$Ambiguous - mergedata$No_Solution
mergedata_per$DiscardedP = 100 - mergedata_per$JoinedP - mergedata_per$AmbiguousP - mergedata_per$No_SolutionP

col_order = c("Joined", "Ambiguous", "No_Solution", "Discarded")
col_order_p = c("JoinedP", "AmbiguousP", "No_SolutionP", "DiscardedP")
 
mergedata_join = mergedata_join[, col_order]
mergedata_per = mergedata_per[, col_order_p]
 
svg(file.path(out_dir, "BBmerge_performance_absolute_report.svg"), width = 12)
 	par(mar = c(14, 4, 6, 2), xpd = TRUE) # increase the bottom margin
 	barplot(
 	t(as.matrix(mergedata_join)),
 	col = c(color_vector, "red"),
 	border = NA,
 	ylab = "Reads Count",
 	xlab = "",
 	main = "BBMerge performance - Absolute",
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
 	legend = c("Joined", "Ambiguous", "No Solution", "Discarded"),
 	fill = c(color_vector, "red"),
 	bty = "n",
 	cex = 1
)
 
dev.off()
 
svg(file.path(out_dir, "BBmerge_performance_percentage_report.svg"), width = 12)
 	par(mar = c(14, 4, 6, 2), xpd = TRUE) # increase the bottom margin
 
 	barplot(
 	t(as.matrix(mergedata_per)),
 	col = c(color_vector, "red"),
 	border = NA,
 	ylab = "Reads Percentage",
 	xlab = "",
 	main = "BBMerge performance - Percentage",
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
 	legend = c("Joined", "Ambiguous", "No Solution", "Discarded"),
 	fill = c(color_vector, "red"),
 	bty = "n",
 	cex = 1
 	)
 
dev.off()
 
mergedata_join <- mergedata_join[order(-mergedata_join[, 4]), ] 
 
svg(file.path(out_dir, "BBmerge_performace_absolute_discarded.svg"), width = 12)
 	par(mar = c(14, 4, 6, 2), xpd = TRUE) # increase the bottom margin
 	barplot(
 	t(as.matrix(mergedata_join)[,4]),
 	col = "red",
 	border = NA,
 	ylab = "Reads Count",
 	xlab = "",
 	main = "BBMerge performance - Discarded",
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
 	legend = c("Discarded"),
 	fill = "red",
 	bty = "n",
 	cex = 1
	)
 
dev.off()
 
