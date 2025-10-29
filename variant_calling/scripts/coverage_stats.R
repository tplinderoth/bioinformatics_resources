#!/usr/bin/env Rscript

# avg_depth.R <sample ID list file> <list of coverage files> <outfile name>
# calculates average depth of coverage

# sample ID and coverage files must be in the same order

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
	cat("\navg_depth.R <sample ID file> <list of coverage files> <outfile name>

---Input descriptions---
Sample ID file: file containing one ID per row
Coverage file list: file listing the names of files produced by SAMtools coverage, one file per row
Outfile name: Name of output file

---Notes---
Sample IDs and coverage files must be listed in the same order in their respective input files.\n\n")
	options(show.error.messages=FALSE)
	stop()
}


idlist <- as.character(args[1])
covlist <- as.character(args[2])
outfile <- as.character(args[3])

id.vec <- read.table(idlist, head=FALSE)$V1
covf.vec <- read.table(covlist, head=FALSE)$V1

out.df <- NULL

for (i in 1:length(id.vec)) {
	id <- id.vec[i]
	cov.df <- read.table(covf.vec[i],head=TRUE,comment.char="")
	cov.df$scaff_len=(cov.df$endpos-cov.df$startpos + 1)
	perc_covered_bp = round(100*(sum(cov.df$covbases)/sum(cov.df$scaff_len)),2)
	wts <- cov.df$scaff_len/sum(cov.df$scaff_len)
	depth_wtavg <- sum(wts * cov.df$meandepth) # weighted average mean depth with weights based on scaffold length
	depth_var = sum(wts*(cov.df$meandepth - depth_wtavg)^2)
	depth_sd = sqrt(depth_var)
	median_avgdepth = median(cov.df$meandepth) # this is not that useful because it doesn't account for differences in scaffold length
	out.df <- rbind(out.df, data.frame(ID = id, PERCENT_COVERED = perc_covered_bp, AVG_DEPTH = round(depth_wtavg,2), SD_DEPTH = round(depth_sd,2), 
	MEDIAN_SCAFFOLD_AVG_DEPTH = round(median_avgdepth,2)))
}

write.table(out.df, file=outfile, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
