#!/usr/bin/env Rscript

# avg_depth.R <sample ID list file> <list of coverage files> <outfile name>
# calculates average depth over covered regions

# sample ID and coverage files must be in the same order

args <- commandArgs(trailingOnly=TRUE)
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
	wts <- cov.df$covbases/sum(cov.df$covbases)
	#wts <- covered.length/sum(covered.length)
	depth_wtavg <- sum(wts * cov.df$meandepth) # coverage-weighted average of mean depths
	depth_var = sum(wts*(cov.df$meandepth - depth_wtavg)^2)
	depth_sd = sqrt(depth_var)
	median_avgdepth = median(cov.df$meandepth) # this is not that useful because it doesn't account for differences in scaffold length or number covered bases
	out.df <- rbind(out.df, data.frame(ID = id, PERCENT_COVERED = perc_covered_bp, COVERED_REGION_AVG_DEPTH = round(depth_wtavg,2), COVERED_REGION_SD_DEPTH = round(depth_sd,2), 
	AVG_DEPTH_MEDIAN = round(median_avgdepth,2)))
}

write.table(out.df, file=outfile, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
