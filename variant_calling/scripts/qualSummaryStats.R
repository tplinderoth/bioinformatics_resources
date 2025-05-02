#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
        cat("\nqualSummaryStats.R <list of files with quality values (one file per row)> <outfile name> [column numbers to analyze]\n\n")
        options(show.error.messages=FALSE)
        stop()
}

flist <- as.character(read.table(args[1], head=FALSE)$V1)
outname <- args[2]
subidx <- NULL
if (length(args) > 1) subidx = as.numeric(args[3:length(args)])

df <- NULL
fheader <- read.table(flist[1],head=TRUE,nrows=1)
classes <- NULL
if (is.null(subidx)) classes <- "character" else classes <- sapply(1:ncol(fheader),function(x,idx){ifelse(x %in% idx,"character","NULL")},idx=subidx)

for (infile in flist) {
	cat(paste0("Merging ", infile, "\n"))
	df <- rbind(df, read.table(infile, head=TRUE, na.strings=".", colClasses=classes))
}

df.stats <- data.frame(stat=colnames(df))
df.stats$mean = sapply(1:ncol(df),function(x,val){mean(as.numeric(val[,x]), na.rm=TRUE)}, val=df)
df.stats$median = sapply(1:ncol(df),function(x,val){median(as.numeric(val[,x]), na.rm=TRUE)}, val=df)
df.stats$var = sapply(1:ncol(df),function(x,val){var(as.numeric(val[,x]), na.rm=TRUE)}, val=df)
df.stats$percentile_1 = sapply(1:ncol(df),function(x,val){quantile(as.numeric(val[,x]), 0.01, na.rm=TRUE)}, val=df)
df.stats$percentile_2 = sapply(1:ncol(df),function(x,val){quantile(as.numeric(val[,x]), 0.02, na.rm=TRUE)}, val=df)
df.stats$percentile_3 = sapply(1:ncol(df),function(x,val){quantile(as.numeric(val[,x]), 0.03, na.rm=TRUE)}, val=df)
df.stats$percentile_5 = sapply(1:ncol(df),function(x,val){quantile(as.numeric(val[,x]), 0.05, na.rm=TRUE)}, val=df)
df.stats$percentile_95 = sapply(1:ncol(df),function(x,val){quantile(as.numeric(val[,x]), 0.95, na.rm=TRUE)}, val=df)
df.stats$percentile_97 = sapply(1:ncol(df),function(x,val){quantile(as.numeric(val[,x]), 0.97, na.rm=TRUE)}, val=df)
df.stats$percentile_98 = sapply(1:ncol(df),function(x,val){quantile(as.numeric(val[,x]), 0.98, na.rm=TRUE)}, val=df)
df.stats$percentile_99 = sapply(1:ncol(df),function(x,val){quantile(as.numeric(val[,x]), 0.99, na.rm=TRUE)}, val=df)

write.table(df.stats,file=outname,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

cat("Finished\n")
