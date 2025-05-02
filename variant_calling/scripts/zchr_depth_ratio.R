# zchr_depth_ratio.R

homolog_file = '/mnt/research/Fitz_Lab/ref/bird/FSJ_V3/FSJV3_convert.short.txt'
homolog = read.table(homolog_file,head=TRUE)
autosome = homolog$FSJV3_SCAFFOLD[which(homolog$ZEBRA_FINCH_CHR != "Z")]
zchr = homolog$FSJV3_SCAFFOLD[which(homolog$ZEBRA_FINCH_CHR == "Z")]

samples = read.table('/mnt/research/Fitz_Lab/projects/mosaic/mosiac_ids.txt',head=FALSE)

df = NULL
for (i in 1:nrow(samples)) {
	id = samples$V1[i]
	coverage = read.table(paste0('/mnt/research/Fitz_Lab/projects/mosaic/map/coverage/',id,'_coverage.txt'),head=FALSE)
	autoset = coverage[which(coverage$V1 %in% autosome),]
	zset = coverage[which(coverage$V1 %in% zchr),]
	len_total = sum(autoset$V3)
	autocov = sum(autoset$V7 * autoset$V3/len_total) # chromosome length weighted average
	zcov = zcov = zset$V7
	dratio = zcov/autocov # ratio of Z chr coverage to autosome coverage
	df = rbind(df, data.frame(SAMPLE = id, AUTOSOME_DEPTH = autocov, Z_DEPTH = zcov, DEPTH_RATIO = dratio))
}

# assign sex based on coverage ratio

df$NGS_SEX = NA
df$NGS_SEX[which(df$DEPTH_RATIO < 0.6)] = "Female"
df$NGS_SEX[which(is.na(df$NGS_SEX) == TRUE)] = "Male"

# write output
# write.table(df,file='/mnt/research/Fitz_Lab/projects/mosaic/map/coverage/ngs_sex_assignment.txt',col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
