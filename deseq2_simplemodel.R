# add libraries
suppressPackageStartupMessages({
  library("DESeq2")
  library(readr)
  library(optparse)
  library("ggplot2")

})

option_list = list(
  make_option(c("-c", "--count"), default=NA, type='character',
              help="count table (tsv)"),
  make_option(c("-m", "--metadata"), default=NA, type='character',
              help="metadata (tsv). Required columns are 'name' and 'condition'. First condition is used as the base level (aka control)."),
  make_option(c("-o", "--output"), default=NA, type='character',
              help="path to output files.")
)
opt = parse_args(OptionParser(option_list=option_list))


# Read count matrix
countdata <- read.table(opt$count, header=TRUE, check.names = FALSE)

# Read smaplesheet , detects conditions, generate the formula
sampleInfo <- read.table(opt$metadata, header=TRUE, check.names = FALSE, stringsAsFactor = F)
sampleInfo$condition <- as.factor(sampleInfo$condition)
sampleInfo$condition <- relevel(sampleInfo$condition, ref = as.character(sampleInfo$condition[1]))

design <- ~ 1 + condition

# build the matrix
dds <- DESeq2::DESeqDataSetFromMatrix(countData = countdata[,sampleInfo$name], colData = sampleInfo, design = design)
dds <- estimateSizeFactors(dds)
dds <- DESeq2::DESeq(dds)

#normalised count
normalized_counts <- counts(dds, normalized=T)
write.table(normalized_counts, file=paste0(opt$output,'normalized_counts.tsv',sep =""), quote=FALSE, sep='\t', row.names=TRUE)

#compute ddr
cond_name <- DESeq2::resultsNames(dds)[2]
ddr <- DESeq2::results(dds, alpha = 0.05, name = cond_name)

# compute shrunk result
ddr_shrunk <- DESeq2::lfcShrink(dds,coef = cond_name,type="apeglm",res=ddr)

# save shrunk result
kd <- unique(sampleInfo$condition)[2]
out_filename <- paste(opt$output, "/", kd, "_vs_ctrl_shrunk.tsv", sep="")
write.table(ddr_shrunk, file=out_filename, quote=FALSE, sep='\t')

# pca
vsd <- DESeq2::vst(dds, blind=FALSE)
out_pca_filename <- paste(opt$output, "/", kd, "_PCA.png", sep="")
ggplt = DESeq2::plotPCA(vsd, intgroup=c("condition", "name"))
ggsave(out_pca_filename, ggplt)
