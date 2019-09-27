library(tximport)
library(readr)

ttg <- read.csv(snakemake@input[[1]])
samples <- read.csv(snakemake@input[[2]], sep="\t")

r <- seq(3, length(snakemake@input)+1)
files <- snakemake@input[r]
files <- unlist(files)
names(files) <- samples[, "sample"]

all(file.exists(files))

txi.kallisto <- tximport(files, type="kallisto", tx2gene = ttg, ignoreAfterBar = TRUE) # txOut means that we are outputting transcripts
# set  if can't match between tx2gene and files
df <- as.data.frame(txi.kallisto$abundance)

df["Gene"] <- rownames(df)

df <- df[,c(ncol(df),1:ncol(df)-1)]

rownames(df) <- NULL

write.table(df, file=snakemake@output[[1]], row.names=FALSE, sep="\t")
