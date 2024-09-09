library(dplyr)
library(ggplot2)

bindev <- read.csv("findBiasedFeatures/processed-data/hpc_bindev_default-sample_svgs-only.csv")
bindev$d.diff_sample = (bindev$dev_default-bindev$dev_sample)/bindev$dev_sample
bindev$r.diff_sample = bindev$rank_sample-bindev$rank_default

mean1 = mean(bindev$d.diff_sample)
sd1 = sd(bindev$d.diff_sample)
bindev$nSD_dev = (bindev$d.diff_sample-mean1)/sd1
summary(bindev$nSD_dev)

mean2 = mean(bindev$r.diff_sample)
sd2 = sd(bindev$r.diff_sample)
bindev$nSD_rank = (bindev$r.diff_sample-mean2)/sd2
summary(bindev$nSD_rank)

sd.interval = 5
bindev$nSD.bin_dev = cut(abs(bindev$nSD_dev), right=FALSE,
                         breaks=seq(0,max(bindev$nSD_dev)+sd.interval, by=sd.interval),
                         include.lowest=TRUE)
bindev$nSD.bin_rank = cut(abs(bindev$nSD_rank), right=FALSE,
                          breaks=seq(0,max(bindev$nSD_rank)+sd.interval, by=sd.interval),
                          include.lowest=TRUE)

bias = filter(bindev, nSD_dev>=10 | nSD_rank>=5)$gene
names(bias) = filter(bindev, nSD_dev>=10 | nSD_rank>=5)$gene_name

saveRDS(bias, "findBiasedFeatures/processed-data/hpc_biased-features.rda")
