
library(optparse)
library(yaml)
library(rtracklayer)

option_list = option_list = list(
  make_option(c("--bigWig"), type = "character", help = "Path to the bigWig file to normalize"),
  make_option(c("--MS_yaml"), type = "character", help = "Path to the yaml file with MS norm factors"),
  make_option(c("--removeSuffix"), type = "character", help = "Suffix to remove from the bw name to match MS yaml file"),
  make_option(c("--exportFolder"), type = "character", help = "Folder to export MS normed bigWig files")
)

opt <- parse_args(OptionParser(usage = "Usage: %prog [options]", option_list = option_list))

yaml_data <- yaml.load_file(opt$MS_yaml)

ms_norm <- function(bw, MS) {
  bw <- import.bw(bw)
  N <- as.numeric(length(bw))
  sum_bw <- as.numeric(sum(bw$score))
  norm_factor <- MS * N / (sum_bw)
  bw$score <- bw$score * norm_factor
  filt <- bw$score > 100
  high_signal_bins <- bw[filt]
  print(paste0("percent > 100 = ", as.integer(length(high_signal_bins) / N * 100), " %"))
  print(quantile(bw$score))
  print(mean(bw$score))
  bw[filt]$score <- rep(100, length(bw[filt]$score)) # Anhadido nuevo
  print("After")
  print(bw[filt])
  # bw <- bw[!filt] # Esto estaba antes, se eliminaban los bins con mas senhal
  print(paste0("Quantiles after filtering:",quantile(bw$score)))
  export_name <- paste0(opt$exportFolder, "/", gsub(opt$removeSuffix, "", basename(opt$bigWig)), ".MSnorm.bigWig")
  print(paste("Exporting name: ", export_name))
  return(export.bw(bw, export_name))
}

print(paste("To search: ", gsub(opt$removeSuffix, "", basename(opt$bigWig))))
ms_norm(bw=opt$bigWig, MS=yaml_data[[gsub(opt$removeSuffix, "", basename(opt$bigWig))]])
