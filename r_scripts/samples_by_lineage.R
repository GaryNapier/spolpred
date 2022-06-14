#!/usr/bin/env Rscript

source("https://raw.githubusercontent.com/GaryNapier/Packages_functions/master/Functions.R")

library(optparse)

TESTING <- 0

if (!TESTING){

  # Arguments ----

  option_list = list(
    make_option(c("-n", "--sample_n"), type="character", default=10,
                help="max number of samples in each sublineage", metavar="character"),
    make_option(c("-l", "--lineage_file"), type="character", default=NULL,
                help="metadata of all samples and their lineage data", metavar="character"),
    make_option(c("-s", "--all_lineages_samples_outfile"), type="character", default=NULL,
                help="", metavar="character"),
    make_option(c("-g", "--grouped_lins_prefix"), type="character", default="results/grouped_samples_",
                help="", metavar="character")
  );
  
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);

  print("ARGUMENTS:")
  print(opt)
  print("---")
  print(str(opt))

  # Setup

  # Files ----

  lineage_file <- opt$lineage_file
  sample_n <- opt$sample_n
  all_lineages_samples_outfile <- opt$all_lineages_samples_outfile
  grouped_lins_prefix <- opt$grouped_lins_prefix

}else{
  
  sample_n <- 10
  
  # in
  setwd("~/Documents/spolpred/")
  results_path <- "results/"
  lineage_file <- paste0(results_path, "lineage_file.csv")
  # out
  all_lineages_samples_file <- paste0(results_path, "all_lineages_samples.txt")
  # all_lineages_metadata_file <- paste0(results_path, "all_lineages_metadata.txt")
  grouped_lins_prefix <- paste0(results_path, "grouped_samples_")
  
}

# Read in data
lineage_data <- read.csv(lineage_file)

# Clean
# lineage_data <- subset(lineage_data, !(main_lin == ""))
lineage_data$main_lin <- ifelse(lineage_data$main_lin == "", "lineage8", lineage_data$main_lin)

# Wrangle

# Group lineages
# 5, 6, 8, 9, animal
# 1, 7
# 2
# 3
# 4

lineage_data$grouped_lin <- ifelse(lineage_data$main_lin %in% c("lineage5", "lineage6", "lineage8", "lineage9", "La1", "La2", "La3"), "L_5_6_8_9_La", 
       ifelse(lineage_data$main_lin %in% c("lineage1", "lineage7"), "L_1_7", lineage_data$main_lin))

lineage_data_split <- split(lineage_data, lineage_data$grouped_lin)

# Go through main lineages and subset by taking x number of samples from each sublin
final_lineage_data <- list()
for(i in seq(lineage_data_split)){
  
  x <- lineage_data_split[[i]]
  x_split <- split(x, x$sublin)
  subset_list <- list()
  for(j in seq(x_split)){
    
    n <- nrow(x_split[[j]])
    
    print(names(lineage_data_split[i]))
    print(paste0("nrow = ", as.character(n)))
    
    if(n < sample_n){
      
      print("adding without sampling")
      
      subset_list[[j]] <- x_split[[j]]
    }else{
      
      print("adding with sampling")
      
      set.seed(123)
      subset_list[[j]] <- x_split[[j]][sample(1:n, sample_n), ]
    }
  }
  final_lineage_data[[i]] <- do.call('rbind', subset_list)

}


# Write groups separately
for(i in seq(final_lineage_data)){
  lin_group <- unique(final_lineage_data[[i]]$grouped_lin)
  file <- paste0(grouped_lins_prefix, lin_group, ".txt")
  write.table(final_lineage_data[[i]]$id, file = file, quote = F, row.names = F, col.names = F)
}

final_lineage_data <- do.call('rbind', final_lineage_data)

# Write all samples...
write.table(final_lineage_data$id, file = all_lineages_samples_file, 
            quote = F, row.names = F, col.names = F)
# # ...and metadata
# write.table(final_lineage_data, file = all_lineages_metadata_file, 
#             quote = F, row.names = F, col.names = F)



















