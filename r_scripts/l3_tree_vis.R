

rm(list=ls()) 

setwd("~/Documents/spolpred/")

library(ape)
library(phytools)
library(ggplot2)
library(ggtree)
library(dplyr)
library(colorspace)
library(scales)

# Functions ----

source("https://raw.githubusercontent.com/GaryNapier/Packages_functions/master/Functions.R")

get_binary_cols <- function(x, col, prefix){
  
  # Data with column of binary data like this:
  # sample      spoligotype
  # ERR038262   1110000000111111111111000000000000011111000
  # ERR038263   1110000000111111111111000000000000011111000
  # SRR2100087  1110000111111111111111000000000000111111111
  # ERR072072   1110000111111111111111000000000000000111111
  
  # Split out the binary col and add in cols like this:
  # sample                                    spoligotype s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 ...etc
  # ERR038262 1110000000111111111111000000000000011111000  1  1  1  0  0  0  0  0  0   0 ...etc
  # ERR038263 1110000000111111111111000000000000011111000  1  1  1  0  0  0  0  0  0   0 ...etc
  # SRR2100087 1110000111111111111111000000000000111111111  1  1  1  0  0  0  0  1  1   1 ...etc
  # ERR072072 1110000111111111111111000000000000000111111  1  1  1  0  0  0  0  1  1   1 ...etc
  # SRR2100127 1110000111111111111111000000000000000111111  1  1  1  0  0  0  0  1  1   1 ...etc
  # SRR2100292 1110000110111111111000000000000000011111111  1  1  1  0  0  0  0  1  1   0 ...etc
  # SRR2100324 1110000111111111111111000000000000011111111  1  1  1  0  0  0  0  1  1   1 ...etc
  # SRR2100410 1100000111111111111000000000000000011111111  1  1  0  0  0  0  0  1  1   1 ...etc
  
  library(stringr)
  
  bin_split <- str_split(x[, col], "")
  bin_mat <- data.frame(do.call("rbind", bin_split))
  names(bin_mat) <- paste0(prefix, 1:ncol(bin_mat))
  cbind(x, bin_mat)
}

ggtree_strip_cols <- function(x, alpha = 0.7){
  
  # Dataframe:
  #            lin_level_1
  # SRR8651556           2
  # SRR8651557           4
  # SRR8651558           2
  # SRR8651559           2
  # SRR8651560           2
  # SRR8651561           4
  
  # Output:
  # 1           2           3           4           5           6           7           9           M.bovis     M.caprae    M.orygis 
  # "#FF0000E6" "#FF8B00E6" "#E8FF00E6" "#5DFF00E6" "#00FF2EE6" "#00FFB9E6" "#00B9FFE6" "#002EFFE6" "#5D00FFE6" "#E800FFE6" "#FF008BE6" 
  
  # One distinct colour per unique value in the dataframe. Colour names are the unique values.
  # NOTE: Need to call the function with the column specified (i.e. read in as a vector):
  # ggtree_strip_cols(df$<col_name>)
  
  cols <- rainbow(length(unique(x)), alpha = alpha)
  names(cols) <- c(sort(unique(x)))
  cols
}

dark_to_light <- function(lin_names, first_col){
  if(length(lin_names) == 1){
    return(first_col)
  }else{
    lighten_to_col <- lighten(first_col, 1-(1/(length(lin_names))))
    fc <- colorRampPalette(c(first_col, lighten_to_col))
    lin_all_cols <- fc(length(lin_names))
    scales::show_col(lin_all_cols)
    names(lin_all_cols) <- lin_names
    return(lin_all_cols)
  }
}


# Paths ----

results_path <- "results/"
newick_path <- "results/newick/"
data_path <- "data/"
itol_templates_path <- "../pipeline/itol_templates/"

# Files ----

# spoligo_lineage_full_file <- paste0(data_path, "spoligo_lineage.full.txt")
spoligo_lineage_full_file <- paste0(data_path, "spoligo_lineage.SNPs.csv")
l3_treefile <- paste0(newick_path, "grouped_samples_lineage3.treefile")
l4_treefile <- paste0(newick_path, "grouped_samples_lineage4.treefile")
l1_7_treefile <- paste0(newick_path, "grouped_samples_L_1_7.treefile")
# l3_treefile <- paste0(newick_path, "lin_3_all_samps.treefile")
itol_binary_template_file <- paste0(itol_templates_path, "itol.binary.txt")
itol_outfile <- paste0(results_path, "itol_all_spol_binary.txt")

# Read in data ----
# spoligo_lineage_full <- read.table(spoligo_lineage_full_file, header = T, sep = "\t", 
#                                    colClasses = c("spoligotype" = "character"))



spoligo_lineage_full <- read.csv(spoligo_lineage_full_file, header = T, 
                                   colClasses = c("spoligotype" = "character"))
l3_tree <- midpoint.root(read.tree(l3_treefile))
l4_tree <- midpoint.root(read.tree(l4_treefile))
l1_7_tree <- midpoint.root(read.tree(l1_7_treefile))
itol_binary_template <- readChar(itol_binary_template_file, nchars = 1e6)


# Clean ----

# l3_samples <- l3_tree$tip.label

# Blank lineages - all 1.2.2.1 apparently
spoligo_lineage_full$lineage <- ifelse(spoligo_lineage_full$lineage == "", "1.2.2.1", spoligo_lineage_full$lineage)

# Remove "lineage"
spoligo_lineage_full$lineage <- gsub("lineage", "", spoligo_lineage_full$lineage)

# Remove and store animal strains
animal <- subset(spoligo_lineage_full, grepl("M", lineage))
spoligo_lineage_full <- subset(spoligo_lineage_full, !(grepl("M", lineage)))

# Expand lineages, rbind animal back in and merge
# exp <- expand_hierarchy(spoligo_lineage_full, "sample", "lineage")
exp <- expand_hierarchy_fill(spoligo_lineage_full, "sample", "lineage")
data <- merge(rbind(spoligo_lineage_full, animal), exp,
              by.x = "sample", by.y = "id", all.x = T, 
              sort = F)

# Fill in cols lin_level_1 and max_lin for animal
data$lin_level_1 <- ifelse(grepl("M", data$lineage), data$lineage, data$lin_level_1)
data$lin_level_2 <- ifelse(grepl("M", data$lineage), data$lineage, data$lin_level_2)
data$lin_level_3 <- ifelse(grepl("M", data$lineage), data$lineage, data$lin_level_3)
data$lin_level_4 <- ifelse(grepl("M", data$lineage), data$lineage, data$lin_level_4)
data$lin_level_5 <- ifelse(grepl("M", data$lineage), data$lineage, data$lin_level_5)
data$max_lin <- ifelse(grepl("M", data$lineage), data$lineage, data$max_lin)

# Remove 0000000000000000000000000000000000000000000 - error apparently
data <- subset(data, !(spoligotype == "0000000000000000000000000000000000000000000"))

# Add rownames for plotting
rownames(data) <- data$sample

# Add binary cols - convert binary spoligotype 00101 etc to columns - for both itol and ggtree
spol_data <- get_binary_cols(data, "spoligotype", "s")

# Clean up trees - only include samples which are in the spol dataset
l3_tree <- ape::keep.tip(l3_tree, intersect(data$sample, l3_tree$tip.label))
l4_tree <- keep.tip(l4_tree, intersect(data$sample, l4_tree$tip.label))
l1_7_tree <- keep.tip(l1_7_tree, intersect(data$sample, l1_7_tree$tip.label))

# ggtree ----

# subset data by relevant cols - i.e. remove lineage and original spoligotype cols for ggplot
lin_lv1_data_ggplot <- select(spol_data, lin_level_1)
lin_lv2_data_ggplot <- select(spol_data, lin_level_2)
spol_data_ggplot <- select(spol_data, -(sample:max_lin))

# Colours

lin_colours_lv_1 <- ggtree_strip_cols(lin_lv1_data_ggplot$lin_level_1)

# Split out data by lineage level 1
spol_data_lv1_split <- split(spol_data, spol_data$lin_level_1)
# Add 'lineage' for easier subsetting
names(spol_data_lv1_split) <- paste0("lineage", names(spol_data_lv1_split))

# Get colours fro each level 2 (4.1, 4.2 etc) - lin 3 and 4 only
l3_lv2_cols <- ggtree_strip_cols(spol_data_lv1_split$lineage3$lin_level_2)
l4_lv2_cols <- ggtree_strip_cols(spol_data_lv1_split$lineage4$lin_level_2)

# Level 2 colours - successively lighten the top lin level colour proportional to the number of sublins
# Only relevant to the L1&7 and L5, 6, etc trees - level 2 for lin 3 and 4 are defined above
lv2_cols_list <- list()
for(lin in names(spol_data_lv1_split)){
  
  # lin_1_7_lv1_cols <- lin_colours_lv_1[c("1", "7")]
  lin_data <- spol_data_lv1_split[[lin]]
  lin_lv1 <- unique(lin_data$lin_level_1)
  lin_top_col <- lin_colours_lv_1[lin_lv1]
  # 7 
  # "#00B9FFB3" 
  lins_lv2 <- sort(unique(lin_data$lin_level_2))
  n_lins_lv2 <- length(lins_lv2)
  
  # Get lin level 2 colours
  # If the first of the level 2 colours is the same as the main lin (top level colour), then keep it the same
  # If not, then start from a lighter colour 
  if(lins_lv2[1] == names(lin_top_col)[1]){
    # 7 
    # "#00B9FF" 
    lv2_cols_list[[lin]] <- dark_to_light(lins_lv2, lin_top_col)
  }else{
    lv2_cols_list[[lin]] <- dark_to_light(lins_lv2, lighten(lin_top_col, 1/n_lins_lv2))
  }
}

# L1 & L7 

lin_1_7_lv2_cols <- c(lv2_cols_list$lineage1, lv2_cols_list$lineage7)

# Plot ----

do_tree <- "1_7"

# Set up ggtree parameters 
width <- 0.05
font_sz <- 3
line_sz <- 0.1
angle <- 30
os <- 0.0010

# L3

  if(do_tree == "3"){
  l3_ggtree <- ggtree(l3_tree, size = line_sz, layout = "circular")
  l3_ggtree <- gheatmap(l3_ggtree, lin_lv2_data_ggplot,
                        color = NA, 
                        width = width,
                        offset = 0,
                        colnames = F)+
    scale_fill_manual(values = l3_lv2_cols, breaks = names(l3_lv2_cols) )+
    labs(fill = "Subineage")
  
  l3_ggtree <- l3_ggtree + ggnewscale::new_scale_fill()
  
  l3_ggtree <- gheatmap(l3_ggtree, spol_data_ggplot,
           offset = os+0.0020,
           color = NA,
           low="white",
           high="black",
           colnames = F) +
    scale_fill_manual(values=c("white", "black"), labels = c("0", "1", "NA"), na.value = "grey")+
    labs(fill = "Spoligotype")+
    ggtitle("Lineage 3")
  
  # ggsave(file = paste0(results_path, "lin3_ggtree.png"), plot = l3_ggtree)
  # ggsave(file = paste0(results_path, "lin3_ggtree.svg"), plot = l3_ggtree)

}else if(do_tree == "4"){

  # L4
  
  # Subset L4, still too big
  set.seed(123)
  l4_tree <- keep.tip(l4_tree, sample(l4_tree$tip.label, floor(length(l4_tree$tip.label)*0.5)))
  
  l4_ggtree <- ggtree(l4_tree, size = line_sz, layout = "circular")
  l4_ggtree <- gheatmap(l4_ggtree, lin_lv2_data_ggplot,
                        color = NA,
                        width = 0.2,
                        offset = 0,
                        colnames = F)+
    scale_fill_manual(values = l4_lv2_cols, breaks = names(l4_lv2_cols) )+
    labs(fill = "Subineage")
  
  l4_ggtree <- l4_ggtree + ggnewscale::new_scale_fill()
  
  l4_ggtree <- gheatmap(l4_ggtree, spol_data_ggplot,
                        offset = os+0.0020,
                        color = NA,
                        low="white",
                        high="black",
                        colnames = F) +
    scale_fill_manual(values=c("white", "black"), labels = c("0", "1", "NA"), na.value = "grey")+
    labs(fill = "Spoligotype")+
    ggtitle("Lineage 4")
  
  # ggsave(file = paste0(results_path, "lin4_ggtree.png"), plot = l4_ggtree)
  # ggsave(file = paste0(results_path, "lin4_ggtree.svg"), plot = l4_ggtree)

}else if(do_tree == "1_7"){

# L1 & L7

  l1_7_ggtree <- ggtree(l1_7_tree, size = line_sz, layout = "circular")
  
  l1_7_ggtree <- gheatmap(l1_7_ggtree, lin_lv1_data_ggplot,
                          color = NA,
                          width = 0.1,
                          offset = 0,
                          colnames = F)+
    scale_fill_manual(values = lin_colours_lv_1[c("1", "7")], 
                      breaks = names(lin_colours_lv_1[c("1", "7")]), 
                      name = "Lineage")
  
  l1_7_ggtree <- l1_7_ggtree + ggnewscale::new_scale_fill()
  
  l1_7_ggtree <- gheatmap(l1_7_ggtree, lin_lv2_data_ggplot,
                          color = NA,
                          width = 0.1,
                          offset = os+0.004,
                          colnames = F)+
    scale_fill_manual(values = lin_1_7_lv2_cols,
                      breaks = names(lin_1_7_lv2_cols), 
                      name = "Sublineage")
  
  l1_7_ggtree <- l1_7_ggtree + ggnewscale::new_scale_fill()
  
  l1_7_ggtree <- gheatmap(l1_7_ggtree, spol_data_ggplot,
                          offset = os+(0.004*3),
                          color = NA,
                          low="white",
                          high="black",
                          colnames = F)+
    scale_fill_manual(values=c("white", "black"), labels = c("0", "1", "NA"), na.value = "grey")+
    labs(fill = "Spoligotype")+
    ggtitle("Lineages 1 & 7")
  
  ggsave(file = paste0(results_path, "lin1_7_ggtree.png"), plot = l1_7_ggtree)
  
  # ggsave(plot = l1_7_ggtree, filename = paste0(results_path, "lin1_7_ggtree.svg"))

}









# itol ----

# Subset for itol
itol_data <- select(spol_data, -(spoligotype:max_lin), sample)
len_spol <- len_str("1110000000111111111111000000000000011111000")

# FIELD_SHAPES,1,1,1
# FIELD_LABELS,f1,f2,f3
# FIELD_COLORS,#000000,#000000,#000000

field_shapes <- paste0("FIELD_SHAPES\t", paste0(paste0(rep("1", len_spol), "\t"), collapse = ""))
field_labels <- paste0("FIELD_LABELS\t", paste0(paste0(paste0("s", 1:len_spol), "\t"), collapse = ""))
field_colours <- paste0("FIELD_COLORS\t", paste0(paste0(rep("#000000", len_spol), "\t"), collapse = ""))

itol_binary_template <- gsub("#FIELD_SHAPES,2,4,5,1", field_shapes, itol_binary_template)
itol_binary_template <- gsub("#FIELD_LABELS,f1,f2,f3,f4", field_labels, itol_binary_template)
itol_binary_template <- gsub("#FIELD_COLORS,#ff0000,#00ff00,#ffff00,#0000ff", field_colours, itol_binary_template)

# Write the template out under the new file name for the study accession
write.table(itol_binary_template, file = itol_outfile, sep="\t",
            row.names=F, col.names=F, quote = F)

# Append the drug resistance data to the template
write.table(itol_data, file = itol_outfile,
            append = T, sep="\t",
            row.names=F, col.names=F, quote = F)



















