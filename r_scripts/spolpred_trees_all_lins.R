

setwd("~/Documents/spolpred/")

library(ggplot2)
library(ggtree)
library(phytools) # Midpoint root
library(ggnewscale) # For adding new heatmaps (annotations)
library(plyr)
library(dplyr)
library(stringr)

# Functions

get_last <- function(x, split_on = ";"){
  unlist(lapply(strsplit(x, split_on), function(x){
    x[length(x)]
  }))
}

len_str <- function(string){
  length(unlist(strsplit(string, split = "")))
}

expand_hierarchy <- function(df, group_by_col_name, hierarchy_to_expand_col_name){
  # Takes df like this:
  #   ID      Group
  # 1 samp_1  4.2.1.1
  # 2 samp_2  1.2.1.2.1
  
  # And makes this:
  #   ID    lin_level_1 lin_level_2 lin_level_3 lin_level_4 lin_level_5   max_lin
  # 1 samp_1          4         4.2       4.2.1     4.2.1.1        <NA>   4.2.1.1
  # 2 samp_2          1         1.2       1.2.1     1.2.1.2   1.2.1.2.1 1.2.1.2.1
  
  split_lins <- str_split(df[[hierarchy_to_expand_col_name]], "\\.")
  max_lin_len <- max(sapply(split_lins, length))
  mat <- matrix(nrow = length(df[[group_by_col_name]]), ncol = max_lin_len+1)
  mat[, 1] <- df[[group_by_col_name]]
  for(i in 1:nrow(mat)){
    for(lin_level in 1:max_lin_len){
      
      len_lin <- length(split_lins[[i]])
      
      if(lin_level > len_lin){
        mat[i, lin_level+1] <- NA
      }else{
        mat[i, lin_level+1] <- paste0(split_lins[[i]][1:lin_level], collapse = ".")
      }
    }
  }
  
  max_lin <- vector()
  for(i in seq(nrow(mat))){
    max_lin[i] <- mat[i, which.max(sapply(mat[i, -1], len_str))+1]
  }
  mat <- data.frame(cbind(mat, max_lin), stringsAsFactors = F)
  names(mat) <- c("ID", paste0("lin_level_", 1:(ncol(mat)-2) ), "max_lin")
  return(mat)
}

# Directories

results_path <- "results/"
newick_path <- "results/newick/"
db_path <- "../pipeline/db/"
# Files ----
# in
all_lineages_file <- paste0(newick_path, "all_lineages.treefile")
metadata_file <- paste0(results_path, "lineage_file.csv")
lin_lookup_file <- paste0(db_path, "lineage_conversions.txt")
# out
all_lineages_outfile <- paste0(results_path, "all_lineages.png")

# Read in data ---- 
all_lineages_tree <- read.tree(all_lineages_file)
metadata <- read.csv(metadata_file)
lin_lookup <- read.table(lin_lookup_file, header = T, sep = "\t")

# Midpoint root ----

all_lineages_tree <- phytools::midpoint.root(all_lineages_tree)
n_samps <- length(all_lineages_tree$tip.label)

# Clean data ----

metadata <- subset(metadata, id %in% all_lineages_tree$tip.label)

# metadata <- subset(metadata, main_lin == "lineage1")
# all_lineages_tree <- keep.tip(all_lineages_tree, metadata$id)

# Clean
# lineage_data <- subset(lineage_data, !(main_lin == ""))
metadata$main_lin <- ifelse(metadata$main_lin == "", "lineage1", metadata$main_lin)
metadata$sublin <- ifelse(metadata$sublin == "", "lineage1", metadata$sublin)

# Lineages - remove 'lineage' and convert to factor
metadata$main_lin <- gsub('lineage', '', metadata$main_lin)
metadata$sublin <- gsub('lineage', '', metadata$sublin)

# Make ID rownames of metadata so heatmap strips work in ggtree
row.names(metadata) <- metadata$id

# Add main lineage to lookup table
lin_lookup$main_lineage <- unlist(lapply(strsplit(lin_lookup$mtbc_lineage, "\\."), function(x){x[1]}))

# Clean up L2
lin_lookup[lin_lookup[, "mtbc_lineage"] == "2.2.2", "rd_number"] <- "None"

# Get the lowest level of the RDs
lin_lookup$rd_lowest_level <- get_last(lin_lookup$rd_number)


# TEST - subset just L2 because layered 
lookup_L2 <- subset(lin_lookup, main_lineage == "2")
meta_L2 <- subset(metadata, main_lin == "2")

# Get the max dept of the sublins
max_depth_lins <- max(unlist(lapply(strsplit(meta_L2$sublin, "\\."), length)))

# Expand out the L2 lineages ready for merging
expanded_lins <- expand_hierarchy(meta_L2, "id", "sublin")


# # Parse out the lin levels into new cols
# meta_L2$lin_level_1 <- substr(meta_L2$sublin, 1, 1)
# meta_L2$lin_level_2 <- substr(meta_L2$sublin, 1, 3)
# meta_L2$lin_level_3 <- substr(meta_L2$sublin, 1, 5)
# meta_L2$lin_level_4 <- substr(meta_L2$sublin, 1, 7)


# Subset the RD-lin lookup to just the lin and lowest level RD
lookup_L2 <- select(lookup_L2, mtbc_lineage, rd_lowest_level)

for(i in seq(max_depth_lins)){
  col <- sprintf("lin_level_%s", i)
  print(col)
  expanded_lins <- merge(expanded_lins, lookup_L2,
                   by.x = col, by.y = "mtbc_lineage",
                   all.x = T,
                   sort = F)
  
  names(expanded_lins)[ncol(expanded_lins)] <- sprintf("rd_level_%s", i)
}





# tableA %>% 
#   inner_join(tableB, by = c("A_id" = "id"), suffix=c("_A","_B"))
# 
# 
# meta_L2 %>% 
#   inner_join(lookup_L2, by = c("lin_level_1" = "mtbc_lineage"), suffix = c("_A", "_B"))


# Get data for separate heatmap strips
lin_data <- dplyr::select(metadata, main_lin)
spol_data <- dplyr::select(metadata, spoligotype)
rd_data <- dplyr::select(metadata, rd)


# MRCA ----

# Get MRCA for each lineage
lin_split <- split(metadata, metadata$main_lin)
lin_mrca <- lapply(lin_split, function(x){
  getMRCA(all_lineages_tree, x$id)
})
lin_mrca_df <- data.frame(node = as.vector(unlist(lin_mrca)), name = names(lin_mrca))

# Get MRCA for each RD
rd_split <- split(metadata, metadata$rd)
rd_mrca <- lapply(rd_split, function(x){
  getMRCA(all_lineages_tree, x$id)
})
rd_mrca_df <- data.frame(node = as.vector(unlist(rd_mrca)), name = names(rd_mrca))
rd_mrca_df <- subset(rd_mrca_df, !(name == "None"))










# Clean up RDs 105 - 150
rd_mrca_df$name <- gsub("RD105;RD207;RD181", "RD181", rd_mrca_df$name)
rd_mrca_df$name <- gsub("RD181;RD142", "RD142", rd_mrca_df$name)
rd_mrca_df$name <- gsub("RD181;RD150", "RD150", rd_mrca_df$name)
rd_mrca_df$name <- gsub(";", "-", rd_mrca_df$name)

# Change col headers to match legends 
colnames(lin_mrca_df) <- c("node", "lineage")

# Colours
alpha <- 0.9
lin_colours <- rainbow(length(unique(metadata$main_lin)), alpha = 1)
spol_colours <- rainbow(length(unique(metadata$spoligotype)), alpha = alpha-0.2)
# rd_colours <- rainbow(length(unique(metadata$rd)), alpha = alpha-0.4)
rd_colours <- rep(c("black", "darkgrey"), length.out = length(rd_df$name))

names(lin_colours) <- c(sort(unique(metadata$main_lin)))
names(spol_colours) <- c(sort(unique(metadata$spoligotype)))
names(rd_colours) <- rd_df$name

# Set up ggtree parameters ----

width <- 0.05
font_sz <- 3
line_sz <- 0.25
angle <- 30

legend_spec <- theme(legend.title = element_text(size = 9),
                     legend.text = element_text(size = 7),
                     legend.key.size = unit(0.3, "cm"))

max_dist <- castor::get_tree_span(all_lineages_tree, as_edge_count=FALSE)$max_distance

# Make tree ----

ggtree_all_lineages <- ggtree(all_lineages_tree, size = line_sz, layout = "rectangular")

lin_tree <- ggtree_all_lineages + 
geom_hilight(data = lin_mrca_df, 
             mapping = aes(node = node, fill = lineage))+
  scale_fill_manual(values = lin_colours)

# # Add lineage data 
# lin_hm <- gheatmap(ggtree_all_lineages, lin_data, 
#                    width = width, 
#                    offset = 0, 
#                    colnames_position = "top",
#                    colnames_angle = angle, 
#                    colnames_offset_y = 1,
#                    hjust = 0,
#                    font.size = font_sz) +
#   # Add the custom colours defined above
#   scale_fill_manual(values = lin_colours, breaks = names(lin_colours) )
#   # Define the legend title
#   # labs(fill = lin_data_lab)

# # Pull the width of the strip from the plot just created in order to set offset for next strips
# ggtree_data <- ggplot2::ggplot_build(lin_hm)
# os <- unique(ggtree_data$data[[3]]$xmax)-unique(ggtree_data$data[[3]]$xmin)

ggtree_data <- ggplot2::ggplot_build(lin_tree)
os <- max(unique(ggtree_data$data[[3]]$xmax)-unique(ggtree_data$data[[3]]$xmin))

# lin_hm + geom_cladelab(data = rd_mrca_df,
lin_tree + geom_cladelab(data = rd_mrca_df,
                         mapping = aes(node = node, label = name), 
                         fontsize = 3, 
                         # offset = os/11,
                         align = F)+
  vexpand(.1)
  # scale_color_manual(values = rd_colours)



# # Need to do this bit of code before adding the next heatmap
# # See - See "7.3.1 Visualize tree with multiple associated matrix" https://yulab-smu.top/treedata-book/chapter7.html
# lin_hm <- lin_hm + ggnewscale::new_scale_fill() 
# 
# 
# # Spoligotype
# spol_hm <- gheatmap(lin_hm, spol_data, 
#                     width = width, 
#                     offset = os, 
#                     colnames_position = "top",
#                     colnames_angle = angle, 
#                     colnames_offset_y = 1,
#                     hjust = 0,
#                     font.size = font_sz) +
#   # Add the custom colours defined above
#   scale_fill_manual(values = spol_colours, breaks = names(spol_colours) )
# 
# spol_hm <- spol_hm + ggnewscale::new_scale_fill() 
# 
# rd_hm <- gheatmap(spol_hm, rd_data, 
#                     width = width, 
#                     offset = os*2, 
#                     colnames_position = "top",
#                     colnames_angle = angle, 
#                     colnames_offset_y = 1,
#                     hjust = 0,
#                     font.size = font_sz) +
#   # Add the custom colours defined above
#   scale_fill_manual(values = rd_colours, breaks = names(rd_colours) )
# 






























