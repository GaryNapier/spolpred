

setwd("~/Documents/spolpred/")

library(ggplot2)
library(ggtree)
library(phytools) # Midpoint root
library(ggnewscale) # For adding new heatmaps (annotations)
library(plyr)
library(dplyr)

# Directories

results_path <- "results/"
newick_path <- "results/newick/"

# Files ----
# in
all_lineages_file <- paste0(newick_path, "all_lineages.treefile")
metadata_file <- paste0(results_path, "lineage_file.csv")
# out
all_lineages_outfile <- paste0(results_path, "all_lineages.png")

# Read in data ---- 
all_lineages_tree <- read.tree(all_lineages_file)
metadata <- read.csv(metadata_file)

# Midpoint root ----
all_lineages_tree <- phytools::midpoint.root(all_lineages_tree)
n_samps <- length(all_lineages_tree$tip.label)

# Clean data ----

metadata <- subset(metadata, id %in% all_lineages_tree$tip.label)


# TEST

# metadata <- subset(metadata, main_lin == "lineage1")
# all_lineages_tree <- keep.tip(all_lineages_tree, metadata$id)




# Clean
# lineage_data <- subset(lineage_data, !(main_lin == ""))
metadata$main_lin <- ifelse(metadata$main_lin == "", "lineage1", metadata$main_lin)
metadata$sublin <- ifelse(metadata$sublin == "", "lineage1", metadata$sublin)

# Lineages - remove 'lineage' and convert to factor
metadata$main_lin <- factor(gsub('lineage', '', metadata$main_lin))
metadata$sublin <- factor(gsub('lineage', '', metadata$sublin))

# Make ID rownames of metadata so heatmap strips work in ggtree
row.names(metadata) <- metadata$id

# Get data for separate heatmap strips
lin_data <- dplyr::select(metadata, main_lin)
spol_data <- dplyr::select(metadata, spoligotype)
rd_data <- dplyr::select(metadata, rd)

# Wrangle data 
rd_split <- split(metadata, metadata$rd)
rd_mrca <- lapply(rd_split, function(x){
  getMRCA(all_lineages, x$id)
})
rd_df <- data.frame(node = as.vector(unlist(rd_mrca)), name = names(rd_mrca))
rd_df <- subset(rd_df, !(name == "None"))

# Change col headers to match legends 
# colnames(lin_data) <- lin_data_lab

# Colours
alpha <- 0.9
lin_colours <- rainbow(length(unique(metadata$main_lin)), alpha = 1)
spol_colours <- rainbow(length(unique(metadata$spoligotype)), alpha = alpha-0.2)
rd_colours <- rainbow(length(unique(metadata$rd)), alpha = alpha-0.4)

names(lin_colours) <- c(sort(unique(metadata$main_lin)))
names(spol_colours) <- c(sort(unique(metadata$spoligotype)))
names(rd_colours) <- c(sort(unique(metadata$rd)))

# Set up ggtree parameters ----

width <- 0.05
font_sz <- 3
line_sz <- 0.25
angle <- 30

y_lim <- c(-10, n_samps + (n_samps * 0.1))
legend_spec <- theme(legend.title = element_text(size = 9),
                     legend.text = element_text(size = 7),
                     legend.key.size = unit(0.3, "cm"))

max_dist <- castor::get_tree_span(all_lineages_tree, as_edge_count=FALSE)$max_distance

# Make tree ----

ggtree_all_lineages <- ggtree(all_lineages_tree, size = line_sz, layout = "rectangular")

# Add lineage data 
lin_hm <- gheatmap(ggtree_all_lineages, lin_data, 
                   width = width, 
                   offset = 0, 
                   colnames_position = "top",
                   colnames_angle = angle, 
                   colnames_offset_y = 1,
                   hjust = 0,
                   font.size = font_sz) +
  # Add the custom colours defined above
  scale_fill_manual(values = lin_colours, breaks = names(lin_colours) ) +
  # Define the legend title
  labs(fill = lin_data_lab)

# # Pull the width of the strip from the plot just created in order to set offset for next strips
ggtree_data <- ggplot2::ggplot_build(lin_hm)
os <- unique(ggtree_data$data[[3]]$xmax)-unique(ggtree_data$data[[3]]$xmin)

rd_tree <- lin_hm + geom_cladelab(data = rd_df, mapping = aes(node = node, label = name, color = name), fontsize = 3, offset = os, align = T)
rd_tree



# + vexpand(.1)











# 
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






























