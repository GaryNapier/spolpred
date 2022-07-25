

setwd("~/Documents/spolpred/")

# Functions ----

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
  
  bin_split <- str_split(x[, col], "")
  bin_mat <- data.frame(do.call("rbind", bin_split))
  names(bin_mat) <- paste0(prefix, 1:ncol(bin_mat))
  cbind(x, bin_mat)
}

# Paths ----

results_path <- "results/"
newick_path <- "results/newick/"
data_path <- "data/"
itol_templates_path <- "../pipeline/itol_templates/"

# Files ----

spoligo_lineage_full_file <- paste0(data_path, "spoligo_lineage.full.txt")
l3_treefile <- paste0(newick_path, "grouped_samples_lineage3.treefile")
# l3_treefile <- paste0(newick_path, "lin_3_all_samps.treefile")
itol_binary_template_file <- paste0(itol_templates_path, "itol.binary.txt")
itol_outfile <- paste0(results_path, "itol_all_spol_binary.txt")

# Read in data ----
spoligo_lineage_full <- read.table(spoligo_lineage_full_file, header = T, sep = "\t", 
                                   colClasses = c("spoligotype" = "character"))
l3_tree <- midpoint.root(read.tree(l3_treefile))
itol_binary_template <- readChar(itol_binary_template_file, nchars = 1e6)


# Clean ----

# l3_samples <- l3_tree$tip.label

# Remove "lineage"
spoligo_lineage_full$lineage <- gsub("lineage", "", spoligo_lineage_full$lineage)

# Remove and store animal strains
animal <- subset(spoligo_lineage_full, grepl("M", lineage))
spoligo_lineage_full <- subset(spoligo_lineage_full, !(grepl("M", lineage)))

# Expand lineages and merge
exp <- expand_hierarchy(spoligo_lineage_full, "sample", "lineage")
data <- merge(rbind(spoligo_lineage_full, animal), exp,
              by.x = "sample", by.y = "id", all.x = T, 
              sort = F)

data$lin_level_1 <- ifelse(grepl("M", data$lineage), data$lineage, data$lin_level_1)
data$max_lin <- ifelse(grepl("M", data$lineage), data$lineage, data$max_lin)

# Remove 0000000000000000000000000000000000000000000 ? 
data <- subset(data, !(spoligotype == "0000000000000000000000000000000000000000000"))

# Remove everything except lin and spoligotype info
data <- dplyr::select(data, sample:lineage, lin_level_1:max_lin)

# Keep only those with 5 or more samples per spoligo
data <- data %>% group_by(spoligotype) %>% filter(n() > 4) %>% data.frame()
# Add rownames for plotting
rownames(data) <- data$sample

# Subset L3 
# l3_spol_data <- subset(data, sample %in% l3_samples)
# Add binary cols
spol_data <- get_binary_cols(data, "spoligotype", "s")
# subset spol data 
spol_data_ggplot <- select(spol_data, s1:s43)


# itol ----

itol_data <- select(spol_data, sample, s1:s43)
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



# Set up ggtree parameters ----

width <- 0.05
font_sz <- 3
line_sz <- 0.25
angle <- 30

tree <- ggtree(l3_tree, size = line_sz, layout = "circular")

gheatmap(tree, spol_data_ggplot,
         # Increase offset
         # offset = offset+200,
         # width = width+0.5,
         # Change color to black
         # color = NULL,
         color="black",
         low="white",
         high="black",
         colnames_position = "top",
         # colnames_angle = angle,
         colnames_offset_y = 1,
         hjust = 0,
         font.size = 2.5) +
  # Define colours
  scale_fill_manual(values=c("white", "black"), labels = c("0", "1", "NA"), na.value = "grey")+
  labs(fill = "Spoligotype")







