
library(dplyr)
library(stringr)
library(ggplot2)
library(grid)
library(gridExtra)

# Functions ----

source("https://raw.githubusercontent.com/GaryNapier/Packages_functions/master/Functions.R")

# Paths ----

setwd("~/Documents/spolpred/r_scripts")
data_path <- "../data/"
results_path <- "../results/"

# Files ----

# in
# spol_data_file <- paste0(data_path, "spoligo_lineage.SNPs.csv")
spol_data_file <- paste0(data_path, "spoligo68.csv")
# animal_file <- paste0(data_path, "animal_spol_data.csv")
animal_file <- paste0(data_path, "animal_spol_data_68.csv")

# out
spol_lin_levels_table_file <- paste0(results_path, "spol_lin_levels_table_68.csv")

# Read in data ----

spol_data <- read.csv(spol_data_file, header = T, 
                      colClasses = c("spoligotype" = "character"))

animal <- read.csv(animal_file, header = T, 
                   colClasses = c("spoligotype" = "character"))

# Clean data ----

# Blank lineages - all 1.2.2.1 apparently
spol_data$lineage <- ifelse(spol_data$lineage == "", "lineage1.2.2.1", spol_data$lineage)
# Remove 0000000000000000000000000000000000000000000 - error apparently
# spol_data <- subset(spol_data, !(spoligotype == "0000000000000000000000000000000000000000000"))
spol_data <- subset(spol_data, 
                    !(spoligotype == "00000000000000000000000000000000000000000000000000000000000000000000"))



# Merge in SIT and family from 43-spacer?




# Clean SIT and family cols
# spol_data$SIT <- ifelse(spol_data$SIT == "", "Unknown", spol_data$SIT)
# spol_data$family <- ifelse(spol_data$family == "", "Unknown", spol_data$family)
# spol_data$family <- ifelse(spol_data$family == "-", "Unknown", spol_data$family)
# spol_data$family <- ifelse(spol_data$family == "Zero-copy", "Unknown", spol_data$family)

# Separate out all the metadata into lookup tables and data tables
# lin_spol_sit_fam <- odr(unique(select(spol_data, lineage, spoligotype, SIT, family)))
# sit_fam_spol <- odr(unique(select(spol_data, SIT, family, spoligotype)))
# sit_fam <- odr(unique(select(spol_data, SIT, family)))

# Separate out the data for just correlating lineage and spol

# Remove "lineage" so that the levels can be split out
spol_data$lineage <- gsub("lineage", "", spol_data$lineage)

# Clean animal spoligotypes (messed up when saving to csv and editing in excel)
animal <- merge(select(animal, -(spoligotype)), 
                select(spol_data, sample, spoligotype), 
                by = "sample", all.x = T, sort = F)
animal <- select(animal, sample, lineage, spoligotype, everything())

# # Just retain sample, lin and spol cols 
# spol_data <- select(spol_data, sample, lineage, spoligotype)
# 
# # Remove and store animal strains - need to process separately
# # Or only analyse non-animal?
# # animal <- subset(spol_data, grepl("La", lineage))
# spol_data <- subset(spol_data, !(grepl("La", lineage)))
# 
# # Expand lineages and merge
# exp <- expand_hierarchy_fill(spol_data, "sample", "lineage")
# # data <- merge(rbind(spol_data, animal), exp,
# #               by.x = "sample", by.y = "id", all.x = T, 
# #               sort = F)
# spol_data <- merge(spol_data, exp,
#                    by.x = "sample", by.y = "id", all.x = T, 
#                    sort = F)

# Just retain sample, lin and spol cols 
spol_data <- select(spol_data, sample, lineage, spoligotype)

# Remove and store animal strains - need to process separately
# Add back in after expanding animal strains manually in excel
# animal <- subset(spol_data, grepl("La", lineage))
# write.csv(animal, file = paste0(data_path, "animal_spol_data.csv"), 
#           quote = F, row.names = F)
spol_data <- subset(spol_data, !(grepl("La", lineage)))

# Expand lineages, rbind animal back in and merge
exp <- expand_hierarchy_fill(spol_data, "sample", "lineage")
# data <- merge(rbind(spol_data, animal), exp,
#               by.x = "sample", by.y = "id", all.x = T, 
#               sort = F)
spol_data <- merge(spol_data, exp,
                   by.x = "sample", by.y = "id", all.x = T, 
                   sort = F)

# Add animal back in after processing manually
spol_data <- rbind(spol_data, animal)

# Keep only those with 5 or more samples per spoligo
spol_data <- spol_data %>% group_by(spoligotype) %>% filter(n() > 4) %>% data.frame()

levels <- c("lin_level_1", "lin_level_2", "lin_level_3", "lin_level_4")

# Loop over levels and calculate top spoligotypes for each lin
spol_lin_levels_table <- list()
for(level in levels){
  # Make freq table for the level
  lv_tab <- tab2df(dplyr::select(spol_data, spoligotype, all_of(level)))
  
  # Calc probs
  lv_col_prob <- col_probs(lv_tab)
  lv_prob <- row_probs(col_probs(lv_tab))
  
  # Save the spoligotype
  spol <- row.names(lv_prob)
  # Save the colnames before adding spoligotype
  col_nms <- colnames(lv_prob)
  lv_tab$spoligotype <- row.names(lv_tab)
  lv_prob$spoligotype <- row.names(lv_prob)
  lv_col_prob$spoligotype <- row.names(lv_col_prob)
  
  
  # Pull out the non-zeros for each lineage and store each lineage (col) as a separate list entry
  # Save full table and a top - X table of results
  lv_list <- list()
  lv_list_top <- list()
  for(col in col_nms){
    df <- lv_prob[, c("spoligotype", col)]
    df <- df[which(df[, col] > 0), ]
    lv_list[[col]] <- sort_df_by_col_name(df, col)
  }
  
  # Merge in sample numbers (freqs) and col probs (freqs as proportion)
  
  for(lin in names(lv_list)){
    
    # Merge in the freqs
    lv_list[[lin]] <- setNames(merge(lv_list[[lin]], lv_tab[, c("spoligotype", lin)], 
                                     by = "spoligotype", 
                                     sort = F), 
                               c("spoligotype", lin, "freq"))
    
    # Merge in the freq as a proportion of the total n samps per lin
    lv_list[[lin]] <- setNames(merge(lv_list[[lin]], lv_col_prob[, c("spoligotype", lin)],
                                     by = "spoligotype",
                                     sort = F), 
                               c("spoligotype", "prob", "freq", "col_prob"))
    
    # names(lv1_list[[lin]]) <- c("spol", "prob", "freq")
    
    # Sort
    lv_list[[lin]] <- dplyr::arrange(lv_list[[lin]], desc(prob), desc(freq))
    
  }
  
  # Add the lineages as a col before r-binding
  
  # Full list
  lv_list <- lapply(seq(lv_list), function(i){
    lin <- names(lv_list[i])
    lv_list[[i]]$lineage <- rep(lin, nrow(lv_list[[i]])); lv_list[[i]]
  })
  lv_df <- do.call("rbind", lv_list)
  
  # Merge the SIT and family cols 
  lv_df <- plyr::join(lv_df, sit_fam_spol, by = "spoligotype", type = "left")
  
  # Tidy
  lv_df$col_prob <- round(lv_df$col_prob*100, 1)
  lv_df <- dplyr::select(lv_df, lineage, spoligotype, SIT, family, prob, freq, col_prob)
  
  # names(lv_df) <- c("lineage", "spoligotype", "SIT", "family", "weighted proportion \nof lineage", "n", "%n of lineage")
  
  spol_lin_levels_table[[level]] <- lv_df
  
}

# Add the level as a column before r-binding
spol_lin_levels_table <- lapply(seq(spol_lin_levels_table), function(i){
  level <- names(spol_lin_levels_table[i])
  spol_lin_levels_table[[i]]$level <- rep(level, nrow(spol_lin_levels_table[[i]])); spol_lin_levels_table[[i]]
})

spol_lin_levels_table <- do.call("rbind", spol_lin_levels_table)
spol_lin_levels_table <- spol_lin_levels_table %>% select(level, everything())
spol_lin_levels_table$level <- gsub("lin_level_", "", spol_lin_levels_table$level)
spol_lin_levels_table$prob <- round(spol_lin_levels_table$prob, 2)

write.csv(spol_lin_levels_table, file = spol_lin_levels_table_file, quote = F, row.names = F)

# Make copy with different headers

spol_lin_levels_table_github <- spol_lin_levels_table

names(spol_lin_levels_table_github) <- c("Level", "Lineage", "Spoligotype", "SIT", "Family", 
                                         "Correlation", "N", "% in level-lineage")

write.csv(spol_lin_levels_table_github, file = "../results/spol_lin_levels_table_68_github.csv", 
          quote = F, row.names = F)



# Plot ----

# Get main lineage
spol_lin_levels_table$main_lin <- unlist(lapply(strsplit(spol_lin_levels_table$lineage, "\\."), 
                                                function(x){x[1]}))
# Clean up animal
spol_lin_levels_table$main_lin <- ifelse(grepl("La", spol_lin_levels_table$main_lin), 
                                         "La", 
                                         spol_lin_levels_table$main_lin)

pc_plot <- ggplot(data = spol_lin_levels_table, 
                  aes(x = main_lin, 
                      y = col_prob, 
                      fill = main_lin))+
  geom_boxplot(outlier.size = 1,
               outlier.stroke = 0.2, 
               size = 0.2)+
  stat_summary(fun = max, geom = "point", shape = 18,
               size = 5, color = "red")+
  xlab("Main lineage")+
  ylab("%")+
  ggtitle("% of samples with a spol. at a lineage-level")+
  facet_wrap(~level)+
  theme_bw()+
  theme(legend.position="none")

cor_plot <- ggplot(data = spol_lin_levels_table, 
                   aes(x = main_lin, 
                       y = prob, 
                       fill = main_lin))+
  geom_boxplot(outlier.size = 1,
               outlier.stroke = 0.2, 
               size = 0.2)+
  stat_summary(fun = max, geom = "point", shape = 18,
               size = 5, color = "red")+
  xlab("Main lineage")+
  ylab("Correlation score")+
  ggtitle("Lineage-level-spoligotype correlation")+
  facet_wrap(~level)+
  theme_bw()+
  theme(legend.position="none")

boxplots <- grid.arrange(pc_plot, cor_plot, 
                         nrow = 2, 
                         top = textGrob("68-spacer spoligotypes",
                                        gp = gpar(fontsize=20)))

ggsave(filename = paste0(results_path, "spol_68_boxplots.png"), boxplots)
ggsave(filename = "spolpred_shiny/www/spol_68_boxplots.png", boxplots)

# rbind_list_with_names <- function(my_list, name_new_col){
#   # rbind a list but add the names of each list element as a new col, supplying a name for that new col
#   
#   # in:
#   # my_list:
#   # $group_1
#     # X  Y
#     # 1 11
#     # 2 12
#     # 3 13
#     # 4 14
#     # 5 15
#   # $group_2
#     # X   Y
#     # 51 101
#     # 52 102
#     # 53 103
#     # 54 104
#     # 55 105
#   
#   # out:
#   # rbind_list_with_names(my_list, "group")
#   # X   Y   group
#   # 1  11   group_1
#   # 2  12   group_1
#   # 3  13   group_1
#   # 4  14   group_1
#   # 5  15   group_1
#   # 51 101  group_2
#   # 52 102  group_2
#   # 53 103  group_2
#   # 54 104  group_2
#   # 55 105  group_2
#   
#   lapply(seq(my_list), function(i){
#     nms <- names(my_list[i])
#     my_list[[i]][name_new_col] <- rep(nms, nrow(my_list[[i]])); my_list[[i]]
#   })
#   my_df <- do.call("rbind", my_list)
#   my_df
# }












