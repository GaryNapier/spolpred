# R
library(sjmisc)
library(dplyr)
library(tibble)

# Functions ----

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
  names(mat) <- c("id", paste0("lin_level_", 1:(ncol(mat)-2) ), "max_lin")
  return(mat)
}

# setwd("/Users/cmeehan2/Desktop/SpoligoLineageComp/Ranalyses")
setwd("~/Documents/spolpred/")

# Paths
data_path <- "data/"

# Files
spol_file <- paste0(data_path, "spoligotypes.csv")

#load the data in

# data <- read.table("spoligo_lineage.collated_noTotals.txt", header = TRUE)
data <- read.csv(spol_file, colClasses = c("octal_spoligotype" = "character", "binary_spoligotype" = "character"))

# Exploratory ----

head(data)
nrow(data)
length(unique(data$id))
unique(data$main_lineage)
unique(data$sublineage)
unique(data$octal_spoligotype)
unique(data$binary_spoligotype)
"0000000000000000000000000000000000000000000" %in% data$binary_spoligotype
subset(data, binary_spoligotype == "0000000000000000000000000000000000000000000")
# How many mixed/blank lins?
sum(grepl(";", data[, "main_lineage"]))
sum(grepl(";", data[, "sublineage"]))
sum(data[, "main_lineage"] == "")
sum(data[, "sublineage"] == "")

# Clean ----

# Removed mixed and blank lins
data <- subset(data, !(main_lineage == "" | grepl(";", main_lineage) | sublineage == "" | grepl(";", sublineage)))

# Remove 'lineage'
data$main_lineage <- gsub("lineage", "", data$main_lineage)
data$sublineage <- gsub("lineage", "", data$sublineage)

# Separate animal (deal with later)
animal_data <- subset(data, grepl("La", main_lineage))
data <- subset(data, !(grepl("La", main_lineage)))

# Expand sublins
expanded <- expand_hierarchy(data, "id", "sublineage")
data <- merge(data, expanded, by = "id", sort = F)


# ----- CONOR NOTES -------

# conda activate r_env
# cd /Users/cmeehan2/Desktop/SpoligoLineageComp/Ranalyses 

# Use the spoligo_lineage.full.txt file and remove all after the . in the lineage. Also put a 's' before the spoligotype to ensure it is not included
# Iteratively do this in other columns so that its sub-lineage (3.1) and then sub-sublineage (3.1.1) etc. if no sub-lineage for that spoligotype pattern,
# retain the higher level lineage definition.
# Results in 4 levels
# Use excel and BBedit to do this

# e.g. on each column
# Find:
# (lineage\d+\.\d+)\..*
# Replace
# \1

# need to fix the M.orygis etc at lineage level 1
# Find:
# lineageM\t(lineageM.[a-z]+)\t
# replace
# \1\t\1\t

#saved as spoligo_lineage_splits.txt

# remove the s0000000000000000000000000000000000000000000 as this is an error
# saved as spoligo_lineage_splits_noErrorSpol.txt 

# Removed the totals from the collated file to use here
# put s before the spoligo names to stop it converting to scientific notation

# ----- CONOR NOTES -------

s <- data[,-1]
rownames(s) <- data[,1]

#Remove s0000000000000000000000000000000000000000000 as this is an error (cannot find spoligo)
# row.names.remove <- c("s0000000000000000000000000000000000000000000")
row.names.remove <- c("0000000000000000000000000000000000000000000")
spol <- s[!(row.names(s) %in% row.names.remove), ]



# STUCK HERE

#percentage how many spoligotypes appear in the data once, 1+, 5+, 10+, 100+
totalsSpols <- as.data.frame(rowSums(spol))

# STUCK HERE





nrow(totalsSpols)										#3146
sum(totalsSpols[,1] == 1)/ nrow(totalsSpols) *100		#58.39
sum(totalsSpols[,1] > 1)/ nrow(totalsSpols) *100		#41.61	
sum(totalsSpols[,1] > 5)/ nrow(totalsSpols) *100		#12.02
sum(totalsSpols[,1] > 10)/ nrow(totalsSpols) *100		#7.41
sum(totalsSpols[,1] > 100)/ nrow(totalsSpols) *100		#1.34

#How many spoligotypes contribute more than 1%, 5% of the data?
dbTotal=colSums(totalsSpols)							#33180
sum(totalsSpols[,1]/ dbTotal *100 > 1)					#10				
sum(totalsSpols[,1]/ dbTotal *100 > 5)					#2
sum(totalsSpols[,1]/ dbTotal *100 > 10)					#2			


#percentage how many lineages appear in the data once, 1+, 5+, 10+, 100+
totalsLineages=as.data.frame(colSums(spol))
nrow(totalsLineages)										#83
sum(totalsLineages[1] == 1)/ nrow(totalsLineages) *100		#0
sum(totalsLineages[,1] > 1)/ nrow(totalsLineages) *100		#100
sum(totalsLineages[,1] > 5)/ nrow(totalsLineages) *100		#95.18
sum(totalsLineages[,1] > 10)/ nrow(totalsLineages) *100		#93.97
sum(totalsLineages[,1] > 100)/ nrow(totalsLineages) *100	#61.45

#How many lineages contribute more than 1%, 5% of the data?
dbTotal=colSums(totalsLineages)								#33180
sum(totalsLineages[,1]/ dbTotal *100 > 1)					#23
sum(totalsLineages[,1]/ dbTotal *100 > 5)					#4
sum(totalsLineages[,1]/ dbTotal *100 > 10)					#1



#count how many spoligotypes are in more than 1 lineage etc
data<- read.table("spoligo_lineage_splits_noErrorSpol.txt",header=TRUE)
s <- data[,-1]
rownames(s) <- data[,1]
lv1Summary <- data %>% count(spoligotype, lv1Lineage)
SpoligoMultiLv1Lineage <- lv1Summary %>%   group_by(spoligotype) %>%   filter(n()>1)
lv2Summary <- data %>% count(spoligotype, lv2Lineage)
SpoligoMultiLv2Lineage <- lv2Summary %>%   group_by(spoligotype) %>%   filter(n()>1)
lv3Summary <- data %>% count(spoligotype, lv3Lineage)
SpoligoMultiLv3Lineage <- lv3Summary %>%   group_by(spoligotype) %>%   filter(n()>1)
lv4Summary <- data %>% count(spoligotype, lv4Lineage)
SpoligoMultiLv4Lineage <- lv4Summary %>%   group_by(spoligotype) %>%   filter(n()>1)

#in more than 2 lineages
SpoligoMultiLv1Lineage2 <- lv1Summary %>%   group_by(spoligotype) %>%   filter(n()>2)
SpoligoMultiLv2Lineage2 <- lv2Summary %>%   group_by(spoligotype) %>%   filter(n()>2)
SpoligoMultiLv3Lineage2 <- lv3Summary %>%   group_by(spoligotype) %>%   filter(n()>2)
SpoligoMultiLv4Lineage2 <- lv4Summary %>%   group_by(spoligotype) %>%   filter(n()>2)

#in more than 3 lineages
SpoligoMultiLv1Lineage3 <- lv1Summary %>%   group_by(spoligotype) %>%   filter(n()>3)
SpoligoMultiLv2Lineage3 <- lv2Summary %>%   group_by(spoligotype) %>%   filter(n()>3)
SpoligoMultiLv3Lineage3 <- lv3Summary %>%   group_by(spoligotype) %>%   filter(n()>3)
SpoligoMultiLv4Lineage3 <- lv4Summary %>%   group_by(spoligotype) %>%   filter(n()>3)

#in more than 4 lineages
SpoligoMultiLv1Lineage4 <- lv1Summary %>%   group_by(spoligotype) %>%   filter(n()>4)
SpoligoMultiLv2Lineage4 <- lv2Summary %>%   group_by(spoligotype) %>%   filter(n()>4)
SpoligoMultiLv3Lineage4 <- lv3Summary %>%   group_by(spoligotype) %>%   filter(n()>4)
SpoligoMultiLv4Lineage4 <- lv4Summary %>%   group_by(spoligotype) %>%   filter(n()>4)


#TO DO
#Keep only those with 5 or more samples per spoligo

SpoligoOver5Samples <- data %>%   group_by(spoligotype) %>%   filter(n()>4)
#29,211 samples
write.table(SpoligoOver5Samples, "/Users/cmeehan2/Desktop/SpoligoLineageComp/Ranalyses/spoligo_lineage_splits_noErrorSpol_over5samples.txt", sep="\t", quote = FALSE, row.names = FALSE)






#Look at correlation between spoligotype and lineage levels


#get the Theil's U correlation between the spoligo and the 4 levels of lineage
# python ~/Dropbox/Scripts/TheilsU_mx.py --matrix spoligo_lineage_splits_noErrorSpol_over5samples.txt


#load the data into Rm on seerver

library(dplyr)

data <- read.table("spoligo_lineage_splits_noErrorSpol_over5samples.txt", header = TRUE)

s <- data[,-1]
rownames(s) <- data[,1]


# Do ANOVA between lineages (USEFUL??)
fm <- aov(as.numeric(spoligotype) ~ lv1Lineage, data=data)
a <- anova(fm)
a
posthoc <- TukeyHSD(x = fm, 'lv1Lineage', conf.level = 0.95)
posthoc
fm <- aov(as.numeric(spoligotype) ~ lv2Lineage, data = data)
a <- anova(fm)
a
fm <- aov(as.numeric(spoligotype) ~ lv3Lineage, data = data)
a <- anova(fm)
a
fm <- aov(as.numeric(spoligotype) ~ lv4Lineage, data = data)
a <- anova(fm)
a


#Build a predictor of lineage from spoligo using Random Forest
#Matthews corrrlation coefficient or f-statistic or MCC?

#random forest with h20
library(h2o)
h2o.init()

dataRF <- h2o.importFile("/mnt/DATA2/conor/spoligoAllLineages/Ranalyses/spoligo_lineage_splits_noErrorSpol_over5samples.txt", header=TRUE)

# Set the predictors and response for Lv1Lineage
predictors <- c("spoligotype")
responselv1 <- "lv1Lineage"

# Split the dataset into a train and valid set: #not used in cross fold validation
data_split <- h2o.splitFrame(data = dataRF, ratios = 0.8)
train <- data_split[[1]]
valid <- data_split[[2]]


#Create predictor
data_drflv1 <- h2o.randomForest(x = predictors, y = responselv1, ntrees = 1000, max_depth = 0, nfolds = 5,
                             min_rows = 1,
                             categorical_encoding = "auto",
                             training_frame = dataRF,
                             balance_classes = TRUE)

#save the model
lv1model_path <- h2o.saveModel(object = data_drflv1, path = getwd(), force = TRUE)
print(lv1model_path)
#returns /mnt/DATA2/conor/spoligoAllLineages/Ranalyses/DRF_model_R_1618824303760_17

# Eval performance:
perflv1 <- h2o.performance(data_drflv1)
perflv1

#main output:
# MSE: (Extract with `h2o.mse`) 0.001598588
# RMSE: (Extract with `h2o.rmse`) 0.03998234
# Logloss: (Extract with `h2o.logloss`) 0.008292407
# Mean Per-Class Error: 0.00162723
# R^2: (Extract with `h2o.r2`) 0.9997604
# Confusion Matrix: Extract with `h2o.confusionMatrix(<model>,train = TRUE)`)
# =========================================================================
# Confusion Matrix: Row labels: Actual class; Column labels: Predicted class
#                 lineage1 lineage2 lineage3 lineage4 lineage5 lineage6 lineage7
# lineage1           14919        7        0       13        0        0        0
# lineage2               0    14914        0        2        0        0        0
# lineage3               0      105    14773        0        0        0        0
# lineage4               2        0        2    14876        0        0        0
# lineage5               0       87        0        0    14795        0        0
# lineage6               0        0        0        0        0    14889        0
# lineage7               0        0        0        0        0        0    14879
# lineageM.bovis         0        0        0        0        0        0        0
# lineageM.orygis        0        0        0        0        0        0        0
# Totals             14921    15113    14775    14891    14795    14889    14879
#                 lineageM.bovis lineageM.orygis  Error            Rate
# lineage1                     0               0 0.0013 =   20 / 14,939
# lineage2                     0               0 0.0001 =    2 / 14,916
# lineage3                     0               0 0.0071 =  105 / 14,878
# lineage4                     0               0 0.0003 =    4 / 14,880
# lineage5                     0               0 0.0058 =   87 / 14,882
# lineage6                     0               0 0.0000 =    0 / 14,889
# lineage7                     0               0 0.0000 =    0 / 14,879
# lineageM.bovis           14872               0 0.0000 =    0 / 14,872
# lineageM.orygis              0           14880 0.0000 =    0 / 14,880
# Totals                   14872           14880 0.0016 = 218 / 134,015

# Hit Ratio Table: Extract with `h2o.hit_ratio_table(<model>,train = TRUE)`
# =======================================================================
# Top-9 Hit Ratios: 
#   k hit_ratio
# 1 1  0.998373
# 2 2  0.999888
# 3 3  1.000000
# 4 4  1.000000
# 5 5  1.000000
# 6 6  1.000000
# 7 7  1.000000
# 8 8  1.000000
# 9 9  1.000000


# Generate predictions on a validation set (if necessary):
predictlv1 <- h2o.predict(data_drflv1, newdata = valid)


#Repeat for level 2

# Set the predictors and response for Lv2Lineage
predictors <- c("spoligotype")
responselv2 <- "lv2Lineage"

# Split the dataset into a train and valid set:
data_split <- h2o.splitFrame(data = dataRF, ratios = 0.8)
train <- data_split[[1]]
valid <- data_split[[2]]

#Create predictor
data_drflv2 <- h2o.randomForest(x = predictors, y = responselv2, ntrees = 1000, max_depth = 0, nfolds = 5,
                             min_rows = 1,
                             categorical_encoding = "auto",
                             training_frame = dataRF,
                             balance_classes = TRUE)
#save the model
lv2model_path <- h2o.saveModel(object = data_drflv2, path = getwd(), force = TRUE)
print(lv2model_path)
#returns /mnt/DATA2/conor/spoligoAllLineages/Ranalyses/DRF_model_R_1618824303760_18


# Eval performance:
perflv2 <- h2o.performance(data_drflv2)
perflv2

#main output:
# MSE: (Extract with `h2o.mse`) 0.1554203
# RMSE: (Extract with `h2o.rmse`) 0.3942338
# Logloss: (Extract with `h2o.logloss`) 0.5240711
# Mean Per-Class Error: 0.1711089
# R^2: (Extract with `h2o.r2`) 0.9964673
# Confusion Matrix: Extract with `h2o.confusionMatrix(<model>,train = TRUE)`)
# =========================================================================
# Confusion Matrix: Row labels: Actual class; Column labels: Predicted class
#            lineage1.1 lineage1.2 lineage1.3 lineage2.1 lineage2.2 lineage3
# lineage1.1       6094         42        182          0          5        0
# lineage1.2        118       6230          0          0          0        0
# lineage1.3        132          0       6216          0          0        0
# lineage2.1          0          0          0       6277          0        0
# lineage2.2          0          0          0          0       6315        0
#            lineage3.1 lineage3.2 lineage4 lineage4.1 lineage4.2 lineage4.3
# lineage1.1          0          0        0          5          0          0
# lineage1.2          0          0        0          0          0          0
# lineage1.3          0          0        0          0          0          0
# lineage2.1          0          0        0         68          0          0
# lineage2.2          0          0        0          0          0          0
#            lineage4.4 lineage4.5 lineage4.6 lineage4.7 lineage4.8 lineage4.9
# lineage1.1          0          0          0          5          0          0
# lineage1.2          0          0          0          0          0          0
# lineage1.3          0          0          0          0          0          0
# lineage2.1          0          0          0          0          0          0
# lineage2.2          0          0          0          0          0          0
#            lineage5 lineage6 lineage7 lineageM.bovis lineageM.orygis  Error
# lineage1.1        0        0        0              0               0 0.0377
# lineage1.2        0        0        0              0               0 0.0186
# lineage1.3        0        0        0              0               0 0.0208
# lineage2.1        0        0        0              0               0 0.0107
# lineage2.2        0        0        0              0               0 0.0000
#                          Rate
# lineage1.1 =      239 / 6,333
# lineage1.2 =      118 / 6,348
# lineage1.3 =      132 / 6,348
# lineage2.1 =       68 / 6,345
# lineage2.2 =        0 / 6,315
# 
# ---
#                 lineage1.1 lineage1.2 lineage1.3 lineage2.1 lineage2.2 lineage3
# lineage5                 0          0          0          0         37        0
# lineage6                 0          0          0          0          0        0
# lineage7                 0          0          0          0          0        0
# lineageM.bovis           0          0          0          0          0        0
# lineageM.orygis          0          0          0          0          0        0
# Totals                6374       6272       6398       6277       6471     8506
#                 lineage3.1 lineage3.2 lineage4 lineage4.1 lineage4.2 lineage4.3
# lineage5                 0          0        0          0          0          0
# lineage6                 0          0        0          0          0          0
# lineage7                 0          0        0          0          0          0
# lineageM.bovis           0          0        0          0          0          0
# lineageM.orygis          0          0        0          0          0          0
# Totals                5662       4767     1994       6040       4324       6905
#                 lineage4.4 lineage4.5 lineage4.6 lineage4.7 lineage4.8
# lineage5                 0          0          0          0          0
# lineage6                 0          0          0          0          0
# lineage7                 0          0          0          0          0
# lineageM.bovis           0          0          0          0          0
# lineageM.orygis          0          0          0          0          0
# Totals                6211       3346       6841       4205      22459
#                 lineage4.9 lineage5 lineage6 lineage7 lineageM.bovis
# lineage5                 0     6312        0        0              0
# lineage6                 0        0     6348        0              0
# lineage7                 0        0        0     6344              0
# lineageM.bovis           0        0        0        0           6349
# lineageM.orygis          0        0        0        0              0
# Totals                1172     6312     6348     6344           6349
#                 lineageM.orygis  Error               Rate
# lineage5                      0 0.0058 =       37 / 6,349
# lineage6                      0 0.0000 =        0 / 6,348
# lineage7                      0 0.0000 =        0 / 6,344
# lineageM.bovis                0 0.0000 =        0 / 6,349
# lineageM.orygis            6348 0.0000 =        0 / 6,348
# Totals                     6348 0.1712 = 24,987 / 145,925
# 
# Hit Ratio Table: Extract with `h2o.hit_ratio_table(<model>,train = TRUE)`
# =======================================================================
# Top-10 Hit Ratios: 
#     k hit_ratio
# 1   1  0.828768
# 2   2  0.891191
# 3   3  0.925400
# 4   4  0.940785
# 5   5  0.959232
# 6   6  0.985609
# 7   7  0.989762
# 8   8  0.998342
# 9   9  0.999959
# 10 10  0.999973



# Generate predictions on a validation set (if necessary):
predictlv2 <- h2o.predict(data_drflv2, newdata = valid)


#Repeat for level 3

# Set the predictors and response for Lv3Lineage
predictors <- c("spoligotype")
responselv3 <- "lv3Lineage"

# Split the dataset into a train and valid set:
data_split <- h2o.splitFrame(data = dataRF, ratios = 0.8)
train <- data_split[[1]]
valid <- data_split[[2]]

#Create predictor
data_drflv3 <- h2o.randomForest(x = predictors, y = responselv3, ntrees = 1000, max_depth = 0, nfolds = 5,
                             min_rows = 1,
                             categorical_encoding = "auto",
                             training_frame = dataRF,
                             balance_classes = TRUE)

#save the model
lv3model_path <- h2o.saveModel(object = data_drflv3, path = getwd(), force = TRUE)
print(lv3model_path)
#returns /mnt/DATA2/conor/spoligoAllLineages/Ranalyses/DRF_model_R_1618824303760_8

# Eval performance:
perflv3 <- h2o.performance(data_drflv3)
perflv3

#main output:
# MSE: (Extract with `h2o.mse`) 0.2201681
# RMSE: (Extract with `h2o.rmse`) 0.4692208
# Logloss: (Extract with `h2o.logloss`) 0.8634433
# Mean Per-Class Error: 0.2434693
# R^2: (Extract with `h2o.r2`) 0.999059
# 
# Top-10 Hit Ratios: 
#     k hit_ratio
# 1   1  0.756458
# 2   2  0.837301
# 3   3  0.873993
# 4   4  0.888551
# 5   5  0.917632
# 6   6  0.934011
# 7   7  0.943905
# 8   8  0.960928
# 9   9  0.980771
# 10 10  0.985279



# Generate predictions on a validation set (if necessary):
predictlv3 <- h2o.predict(data_drflv3, newdata = valid)


#Repeat for level 4

# Set the predictors and response for Lv3Lineage
predictors <- c("spoligotype")
responselv4 <- "lv4Lineage"

# Split the dataset into a train and valid set:
data_split <- h2o.splitFrame(data = dataRF, ratios = 0.8)
train <- data_split[[1]]
valid <- data_split[[2]]

#Create predictor
data_drflv4 <- h2o.randomForest(x = predictors, y = responselv4, ntrees = 1000, max_depth = 0, nfolds = 5,
                             min_rows = 1,
                             categorical_encoding = "auto",
                             training_frame = dataRF,
                             balance_classes = TRUE)

#save the model
lv4model_path <- h2o.saveModel(object = data_drflv4, path = getwd(), force = TRUE)
print(lv4model_path)
#returns /mnt/DATA2/conor/spoligoAllLineages/Ranalyses/DRF_model_R_1618824303760_11


# Eval performance:
perflv4 <- h2o.performance(data_drflv4)
perflv4

#main output:
# MSE: (Extract with `h2o.mse`) 0.2611564
# RMSE: (Extract with `h2o.rmse`) 0.5110346
# Logloss: (Extract with `h2o.logloss`) 1.055343
# Mean Per-Class Error: 0.290533
# R^2: (Extract with `h2o.r2`) 0.9994977
# Top-10 Hit Ratios: 
#     k hit_ratio
# 1   1  0.709479
# 2   2  0.805775
# 3   3  0.855677
# 4   4  0.891904
# 5   5  0.908387
# 6   6  0.921357
# 7   7  0.936508
# 8   8  0.948893
# 9   9  0.954528
# 10 10  0.966473

# Generate predictions on a validation set (if necessary):
predictlv4 <- h2o.predict(data_drflv4, newdata = valid)

