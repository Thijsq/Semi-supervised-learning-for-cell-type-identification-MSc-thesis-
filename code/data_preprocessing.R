####################
# Author: Thijs Quast
# Contact: thijsquast@gmail.com
####################

# Load Macosko dataset .rds file
library('SingleCellExperiment')
macosko = readRDS('macosko.rds')

# Obtain raw count matrix:
macosko_counts = counts(macosko)

# Transform raw count matrix to R dataframe
macosko_counts = as.data.frame(macosko_counts)

# Load csv file with cell labels
macosko_cell_types = read.csv2("Macosko_cell_types.csv", header = FALSE)

# Explore dataset and basic manipulations
summary(macosko_cell_types)
macosko_cell_types = apply(macosko_cell_types, 1, as.character)
colnames(macosko_counts) = macosko_cell_types

# Extract overlapping cell types with other datset
# Store each label as a number as well, these numbers can later be used for several Python implementations
# which require numbers instead of characters

correct_cell_types_macosko = colnames(macosko_counts) == "amacrine"  | colnames(macosko_counts) == "bipolar" | colnames(macosko_counts) == "cones" | colnames(macosko_counts) == "muller" | colnames(macosko_counts) == "rods" 
macosko_labels = macosko_cell_types[correct_cell_types_macosko]
macosko_labels = as.data.frame(macosko_labels)
macosko_labels$nr = NA
macosko_labels$nr[macosko_labels$macosko_labels == "amacrine"] = 1
macosko_labels$nr[macosko_labels$macosko_labels == "bipolar"] = 2
macosko_labels$nr[macosko_labels$macosko_labels == "cones"] = 3
macosko_labels$nr[macosko_labels$macosko_labels == "muller"] = 4
macosko_labels$nr[macosko_labels$macosko_labels == "rods"] = 5
# To save file, uncomment line below
#write.csv(macosko_labels, "macosko_labels.csv")

macosko_counts = macosko_counts[,correct_cell_types_macosko]
# To save file, uncomment line below
#write.csv(macosko_counts, "macosko_transposed.csv")


# Repeat all steps above for other dataset
# Load shekhar file
shekhar = readRDS('shekhar.rds')

# get counts:
shekhar_counts = counts(shekhar)
shekhar_counts = as.data.frame(shekhar_counts)
shekhar_cell_types = read.csv2("Shekhar_cell_types.csv", header = FALSE)
summary(shekhar_cell_types)
shekhar_cell_types = apply(shekhar_cell_types, 1, as.character)
colnames(shekhar_counts) = shekhar_cell_types
correct_cell_types_shekhar = colnames(shekhar_counts) == "amacrine"  | colnames(shekhar_counts) == "bipolar" | colnames(shekhar_counts) == "cones" | colnames(shekhar_counts) == "muller" | colnames(shekhar_counts) == "rods"

shekhar_labels = shekhar_cell_types[correct_cell_types_shekhar]
shekhar_labels = as.data.frame(shekhar_labels)
shekhar_labels$nr = NA
shekhar_labels$nr[shekhar_labels$shekhar_labels == "amacrine"] = 1
shekhar_labels$nr[shekhar_labels$shekhar_labels == "bipolar"] = 2
shekhar_labels$nr[shekhar_labels$shekhar_labels == "cones"] = 3
shekhar_labels$nr[shekhar_labels$shekhar_labels == "muller"] = 4
shekhar_labels$nr[shekhar_labels$shekhar_labels == "rods"] = 5
#write.csv(shekhar_labels, "shekhar_labels.csv")


shekhar_counts = shekhar_counts[, correct_cell_types_shekhar]
#write.csv(shekhar_counts, "shekhar_transposed.csv")



########### shekhar downsampling on bipolar cells ##########

shekhar = readRDS('shekhar.rds')

# get counts:
shekhar_counts = counts(shekhar)
shekhar_counts = as.data.frame(shekhar_counts)
shekhar_cell_types = read.csv2("Shekhar_cell_types.csv", header = FALSE)
summary(shekhar_cell_types)
shekhar_cell_types = apply(shekhar_cell_types, 1, as.character)
colnames(shekhar_counts) = shekhar_cell_types
correct_cell_types_shekhar = colnames(shekhar_counts) == "amacrine"  | colnames(shekhar_counts) == "cones" | colnames(shekhar_counts) == "muller" | colnames(shekhar_counts) == "rods"

# Obtain only bipolar cells
cells_bipolar = colnames(shekhar_counts) == "bipolar"  

# Extract 2945 bipolar cells
first_bipolar = which(cells_bipolar == TRUE)
first_bipolar = first_bipolar[1:2945]

correct_cell_types_shekhar = c(which(correct_cell_types_shekhar == TRUE), first_bipolar)

shekhar_labels = shekhar_cell_types[correct_cell_types_shekhar]
shekhar_labels = as.data.frame(shekhar_labels)
shekhar_labels$nr = NA
shekhar_labels$nr[shekhar_labels$shekhar_labels == "amacrine"] = 1
shekhar_labels$nr[shekhar_labels$shekhar_labels == "bipolar"] = 2
shekhar_labels$nr[shekhar_labels$shekhar_labels == "cones"] = 3
shekhar_labels$nr[shekhar_labels$shekhar_labels == "muller"] = 4
shekhar_labels$nr[shekhar_labels$shekhar_labels == "rods"] = 5
#write.csv(shekhar_labels, "shekhar_labels_downsampled.csv")



#### transposing shekhar dataframe to conventional statistics dataset setup ######
shekhar_normal <- t(shekhar_counts)
shekhar_normal <- as.data.frame(shekhar_normal)
shekhar_normal <- shekhar_normal[correct_cell_types_shekhar,]
shekhar_normal$cell_type <- shekhar_labels[,1]
rownames(shekhar_normal) <- c(1:nrow(shekhar_normal))
#write.csv(shekhar_normal, "shekhar_normal_downsampled.csv")

#### transposing macosko dataframe to conventional statistics dataset setup ####
macosko_normal <- t(macosko_counts)
macosko_normal <- as.data.frame(macosko_normal)
macosko_normal$cell_type <- macosko_labels[,1]
rownames(macosko_normal) <- c(1:nrow(macosko_normal))
#write.csv(macosko_normal, "macosko_normal.csv")


# extracting common genes for both datasets #
common_genes <- intersect(names(shekhar_normal[,-ncol(shekhar_normal)]), names(macosko_normal[,-ncol(macosko_normal)]))
shekhar_normal <- shekhar_normal[, c(common_genes, "cell_type")]

macosko_normal <- macosko_normal[, c(common_genes, "cell_type")]

#write.csv(shekhar_normal, "shekhar_normal_downsampled.csv", row.names = FALSE)
#write.csv(macosko_normal, "macosko_normal.csv", row.names = FALSE)

library(readr)
### Creating full dataset of shekhar_normal_downsampled and macosko_normal, cells x genes ###
shekhar_normal_downsampled <- read_csv("shekhar_normal_downsampled.csv")
macosko_normal <- read_csv("macosko_normal.csv")

full_normal <- as.data.frame(rbind(shekhar_normal_downsampled, macosko_normal))
full_normal <- full_normal[,-ncol(full_normal)]
#write.csv(full_normal, "full_normal_shekhar_downsampled.csv", row.names = FALSE)
full_normal <- read_csv("/data/thijs/thesis/full_normal_shekhar_downsampled.csv")
full_normal_transposed <- as.data.frame(t(full_normal))

#write.csv(full_normal_transposed, "full_normal_transposed_shekhar_downsampled.csv")

### Downsampling overrepresented classes in test data, in case necessary ###

library('SingleCellExperiment')
macosko = readRDS('macosko.rds')

# get counts:
macosko_counts = counts(macosko)
macosko_counts = as.data.frame(macosko_counts)
macosko_cell_types = read.csv2("Macosko_cell_types.csv", header = FALSE)
summary(macosko_cell_types)
macosko_cell_types = apply(macosko_cell_types, 1, as.character)
colnames(macosko_counts) = macosko_cell_types
correct_cell_types_macosko = colnames(macosko_counts) == "amacrine"  | colnames(macosko_counts) == "bipolar" | colnames(macosko_counts) == "muller"
cells_cones = colnames(macosko_counts) == "cones"  
cells_rods = colnames(macosko_counts) == "rods"  

first_cones = which(cells_cones == TRUE)
first_cones = first_cones[1:480]

first_rods = which(cells_rods == TRUE)
first_rods = first_rods[1:910]

correct_cell_types_macosko = c(which(correct_cell_types_macosko == TRUE), first_cones, first_rods)

macosko_labels = macosko_cell_types[correct_cell_types_macosko]
macosko_labels = as.data.frame(macosko_labels)
macosko_labels$nr = NA
macosko_labels$nr[macosko_labels$macosko_labels == "amacrine"] = 1
macosko_labels$nr[macosko_labels$macosko_labels == "bipolar"] = 2
macosko_labels$nr[macosko_labels$macosko_labels == "cones"] = 3
macosko_labels$nr[macosko_labels$macosko_labels == "muller"] = 4
macosko_labels$nr[macosko_labels$macosko_labels == "rods"] = 5
#write.csv(macosko_labels, "macosko_labels_downsampled.csv")

macosko_counts = macosko_counts[,correct_cell_types_macosko]
#write.csv(macosko_counts, "macosko_transposed_downsampled.csv")


# creating downsampled macosko dataset in conventional statistics setup # 
macosko_normal_downsampled <- t(shekhar_counts)
shekhar_normal <- as.data.frame(shekhar_normal)
shekhar_normal$cell_type <- shekhar_labels[,1]
rownames(shekhar_normal) <- c(1:nrow(shekhar_normal))
#write.csv(shekhar_normal, "shekhar_normal_downsampled.csv")