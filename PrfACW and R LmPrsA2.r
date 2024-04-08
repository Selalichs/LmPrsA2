
#Remove anything from work space

rm(list = ls())


install.packages("DEP")
install.packages("dplyr")
#if the above did not work, use the biomanager-bioconducter install

#load packages
library(DEP)
library(dplyr)

data <- read.csv("C:/Users/selal/Desktop/PrsA2/count.csv", header=TRUE, stringsAsFactors=FALSE)
UbiLength_ExpDesign <- read.csv("C:/Users/selal/Desktop/PrsA2/col.csv", header=TRUE, stringsAsFactors=FALSE)


data <- count

UbiLength_ExpDesign <- col


# We filter for contaminant proteins and decoy database hits, which are indicated by "+" in the columns "Potential.contaminants" and "Reverse", respectively. Not for my Data
data <- filter(data, Reverse != "+", Potential.contaminant != "+") 

dim(data)

#The data.frame has the following column names:
colnames(data)

# Are there any duplicated gene names?
data$Gene.names %>% duplicated() %>% any()

# Make a table of duplicated gene names (Dont run)
data %>% group_by(Gene.names) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)

# Make unique names using the annotation in the "Gene.names" column as primary names and the annotation in "Protein.IDs" as name for those that do not have an gene name.
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")

# Are there any duplicated names?
data$name %>% duplicated() %>% any()

LFQ_columns <- grep("LFQ.", colnames(data_unique)) # get LFQ column numbers
experimental_design <- UbiLength_ExpDesign
data_se <- make_se(data_unique, LFQ_columns, experimental_design)

# Let's have a look at the SummarizedExperiment object
data_se

# Plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)

# Filter for proteins that are identified in all replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 0)

# Less stringent filtering:
# Filter for proteins that are identified in 2 out of 3 replicates of at least one condition
data_filt2 <- filter_missval(data_se, thr = 1)

# Plot a barplot of the number of identified proteins per samples
plot_numbers(data_filt)

# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_filt)

# Normalize the data
data_norm <- normalize_vsn(data_filt)

# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm)

# Plot a heatmap of proteins with missing values
plot_missval(data_filt)

# Plot intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(data_filt)

# All possible imputation methods are printed in an error, if an invalid function name is given.
impute(data_norm, fun = "")

# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)

# Impute missing data using random draws from a manually defined left-shifted Gaussian distribution (for MNAR)
data_imp_man <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)

# Impute missing data using the k-nearest neighbour approach (for MAR)
data_imp_knn <- impute(data_norm, fun = "knn", rowmax = 0.9)

# Plot intensity distributions before and after imputation
plot_imputation(data_norm, data_imp)

# Differential enrichment analysis  based on linear models and empherical Bayes statistics

# Test every sample versus control/Change this for release accordingly
data_diff <- test_diff(data_imp, type = "control", control = "prfacw")

# Test all possible comparisons of samples
data_diff_all_contrasts <- test_diff(data_imp, type = "all")

# Denote significant proteins based on user defined cutoffs/Note: cutoff was determined for CW and Release for PrfA background
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(0.5))

# Plot the first and second principal components
plot_pca(dep, x = 1, y = 2, n = 225, point_size = 4)
SaveFigure(plot_pca, "pca", height = 8, width = 8, res = 300)

# Plot the Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = -1, upper = 1, pal = "Reds")

SaveFigure(plot_cor, "correlation_2", height = 8, width = 8, res = 300)

# Plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = TRUE,
             indicate = c("condition", "replicate"))


# Plot a heatmap of all significant proteins (rows) and the tested contrasts (columns)
plot_heatmap(dep, type = "contrast", kmeans = TRUE, 
             k = 10, col_limit = 10, show_row_names = TRUE)
SaveFigure(heatmap_comp, "heatmap_comparisons", height = 12, width = 6, res = 300)



# Plot a volcano plot for the contrast 
plot_volcano(dep, contrast = "prfacwp2_vs_prfacw", label_size = 2, add_names = TRUE)
SaveFigure(volcano, "volcano", height = 8, width = 8, res = 300)

# Plot a barplot for prsA2 and hly
plot_single(dep, proteins = c("prsA2", "hly"))
SaveFigure(plot_single, "prsA2_hly_comparisons", height = 8, width = 12, res = 300)


plot_single(dep, proteins = c("hly", "actA"))
SaveFigure(plot_single, "prsA2_actA_comparisons", height = 8, width = 12, res = 300)


# Plot a barplot for the protein hly with the data centered
plot_single(dep, proteins = "hly", type = "centered")
SaveFigure(plot_single, "hly alone", height = 8, width = 12, res = 300)

plot_single_centered <- plot_single(dep, proteins = c("hly", "prsA2"), type = "centered") #produces centered log2-intensity
SaveFigure(plot_single_centered, "hly_prsA2_log2", height = 8, width = 12, res = 300)


# Plot a frequency plot of significant proteins for the different conditions
plot_cond(dep)

# Generate a results table
data_results <- get_results(dep)

# Number of significant proteins
data_results %>% filter(significant) %>% nrow()

# Column names of the results table
colnames(data_results)

# Save analyzed data
save(data_se, data_norm, data_imp, data_diff, dep, file = "Pcwdata.RData")
# These data can be loaded in future R sessions using this command
load("Pcwdata.RData")