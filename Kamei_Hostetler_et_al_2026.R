#R Script associated with Kamei & Hostetler et al., 2026####
library(agricolae)
library(car)
library(cowplot)
library(DESeq2)
library(dplyr)
library(emmeans)
library(ggplot2)
library(gridExtra)
library(lme4)
library(multcomp)
library(multcompView) 
library(pheatmap)
library(purrr)
library(rcompanion)
library(stringr)
library(tidyr)
library(readr)
R.version

cat("\014")
rm(list=ls()) 
ls() 
setwd("~/Desktop/Kamei_Hostetler_et_al_2026//")
getwd()

# SBP transcription factors are candidate regulators of brace root initiation ####
## FigureS1C ####
cat("\014")
rm(list=ls()) 
primordia <- read.csv("~/Desktop/Kamei_Hostetler_et_al_2026/Data/Phenotype/BR_primordia.csv")
primordia <- primordia %>%
  mutate(
    stage = factor( # convert stage to factor
      stage,
      levels = c("v4","v6","v8","v10"), #relabel stage from lowercase to uppercase
      labels = c("V4","V6","V8","V10")),
    node = factor(node, levels = c("1","2","3")),
    plant = factor(plant))
primordia$plantID <- interaction(      # Assign unique plant ID within stage
  primordia$stage,
  primordia$plant)
model <- lmer(BR_primordia ~ stage * node + (1 | plantID), data = primordia) #Fit linear mixed effects model (lme4 package); measurements are not independent
anova(model)
emm <- emmeans(model, ~ node | stage) #estimates least-square means 
cld_df <- cld( #assigns letter displays to summarize pairwise comparisons 
  emm,
  adjust = "sidak",
  Letters = letters)
hsd_nodes_by_stage <- cld_df %>%
  dplyr::select(stage, node, .group) %>% #selects columns needed
  dplyr::mutate(.group = gsub(" ", "", .group)) %>% #removes extra spaces 
  dplyr::rename(groups = .group) #renames .group to group 
label_positions <- primordia %>% #determines where to plot letters in ggplot
  group_by(stage, node) %>%
  summarise(
    y = max(BR_primordia, na.rm = TRUE) + 0.5,
    .groups = "drop")
plot_letters <- left_join( #combine hsd and position data frames
  hsd_nodes_by_stage,
  label_positions,
  by = c("stage", "node"))
plot_letters$node <- factor( 
  plot_letters$node,
  levels = c("1","2","3"))
#pdf(file="/Users/ashley/Desktop/Kamei_Hostetler_et_al_2026/Manuscript/FigureS1C.pdf", width=6, height=4)
ggplot(primordia, aes(x = stage, y = BR_primordia, fill = node)) +
  geom_boxplot(alpha = 0.85, position = position_dodge(width = 0.75)) +
  geom_text(data = plot_letters, 
            aes(x = stage, y = y, label = groups, group = node),
            position = position_dodge(width = 0.75), vjust = 0, size = 4) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 1)) +
  scale_fill_manual(
    values = c("pink4", "slategray3", "burlywood"),
    breaks = c("3", "2", "1"),
    labels = c(
      "3" = "3 - Highest Node (Induction Stage)",
      "2" = "2 - Middle Node (Initiation Stage)",
      "1" = "1 - Lowest Node (Emergence Stage)")) +
  labs(
    x = "Vegetative Leaf Stage",
    y = "Brace Root Primordia (count)",
    fill = "Aboveground Node (BR Stage)") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold"))
#dev.off()

## B73 Transcriptomics ####
cat("\014")
rm(list = setdiff(ls(), c()))
annotation_data <- read.csv("Data/RNAseq/Zm-B73-REFERENCE-NAM-5.0-GMs-GOTerms.csv")
annotation_data <- annotation_data %>%
  distinct(gene_model, .keep_all = TRUE)
raw_counts <- read.csv("~/Desktop/Kamei_Hostetler_et_al_2026/Data/RNAseq/WholeNode/B73_RawCounts_v5.csv")

merged_data <- merge( #merge annotation file and RNA seq file
  raw_counts,
  annotation_data[, c("gene_model", "locus_symbol", "locus_name", "gene_products")],
  by.x = "Geneid",
  by.y = "gene_model",
  all.x = TRUE) %>%
  distinct(Geneid, .keep_all = TRUE)

### Prepare count matrix
count_data <- as.data.frame(raw_counts[, -1]) #remove the first column from raw counts
rownames(count_data) <- raw_counts$Geneid #assign the info from column 1 (above) to row names 

### Define sample sets
sample_sets <- list(
  NoPrimordia = c("No primordia4", "No primordia3", "No primordia2"),
  Primordia   = c("Primordia4", "Primordia3", "Primordia2"),
  Emergence   = c("Emergence4", "Emergence3", "Emergence2"))

### Filter low-expression genes
keep_genes <- rowSums(count_data > 1) >= 3 #must be greater to 1 in at least three samples
filtered_counts <- count_data[keep_genes, ]

### Create sample metadata
sample_names <- colnames(filtered_counts) #pull the column names from the data frame
condition_labels <- gsub("[0-9]", "", sample_names) #removes the number (biological replicate info) fro the sample names
condition_labels <- trimws(condition_labels)
condition_labels <- gsub("\\.", "_", condition_labels) #change No.primordia to No_primordia (change . to _)
stage_labels <- gsub("[^0-9]", "", sample_names) #pull the biological replicate info 
sample_metadata <- data.frame( #create a data frame with the condition and stage
  sampleName = sample_names,
  condition  = factor(condition_labels,
                      levels = c("No_primordia", "Primordia", "Emergence")),
  stage      = factor(stage_labels))
rownames(sample_metadata) <- sample_metadata$sampleName

### Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = filtered_counts,
  colData   = sample_metadata,
  design    = ~ condition)

### Run DESeq2
dds <- DESeq(dds)
vsd <- vst(dds, blind = TRUE) #variance-stablized matrix also normalized for library size; used for PCA, clustering, heatmpas, z-scores 
mat <- assay(vsd) #extract vsd matrix (matrix)
res <- results(dds) #res = differential expression statistics; based on normalized countes ; used for DE gene lists

### Extract DESeq2 normalized counts
norm_counts <- counts(dds, normalized = TRUE) #returns counts that have been normalized for library size
norm_counts_df <- as.data.frame(norm_counts) #norm_counts = raw counts adjusted for library size differences; used for exporting expression tables
norm_counts_df$Geneid <- rownames(norm_counts_df) #creates a column called Geneid from the rownames

### Rename normalized columns for clarity
colnames(norm_counts_df)[colnames(norm_counts_df) != "Geneid"] <-
  paste0(colnames(norm_counts_df)[colnames(norm_counts_df) != "Geneid"], "_norm")

## Differential Expression
contrast_list <- list(
  NP_vs_P = c("condition","No_primordia","Primordia"),
  P_vs_NP = c("condition","Primordia","No_primordia"),
  
  E_vs_P = c("condition","Emergence","Primordia"),
  P_vs_E = c("condition","Primordia","Emergence")
)
output_dir <- "Data/DE_results/B73_V6/"
dir.create(output_dir, showWarnings = FALSE)

for (contrast_name in names(contrast_list)) {
  res <- results(dds, contrast = contrast_list[[contrast_name]])
  res_df <- as.data.frame(res)
  res_df$Geneid <- rownames(res_df)
  
  res_final <- res_df %>%
    left_join(
      merged_data[, c("Geneid","locus_symbol","locus_name","gene_products")],
      by = "Geneid"
    ) %>%
    left_join(raw_counts, by = "Geneid") %>%
    left_join(norm_counts_df, by = "Geneid")
  
  write.csv(
    res_final,
    file.path(output_dir, paste0(contrast_name, "_results.csv")),
    row.names = FALSE
  )
  
  cat("Saved:", contrast_name, "\n")
}

### Table S3-S4: Candidate Genes ####
cat("\014")
rm(list = setdiff(ls(), c("res","norm_counts","vsd","mat","sample_metadata")))
# Read the datasets
Initiation_Candidates <-  read.csv("Data/DE_results/B73_V6/P_vs_NP_results.csv")
Emergence_Candidates <- read.csv("Data/DE_results/B73_V6/E_vs_P_results.csv")
All_data <- read.csv("~/Desktop/Kamei_Hostetler_et_al_2026/Data/RNAseq/WholeNode/B73_V5_gene_list.csv")
Initiation_Candidates <- Initiation_Candidates %>%
  filter((log2FoldChange > 1 | log2FoldChange < -1) & padj < 0.01 & baseMean > 1)
Emergence_Candidates <- Emergence_Candidates %>%
  filter((log2FoldChange > 1 | log2FoldChange < -1) & padj < 0.01 & baseMean > 1)
Initiation_Candidates <- Initiation_Candidates %>%
  left_join(
    All_data %>% dplyr::select(gene.ID, family),
    by = c("Geneid" = "gene.ID"))
Emergence_Candidates <- Emergence_Candidates %>%
  left_join(
    All_data %>% dplyr::select(gene.ID, family),
    by = c("Geneid" = "gene.ID"))
#write.csv(Initiation_Candidates, "~/Desktop/Kamei_Hostetler_et_al_2026/Manuscript/TableS3.csv", row.names = FALSE)
#write.csv(Emergence_Candidates, "~/Desktop/Kamei_Hostetler_et_al_2026/Manuscript/TableS4.csv", row.names = FALSE)

## Figure S2 - PCA Plot ####
pca_coords <- prcomp(t(mat))$x #transpose expression mattrix
pca <- prcomp(t(mat)) #run PCA

percent_var <- (pca$sdev^2 / sum(pca$sdev^2)) * 100 #calculate variance explained
pc1_var <- round(percent_var[1], 1)
pc2_var <- round(percent_var[2], 1)
pca_df <- data.frame(
  PC1 = pca_coords[,1], #extract sample coordinates
  PC2 = pca_coords[,2], #extract sample coordinates
  condition = sample_metadata$condition)
pca_df$condition <- factor(
  pca_df$condition,
  levels = c("Emergence","Primordia","No_primordia"),
  labels = c("Emergence","Primordia","No Primordia"))
#pdf(file="/Users/ashley/Desktop/Kamei_Hostetler_et_al_2026/Manuscript/FigS2.pdf", width=6, height=4)
ggplot(pca_df, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  scale_color_manual(values = c("pink4", "slategray3", "burlywood"),
                     breaks = c("No Primordia", "Primordia", "Emergence"),
                     labels = c(
                       "No Primordia" = "3 - Highest Node",
                       "Primordia" = "2 - Middle Node",
                       "Emergence" = "1 - Lowest Node")) +
  labs(x = "PC1 (56.3 PVE)", y = "PC2 (24.6 PVE)", color = "Aboveground Node") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size=10),
    axis.text.y = element_text(size=10),
    axis.text = element_text(size=10),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))
#dev.off()

## Figure 1B - Heatmap ####
#mat = variance-stablized matrix also normalized for library size
stage_means <- data.frame( #averaging replicates by stage
  "No primordia" = rowMeans(mat[, c("No.primordia4", "No.primordia3", "No.primordia2")]),
  "Primordia" = rowMeans(mat[, c("Primordia4", "Primordia3", "Primordia2")]),
  "Emergence" = rowMeans(mat[, c("Emergence4", "Emergence3", "Emergence2")]))
stage_means_z <- t(scale(t(stage_means))) #z-scores; relative expression patterns across stages (mean expression becomes 0 and SD becomes 1)
# Remove rows with NaN or Inf values in the normalized data
stage_means_z[is.infinite(stage_means_z)] <- 0  # Replace Inf with 0; genes with 0 variance across stages
stage_means_z[is.nan(stage_means_z)] <- 0       # Replace NaN with 0; genes with 0 variance across stages
# Calculate the maximum absolute Z-score value
max_abs_z <- max(abs(stage_means_z), na.rm = TRUE)
# Define breaks to have a midpoint at 0
breaks <- seq(-max_abs_z, max_abs_z, length.out = 101)
z_score_labels <- c("Min", "", "Max")
z_score_breaks <- c(-max_abs_z, 0, max_abs_z)
new_labels <- c(
  "Emergence" = "Emerging",
  "Primordia" = "Primordia",
  "No.primordia" = "No Primordia")
pheatmap(
  mat = stage_means_z,
  scale = "none",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  fontsize_row = 8,
  fontsize_col = 12,
  cluster_cols = FALSE,
  cluster_rows = TRUE, #genes with similar stage-specific behavior clsuter together
  angle_col = 45,
  show_rownames = FALSE,
  show_colnames = TRUE,
  cellwidth = 25,
  breaks = breaks,
  labels_col = new_labels[colnames(stage_means_z)],
  legend_labels = z_score_labels,
  legend_breaks = z_score_breaks,
  color_col = "black",
  fontface_col = "bold")
SBP <- Initiation_Candidates %>% #Filter for SBP gene family
  filter(str_detect(gene_products, "SBP")) %>%
  filter(baseMean > 10)
sbp_geneids <- SBP$Geneid
SBP <- SBP %>%
  mutate(
    No.primordia = rowMeans( #Takes the average of the three No primordia columns
      dplyr::select(
        ., #means it is selecting the current data frame (SBP %>% above)
        No.primordia4_norm,
        No.primordia3_norm,
        No.primordia2_norm),
      na.rm = TRUE), #ignores missing values when calculating mean
    Primordia = rowMeans(
      dplyr::select(
        .,
        Primordia4_norm,
        Primordia3_norm,
        Primordia2_norm),
      na.rm = TRUE),
    Emergence = rowMeans(
      dplyr::select(
        .,
        Emergence4_norm,
        Emergence3_norm,
        Emergence2_norm),
      na.rm = TRUE))
SBP_long <- SBP %>% #take wide format and make it long
  dplyr::select(
    locus_symbol,
    No.primordia,
    Primordia,
    Emergence
  ) %>%
  pivot_longer(
    cols = c(No.primordia, Primordia, Emergence),
    names_to = "Stage",
    values_to = "Mean_normalized_expression"
  )
gene_labels <- SBP$locus_symbol #character list of each locus symbol 
names(gene_labels) <- SBP$Geneid  #assigns gene ID to locus symbol 
row_labels <- ifelse(
  rownames(stage_means_z) %in% sbp_geneids, #Assign a rowname to zscore dataframe if the gene id is found in the SBP gene id 
  gene_labels[rownames(stage_means_z)],"")
#pdf(file="/Users/ashley/Desktop/Kamei_Hostetler_et_al_2026/Manuscript/Figure1B.pdf", width=4, height=6)
pheatmap(
  mat = stage_means_z,
  scale = "none",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  fontsize_row = 12,
  fontsize_col = 12,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  angle_col = 45,
  show_rownames = TRUE,
  show_colnames = TRUE,
  cellwidth = 25,
  breaks = breaks,
  labels_row = row_labels,
  labels_col = new_labels[colnames(stage_means_z)],
  legend_labels = z_score_labels,
  legend_breaks = z_score_breaks,
  color_col = "black",
  fontface_col = "bold")
#dev.off()
keep <- row_labels != "" #only keep rows that have labels with them
mat_small <- stage_means_z[keep, ] #subset expression matrix to only the rows we kept
labels_small <- row_labels[keep]
#pdf(file="/Users/ashley/Desktop/Kamei_Hostetler_et_al_2026/Manuscript/Figure1C.pdf", width=3, height=2)
pheatmap(
  mat_small,
  scale = "none",
  color = colorRampPalette(c("blue","white","red"))(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  labels_row = labels_small,
  labels_col = new_labels[colnames(mat_small)],
  breaks = breaks,
  legend_labels = z_score_labels,
  legend_breaks = z_score_breaks
)
#dev.off()

#Brace root initiation is increased in sbp mutants #####
##Figure 2 ####
cat("\014")
rm(list=ls()) 
setwd("~/Desktop/Kamei_Hostetler_et_al_2026/")
data <- read.csv("Data/Phenotype/2024_Field_ManualPhenotyping.csv")
data[] <- lapply(data, trimws)
data$Total_BR_Nodes_count = as.numeric(data$Total_BR_Nodes_count)
summary_table <- data %>% #Identify the mean and standard deviation for each genotype for BR Node (count)
  group_by(Genotype) %>%
  summarise(
    n = n(),
    Total_BR_Nodes_count_mean = mean(Total_BR_Nodes_count, na.rm = TRUE),
    Total_BR_Nodes_countR_sd   = sd(Total_BR_Nodes_count, na.rm = TRUE)
  ) %>%
  ungroup()
summary_table
data$Genotype <- factor(data$Genotype, levels = c("W22", "tsh4", "ub3", "ub2;ub3", "ub2;ub3;tsh4"))
methods <- data %>%
  group_by(Genotype, PlotNumber) %>%
  summarise(
    n_plants = n_distinct(PlantNumber), #Identify how many plants for each Plot (and Genotype) combo
    .groups = "drop" #Remove grouping (Genotype and Plot) after summarizing
  ) %>%
  group_by(Genotype) %>%  #Now group summarized data by Genotype
  summarise(
    n_plots = n(), #Identify how many plots for each genotype
    min_plants = min(n_plants), #Identify the least number of plants observed for any plot
    max_plants = max(n_plants), #Identify the most number of plants observed for any plot
    mean_plants = mean(n_plants), #Identify the average number of plants observed for any plot
    total_plants = sum(n_plants), #Total number of plants measured for this genotype 
    .groups = "drop") #Remove grouping
methods

###Figure 2B ####
attach(data)
lm <- lm(Total_BR_Nodes_count ~ Genotype) #Run a one-way linear model to test whetehre genotype affects the total number of brace root nodes
anova(lm) #Produce an ANOVA table to test the null hypothesis (H0: All genotype means are equal)
lm_aov=aov(lm)  #Convert to an ANOVA object
resid = residuals(object = lm_aov) #Extract model residuals
shapiro.test(x=resid) #Test whether residuals are normally distributed 
leveneTest(Total_BR_Nodes_count ~ Genotype) #Test for equal variances among genotypes
tuk <- TukeyHSD(lm_aov) #Perform Tukeys HSD test 
tuk_df <- as.data.frame(tuk$Genotype) 
tuk_df$comparison <- rownames(tuk_df) #Save genotype comparisons as a column 
tuk_letters <- multcompLetters4(lm_aov, tuk) #Generate letter displayes 
letters_df <- data.frame( #Convert the letter display into a tidy data frame 
  Genotype = names(tuk_letters$Genotype$Letters),
  Letters = tuk_letters$Genotype$Letters)
detach(data)

positions <- data %>%
  group_by(Genotype) %>% #Group each observation by genotype
  summarize(y = max(Total_BR_Nodes_count, na.rm = TRUE)) #For each genotype, identify the largest observed value of Total BR Nodes
letters_df <- left_join(letters_df, positions, by = "Genotype") #Merge Tukey significance letters with the y-axis positions 

Figure2B = ggplot(data,aes(x = Genotype, y = Total_BR_Nodes_count)) +
  geom_boxplot(fill = "gray1", outliers = TRUE) +
  labs(x = NULL, y = "Total BR Nodes (count)") +
  scale_y_continuous(limits = c(0, 10),breaks = seq(0, 10, by = 2))+
  geom_text(data = letters_df, aes(y = y + 1, label = Letters), size = 4) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=10),
    axis.ticks.x = element_blank(), 
    axis.text = element_text(size=10),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))
Figure2B

rm(list = setdiff(ls(), c("Figure2B")))

data <- read.csv("~/Desktop/Kamei_Hostetler_et_al_2026/Data/Phenotype/2024_Field_ImageAnalysis.csv")
data$CollectionDate = as.Date(data$CollectionDate,
                              format = "%m/%d/%Y")
data$PlantingDate = as.Date(data$PlantingDate,
                            format = "%m/%d/%Y")
data$DAP = as.numeric(data$CollectionDate - data$PlantingDate)
data$Genotype <- factor(data$Genotype, levels = c("W22", "tsh4", "ub3", "ub2;ub3", "ub2;ub3;tsh4"))
data[] <- lapply(data, trimws)
str(data)
colnames(data)
i = c(9:12)
data[,i] = apply(data[ ,i], 2,
                 function(x) as.numeric(x))
data$Total_BR_count = data$Total_BR_count*2
data$Tier1_count = data$Tier1_count*2
data$Tier2_count = data$Tier2_count*2
summary_table <- data %>%
  group_by(Genotype, CollectionDate, Node) %>% #Group by Genotype and Collection Date/Node
  summarise(
    n = n(),
    
    Total_BR_count_mean = mean(Total_BR_count, na.rm = TRUE), #summarise the mean BR number
    Total_BR_count_sd   = sd(Total_BR_count, na.rm = TRUE), #summarize the standard deviation of BR number
    
    Tier1_count_mean = mean(Tier1_count, na.rm = TRUE), #summarise the mean BR number on Tier 1
    Tier1_count_sd   = sd(Tier1_count, na.rm = TRUE), #summarize the standard deviation of BR number on Tier 1
    
    Tier2_count_mean = mean(Tier2_count, na.rm = TRUE), #summarise the mean BR number on Tier 2
    Tier2_count_sd   = sd(Tier2_count, na.rm = TRUE)) %>% #summarize the standard deviation of BR number on Tier 2
  ungroup()
summary_table
methods<- data %>%
  group_by(Genotype, DAP, Node, PlotNumber) %>% #Group by Genotype, plot, Node, DAP
  summarise(
    n_plants = n_distinct(PlantNumber),
    .groups = "drop"
  ) %>%
  group_by(Genotype, DAP, Node) %>% 
  summarise(
    n_plots = n(), #Identify how many plots for each genotype, DAP, Node combo
    min_plants = min(n_plants), #Identify the minimum number of plants in any plot for each genotype, DAP, Node combo
    max_plants = max(n_plants), #Identify the maximum number of plants in any plot for each genotype, DAP, Node combo
    mean_plants = mean(n_plants), #Identify the mean number of plants in any plot for each genotype, DAP, Node combo
    sd_plants = sd(n_plants), #Identify the standard deviation for number of plants in any plot for each genotype, DAP, Node combo
    total_plants = sum(n_plants), #Identify the total number of plants in all plots for each genotype, DAP, Node combo
    .groups = "drop")
methods
data_long <- data %>%
  pivot_longer(
    cols = c(Tier1_count, Tier2_count), #columns to combine into a single measurement column
    names_to = "Tier", #Name of the new column indicating the original source column 
    values_to = "BR_count") #Name of the new column containing the brace root counts
data_long$Genotype <- factor(data_long$Genotype, levels = c("W22", "tsh4", "ub3", "ub2;ub3", "ub2;ub3;tsh4")) #Order genotypes this way for plotting
str(data_long)
data_long$Node = as.character(data_long$Node)
###Figure 2C ####
nodes <- sort(unique(data_long$Node)) #Create a sorted list of all unique node values in the dataset
letters_df <- bind_rows( ## Perform a separate one-way ANOVA and Tukey HSD test for each node, then collect the significance letters into a single data frame.
  lapply(nodes, function(n){ #apply the following analysis to each node in the "nodes" vector (above)
    
    dat <- subset(data_long, Node == n) #Subset the data to only include observations from the current node 
    
    fit <- aov(Total_BR_count ~ Genotype, data = dat) #Fit a one-way ANOVA testing whether genotype affects the total brace root count at this node
    
    tuk <- TukeyHSD(fit) #Perform Tukey HSD to compre genotype pairs
    
    letters <- multcompLetters4(fit, tuk) #Generate compact letter dispaly for Tukey HSD results 
    
    data.frame( #Store results in a data frame 
      Node = n,
      Genotype = names(letters$Genotype$Letters),
      Letters = letters$Genotype$Letters)}))
positions <- data_long %>%
  group_by(Node, Genotype) %>% #Group by Node and Genotype
  summarize( 
    y = max(Total_BR_count, na.rm = TRUE), #summarize the largest value for total BR count
    .groups = "drop") #remove grouping after summarizing 
letters_df <- left_join( #Join the letters data frame and the positions data frame by Node and Genotype
  letters_df,
  positions,
  by = c("Node", "Genotype"))


Figure2C = ggplot(data_long,aes(x = Genotype, y = Total_BR_count, fill = Node)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  labs(x = NULL, y = "Total BR (count)") +
  scale_fill_manual(name = "Node",values = c("3" = "gray85","4" = "gray60"),
                    labels = c(
                      "3" = "Node 3",
                      "4" = "Node 4"))  +
  scale_y_continuous(limits = c(0, NA))+
  geom_text(data = letters_df, aes(x = Genotype, y = y + 1, label = Letters, group = factor(Node)),
            position = position_dodge(width = 0.75),inherit.aes = FALSE,size = 3,angle = 0) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.ticks.x = element_blank(), 
        axis.text = element_text(size=10),
        axis.title = element_text(face="bold", size=10),
        axis.title.x = element_text(face="bold", size=10),
        axis.title.y = element_text(face="bold", size=10))

Figure2C

data$Node = as.character(data$Node) #Make Node a charachter instead of number
data$Genotype <- factor(data$Genotype, levels = c("W22", "tsh4", "ub3", "ub2;ub3", "ub2;ub3;tsh4")) #Factor genotypes in this order for plotting

###Figure 2D ####
nodes <- sort(unique(data$Node))
letters_df_t1 <- bind_rows(
  lapply(nodes, function(n){
    
    dat <- subset(data, Node == n)
    
    fit <- aov(Tier1_count ~ Genotype, data = dat)
    
    tuk <- TukeyHSD(fit)
    
    letters <- multcompLetters4(fit, tuk)
    
    data.frame(
      Node = n,
      Genotype = names(letters$Genotype$Letters),
      Letters = letters$Genotype$Letters)}))
positions_t1 <- data %>%
  group_by(Node, Genotype) %>%
  summarize(
    y = max(Tier1_count, na.rm = TRUE),
    .groups = "drop")
letters_df_t1 <- left_join(
  letters_df_t1,
  positions_t1,
  by = c("Node", "Genotype"))

Figure2D = ggplot(data,aes(x = Genotype, y = Tier1_count, fill = Node)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  labs(x = NULL, y = "Tier 1 BR (count)")+
  scale_fill_manual(
    name = "Node",values = c("3" = "gray85","4" = "gray60"),
    labels = c(
      "3" = "Node 3",
      "4" = "Node 4"))  +
  scale_y_continuous(limits = c(0, NA))+
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.text = element_text(size=10),
        axis.title = element_text(face="bold", size=10),
        axis.title.x = element_text(face="bold", size=10),
        axis.title.y = element_text(face="bold", size=10))
Figure2D

####Figure 2E ####
nodes <- sort(unique(data$Node))
letters_df_t2 <- bind_rows(
  lapply(nodes, function(n){
    
    dat <- subset(data, Node == n)
    
    fit <- aov(Tier2_count ~ Genotype, data = dat)
    
    tuk <- TukeyHSD(fit)
    
    letters <- multcompLetters4(fit, tuk)
    
    data.frame(
      Node = n,
      Genotype = names(letters$Genotype$Letters),
      Letters = letters$Genotype$Letters)}))
positions_t2 <- data %>%
  group_by(Node, Genotype) %>%
  summarize(
    y = max(Tier2_count, na.rm = TRUE),
    .groups = "drop")
letters_df_t2 <- left_join(
  letters_df_t2,
  positions_t2,
  by = c("Node", "Genotype"))

Figure2E = ggplot(data,aes(x = Genotype, y = Tier2_count, fill = Node)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  labs(x = NULL, y = "Tier 2 BR (count)")+
  scale_fill_manual(
    name = "Node",values = c("3" = "gray85","4" = "gray60"),
    labels = c(
      "3" = "Node 3",
      "4" = "Node 4"))  +
  geom_text(data = letters_df_t2,
            aes(x = Genotype,y = y + 1,label = Letters,group = factor(Node)),position = position_dodge(width = 0.75),inherit.aes = FALSE,size = 3,angle = 0) +
  scale_y_continuous(limits = c(0, NA))+
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.text = element_text(size=10),
        axis.title = element_text(face="bold", size=10),
        axis.title.x = element_text(face="bold", size=10),
        axis.title.y = element_text(face="bold", size=10))
Figure2E

#pdf(file="/Users/ashley/Desktop/Kamei_Hostetler_et_al_2026/Manuscript/Figure2.pdf", width=6, height=4)
ggdraw() +
  draw_plot(Figure2D, x = 0, y = 0, width = 0.50, height = 0.5) +
  draw_plot(Figure2E, x = 0.50, y = 0, width = 0.50, height = 0.5) +
  draw_plot(Figure2B, x = 0.0, y = 0.5, width = 0.50, height = 0.5) +
  draw_plot(Figure2C, x = 0.50, y = 0.5, width = 0.50, height = 0.5) +
  draw_plot_label(label = c("B", "C", "D", "E"), 
                  size = 12,
                  x = c(0,0.5,0,0.5), 
                  y = c(1,1,0.55,0.55))
#dev.off()

#SBP TFs limit the size and competency of the stem pulvinus  ####
##Figure 3 #####
cat("\014")
rm(list=ls()) 
setwd("~/Desktop/Kamei_Hostetler_et_al_2026/")
data <- read.csv("~/Desktop/Kamei_Hostetler_et_al_2026/Data/Phenotype/2024_Field_ImageAnalysis.csv")
data$Genotype <- factor(data$Genotype, levels = c("W22", "tsh4", "ub3", "ub2;ub3", "ub2;ub3;tsh4")) #Order these genotypes for plotting
data[] <- lapply(data, trimws)
str(data)
colnames(data)
i = c(10:12)
data[,i] = apply(data[ ,i], 2,
                 function(x) as.numeric(x)) 
data$Total_BR_count = data$Total_BR_count*2 #Doubling the number of nodes because previous work has shown even distributions so side A = side B 
data$Tier1_count = data$Tier1_count*2 #Doubling the number of nodes because previous work has shown even distributions so side A = side B 
data$Tier2_count = data$Tier2_count*2 #Doubling the number of nodes because previous work has shown even distributions so side A = side B 

methods <- data %>%
  group_by(Genotype, Node, PlotNumber) %>%
  summarise(
    n_plants = n_distinct(PlantNumber),
    .groups = "drop"
  ) %>%
  group_by(Genotype, Node) %>%
  summarise(
    n_plots = n(),
    min_plants = min(n_plants),
    max_plants = max(n_plants),
    mean_plants = mean(n_plants),
    sd_plants = sd(n_plants),
    total_plants = sum(n_plants),
    .groups = "drop")
methods
data_long <- data %>%
  pivot_longer(
    cols = c(Tier1_count, Tier2_count),
    names_to = "Tier",
    values_to = "Value")
data_long$Genotype <- factor(data_long$Genotype, levels = c("W22", "tsh4", "ub3", "ub2;ub3", "ub2;ub3;tsh4"))
data = subset(data_long, Node == 4)
data = subset(data, Tier == "Tier2_count")
data$ID = paste(data$PlotNumber, data$PlantNumber, sep ="_") #Create a unique ID of Plot + Plant Number for easy merging
pheno = data

data <- read.csv("~/Desktop/Kamei_Hostetler_et_al_2026/Data/Phenotype/2024_Pulvinus.csv")
data$Genotype <- factor(data$Genotype, levels = c("W22", "tsh4", "ub3", "ub2;ub3", "ub2;ub3;tsh4"))
data[] <- lapply(data, trimws)
str(data)
i = c(8:10)
data[,i] = apply(data[ ,i], 2,
                 function(x) as.numeric(x))
data$Genotype <- factor(data$Genotype, levels = c("W22", "tsh4", "ub3", "ub2;ub3", "ub2;ub3;tsh4"))
data$ID = paste(data$Plot, data$Plant, sep ="_") #Create a unique ID of Plot + Plant Number for easy merging 
data = data[,c(11,1,10)]
data$pulvinus_cm = data$average * 2.54 #convert inches to cm
pheno = as.data.frame(pheno)
pheno = pheno[,c(14,5,13)]
colnames(pheno)[3] = "Tier2_count" #rename column 3 
data1 = merge(data, pheno, by = "ID")
data1 = data1[,c(1,2,4,6)]
colnames(data1)[2] = "Genotype"
data1 = subset(data1, Genotype != "W22") #Remove W22 data because W22 lacks a Tier 2

model <- lm(Tier2_count ~ pulvinus_cm, data = data1) #Fi a linear regression model to evaluate the relationship between pulvinus size and number of brace roots on tier 2
r2 <- summary(model)$r.squared #Extract the R2 value from model summary; R2 is the proportion of varitation in Tier2 brace root count by the size of pulvinus 
r2
###Figure 3B ####
Figure3B = ggplot(data1, aes(x=pulvinus_cm, y = Tier2_count))+
  geom_point(size = 2) +
  geom_smooth(method = "lm", se= FALSE)+
  labs(x = NULL, y = "Tier 2 BR (count)")+
  scale_y_continuous(limits = c(0, 15),breaks = seq(0, 15, by = 2))+
  theme_bw() +
  theme(
        axis.text.x = element_blank(),
        #axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.text = element_text(size=10),
        axis.title = element_text(face="bold", size=10),
        axis.title.x = element_text(face="bold", size=10),
        axis.title.y = element_text(face="bold", size=10))
Figure3B

###Figure 3A ####
attach(data)
lm <- lm(pulvinus_cm ~ Genotype) #Generate a linear model to test whether mean pulvinus height differs among genotypes
anova(lm) #convert linear model to an ANOVA object 
lm_aov=aov(lm) 
resid = residuals(object = lm_aov) #extract residuals 
shapiro.test(x=resid) #test whether residuals are normally distributed 
leveneTest(pulvinus_cm ~ Genotype) #test if there are equal variances among genotypes
tuk <- TukeyHSD(lm_aov) #Perform a Tukey HSD test to compare all genotype pairs while controlling for multiple comparisons
tuk_df <- as.data.frame(tuk$Genotype) 
tuk_df$comparison <- rownames(tuk_df) #Move genotype comparisons from row names to a separte column 
tuk_letters <- multcompLetters4(lm_aov, tuk) #Generate compact letter displays from Tukey HSD results
letters_df <- data.frame( 
  Genotype = names(tuk_letters$Genotype$Letters),
  Letters = tuk_letters$Genotype$Letters)
detach(data)

positions <- data %>%
  group_by(Genotype) %>%
  summarize(y = max(pulvinus_cm, na.rm = TRUE))
letters_df <- left_join(letters_df, positions, by = "Genotype")
Figure3A = ggplot(data, aes(x = Genotype, y = pulvinus_cm)) +
  geom_boxplot() +
  labs(x = NULL, y = "Pulvinus Thickness (cm)") +
  scale_y_continuous(limits = c(0, NA))+
  geom_text(data = letters_df, aes(y = y+0.1, label = Letters), size = 4) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    #axis.text.x = element_text(size=10, angle = 60, hjust=1),
    axis.text.y = element_text(size=10),
    axis.text = element_text(size=10),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))
Figure3A
Figure3B

#pdf(file="/Users/ashley/Desktop/Kamei_Hostetler_et_al_2026/Manuscript/Figure3.pdf", width=6, height=4)
ggdraw() +
  draw_plot(Figure3A, x = 0, y = 0, width = 0.5, height = 0.5) +
  draw_plot(Figure3B, x = 0.5, y = 0, width = 0.5, height = 0.5)
#dev.off()

#SBP regulatory networks during brace root initiation are distinct from those during ear development ####
## Table S5 ####
rm(list = setdiff(ls(), c()))
cat("\014")
data <- read.csv("Data/Grassius_UB3-DAP-seq.csv")
data = data[,c(3,4)]
data <- unique(data)
data2 <- read.csv("Data/TSH4_ChIP.csv")
data2 = as.data.frame(data2[,c(4)])
colnames(data2)[1] = "v4"
data2 <- unique(data2)
merged_data <- merge( #merged data is the overlap of the two datasets (CHIP and DAP seq); only what is shared is retained
  data, #UB3 DAP seq
  data2, #TSH4 CHIP seq
  by.x = "Target.Gene", 
  by.y = "v4")
data = merged_data #data is shared CHIP and DAP seq targets (ears); these are B73 v4 gene IDs
rm(list = setdiff(ls(), c("data")))

annot <- read.csv("Data/RNAseq/B73_v5_B73_v4.csv")
annot = annot[,c(1,2)]
data <- merge( #merge & retain only genes that have a V5 gene ID
  data,
  annot,
  by.x = "Target.Gene",
  by.y = "B73_v4")

data1 <- read.csv("Data/DE_results/B73_V6/P_vs_NP_results.csv")
data1 <- data1 %>% #retain initiation targets that are significant
  filter((log2FoldChange > 1 | log2FoldChange < -1) & padj < 0.01 & baseMean > 1)
data <- merge( #Merge and retain data frame (CHIPseq-ear, DAPseq-ear, & v5 gene id) w/ significant initiation candidates
  data,
  data1,
  by.x = "B73_v5",
  by.y = "Geneid")
table(ifelse(data$log2FoldChange > 0, "Upregulated", "Downregulated"))
upregulated <- data %>%
  filter(log2FoldChange > 0)
downregulated <- data %>%
  filter(log2FoldChange < 0)
upregulated$Regulation = "Upregulated"
downregulated$Regulation = "Downregulated"
TableS5 = rbind(upregulated, downregulated)
#write.csv(TableS5, file = "Manuscript/TableS5.csv", row.names = FALSE)
targets = data[,c(1,2,3,5,9,10,11)] #targets data frame is CHIPseq-ear w/ B73 V5 gene ID, DAPseq-ear w/ B73 V5 gene ID, & initiation candidates

rm(list = setdiff(ls(), c("targets")))
cat("\014")
## Table S6-S11 ####
# Gene Conversion Files
B73_v5_v3 <- read.csv("Data/RNAseq/B73_v5_B73_v3.csv")
B73_v5_v3 = B73_v5_v3[,c(1,2)]
W22_B73_v5 <- read.csv("Data/RNAseq/W22_B73_v5.csv")
W22_B73_v5 = W22_B73_v5[,c(1,2)]
B73_v5_B73_v4 <- read.csv("Data/RNAseq/B73_v5_B73_v4.csv")
B73_v5_B73_v4 = B73_v5_B73_v4[,c(1,2)]
data <- merge(B73_v5_v3, B73_v5_B73_v4, by.x = "v5", by.y = "B73_v5", all = TRUE)
data <- merge(data, W22_B73_v5, by.x = "v5", by.y = "B73_v5", all = TRUE) #data is a dataframe that contains all the matched gene ids for B73 V3, B73 V4, B73 V5, and W22


annotation_data <- read.csv("Data/RNAseq/Zm-B73-REFERENCE-NAM-5.0-GMs-GOTerms.csv", header = TRUE)
annotation_data <- annotation_data %>% #merge Gene ID file with annotation file
  left_join(data, by = c("gene_model" = "v5"))
annotation_data <- annotation_data %>% distinct(gene_model, .keep_all = TRUE) #removes duplicate gene annotation entries based on gene_model (B73 V5)

# Load raw count data
raw_counts <- read.csv("Data/RNAseq/Mutants/mutantRNAseq_rawcount.csv", header = TRUE) 
merged_data <- merge(raw_counts, annotation_data[, c("W22","gene_model", "B73_v4", "v3", "locus_symbol", "locus_name", "gene_products")], 
                     by.x = "Geneid", by.y = "W22", all.x = TRUE) # Merge raw data with annotation
merged_data <- merged_data %>% #merged data is mutant RNAseq + gene conversion file + annotation file 
  dplyr::rename(B73_v5 = gene_model) #rename B73_v5 to gene_modek
merged_data <- merged_data %>% distinct(Geneid, .keep_all = TRUE) #keeps only one occurance of W22

# Prepare count data for DESeq2
gene_ids <- raw_counts$Geneid #extract W22 Gene IDs from RNAseq data
count_data <- raw_counts[, -1] #remove GeneID (W22) column from data frame
rownames(count_data) <- gene_ids #Assign GeneID (W22) to rownames 

# Define sample sets for different genotypes and conditions
sample_sets <- list(
  tsh4_set1 = c("tsh4_P1S1", "tsh4_P5S1", "tsh4_P4S1", "tsh4_P3S1"),
  tsh4_set2 = c("tsh4_P5S2", "tsh4_P4S2", "tsh4_P3S2", "tsh4_P1S2"),
  ub2ub3_set1 = c("ub2.ub3_P1S1", "ub2.ub3_P5S1", "ub2.ub3_P4S1", "ub2.ub3_P3S1"),
  ub2ub3_set2 = c("ub2.ub3_P5S2", "ub2.ub3_P4S2", "ub2.ub3_P3S2", "ub2.ub3_P1S2"),
  W22_set1 = c("W22_P6S1", "W22_P3S1", "W22_P2S1", "W22_P1S1"),
  W22_set2 = c("W22_P6S2", "W22_P3S2", "W22_P2S2", "W22_P1S2"))

# Filter genes: retain genes expressed with CPM > 1 in at least 4 samples per set
keep_genes <- rowSums(count_data[, sample_sets$tsh4_set1, drop = FALSE] > 1) >= 4 |
  rowSums(count_data[, sample_sets$tsh4_set2, drop = FALSE] > 1) >= 4 |
  rowSums(count_data[, sample_sets$ub2ub3_set1, drop = FALSE] > 1) >= 4 |
  rowSums(count_data[, sample_sets$ub2ub3_set2, drop = FALSE] > 1) >= 4 |
  rowSums(count_data[, sample_sets$W22_set1, drop = FALSE] > 1) >= 4 |
  rowSums(count_data[, sample_sets$W22_set2, drop = FALSE] > 1) >= 4

filtered_counts <- count_data[keep_genes, ] #count_data is mutant RNAseq raw counts with genes assigned to rownames; subsetting data to only those that have more than 1 count in all 4 samples
filtered_gene_ids <- gene_ids[keep_genes] #list of genes with more than 1 count in all 4 samples

# Create metadata for samples
sample_names <- colnames(filtered_counts) #extract column names (which are samples) from the dataset
genotypes <- gsub("_.*", "", sample_names) #remove everything after the _ and replace with nothing (leaves just genotypes and not biological reps or tier info)
sets <- ifelse(grepl("S1$", sample_names), "S1", "S2") #grepl searches each element of sample name for a pattern; S1$ means that S1 at the end of the string (Does this sample name end with S1), if TRUE it assings S1, if FALSE it assigns S2
conditions <- factor(paste(genotypes, sets, sep = "_")) #this combines the genotype and sets indormation into a condition and converts it to a factor 

sample_metadata <- data.frame(sampleName = sample_names, condition = conditions)
rownames(sample_metadata) <- sample_metadata$sampleName

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = filtered_counts,
                              colData = sample_metadata,
                              design = ~ condition)

# Normalize and perform differential expression analysis
dds <- DESeq(dds)

# Rename normalized count columns so they are clearly labeled
norm_counts <- as.data.frame(counts(dds, normalized = TRUE)) #extracts normalized counts from dds object and converts them to a df
colnames(norm_counts) <- paste0("Normalized_", colnames(norm_counts)) #adds the prefix Normalized_ to each column name 
norm_counts$Geneid <- rownames(norm_counts) #This line creates a new column called Geneid and fills it with the current row names of the norm_counts data frame.

# Define contrasts
contrast_list <- list(
  tsh4_vs_W22_S1 = c("condition", "tsh4_S1", "W22_S1"), #Notion: factor, numerator, denominator
  tsh4_vs_W22_S2 = c("condition", "tsh4_S2", "W22_S1"),
  ub2ub3_vs_W22_S1 = c("condition", "ub2.ub3_S1", "W22_S1"),
  ub2ub3_vs_W22_S2 = c("condition", "ub2.ub3_S2", "W22_S1"))

# Create output directory
output_dir <- "Data/DE_results/MutantResults" #folder where all results will be saved
dir.create(output_dir, showWarnings = FALSE) #creates this folder if it doesnt exist, showWarnings = FALSE : prevents R from printing a warning if this folder already exists


fdr_threshold <- 0.01 
logFC_threshold <- 1 #absolute log fold change has to be greater than 1 (2 fold increase or decrease)

for (contrast_name in names(contrast_list)) { #loop though every contrast for each comparison in contrast_list
  contrast <- contrast_list[[contrast_name]] #retrieve the current contrast in the loop
  res <- results(dds, contrast = contrast) #run DEseq2
  
  res_df <- as.data.frame(res) #convert to a data frame
  res_df$Geneid <- rownames(res_df) #save gene IDs as a column
  
  res_with_annotation <- merge( #add gene annotaton information
    res_df,
    merged_data[, c("Geneid","B73_v5","B73_v4","v3", #merged data contains gene conversion file + annotation file 
                    "locus_symbol","locus_name","gene_products")],
    by = "Geneid", #merge by Geneid (W22)
    all.x = TRUE #keep all res_df rows (DEseq2)
  )
  
  res_with_annotation <- merge( #Add normalized counts for every sample
    res_with_annotation,
    norm_counts,
    by = "Geneid",
    all.x = TRUE
  )
  
  res_with_annotation <- res_with_annotation[!duplicated(res_with_annotation$Geneid), ] #remove duplicate genes based on W22 gene ID
  
  # Export all genes
  write.csv( 
    res_with_annotation,
    file.path(output_dir, paste0(contrast_name, "_all_genes.csv")),
    row.names = FALSE
  )
  
  # Keep only significant DEGs
  sig_res <- res_with_annotation %>%
    filter(
      !is.na(padj),
      padj < fdr_threshold,
      abs(log2FoldChange) > logFC_threshold
    )
  
  # Export significant DEGs only
  write.csv(
    sig_res,
    file.path(output_dir, paste0(contrast_name, "_significant_DEGs.csv")),
    row.names = FALSE
  )
}

#Set filtering thresholds
fdr_threshold <- 0.01
logFC_threshold <- 1

#Create a list of files
files <- list(
  list(file = "Data/DE_results/MutantResults/tsh4_vs_W22_S1_all_genes.csv",
       title = "tsh4_vs_W22_S1"),
  list(file = "Data/DE_results/MutantResults/tsh4_vs_W22_S2_all_genes.csv",
       title = "tsh4_vs_W22_S2"),
  list(file = "Data/DE_results/MutantResults/ub2ub3_vs_W22_S1_all_genes.csv",
       title = "ub2ub3_vs_W22_S1"),
  list(file = "Data/DE_results/MutantResults/ub2ub3_vs_W22_S2_all_genes.csv",
       title = "ub2ub3_vs_W22_S2"))

results <- lapply(files, function(x) { #loops over every element of files and returns a list 
  
  read_csv(x$file, show_col_types = FALSE) %>% #read each csv in the iteration 
    filter( #filter for significant genes
      padj < fdr_threshold,
      abs(log2FoldChange) > logFC_threshold
    ) %>%
    mutate( #add regulation column information
      Regulation = if_else(
        log2FoldChange > 0,
        "Upregulated",
        "Downregulated"
      ),
      Comparison = x$title
    )
  
})

names(results) <- sapply(files, `[[`, "title") #Assign names to each element in my results list

#Extract each comparison results 
tsh4_vs_W22_S1 <- results[["tsh4_vs_W22_S1"]]
tsh4_vs_W22_S2 <- results[["tsh4_vs_W22_S2"]]
ub2ub3_vs_W22_S1 <- results[["ub2ub3_vs_W22_S1"]]
ub2ub3_vs_W22_S2 <- results[["ub2ub3_vs_W22_S2"]]

##Table S9 ####
shared_Geneid <- Reduce( #identify the DEGs from each list that are shared among all four lists
  intersect,
  list(
    tsh4_vs_W22_S1$Geneid,
    tsh4_vs_W22_S2$Geneid,
    ub2ub3_vs_W22_S1$Geneid,
    ub2ub3_vs_W22_S2$Geneid))
length(shared_Geneid)
shared_Geneid <- data.frame(W22 = shared_Geneid) #List of genes (W22) shared in all 4 DE sets

colnames(targets)[4] = "Initiation_log2FoldChange"
colnames(targets)[5] = "Initiation_FDR"

shared_Geneid <- shared_Geneid %>% #Add B73 V5 gene conversion for W22 to the shared_Geneid data frame (shared among four DEG lists)
  left_join(W22_B73_v5, by = c("W22" = "W22"))
X = unique(shared_Geneid$B73_v5) #X = number of B73 genes in the shared DEG list

genes_of_interest <- shared_Geneid$W22 #W22 genes that are shared among all 4 data DEG datasets

comparison_tables <- lapply(names(results), function(comp){ #creates a list of summary tables (one for each DEseq comparison)
  
  results[[comp]] %>%
    dplyr::filter(Geneid %in% genes_of_interest) %>%
    dplyr::select(Geneid, log2FoldChange, padj, Regulation) %>%
    dplyr::rename_with(
      ~ paste0(comp, "_", .x),
      c(log2FoldChange, padj, Regulation)
    )
  
})
merged_genes <- reduce(comparison_tables, full_join, by = "Geneid") #merge esch list by GeneID
head(shared_Geneid)
final_table <- shared_Geneid %>% #adding B73 conversion to the merged_genes data frame
  left_join( #shared_Geneid is the shared DEG among 4 DESeq 
    merged_genes, #merged genes is the FC, P value, and Regulation info for shared_Geneid 
    by = c("W22" = "Geneid")
  )
X = unique(final_table$B73_v5) #X is th enumber of unique B73 V5 genes identified
final_table <- final_table %>%
  dplyr::select(-W22) %>% #Remove W22 and keep only one occurance of B73 V5
  dplyr::distinct(B73_v5, .keep_all = TRUE)

# Load annotation data
annotation_data <- read.csv("Data/RNAseq/Zm-B73-REFERENCE-NAM-5.0-GMs-GOTerms.csv", header = TRUE)
final_table1 <- final_table %>% #Add annotation information to the final table (which is the 4 DEG lists)
  left_join(annotation_data,
            by = c("B73_v5" = "gene_model"))
final_table1 <- final_table1 %>%
  filter(!is.na(B73_v5)) %>% #remove rows that have an NA in the B73 V5 column
  distinct(B73_v5, .keep_all = TRUE) #Keep only one row per uniwu B73 V5 column 
#write.csv(final_table1, file = "Manuscript/TableS10.csv", row.names = FALSE)
TableS10 = final_table

##Table S11 ####
X = unique(targets$B73_v5) #X is the number of unique B73 V5 gene IDs in the targets data frame (intersection of CHIP, DAP, and Initiation data)
shared_Geneid = merge(shared_Geneid, targets, by = "B73_v5")  #Merges the 4 mutant RNA seq data set with the targets data set
genes_of_interest <- shared_Geneid$W22 #identify the genes that were shared among all data sets (W22)
comparison_tables <- lapply(names(results), function(comp){ #loop through each comparison
  
  results[[comp]] %>% #Access each comparison in the results list (mutant DESeq) 
    dplyr::filter(Geneid %in% genes_of_interest) %>% #keep only genes of interest (shared in all datasets)
    dplyr::select(Geneid, log2FoldChange, padj, Regulation) %>% #Keep only these columns: Geneid, log2FoldChange, padj, Regulation
    dplyr::rename_with( #Rename the columns to include which comparison it is coming from 
      ~ paste0(comp, "_", .x),
      c(log2FoldChange, padj, Regulation)
    )
  
})

merged_genes <- reduce(comparison_tables, full_join, by = "Geneid") #turns the individual comparison tables into one large table matched by Gene ID
final_table <- shared_Geneid %>% #shared_Geneid (all datasets merged) with only gene info
  left_join(
    merged_genes, #merged_genes contains log2FC, Padj, Regulation info
    by = c("W22" = "Geneid")
  )
#write.csv(final_table, file = "Manuscript/TableS11.csv", row.names = FALSE)

## Table S12 ####
data1 <- read.csv("Data/DE_results/B73_V6/P_vs_NP_results.csv")
data1 <- data1 %>% #Selecting only significant genes for initiation candidates
  filter((log2FoldChange > 1 | log2FoldChange < -1) & padj < 0.01 & baseMean > 1)
TableS10 <- merge(  #Merging/Overlapping initiation candidates with mutant DEGs
  data1, #Inititation Candidates
  TableS10, #TableS10 is mutant DEG
  by.x = "Geneid",
  by.y = "B73_v5")
X = unique(TableS10$Geneid) #X is the number of unique B73 V5 genes from the overlap of Initiation candidates + mutant RNAseq
table(ifelse(TableS10$log2FoldChange > 0, "Upregulated", "Downregulated")) #Identify how many genes are up and downregulated
TableS10 <- TableS10 %>% #Assigns regulation to initiation column
  mutate(
    Initiation_Regulation = case_when(
      log2FoldChange > 0 ~ "Upregulated",
      log2FoldChange < 0 ~ "Downregulated",
      TRUE ~ "No Change"))
TableS10 <- TableS10 %>% #Keep only genes where the direction of regulation matches across all four mutant comparisons + initiation comparisons
  filter(
    tsh4_vs_W22_S1_Regulation == Initiation_Regulation &
      tsh4_vs_W22_S2_Regulation == Initiation_Regulation &
      ub2ub3_vs_W22_S1_Regulation == Initiation_Regulation &
      ub2ub3_vs_W22_S2_Regulation == Initiation_Regulation)
TableS10 %>%
  summarize( #creates a new summary table 
    Up_all = sum( #Counts all the genes that were upregulated in all the comparisons below
      Initiation_Regulation == "Upregulated" &
        tsh4_vs_W22_S1_Regulation == "Upregulated" &
        tsh4_vs_W22_S2_Regulation == "Upregulated" &
        ub2ub3_vs_W22_S1_Regulation == "Upregulated" &
        ub2ub3_vs_W22_S2_Regulation == "Upregulated"
    ),
    Down_all = sum( #counts all the genes that were downregulated in the comparisons below
      Initiation_Regulation == "Downregulated" &
        tsh4_vs_W22_S1_Regulation == "Downregulated" &
        tsh4_vs_W22_S2_Regulation == "Downregulated" &
        ub2ub3_vs_W22_S1_Regulation == "Downregulated" &
        ub2ub3_vs_W22_S2_Regulation == "Downregulated"
    )
  )
#write.csv(TableS10, file = "Manuscript/TableS12.csv", row.names = FALSE)


#SBP TFs inhibit brace root initiation by inhibiting cell division by limiting ethylene####
## Figure 4 ####
cat("\014")
rm(list=ls()) 
setwd("~/Desktop/Kamei_Hostetler_et_al_2026/Data/Phenotype/")
data <- read.csv("2025_GH_HormoneData.csv")
str(data)
colnames(data)
i = c(9:16)
data[,i] = apply(data[ ,i], 2,
                 function(x) as.numeric(x))
data$Genotype.1 <- factor(data$Genotype.1, levels = c("W22", "tsh4", "ub2;ub3"))
# Hormone columns
colnames(data)
hormones <- colnames(data)[9:16]
# Convert to long format
data_long <- data %>%
  pivot_longer(
    cols = all_of(hormones),
    names_to = "Hormone",
    values_to = "Value"
  )
letters_df <- map_df(hormones, function(horm){ #loops through every hormone 
  
  tmp <- data[, c("Genotype.1", horm)] #extracts genotype and hormone information
  colnames(tmp) <- c("Genotype", "Value") #Renames columns
  
  mod <- aov(Value ~ Genotype, data = tmp) #Runs an ANOVA on the new dataframe 
  
  hsd <- HSD.test(mod, "Genotype", group = TRUE) #Performs Tukey HSD test and assigns letters (group = TRUE)
  
  groups <- as.data.frame(hsd$groups) #stores as a dataframe
  groups$Genotype <- rownames(groups) 
  
  data.frame(
    Genotype = groups$Genotype,
    Letters = groups$groups,
    Hormone = horm
  )
})
letter_pos <- data_long %>% #Calculate letter positions for plotting based on the maximum value within a group 
  group_by(Hormone, Genotype.1) %>%
  summarise(
    y = max(Value, na.rm = TRUE) * 1.1,
    .groups = "drop"
  )
letters_df <- left_join( #Merge letter information and position information into a single dataframe
  letters_df,
  letter_pos,
  by = c("Hormone",
         "Genotype" = "Genotype.1")
)
#pdf(file="/Users/ashley/Desktop/Kamei_Hostetler_et_al_2026/Manuscript/Figure4.pdf", width=6, height=4)
ggplot(data_long,aes(x = Genotype.1,y = Value)) +
  geom_boxplot() +
  geom_text(data = letters_df,aes(x = Genotype,y = y,label = Letters),inherit.aes = FALSE,size = 4) +
  facet_wrap(~ Hormone,nrow = 2,ncol = 4,scales = "free_y") +
  scale_y_continuous(limits = c(0, NA))+
  labs(x = NULL, y = "ng/g FW")+
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size=10, angle = 45, hjust=1),
        axis.text.y = element_text(size=10),
        axis.text = element_text(size=10),
        axis.title = element_text(face="bold", size=10),
        axis.title.x = element_text(face="bold", size=10),
        axis.title.y = element_text(face="bold", size=10))
#dev.off()

#Gibberellin promotes brace root growth after initiation####
## Table S13 ####
cat("\014")
rm(list=ls()) 
setwd("~/Desktop/Kamei_Hostetler_et_al_2026/")
# Gene Conversion Files
B73_v5_v3 <- read.csv("Data/RNAseq/B73_v5_B73_v3.csv")
B73_v5_v3 = B73_v5_v3[,c(1,2)]
W22_B73_v5 <- read.csv("Data/RNAseq/W22_B73_v5.csv")
W22_B73_v5 = W22_B73_v5[,c(1,2)]
B73_v5_B73_v4 <- read.csv("Data/RNAseq/B73_v5_B73_v4.csv")
B73_v5_B73_v4 = B73_v5_B73_v4[,c(1,2)]
data <- merge(B73_v5_v3, B73_v5_B73_v4, by.x = "v5", by.y = "B73_v5", all = TRUE)
data <- merge(data, W22_B73_v5, by.x = "v5", by.y = "B73_v5", all = TRUE)

# Load annotation data
annotation_data <- read.csv("Data/RNAseq/Zm-B73-REFERENCE-NAM-5.0-GMs-GOTerms.csv", header = TRUE)
annotation_data <- annotation_data %>% #merge annotation information with gene conversion file
  left_join(data, by = c("gene_model" = "v5"))
annotation_data <- annotation_data %>% distinct(gene_model, .keep_all = TRUE) #keep only a single instance of B73 V5 

# Load raw count data
raw_counts <- read.csv("Data/RNAseq/Mutants/mutantRNAseq_rawcount.csv", header = TRUE) #Raw RNAseq mutantdata 

# Merge raw data with annotation
merged_data <- merge(raw_counts, annotation_data[, c("W22","gene_model", "B73_v4", "v3", "locus_symbol", "locus_name", "gene_products")], 
                     by.x = "Geneid", by.y = "W22", all.x = TRUE)
merged_data <- merged_data %>%
  dplyr::rename(B73_v5 = gene_model)
merged_data <- merged_data %>% distinct(Geneid, .keep_all = TRUE)

# Prepare count data for DESeq2
gene_ids <- raw_counts$Geneid
count_data <- raw_counts[, -1]
rownames(count_data) <- gene_ids

# Define sample sets for different genotypes and conditions
sample_sets <- list(
  tsh4_set1 = c("tsh4_P1S1", "tsh4_P5S1", "tsh4_P4S1", "tsh4_P3S1"),
  ub2ub3_set1 = c("ub2.ub3_P1S1", "ub2.ub3_P5S1", "ub2.ub3_P4S1", "ub2.ub3_P3S1"))
sample_sets

# Filter genes: retain genes expressed with CPM > 1 in at least 4 samples per set
keep_genes <- rowSums(count_data[, sample_sets$tsh4_set1, drop = FALSE] > 1) >= 4 |
  rowSums(count_data[, sample_sets$ub2ub3_set1, drop = FALSE] > 1) >= 4
filtered_counts <- count_data[keep_genes, ]
filtered_gene_ids <- gene_ids[keep_genes]

# Create metadata for samples
sample_names <- colnames(filtered_counts)
genotypes <- gsub("_.*", "", sample_names)
sets <- ifelse(grepl("S1$", sample_names), "S1", "S2")
conditions <- factor(paste(genotypes, sets, sep = "_"))

sample_metadata <- data.frame(sampleName = sample_names, condition = conditions)
rownames(sample_metadata) <- sample_metadata$sampleName

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = filtered_counts,
                              colData = sample_metadata,
                              design = ~ condition)

# Normalize and perform differential expression analysis
dds <- DESeq(dds)

# Rename normalized count columns so they are clearly labeled
norm_counts <- as.data.frame(counts(dds, normalized = TRUE))
colnames(norm_counts) <- paste0("Normalized_", colnames(norm_counts))
norm_counts$Geneid <- rownames(norm_counts)

# Define contrasts
contrast_list <- list(
  tsh4_vs_ub2ub3_S1 = c("condition", "tsh4_S1", "ub2.ub3_S1"))

# Create output directory
output_dir <- "Data/DE_results/MutantResults"
dir.create(output_dir, showWarnings = FALSE)


fdr_threshold <- 0.01
logFC_threshold <- 1

for (contrast_name in names(contrast_list)) {
  contrast <- contrast_list[[contrast_name]]
  res <- results(dds, contrast = contrast)
  
  res_df <- as.data.frame(res)
  res_df$Geneid <- rownames(res_df)
  
  res_with_annotation <- merge(
    res_df,
    merged_data[, c("Geneid","B73_v5","B73_v4","v3",
                    "locus_symbol","locus_name","gene_products")],
    by = "Geneid",
    all.x = TRUE
  )
  
  res_with_annotation <- merge(
    res_with_annotation,
    norm_counts,
    by = "Geneid",
    all.x = TRUE
  )
  
  res_with_annotation <- res_with_annotation[!duplicated(res_with_annotation$Geneid), ]
  
  # Export all genes
  write.csv(
    res_with_annotation,
    file.path(output_dir, paste0(contrast_name, "_all_genes.csv")),
    row.names = FALSE
  )
  
  # Keep only significant DEGs
  sig_res <- res_with_annotation %>%
    filter(
      !is.na(padj),
      padj < fdr_threshold,
      abs(log2FoldChange) > logFC_threshold
    )
  
  # Export significant DEGs only
  write.csv(
    sig_res,
    file.path(output_dir, paste0(contrast_name, "_significant_DEGs.csv")),
    row.names = FALSE
  )
}

#Generate same info from volcano plot
fdr_threshold <- 0.01
logFC_threshold <- 1

files <- list(
  list(file = "Data/DE_results/MutantResults/tsh4_vs_ub2ub3_S1_all_genes.csv",
       title = "tsh4_vs_ub2ub3_S1"))
results <- lapply(files, function(x) {
  
  read_csv(x$file, show_col_types = FALSE) %>%
    filter(
      padj < fdr_threshold,
      abs(log2FoldChange) > logFC_threshold
    ) %>%
    mutate(
      Regulation = if_else(
        log2FoldChange > 0,
        "Upregulated",
        "Downregulated"
      ),
      Comparison = x$title
    )
  
})
names(results) <- sapply(files, `[[`, "title")
tsh4_vs_ub2ub3_S1 <- results[["tsh4_vs_ub2ub3_S1"]]
data = results$tsh4_vs_ub2ub3_S1
data1 <- read.csv("Data/DE_results/B73_V6/P_vs_NP_results.csv")
data1 <- data1 %>%
  filter((log2FoldChange > 1 | log2FoldChange < -1) & padj < 0.01 & baseMean > 1)
data1 = data1[,c(7,2,6)]
colnames(data1)[2] = "log2FoldChange_Initiation"
colnames(data1)[3] = "padj_Initiation"
data <- merge( #Merge mutant RNAseq (significant DEGs) with Initiation candidates (significant DEGs)
  data,
  data1,
  by.x = "B73_v5",
  by.y = "Geneid")
X = unique(data$B73_v5) #Number of unique B73 V5 DEGs shared in mutants and initiation candidates
table(ifelse(data$log2FoldChange_Initiation > 0, "Upregulated", "Downregulated")) #Regulation information for initiation candidates
table(ifelse(data$log2FoldChange > 0, "Upregulated", "Downregulated")) #Regulation information for mutant comparisons 
X = unique(data$B73_v5) #Number of DEGs shared between tsh4/ub2ub3S1 + initiation 
#write.csv(data, file = "Manuscript/TableS13.csv", row.names = FALSE)

##Fig 5A ####
cat("\014")
rm(list=ls()) 
data = read.csv("Data/Phenotype/2023_FieldDataBRPhenos.csv", header = TRUE, na.strings = "NA")
data <- data %>%
  dplyr::rename(PlantID = Plot) %>% #PlantID becomes Plot
  tidyr::separate( #Separate PlantID into Plot and Plant component
    PlantID, 
    into = c("Plot", "Plant"),
    sep = "-",
    remove = FALSE,
    convert = TRUE)
data <- data %>% #Rename ub2,ub3 to ub2;ub3
  mutate(
    Genotype = dplyr::recode(
      Genotype,
      "ub2,3" = "ub2;ub3"))
methods <- data %>%
  group_by(Genotype, Plot) %>%
  summarise(
    n_plants = n_distinct(Plant),
    .groups = "drop"
  ) %>%
  group_by(Genotype) %>%
  summarise(
    n_plots = n(),
    min_plants = min(n_plants),
    max_plants = max(n_plants),
    mean_plants = mean(n_plants),
    total_plants = sum(n_plants),
    .groups = "drop")
data <- data %>% #Calculate hte average brace root diameter across the 4 diameters measured
  mutate(
    Average_BR_Diameter = rowMeans(across(BR_diameter_1_mm:BR_diameter_4_mm), na.rm = TRUE))
data$Genotype <- factor(data$Genotype, levels = c("tsh4", "ub2;ub3"))

attach(data)
head(data)
lm <- lm(Average_BR_Diameter ~ Genotype)
anova(lm)
lm_aov=aov(lm) 
resid = residuals(object = lm_aov)
shapiro.test(x=resid)
leveneTest(Average_BR_Diameter ~ Genotype)
tuk <- TukeyHSD(lm_aov)
tuk_df <- as.data.frame(tuk$Genotype)
tuk_df$comparison <- rownames(tuk_df)
tuk_letters <- multcompLetters4(lm_aov, tuk)
letters_df <- data.frame(
  Genotype = names(tuk_letters$Genotype$Letters),
  Letters = tuk_letters$Genotype$Letters)
letters_df
detach(data)
positions <- data %>%
  group_by(Genotype) %>%
  summarize(y = max(Average_BR_Diameter, na.rm = TRUE))
letters_df <- left_join(letters_df, positions, by = "Genotype")

Figure5A = ggplot(data, aes(x = Genotype, y = Average_BR_Diameter)) +
  geom_boxplot() +
  labs(x = NULL, y = "BR Diameter (mm)") +
  scale_y_continuous(limits = c(0, NA))+
  geom_text(data = letters_df, aes(y = y+1, label = Letters), size = 4) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size=10, angle = 60, hjust=1),
    axis.text.y = element_text(size=10),
    axis.text = element_text(size=10),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))
Figure5A
data_long <- data %>%
  pivot_longer(
    cols = c(Number_BR_W2, Number_BR_W3),
    names_to = "Node2",
    values_to = "BR_count")

attach(data_long)
head(data_long)
lm <- lm(BR_count ~ Genotype*Node2)
anova(lm)
lm_aov=aov(lm) 
resid = residuals(object = lm_aov)
shapiro.test(x=resid)
leveneTest(BR_count ~ Genotype*Node2)
tuk <- TukeyHSD(lm_aov)
tuk_df <- as.data.frame(tuk$Genotype)
tuk_df$comparison <- rownames(tuk_df)
tuk_letters <- multcompLetters4(lm_aov, tuk)
letters_df <- data.frame(
  Genotype = names(tuk_letters$Genotype$Letters),
  Letters = tuk_letters$Genotype$Letters)
letters_df
detach(data_long)

positions <- data_long %>%
  group_by(Genotype) %>%
  summarize(y = max(BR_count, na.rm = TRUE))
letters_df <- left_join(letters_df, positions, by = "Genotype")
data_long = as.data.frame(data_long)
Figure5B = ggplot(data_long, aes(x = Genotype, y = BR_count)) +
  geom_boxplot() +
  labs(x = NULL, y = "BRs per Node (count)") +
  scale_y_continuous(limits = c(0, NA))+
  geom_text(data = letters_df, aes(y = y+1, label = Letters), size = 4) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size=10, angle = 60, hjust=1),
    axis.text.y = element_text(size=10),
    axis.text = element_text(size=10),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))
Figure5B
#pdf(file="/Users/ashley/Desktop/Kamei_Hostetler_et_al_2026/Manuscript/Figure5.pdf", width=6, height=4)
ggdraw() +
  draw_plot(Figure5A, x = 0, y = 0, width = 0.33, height = 0.5) +
  draw_plot(Figure5B, x = 0.33, y = 0, width = 0.33, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C"), 
                  size = 12,
                  x = c(0, 0.33, 0.66), 
                  y = c(0.51,0.51, 0.51))
#dev.off()

##Fig S5#####
cat("\014")
data <- read.csv("Data/Phenotype/2026_d3_field.csv")
colnames(data)
str(data)
i = c(5:7)
data[,i] = apply(data[ ,i], 2,
                 function(x) as.numeric(x))
i = c(1:4,8)
data[,i] = apply(data[ ,i], 2,
                 function(x) as.character(x))
data <- data %>%
  mutate(
    Genotype = dplyr::recode(
      Genotype,
      "MT" = "d3"))

data$BR_Width_mm <- data$BR.Width..px. * (10 / 50) #convert px to mm according to scale
data$Stem_Width_cm <- data$Stem.Width..px. / 50

attach(data)
head(data)
lm <- lm(BR_Width_mm ~ Genotype)
anova(lm)
lm_aov=aov(lm) 
resid = residuals(object = lm_aov)
shapiro.test(x=resid)
leveneTest(BR_Width_mm ~ Genotype)
detach(data)

data$Genotype <- factor(data$Genotype, levels = c("WT","d3"))
PlotS5A = ggplot(data, aes(x = Genotype, y = BR_Width_mm)) +
  geom_boxplot() +
  labs(x = NULL, y = "BR Diameter (mm)") +
  scale_y_continuous(limits = c(0, 22), breaks = seq(0, 22, 2)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size=10),
    axis.text.y = element_text(size=10),
    axis.text = element_text(size=10),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))
PlotS5A

attach(data)
head(data)
lm <- lm(Number..x2. ~ Genotype)
anova(lm)
lm_aov=aov(lm) 
resid = residuals(object = lm_aov)
shapiro.test(x=resid)
leveneTest(Number..x2. ~ Genotype)
detach(data)

PlotS5B = ggplot(data, aes(x = Genotype, y = Number..x2.)) +
  geom_boxplot() +
  labs(x = NULL, y = "BRs per Node (count)") +
  scale_y_continuous(limits = c(0, 18), breaks = seq(0, 18, 2)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size=10),
    axis.text.y = element_text(size=10),
    axis.text = element_text(size=10),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))
PlotS5B
attach(data)
head(data)
lm <- lm(Stem_Width_cm ~ Genotype)
anova(lm)
lm_aov=aov(lm) 
resid = residuals(object = lm_aov)
shapiro.test(x=resid)
leveneTest(Stem_Width_cm ~ Genotype)
detach(data)

PlotS5C = ggplot(data, aes(x = Genotype, y = Stem_Width_cm)) +
  geom_boxplot() +
  labs(x = NULL, y = "Stem Diameter (cm)") +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 1)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size=10),
    axis.text.y = element_text(size=10),
    axis.text = element_text(size=10),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))
PlotS5C

#pdf(file="/Users/ashley/Desktop/Kamei_Hostetler_et_al_2026/Manuscript/FigureS5.pdf", width=6, height=4)
ggdraw() +
  draw_plot(PlotS5A, x = 0, y = 0, width = 0.33, height = 0.5) +
  draw_plot(PlotS5B, x = 0.33, y = 0, width = 0.33, height = 0.5) +
  draw_plot(PlotS5C, x = 0.66, y = 0, width = 0.33, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C"), 
                  size = 12,
                  x = c(0,0.33,0.66), 
                  y = c(0.51,0.51,0.51))
#dev.off()



##Fig S6 ####
cat("\014")
rm(list=ls()) 
data <- read.csv("Data/Phenotype/2023_BRPhenotype.csv")
data$Accession = dplyr::recode( #Rename accessions for clarity/plotting
  data$Accession,
  "ub2ub3tsh4_W22" = "ub2;ub3;tsh4 (W22)",
  "ub2ub3tsh4_B73"      = "ub2ub3tsh4 (B73)")
attach(data)
lm <- lm(Number.of.nodes.with.emerge.BR ~ Accession)
anova(lm)
lm_aov=aov(lm) 
resid = residuals(object = lm_aov)
shapiro.test(x=resid)
leveneTest(Number.of.nodes.with.emerge.BR ~ Accession)
tuk <- TukeyHSD(lm_aov)
tuk_df <- as.data.frame(tuk$Accession)
tuk_df$comparison <- rownames(tuk_df)
tuk_letters <- multcompLetters4(lm_aov, tuk)
letters_df <- data.frame(
  Accession = names(tuk_letters$Accession$Letters),
  Letters = tuk_letters$Accession$Letters)
letters_df
detach(data)

positions <- data %>%
  group_by(Accession) %>%
  summarize(y = max(Number.of.nodes.with.emerge.BR, na.rm = TRUE))
letters_df <- left_join(letters_df, positions, by = "Accession")
data$Accession <- factor(data$Accession, levels = c("ub2ub3tsh4 (B73)", "B73", "ub2;ub3;tsh4 (W22)", "W22"))
#pdf(file="/Users/ashley/Desktop/Kamei_Hostetler_et_al_2026/Manuscript/FigureS6.pdf", width=6, height=4)
ggplot(data, aes(x = Accession, y = Number.of.nodes.with.emerge.BR)) +
  geom_boxplot() +
  labs(x = NULL, y = "Total BR Nodes (count)") +
  scale_y_continuous(limits = c(0, 8),breaks = seq(0, 8, by = 2))+
  geom_text(data = letters_df, aes(y = y + 1, label = Letters), size = 4) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size=10, angle = 60, hjust=1),
    axis.text.y = element_text(size=10),
    axis.text = element_text(size=10),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))
#dev.off()

#A balance between ub3 and tsh4 expression underlies natural variation in brace root size####
##Figure S7####
cat("\014")
rm(list=ls()) 
data = read.csv("Data/Phenotype/2019_RootPhoto.csv", header = TRUE, na.strings = "NA")
ordered_genotypes = with(data, reorder(Genotype, stalk_width_cm, function(x) mean(x, na.rm = TRUE))) #Order genotypes based on hte mean stalk width
genotype_order = levels(ordered_genotypes)
genotype_order = as.list(genotype_order)
data$Genotype_ordered = ordered_genotypes #Add a new column with this information 
summary(data$stalk_width_cm)
attach(data)
lm_x = lm(stalk_width_cm ~ Genotype)
anova(lm_x)
lm_x_aov=aov(lm_x) 
X = HSD.test(lm_x_aov, trt = c("Genotype"), console = TRUE)
Means = X$means
par(mfrow=c(2,2))
plot(lm_x)
par(mfrow=c(2,1))
plot(data$stalk_width_cm)
hist(data$stalk_width_cm)
resid = residuals(object = lm_x_aov)
shapiro.test(x=resid)
#transformation data to meet assumptions
par(mfrow=c(2,2))
data$tukey = transformTukey(
  data$stalk_width_cm,
  start = -10,
  end = 10,
  int = 0.025,
  plotit = TRUE,
  verbose = FALSE,
  quiet = FALSE,
  statistic = 1,
  returnLambda = FALSE
)
hist(data$tukey)
lm_x2 = lm(data$tukey ~ Genotype)
Z = anova(lm_x2)
Z = as.data.frame(Z)
lm_x2_aov2=aov(lm_x2) 
Y = HSD.test(lm_x2_aov2, trt = c("Genotype"), console = TRUE)
Groups = Y$groups
data$tukey = NULL
detach(data)


groups_df = data.frame( #Extract HSD Results 
  Genotype = rownames(X$groups),
  HSD_group = X$groups$groups,
  mean_stalk = X$groups[,1],
  stringsAsFactors = FALSE)
data2 <- left_join( #Merge HSD results with data
  data,
  dplyr::select(groups_df, Genotype, HSD_group, mean_stalk),
  by = "Genotype")
order_levels <- groups_df %>% #order genotypes based on the mean stalk width
  arrange(desc(mean_stalk)) %>%
  pull(Genotype)
data2$Genotype = factor(data2$Genotype, levels = order_levels)
head(data2)
summary_df = data2 %>% #Calculate the means talk width, root width, and number of roots for each geno
  group_by(Genotype, HSD_group) %>%
  summarise(
    mean_stalk   = mean(stalk_width_cm, na.rm = TRUE),
    mean_root_w  = mean(single_root_width_cm, na.rm = TRUE),
    mean_w1      = mean(w1, na.rm = TRUE),
    n = n()
  ) %>%
  ungroup()
summary_df
#pdf(file="/Users/ashley/Desktop/Kamei_Hostetler_et_al_2026/Manuscript/FigS7.pdf", width=8, height=6)
ggplot(summary_df, aes(x = mean_w1, y = mean_root_w, color = HSD_group, label = Genotype)) +
  geom_point(size = 3) +
  geom_text(hjust = 1.2, vjust = 0.5, size = 3, check_overlap = FALSE, show.legend = FALSE) +
  scale_color_manual(values = c("gold","steelblue1","steelblue4","pink2","darkmagenta"))+
  expand_limits(y = 0, x=0)+
  labs(
    x = "Number of Roots (Node 1)",
    y = "Single Root Width (cm)",
    color = "HSD group (stalk width)")+
  theme_bw() +
  theme(legend.position = "top",
    axis.text.x = element_text(size=10,),
    axis.text.y = element_text(size=10),
    axis.text = element_text(size=10),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))

#dev.off()

##Fig 6B ####
df2 <- data %>% #Select only CML333 and PHHB4 into a new data frame
  filter(Genotype %in% c("CML333", "PHHB4"))
rm(list = setdiff(ls(), c("df2")))
df2 = df2[,c(2,6,8)]
df2 <- df2 %>%
  mutate(RootWidth_mm = single_root_width_cm * 10) #Convert cm to mm
df2 = df2[,c(1,2,4)]

data = read.csv("Data/Phenotype/2023_FieldDataBRPhenos.csv", header = TRUE, na.strings = "NA")
data <- data %>%
  mutate(
    Genotype = dplyr::recode( #Rename ub2,ub3 to ub2;ub3 
      Genotype,
      "ub2,3" = "ub2;ub3"))
data <- data %>%
  mutate(
    Average_BR_Diameter = rowMeans(across(BR_diameter_1_mm:BR_diameter_4_mm), na.rm = TRUE))
data$Genotype <- factor(data$Genotype, levels = c("tsh4", "ub2;ub3"))

data1 = data[,c(7,16,13)]
colnames(data1)[2] = "RootWidth_W2_mm"
colnames(data1)[3] = "BraceRootNumber_W2"

unique(df2$Genotype)
head(df2)
df2$Genotype <- factor(df2$Genotype, levels = c("PHHB4", "CML333"))
Fig6B = ggplot(df2, aes(x = w1, y = RootWidth_mm, shape = Genotype, color = Genotype)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c("PHHB4" = 16, "CML333" = 17))+
  scale_color_manual(values = c("gray1","gray70"))+
  expand_limits(y = 0, x=0)+
  labs(
    x = "Number of Brace Roots (count)",
    y = "Brace Root Width (mm)",
    color = "Genotype")+
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        #axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.text = element_text(size=10),
        axis.title = element_text(face="bold", size=10),
        #axis.title.x = element_text(face="bold", size=10),
        axis.title.x= element_blank(),
        axis.title.y = element_text(face="bold", size=10))
Fig6B
##Fig 6A####
data1$Genotype <- factor(data1$Genotype, levels = c("tsh4", "ub2;ub3"))
Fig6A = ggplot(data1, aes(x = BraceRootNumber_W2, y = RootWidth_W2_mm, shape = Genotype, color = Genotype)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c("tsh4" = 15, "ub2;ub3" = 8))+
  scale_color_manual(values = c("gray1","gray70"))+
  expand_limits(y = 0, x=0)+
  labs(
    x = "Number of Brace Roots (count)",
    y = "Brace Root Width (mm)",
    color = "Genotype")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        #axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.text = element_text(size=10),
        axis.title = element_text(face="bold", size=10),
        #axis.title.x = element_text(face="bold", size=10),
        axis.title.x= element_blank(),
        axis.title.y = element_text(face="bold", size=10))

Fig6A

##Fig 6C ####
cat("\014")
rm(list = setdiff(ls(), c("Fig6A","Fig6B")))
dat <- read.csv("Data/qPCR_SBP_CML333&PHHB4.csv")
unique(dat$Sample.Name)

df <- dat %>% #Average technical replicates
  group_by(Sample.Name, Target.Name) %>%
  summarise(
    Ct = mean(CT, na.rm = TRUE),
    Ct_sd = sd(CT, na.rm = TRUE),
    .groups = "drop")
df1 <- df %>% #Convert to wide format
  dplyr::select(Sample.Name, Target.Name, Ct) %>%
  pivot_wider(
    names_from = Target.Name,
    values_from = Ct
  )
df1 <- df1 %>% 
  mutate(
    dCt_TSH4 = TSH4 - ACT, #Calculates dCT (when did target gene reach detection threshold relative to housekeeping gene)
    dCt_UB3  = UB3  - ACT, #2^X = Y where X is the cycle at which target gene reaches threshold and Y is Yfold less abundant than HK gene
    
    Q_UB3  = 2^(-dCt_UB3),#Converts to normalized expression (Q = quantity of UB3 ; expression) 
    Q_TSH4 = 2^(-dCt_TSH4), #negative sign because higher Ct means less RNA; without negative sign, larger Ct would give larger expression values 
                          
    
    UB3_TSH4_ratio = Q_UB3 / Q_TSH4) #Calculate UB3 relative to TSH4 within each sample; both genes are already normalized by ACT
                                      #How much UB3 transcript is present for every unit of TSH4 transcript?
                                      #UB3/TSH4 = 0.20 means UB3 is expressed at 20% of the level of TSH4 
df1 <- df1 %>%
  mutate(
    Genotype = case_when( #Assign genotype information based on what it was called when qPCR was set up
      grepl("^C", Sample.Name) ~ "CML333",
      grepl("^P", Sample.Name) ~ "PHHB4",
      TRUE ~ NA_character_
    )
  )
summary_table <- df1 %>%
  group_by(Genotype) %>%
  summarise(
    mean_Q_TSH4 = mean(Q_TSH4, na.rm = TRUE),
    sd_Q_TSH4   = sd(Q_TSH4, na.rm = TRUE),
    
    mean_Q_UB3 = mean(Q_UB3, na.rm = TRUE),
    sd_Q_UB3   = sd(Q_UB3, na.rm = TRUE),
    
    mean_ratio = mean(UB3_TSH4_ratio, na.rm = TRUE),
    sd_ratio   = sd(UB3_TSH4_ratio, na.rm = TRUE)
  )

(1.20 - 0.814)/0.814

Fig6C = ggplot(df1, aes(x = Genotype, y = UB3_TSH4_ratio)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  annotate( # TSH4-enriched region
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = 0,
    ymax = 1,
    alpha = 0.05) +
  annotate(   # UB3-enriched region
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = 1,
    ymax = 2,
    alpha = 0.3) +
  labs(
    x = "Number of Brace Roots (count)",
    y = "UB3/TSH4 Expression Ratio")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        #axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.text = element_text(size=10),
        axis.title = element_text(face="bold", size=10),
        #axis.title.x = element_text(face="bold", size=10),
        axis.title.x= element_blank(),
        axis.title.y = element_text(face="bold", size=10))
#pdf(file="/Users/ashley/Desktop/Kamei_Hostetler_et_al_2026/Manuscript/Figure6.pdf", width=6, height=4)
ggdraw() +
  draw_plot(Fig6A, x = 0, y = 0, width = 0.33, height = 0.5) +
  draw_plot(Fig6B, x = 0.33, y = 0, width = 0.33, height = 0.5) +
  draw_plot(Fig6C, x = 0.66, y = 0, width = 0.33, height = 0.5)+
  draw_plot_label(label = c("A", "B", "C"), 
                  size = 12,
                  x = c(0,0.33,0.66), 
                  y = c(0.5,0.5,0.5))
#dev.off()

#Table S2####
cat("\014")
rm(list=ls()) 
data = read.csv("Data/DE_results/B73_V6/P_vs_NP_results.csv", header = TRUE, na.strings = "NA")
data = data[,c(7,8,9,10,20:28)]
df = data
rm(list = setdiff(ls(), c("df")))
data = read.csv("Data/DE_results/MutantResults/ub2ub3_vs_W22_S1_all_genes.csv", header = TRUE, na.strings = "NA")
data = data[,c(8,14:30,32,34,36)]
df1 = data
rm(list = setdiff(ls(), c("df","df1")))
colnames(df1)[1] = "Geneid"
sum(duplicated(df1$Geneid))
df1[df1$Geneid %in% df1$Geneid[duplicated(df1$Geneid)], ]
df1 <- df1[!is.na(df1$Geneid) & !duplicated(df1$Geneid), ]
data <- merge(df, df1, by = "Geneid", all = TRUE)
#write.csv(data, "~/Desktop/Kamei_Hostetler_et_al_2026/Manuscript/TableS2.csv", row.names = FALSE)
