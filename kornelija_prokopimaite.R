####################################
### Installing required packages ### 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")  
BiocManager::install("minfi")
BiocManager::install("IlluminaHumanMethylation450kmanifest")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

install.packages("factoextra")
install.packages("qqman")
install.packages("gplots")
install.packages("future.apply")
install.packages("viridis")
install.packages("pheatmap") 


### 1. Load raw data
rm(list = ls())
library(minfi)
setwd('/Users/kornelijaprok/Desktop/DNA_RNA/fluorescence_info')
baseDir <- getwd()
sampleSheetFile <- file.path(baseDir, "Samplesheet_report_2024.csv")
targets <- read.metharray.sheet(baseDir, pattern = basename(sampleSheetFile))
rownames(targets) <- basename(targets$Basename)
RGset <- read.metharray.exp(targets = targets)
save(RGset, file = "RGset.RData")

# Load and clean Illumina 450k Manifest
manifest <- read.csv("/Users/kornelijaprok/Desktop/DNA_RNA/humanmethylation450_15017482_v1-2.csv", 
                     header = TRUE, skip = 7)
Illumina450Manifest_clean <- manifest[, c("AddressA_ID", "Infinium_Design_Type", "IlmnID")]
save(Illumina450Manifest_clean, file = "Illumina450Manifest_clean.RData")


### 2. Create R/G dataframes
Red <- data.frame(getRed(RGset))
head(Red)
Green <- data.frame(getGreen(RGset))
head(Green)


### 3. Get probe info 
load("Illumina450Manifest_clean.RData")
probe_address <- "42639338"

probe_red <- Red[rownames(Red) == probe_address, ]
probe_green <- Green[rownames(Green) == probe_address, ]

type_val <- Illumina450Manifest_clean$Infinium_Design_Type[
  Illumina450Manifest_clean$AddressA_ID == probe_address]

table_address <- data.frame(
  Sample = colnames(probe_red),
  Red_fluor = unlist(probe_red, use.names = FALSE),
  Green_fluor = unlist(probe_green, use.names = FALSE),
  Type = type_val,
  stringsAsFactors = FALSE
)

print(table_address)


### 4. Create the object MSet.raw 

MSet.raw <- preprocessRaw(RGset)
MSet.raw


### 5. Quality check

qc <- getQC(MSet.raw)
dev.new()      
plotQC(qc)

dev.new()
controlStripPlot(RGset, controls="NEGATIVE")

sample_names <- colnames(RGset)

detP <- detectionP(RGset)
failed <- detP > 0.01
num_failed <- colSums(failed)

df_failed <- data.frame(
  Sample = sample_names,
  Failed_Positions = num_failed,
  stringsAsFactors = FALSE
)

rownames(df_failed) <- NULL
print(df_failed)


### 6. Beta and M values
# Subsetting MSet.raw for WT and MUT groups
wt <- basename(targets[targets$Group=="WT", "Basename"])
mut <- basename(targets[targets$Group=="MUT", "Basename"])
wt_set <- MSet.raw[,colnames(MSet.raw) %in% wt]
mut_set <- MSet.raw[,colnames(MSet.raw) %in% mut]

wtBeta <- getBeta(wt_set)
wtM <- getM(wt_set)
mutBeta <- getBeta(mut_set)
mutM <- getM(mut_set)

mean_wtBeta <- apply(wtBeta, MARGIN=1, mean, na.rm=TRUE)
mean_mutBeta <- apply(mutBeta, MARGIN=1, mean, na.rm=TRUE)
mean_wtM <- apply(wtM, MARGIN=1, mean, na.rm=TRUE)
mean_mutM <- apply(mutM, MARGIN=1, mean, na.rm=TRUE)

d_mean_wtBeta <- density(mean_wtBeta, na.rm=TRUE)
d_mean_mutBeta <- density(mean_mutBeta, na.rm=TRUE)
d_mean_wtM <- density(mean_wtM, na.rm=TRUE)
d_mean_mutM <- density(mean_mutM, na.rm=TRUE)
par(mfrow=c(1,2)) 

# Plotting density of mean Beta values
dev.new()
plot(d_mean_wtBeta, main="Density of Beta Values", col="#1E90FF", lwd=2.5, xlab="Beta Values", ylab="Density") # Blue for WT
lines(d_mean_mutBeta, col="#FF4500", lwd=2.5) # Orange-Red for MUT
legend('topright', legend=c("WT", "MUT"), fill=c("#1E90FF", "#FF4500"), cex=1)

# Plotting density of mean M values
de
plot(d_mean_wtM, main="Density of M Values", col="#1E90FF", lwd=2.5, xlab="M Values", ylab="Density") # Blue for WT
lines(d_mean_mutM, col="#FF4500", lwd=2.5) # Orange-Red for MUT
legend('topright', legend=c("WT", "MUT"), fill=c("#1E90FF", "#FF4500"), cex=1)



### 7. Normalisation

dfI <- subset(Illumina450Manifest_clean, Infinium_Design_Type == "I")
dfII <- subset(Illumina450Manifest_clean, Infinium_Design_Type == "II")

beta <- getBeta(MSet.raw)
beta_I <- beta[rownames(beta) %in% dfI$IlmnID, ]
beta_II <- beta[rownames(beta) %in% dfII$IlmnID, ]

#Calculating raw means and standard deviations for each probe type:
mean_beta_I <- apply(beta_I, 1, mean, na.rm=TRUE)
mean_beta_II <- apply(beta_II, 1, mean, na.rm=TRUE)

sd_beta_I <- apply(beta_I, 1, sd, na.rm=TRUE)
sd_beta_II <- apply(beta_II, 1, sd, na.rm=TRUE)

#Calculating density estimates for means and SDs:
d_mean_beta_I <- density(mean_beta_I, na.rm=TRUE)
d_mean_beta_II <- density(mean_beta_II, na.rm=TRUE)

d_sd_beta_I <- density(sd_beta_I, na.rm=TRUE)
d_sd_beta_II <- density(sd_beta_II, na.rm=TRUE)

#Normalising the data:
MSet.norm <- preprocessNoob(RGset)

#Extracting normalised beta values:
beta_norm <- getBeta(MSet.norm)
beta_norm_I <- beta_norm[rownames(beta_norm) %in% dfI$IlmnID, ]
beta_norm_II <- beta_norm[rownames(beta_norm) %in% dfII$IlmnID, ]

#Calculating mean and SD densities for normalised data:
mean_beta_norm_I <- apply(beta_norm_I, 1, mean, na.rm=TRUE)
mean_beta_norm_II <- apply(beta_norm_II, 1, mean, na.rm=TRUE)

sd_beta_norm_I <- apply(beta_norm_I, 1, sd, na.rm=TRUE)
sd_beta_norm_II <- apply(beta_norm_II, 1, sd, na.rm=TRUE)

d_mean_beta_norm_I <- density(mean_beta_norm_I, na.rm=TRUE)
d_mean_beta_norm_II <- density(mean_beta_norm_II, na.rm=TRUE)

d_sd_beta_norm_I <- density(sd_beta_norm_I, na.rm=TRUE)
d_sd_beta_norm_II <- density(sd_beta_norm_II, na.rm=TRUE)

#Plotting the comparison:
dev.new()
par(mfrow=c(2,3))

group_colors <- c("WT" = "#1E90FF", "MUT" = "#FF4500")
boxplot_colors <- group_colors[as.character(targets$Group)]

plot(d_mean_beta_I, col="#1E90FF", main="Raw Beta Mean (Type I)", xlim=c(0,1), ylim=c(0,5)) # Blue for Type I
lines(d_mean_beta_II, col="#FF4500") # Orange-Red for Type II
legend("topright", legend=c("Type I","Type II"), col=c("#1E90FF","#FF4500"), lty=1, cex=0.8)

plot(d_sd_beta_I, col="#1E90FF", main="Raw Beta SD (Type I)", xlim=c(0,0.6), ylim=c(0,60)) # Blue for Type I
lines(d_sd_beta_II, col="#FF4500") # Orange-Red for Type II
legend("topright", legend=c("Type I","Type II"), col=c("#1E90FF","#FF4500"), lty=1, cex=0.8)

boxplot(beta, col=boxplot_colors, main="Raw Beta Values")

plot(d_mean_beta_norm_I, col="#1E90FF", main="Normalized Beta Mean (Type I)", xlim=c(0,1), ylim=c(0,5)) # Blue for Type I
lines(d_mean_beta_norm_II, col="#FF4500") # Orange-Red for Type II
legend("topright", legend=c("Type I","Type II"), col=c("#1E90FF","#FF4500"), lty=1, cex=0.8)

plot(d_sd_beta_norm_I, col="#1E90FF", main="Normalized Beta SD (Type I)", xlim=c(0,0.6), ylim=c(0,60)) # Blue for Type I
lines(d_sd_beta_norm_II, col="#FF4500") # Orange-Red for Type II
legend("topright", legend=c("Type I","Type II"), col=c("#1E90FF","#FF4500"), lty=1, cex=0.8)

boxplot(beta_norm, col=boxplot_colors, main="Normalized Beta Values")

par(mfrow=c(1,1)) 


### 8. Perform PCA

pca_results <- prcomp(t(beta_norm), scale = TRUE) 

dev.new()
plot(pca_results)

palette(c("#FAC42A", "#AB2494"))

# Plotting PCA by Group
dev.new()
plot(pca_results$x[, 1], pca_results$x[, 2], cex = 1, pch = 19, col = targets$Group,
     xlab = "PC1", ylab = "PC2", main = "PCA - Groups")
text(pca_results$x[, 1], pca_results$x[, 2], labels = rownames(pca_results$x), cex = 0.4, pos = 3)
legend("bottomright", legend = levels(targets$Group), col = 1:length(levels(targets$Group)), pch = 19)

# Plotting PCA by Sex
dev.new()
palette(c("#C364CA", "#8BC4F9"))
plot(pca_results$x[, 1], pca_results$x[, 2], cex = 1, pch = 19, col = targets$Sex,
     xlab = "PC1", ylab = "PC2", main = "PCA - Sex")
text(pca_results$x[, 1], pca_results$x[, 2], labels = rownames(pca_results$x), cex = 0.4, pos = 3)
legend("bottomright", legend = levels(targets$Sex), col = 1:length(levels(targets$Sex)), pch = 19)

# Plotting PCA by Batch
dev.new()
palette(c("#9e0059", "yellow", "#2c7bb6", "#d7191c", "#fdae61"))  # More colors if needed
plot(pca_results$x[, 1], pca_results$x[, 2], cex = 1, pch = 19, col = targets$Batch,
     xlab = "PC1", ylab = "PC2", main = "PCA - Batch")
text(pca_results$x[, 1], pca_results$x[, 2], labels = rownames(pca_results$x), cex = 0.4, pos = 3)
legend("bottomright", legend = levels(targets$Batch), col = 1:length(levels(targets$Batch)), pch = 19)


### 9. T-test

targets$Group <- as.factor(targets$Group)

# Subsetting beta values for WT and MUT groups
wt_samples <- targets$Basename[targets$Group == "WT"]
mut_samples <- targets$Basename[targets$Group == "MUT"]

wt_beta <- beta_norm[, colnames(beta_norm) %in% wt_samples]
mut_beta <- beta_norm[, colnames(beta_norm) %in% mut_samples]

# Performing t-test for each probe
t_test_results <- apply(beta_norm, 1, function(probe_values) {
  t.test(probe_values[targets$Group == "WT"],
         probe_values[targets$Group == "MUT"])$p.value
})

adjusted_p_values <- p.adjust(t_test_results, method = "BH")

hist(t_test_results,
     breaks = 50,
     col = "lightblue",
     main = "Histogram of P-Values",
     xlab = "P-Value",
     ylab = "Frequency")
abline(v = 0.05, col = "red", lwd = 2, lty = 2)

differential_methylation_results <- data.frame(
  Probe = rownames(beta_norm),
  P_Value = t_test_results,
  Adjusted_P_Value = adjusted_p_values,
  stringsAsFactors = FALSE
)

### 10. Multiple test correction

bonferroni_p_values <- p.adjust(t_test_results, method = "bonferroni")
bh_p_values <- p.adjust(t_test_results, method = "BH")

nominal_significant_probes <- sum(t_test_results < 0.05)
bonferroni_significant_probes <- sum(bonferroni_p_values < 0.05)
bh_significant_probes <- sum(bh_p_values < 0.05)

cat("Nominal Significant Probes:", nominal_significant_probes, "\n")
cat("Bonferroni Significant Probes:", bonferroni_significant_probes, "\n")
cat("BH Significant Probes:", bh_significant_probes, "\n")

significance_counts <- c(nominal_significant_probes, bonferroni_significant_probes, bh_significant_probes)
labels <- c("Nominal", "Bonferroni", "BH-adjusted")

barplot(
  height = significance_counts,
  names.arg = labels,
  col = c("#4CAF50", "#2196F3", "#FFC107"),
  main = "Significant Probes at Different Correction Thresholds",
  ylab = "Number of Significant Probes",
  xlab = "Correction Method",
  las = 1
)



### 11. Volcano and Manhattan plots

targets$Basename <- basename(targets$Basename)

wt_samples <- targets$Basename[targets$Group == "WT"]
mut_samples <- targets$Basename[targets$Group == "MUT"]

wt_beta <- beta_norm[, colnames(beta_norm) %in% wt_samples, drop=FALSE]
mut_beta <- beta_norm[, colnames(beta_norm) %in% mut_samples, drop=FALSE]

mean_difference <- rowMeans(wt_beta, na.rm=TRUE) - rowMeans(mut_beta, na.rm=TRUE)
differential_methylation_results$Mean_Difference <- mean_difference

# Volcano plot
plot(
  x = differential_methylation_results$Mean_Difference,
  y = -log10(differential_methylation_results$P_Value),
  pch = 20,
  col = ifelse(differential_methylation_results$Adjusted_P_Value < 0.05, "red", "gray"),
  xlab = "Mean Difference (WT - MUT)",
  ylab = "-log10(P-Value)",
  main = "Volcano Plot"
)
abline(h = -log10(0.05), col = "blue", lty = 2)

#Manhattan plot
library(qqman)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annotation <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

annotated_results <- merge(
  differential_methylation_results,
  annotation[, c("Name", "chr", "pos")],
  by.x = "Probe",
  by.y = "Name",
  all.x = TRUE
)

annotated_results$CHR <- gsub("^chr", "", annotated_results$chr)
annotated_results$CHR[annotated_results$CHR == "X"] <- "23"
annotated_results$CHR[annotated_results$CHR == "Y"] <- "24"
annotated_results$CHR[annotated_results$CHR %in% c("MT", "M")] <- "25"
annotated_results$CHR <- as.numeric(annotated_results$CHR)

# Filtering rows with valid CHR, pos, and P-value
valid <- !is.na(annotated_results$CHR) & !is.na(annotated_results$pos) & !is.na(annotated_results$P_Value) &
  annotated_results$P_Value > 0 & annotated_results$P_Value <= 1

annotated_results_clean <- annotated_results[valid, ]

manhattan_data <- data.frame(
  SNP = annotated_results_clean$Probe,
  CHR = annotated_results_clean$CHR,
  BP = annotated_results_clean$pos,
  P = annotated_results_clean$P_Value
)

# Removing duplicated SNPs
manhattan_data <- manhattan_data[!duplicated(manhattan_data$SNP), ]

dev.new()
if(nrow(manhattan_data) > 0) {
  manhattan(manhattan_data, main = "Manhattan Plot of Differential Methylation")
} else {
  cat("No valid data available to plot Manhattan plot.\n")
}

### 12. Heatmap

library(pheatmap)

top_probes <- head(differential_methylation_results[order(differential_methylation_results$P_Value), ], 100)$Probe
heatmap_data <- beta_norm[top_probes, , drop = FALSE]
heatmap_data_scaled <- t(scale(t(heatmap_data)))
annotation_col <- data.frame(Group = targets$Group)
rownames(annotation_col) <- targets$Basename

dev.new()  
pheatmap(
  heatmap_data_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  annotation_col = annotation_col,
  main = "Heatmap of Top 100 Differentially Methylated Probes"
)
