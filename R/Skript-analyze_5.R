require(openxlsx, quietly = TRUE)
require(lme4, quietly = TRUE)
require(pheatmap, quietly = TRUE)

ARGS <- commandArgs(trailingOnly = TRUE)
isExpr <- which(ARGS == "-e") + 1
isAnno <- which(ARGS == "-a") + 1
isGroup <- which(ARGS == "-g") + 1
isCov <- which(ARGS == "-c") + 1
isQuan <- which(ARGS == "-q") + 1

## import
cat("Step 01: Importing data...\n")
MAT <- read.xlsx(ARGS[isExpr], rowNames = TRUE)
ANNO <- read.xlsx(ARGS[isAnno])
GROUP <- read.xlsx(ARGS[isGroup])
if (length(isCov) != 0) COVAR <- read.xlsx(ARGS[isCov])
if (length(isQuan) != 0) QUAN <- as.numeric(ARGS[isQuan]) else QUAN <- 0.1
NC <- ncol(MAT)

## 0 => 0.01 for avoiding log-errors
MAT[MAT == 0] <- 0.01; MAT[is.na(MAT)] <- 0.01
MAT <- log2(MAT)

## removing background transcripts < 0.1 quantile in at least 50% of samples
cat("Step 02: Removing transcripts near backgound (at least 50% with mean < 0.1 quantile) ...\n")
MEAN <- rowMeans(MAT, na.rm = TRUE)
QUANTILE <- quantile(MEAN, QUAN)
SUM <- apply(MAT, 1, function(x) sum(x < QUANTILE, na.rm = TRUE))
sel <- which(SUM > 0.5 * NC)
MAT <- MAT[-sel, ]
cat(" => Removed", length(sel), "transcripts\n")

## check ID's
if (!all.equal(colnames(MAT), GROUP[, 1])) stop("Sample names in GROUP and EXPRESSION must be identical!")

## set group colors, transparent versions of "cadetblue" and "coral3"
COL <- c("#5F9EA0AA", "#CD5B45AA")[GROUP[, 2] + 1]

## initialize pdf
pdf(file = "Figures.pdf", paper = "a4r", width = 12, height = 10)

## Figure 1: Raw expression boxplot
cat("Step 03: Creating boxplot of Read Counts...\n")
par(mar = c(6, 5, 3, 1))
boxplot(MAT, outline = FALSE, las = 2, ylab = "log2(Counts)", col = COL, border = COL, cex.lab = 1.5,
        main = "Figure 1: Boxplot of log2-transformed Read Counts", cex.main = 1.5, boxwex = 0.5)

## Figure 2: PCA based on all transcripts
cat("Step 04: Creating PCA based on all transcripts...\n")
PCA <- prcomp(t(MAT))
PCs <- PCA$x[, 1:4]
pairs(PCs, col = COL, pch = 16, cex = 2, cex.labels = 2, upper.panel = NULL, 
      main = "Figure 2: PCA based on all Reads, 4 components", cex.main = 1.5)
OUT1 <- cbind(Sample = colnames(MAT), PCs)

## Figure 3: PVCA analysis
cat("Step 05: Calculating variance contribution of the covariates...\n")
if (exists("COVAR")) {
  source("PVCA.R")
  PVCA <- pvca(MAT, COVAR) * 100
  barplot(PVCA, ylim = c(0, 100), ylab = "Contribution [%]", cex.lab = 1.5, 
          main = "Figure 3: Variance contribution of the covariates", cex.main = 1.5, cex.names = 1.5)
  OUT2 <- t(PVCA)
} else OUT2 <- NA

## filter out top 200 variable transcripts
cat("Step 06: Filtering top 2000 variable transcripts...\n")
VAR <- apply(MAT, 1, function(x) var(x, na.rm = TRUE))
MAT2 <- MAT; MAT2 <- cbind(MAT, VAR); MAT2 <- MAT2[order(MAT2$VAR, decreasing = TRUE), ]
MAT3 <- MAT2[1:2000, 1:NC]

## Figure 4: PCA based on top 2000 variable transcripts
cat("Step 07: Creating PCA based on top 2000 variable transcripts...\n")
PCA <- prcomp(t(MAT3))
PCs <- PCA$x[, 1:4]
pairs(PCs, col = COL, pch = 16, cex = 2, cex.labels = 2, upper.panel = NULL, 
      main = "Figure 4: PCA based on top 2000 variable transcripts, 4 components", cex.main = 1.5)
OUT3 <- cbind(Sample = colnames(MAT), PCs)

## Figure 5: Heatmap based on top 200 variable transcripts
cat("Step 08: Creating Heatmap based on top 200 variable transcripts...\n")
aC <- GROUP[, 2, drop = FALSE]; rownames(aC) <- GROUP[, 1]
cC <- list(Group = unique(COL))
pheatmap(MAT3[1:200, ], show_rownames = FALSE, clustering_distance_rows = "euclidean", clustering_distance_cols = "correlation",
         annotation_col = aC, annotation_colors = cC, main = "Figure 5: Heatmap based on top 200 variable transcripts")
Sys.sleep(10)

## Figure 6: Profile plot of top 10 variable transcripts
cat("Step 09: Creating profile plot of top 10 variable transcripts...\n")
par(mfrow = c(10, 1)); par(mar = c(0.5, 4, 0, 0))
for (i in 1:10) {
  plot(as.numeric(MAT3[i, ]), col = 1, pch = 16, type = "l", xaxt = if (i == 10) NULL else "n", ylab = "log2(Counts)")
  points(as.numeric(MAT3[i, ]), col = COL, pch = 16, xaxt = "n", cex = 2)
}
mtext("Figure 6: Profile plot of top 10 variable transcripts", outer = TRUE,  cex = 1.5, line = -2)

## filter using linear model with covariates
cat("Step 10: Calculating (adjusted) linear model...\n")
if (exists("COVAR")) DAT <- data.frame(y = as.numeric(MAT[1, ]), group = GROUP[, 2], COVAR = COVAR) else DAT <- data.frame(y = as.numeric(MAT[1, ]), group = GROUP[, 2])
statMat <- matrix(NA_real_, nrow = nrow(MAT), ncol = 4)
colnames(statMat) <- c("log(Estimate)", "SE", "t", "P")
for (i in 1:nrow(MAT)) {
  DAT[, 1] <- as.numeric(MAT[i, ])
  LM <- lm(y ~ ., data = DAT)
  statMat[i, ] <- summary(LM)$coefficients[2, ]
}
MEAN <- rowMeans(MAT, na.rm = TRUE)

## Figure 7: MA plot
par(mfrow = c(1, 1));  par(mar = c(5.1, 4.1, 4.1, 2.1))
cat("Step 11: Creating MA plot...\n")
plot(MEAN, statMat[, 1], col = ifelse(statMat[, 1] > 0.585 | statMat[, 1] < -0.585, "#8B000088", "#00000088"), pch = 16,
     xlab = "Mean log(Counts)", ylab = "log(Ratio)", main = "Figure 7: MA plot", cex.main = 1.5)

## Figure 8: Volcano plot
Padj <- p.adjust(statMat[, 4], "fdr")
cat("Step 12: Creating Volcano plot...\n")
plot(statMat[, 1], -log10(Padj), col = ifelse((statMat[, 1] > 0.585 | statMat[, 1] < -0.585) & -log10(Padj) > 1.3, "#8B000088", "#00000088"),
     pch = 16, xlab = "log(Ratio)", ylab = "-log10(P adjusted)", main = "Figure 8: Volcano plot", cex.main = 1.5)
abline(v = c(-0.585, 0.585), lwd = 2, lty = 2); abline(h = 1.3, lwd = 2, lty = 2)

cat("Step 13: Exporting complete result matrix sorted by adjusted P-value...\n")
m <- match(rownames(MAT), ANNO[, 1])
OUT4 <- cbind(ANNO[m, ], MAT, statMat, Var = VAR,  Padj, Ratio = 2^statMat[, 1])
OUT4 <- OUT4[order(OUT4$Padj), ]

dev.off()

## export all data
write.xlsx(list(PCA_all = OUT1, PVCA = OUT2, PCA_2000 = OUT3, RESULTS = OUT4), "Results.xlsx")



