library(GEOquery)
library(limma)
library(ggplot2)

source("../R/common.R")

gse <- getGEO("GSE62061")[[1]];

x <- exprs(gse);
hist(x, breaks=100);

min(x)
x <- x - min(x);

x <- log2(x + 1);
hist(x, breaks=100);
min(x)
max(x)

expr.clip <- quantile(x, c(0.005));
x <- x - expr.clip;
x[x < 0] <- 0;
hist(x, breaks=100);

x.f <- filter_undetected(x);

dim(x)
dim(x.f)

# prepare phenotype data

pheno <- data.frame(
	cell_line = rep(NA, ncol(gse)),
	group = rep(NA, ncol(gse))
);
description <- pData(phenoData(gse))[, c("cell line:ch1")];
pheno$group[grep("sensitive", description)] <- "erlotinib-sensitive";
pheno$group[grep("resistant", description)] <- "erlotinib-resistant";
pheno$cell_line[grep("SCC-25", description)] <- "SCC-25";
pheno$cell_line[grep("Cal-27", description)] <- "Cal-27";
pheno$cell_line[grep("FaDu", description)] <- "FaDu";
pheno$cell_line[grep("SQ20B", description)] <- "SQ20B";

pheno$group <- factor(pheno$group,
	levels = c("erlotinib-sensitive", "erlotinib-resistant")
);


# differential gene expression analysis

design <- model.matrix(~ group + cell_line, data = pheno);

fit <- lmFit(x, design);
fit <- eBayes(fit);
res <- topTable(fit, coef="grouperlotinib-resistant", number=Inf);
head(res)

idx <- match(rownames(res), pData(featureData(gse))$ID);
res$gene <- pData(featureData(gse))[idx, "Symbol"];

# target false discovery rate
fdr <- 0.01;
# number of signficant results
sum(res$adj.P.Val < fdr);
# expected number of false positives
sum(res$adj.P.Val < fdr) * fdr;

genes <- "AXL";

res[res$gene %in% genes, ]

pd <- pData(featureData(gse));
probes <- pd$ID[pd$Symbol %in% genes];
x.s <- t(x.f[rownames(x.f) %in% probes, , drop=FALSE]);
d <- data.frame(
	pheno,
	x.s
);

ggplot(d, aes(x = group, y = ILMN_1701877)) +
	theme_bw() +
	geom_jitter(width=0.1)

ggplot(d, aes(x = group, y = ILMN_1701877, colour=cell_line)) +
	theme_bw() +
	geom_jitter(width=0.1)


# ---

# adding group:cell_line interaction term allow us to find
# gene expression differences that are different across cell lines
design <- model.matrix(~ group + cell_line + group:cell_line, data = pheno);

fit <- lmFit(x, design);
fit <- eBayes(fit);
res <- topTable(fit, coef="grouperlotinib-resistant", number=Inf);
head(res)

colnames(coef(fit))

topTable(fit, coef="grouperlotinib-resistant:cell_lineSCC-25", number=Inf)["ILMN_2188862", ]

idx <- match(rownames(res), pData(featureData(gse))$ID);
res$gene <- pData(featureData(gse))[idx, "Symbol"];

res[res$adj.P.Val < fdr, ]

genes <- "GDF15";

res[res$gene %in% genes, ]

pd <- pData(featureData(gse));
probes <- pd$ID[pd$Symbol %in% genes];
x.s <- t(x.f[rownames(x.f) %in% probes, , drop=FALSE]);
d <- data.frame(
	pheno,
	x.s
);

# in one cell line, GDF15 expression is reduced in the resistance group
ggplot(d, aes(x = group, y = ILMN_2188862, colour=cell_line)) +
	theme_bw() +
	geom_jitter(width=0.1)

