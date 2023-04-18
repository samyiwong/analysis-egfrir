library(GEOquery)
library(limma)
library(ggplot2)
library(io)

# platform: GPL10558
#           Illumina HumanHT-12 V4.0 expression beadchip
library(illuminaHumanv4.db)
db <- illuminaHumanv4.db;

source("../R/common.R")

accession <- "GSE62061";
out.fn <- filename(accession);

gse <- getGEO(accession)[[1]];

x <- exprs(gse);
hist(x, breaks=100);

min(x)
# shift towards 0 because x contains negative values
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
	row.names = colnames(gse),
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


# differential gene expression analysis for each cell_line

pheno.g <- split(pheno, pheno$cell_line);
x.g <- lapply(pheno.g, function(pheno) x[, rownames(pheno)]);

res.g <- mapply(
	function(x, pheno) {
		diff_expr_erlotinib_resistant(x, pheno, db)
	},
	x.g, pheno.g,
	SIMPLIFY=FALSE
);

probes <- rownames(x.f);
features <- data.frame(
	row.names = probes,
	ensembl = mapIds(db, keys=probes, keytype="PROBEID", column="ENSEMBL"),
	gene = mapIds(db, keys=probes, keytype="PROBEID", column="SYMBOL")
);

gse.g <- mapply(
	function(x, pheno) {
		ExpressionSet(x,
			phenoData = AnnotatedDataFrame(pheno),
			featureData = AnnotatedDataFrame(features),
			annotation = annotation(gse)
		)
	},
	x.g, pheno.g,
	SIMPLIFY = FALSE
);

qwrite(res.g, insert(out.fn, tag="limma", ext="rds"));
qwrite(gse.g, insert(out.fn, tag="eset", ext="rds"));

# ---

# integrative differential analysis across cell lines

design <- model.matrix(~ cell_line + group, data = pheno);

fit <- lmFit(x.f, design);
fit <- eBayes(fit);
res <- topTable(fit, coef="grouperlotinib-resistant", number=Inf);

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

