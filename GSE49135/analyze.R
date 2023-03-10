library(GEOquery)
library(limma)

source("../R/common.R");

gse <- getGEO("GSE49135")[[1]];

x <- exprs(gse);
hist(x, breaks=100);
summary(x)
min(x)

# left-shift towards 0
x <- x - min(x);
hist(x, breaks=100);


x.f <- filter_undetected(x);

dim(x)
dim(x.f)

hist(x.f, breaks=100);
boxplot(x.f)
max(x.f)

pheno <- pData(phenoData(gse))[, c("cell line:ch1", "phenotype:ch1")];
pheno <- clean_colnames(pheno);
pheno$phenotype <- factor(pheno$phenotype,
	levels = c("Erlotinib sensitive", "Erlotinib resistant"),
	labels = c("erlotinib-sensitive", "erlotinib-resistant")
);

design <- model.matrix( ~ phenotype, data = pheno);
fit <- lmFit(x.f, design);
fit <- eBayes(fit);
res <- topTable(fit, coef="phenotypeerlotinib-resistant", number=Inf);

idx <- match(rownames(res), pData(featureData(gse))$ID);
res$gene <- pData(featureData(gse))[idx, "Symbol"];

head(res, 100)
pheno
x.f["ILMN_1666733", ]
x.f["ILMN_2181593", ]
x.f[c("ILMN_2192072", "ILMN_1685403"), ]

