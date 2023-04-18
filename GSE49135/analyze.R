library(GEOquery)
library(limma)
library(io)

# platform: GPL10558
#           Illumina HumanHT-12 V4.0 expression beadchip
library(illuminaHumanv4.db)
db <- illuminaHumanv4.db;

source("../R/common.R");

accession <- "GSE49135";
out.fn <- filename(accession);

gse <- getGEO(accession)[[1]];

x <- exprs(gse);

# examine distribution of expression data
hist(x, breaks=100);
summary(x)
min(x)

# left-shift towards 0 because x has negative values
x <- x - min(x);
hist(x, breaks=100);

x.f <- filter_undetected(x);

# skip filtering undetected genes
#x.f <- x;

dim(x)
dim(x.f)

# examine distribution of pre-processed data
hist(x.f, breaks=100);
boxplot(x.f)
max(x.f)

pheno <- pData(phenoData(gse))[, c("cell line:ch1", "phenotype:ch1")];
pheno <- clean_colnames(pheno);
pheno$group <- factor(pheno$phenotype,
	levels = c("Erlotinib sensitive", "Erlotinib resistant"),
	labels = c("erlotinib-sensitive", "erlotinib-resistant")
);
pheno$phenotype <- NULL;

# this study only has one cell line
stopifnot(length(unique(pheno$cell_line)) == 1)
cell.line <- pheno$cell_line[1];

res <- diff_expr_erlotinib_resistant(x.f, pheno, db);
res.g <- list(res);
names(res.g) <- cell.line;

# construct feature data
probes <- rownames(x.f);
features <- data.frame(
	row.names = probes,
	ensembl = mapIds(db, keys=probes, keytype="PROBEID", column="ENSEMBL"),
	gene = mapIds(db, keys=probes, keytype="PROBEID", column="SYMBOL")
);

gse2 <- ExpressionSet(x.f,
	phenoData = AnnotatedDataFrame(pheno),
	featureData = AnnotatedDataFrame(features),
	experimentData = experimentData(gse),
	annotation = annotation(gse)
);

gse.g <- list(gse2);
names(gse.g) <- cell.line;

qwrite(res.g, insert(out.fn, tag="limma", ext="rds"));
qwrite(gse.g, insert(out.fn, tag="eset", ext="rds"));

