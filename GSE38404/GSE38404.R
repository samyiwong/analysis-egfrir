#GSE38404

library(GEOquery)
library(limma)
library(ggplot2)
library(dplyr)

gse_1 <- getGEO("GSE38404")[[1]];

x <- exprs(gse_1);
hist(x, breaks=100);
summary(x)
min(x)

x <- x - min(x);
hist(x, breaks=100);

x <- log2(x+1)
hist(x, breaks=50)

x.f <- filter_undetected(x);

dim(x)
dim(x.f)

hist(x.f, breaks=100);
boxplot(x.f)
max(x.f)

View(pData(phenoData(gse_1)))

pheno_1 <- pData(phenoData(gse_1))[, c("title", "pf299804:ch1")];
#pheno_1 <- clean_colnames(pheno);
colnames(pheno_1)[2]  <- "phenotype"
pheno_1$phenotype <- factor(pheno_1$phenotype,
                          levels = c("Sensitive", "Resistant"),
                          labels = c("gefitinib-sensitive", "gefitinib-resistant")
);

design_1 <- model.matrix( ~ phenotype, data = pheno_1);
fit_1 <- lmFit(x.f, design_1);
fit_1 <- eBayes(fit_1);
res_1 <- topTable(fit_1, coef="phenotypegefitinib-resistant", number=Inf);

rownames(res_1)
View(pData(featureData(gse_1)))
idx_1 <- match(rownames(res_1), pData(featureData(gse_1))$ID);
res_1$gene <- pData(featureData(gse_1))[idx_1, "Symbol"];

y <- pData(featureData(gse_1))[,c(1,11,12)]     
colnames(y)[2]  <- "gene_symbol"
colnames(y)[3]  <- "entrez_gene_id"
#no,this is wrong. res_1<-cbind(res_1, y$gene_symbol,y$entrez_gene_id)
#no, x1 <- cbind(y,x)
#no, match(rownames(res_1), )

ID<-rownames(res_1)
res_1<- cbind(ID, res_1)
res_1_1 <- inner_join(res_1, y, by="ID")

rank(res_1_1$t)

