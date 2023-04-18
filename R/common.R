filter_undetected <- function(x, expr.cut=0.5, prop.cut = 0.2) {
	p.detected <- apply(x, 1, function(z) mean(z > expr.cut));
	idx <- p.detected > prop.cut;
	x.f <- x[idx, ];
	attr(x.f, "p.detected") <- p.detected;
	attr(x.f, "selected") <- idx;
	x.f
}

clean_colnames <- function(d) {
	colnames(d) <- gsub(" ", "_",
		sub(":ch1", "", colnames(d))
	);
	d
}

# annotate limma results
# each row is a probe
annotate_probes <- function(res, db) {
	# map probe id to gene symbol
	res$gene <- mapIds(
		db, keys=rownames(res), keytype="PROBEID",
		column="SYMBOL"
	);

	# map probe id to ensembl id
	res$ensembl <- mapIds(
		db, keys=rownames(res), keytype="PROBEID",
		column="ENSEMBL"
	);

	res
}

# perform differential gene expression 
diff_expr_erlotinib_resistant <- function(x, pheno, db) {
	# ensure that samples are in the same order
	stopifnot(colnames(x) == rownames(pheno))

	design <- model.matrix(~ group, data = pheno);

	fit <- lmFit(x, design);
	fit <- eBayes(fit);
	res <- topTable(fit, coef="grouperlotinib-resistant", number=Inf);
	annotate_probes(res, db)
}

