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
