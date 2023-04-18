library(io)

res.g1 <- qread("../GSE49135/GSE49135_limma.rds");
res.g2 <- qread("../GSE62061/GSE62061_limma.rds");

# match probes so that they are in the same order across studies
probes <- intersect(rownames(res.g1[[1]]), rownames(res.g2[[1]]));
res.g1.m <- lapply(res.g1, function(res) res[probes, ]);
res.g2.m <- lapply(res.g2, function(res) res[probes, ]);

# check that the probes are in the same order
for (i in 1:length(res.g1.m)) {
	for (j in 1:length(res.g2.m)) {
		stopifnot(rownames(res.g1.m[[i]]) == rownames(res.g2.m[[j]]))
	}
}

# ---

smooth_scatter <- function(x, y, ...) {
	xlim <- range(x);
	ylim <- range(y);
	lim <- c(min(xlim[1], ylim[1]), max(xlim[2], ylim[2]));
	smoothScatter(x, y, xlim=lim, ylim=lim, ...);
	abline(a=0, b=1, col="grey30");
}

# ideal correlation for the same cell line (i.e. with itself)
smooth_scatter(res.g1.m[[1]]$t, res.g1.m[[1]]$t);

# cell line 1 in study 1 vs. cell line 1 in study 2
smooth_scatter(res.g1.m[[1]]$t, res.g2.m[[1]]$t);

# cell line 1 vs. 2 in study 2
smooth_scatter(res.g2.m[[1]]$t, res.g2.m[[2]]$t);

# cell line 1 vs. 3 in study 2
smooth_scatter(res.g2.m[[1]]$t, res.g2.m[[3]]$t);

# cell line 1 vs. 4 in study 2
smooth_scatter(res.g2.m[[1]]$t, res.g2.m[[4]]$t);

# ---

# add rank of t statistic to results
add_rank <- function(res) {
	res$rank.t <- rank(res$t);
	res
}

res.g1.m <- lapply(res.g1.m, add_rank);
res.g2.m <- lapply(res.g2.m, add_rank);

# between-study comparisons

# cell line 1 in study 1 vs. cell line 1 in study 2
smooth_scatter(res.g1.m[[1]]$rank.t, res.g2.m[[1]]$rank.t);

# cell line 1 in study 1 vs. cell line 2 in study 2
smooth_scatter(res.g1.m[[1]]$rank.t, res.g2.m[[2]]$rank.t);

# cell line 1 in study 1 vs. cell line 2 in study 2
smooth_scatter(res.g1.m[[1]]$rank.t, res.g2.m[[3]]$rank.t);

# cell line 1 in study 1 vs. cell line 2 in study 2
smooth_scatter(res.g1.m[[1]]$rank.t, res.g2.m[[4]]$rank.t);

# within-study comparisons

# cell line 1 vs. 2 in study 2
smooth_scatter(res.g2.m[[1]]$rank.t, res.g2.m[[2]]$rank.t);

# cell line 1 vs. 3 in study 2
smooth_scatter(res.g2.m[[1]]$rank.t, res.g2.m[[3]]$rank.t);

# cell line 1 vs. 4 in study 2
smooth_scatter(res.g2.m[[1]]$rank.t, res.g2.m[[4]]$rank.t);

