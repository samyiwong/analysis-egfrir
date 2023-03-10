library(GEOquery)

gse <- getGEO("GSE83666");
gse <- gse[[1]];

x <- exprs(gse);
hist(x, breaks=100);

x <- log2(x + 1);
hist(x, breaks=100);
summary(x)

# minimum value appears to be 4.392317 (atypical)
# left shift x towards 0
x <- x - min(x); 

expr.cut <- 0.5;
prop.cut <- 0.2;

p.detected <- apply(x, 1, function(z) mean(z > expr.cut));
summary(p.detected)
hist(p.detected, breaks=50)

x.f <- x[p.detected > prop.cut, ];

dim(x)
dim(x.f)

hist(x.f, breaks=100)

# every few genes are expressed: poor hybridization?
