library(GEOquery)

source("../R/common.R");

gse <- getGEO("GSE83666")[[1]];

x <- exprs(gse);
hist(x, breaks=100);

x <- log2(x + 1);
hist(x, breaks=100);
summary(x)

# minimum value appears to be 4.392317 (atypical)
# left shift x towards 0
x <- x - min(x); 

x.f <- filter_undetected(x);

dim(x)
dim(x.f)

hist(x.f, breaks=100)

# every few genes are expressed: poor hybridization?

# ---

supp <- getGEOSuppFiles("GSE83666");
