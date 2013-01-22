library(Rcpp)
library(inline)

# load source code and compile
src <- paste(readLines("baydar.cpp"), collapse="\n")
baydar <- cfunction(signature(d="matrix", o="numeric", cntl="list",
	f="function"), body=src, Rcpp=TRUE, includes=c("#include <gsl/gsl_randist.h>", 
	"#include <gsl/gsl_rng.h>", "#include <gsl/gsl_cdf.h>"), libargs="-lgsl -lgslcblas")
print("Function compiled")

# generate test data
testData <- matrix(nrow=200, ncol=50, data=rnorm(50*200), byrow=TRUE)
testData[11,] <- c(rnorm(25, 1), rnorm(25, -1))
groups <- c(rep(0, 25), rep(1, 25))

# define function that returns statistic and apply to test data
funTStat <- function(z) t.test(x=z[groups == 0], y=z[groups == 1])$statistic
testStats <- apply(testData, 1, funTStat)

# calculate p-values analytically
preal <- apply(testData, 1, function(z) t.test(x=z[groups == 0], 
	y=z[groups == 1])$p.value)

# calculate bootstrap p-values
print("Starting BayDAR")
pboot <- baydar(d=testData, o=testStats, cntl=list(B0=100, B=500, K=100,
	p0=0.05, tail=1), f=funTStat)

