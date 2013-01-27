Bayesian Allocation for Differential Resampling
======
# Motivation
In genomic studies, one is often faced with mutiple-hypotheses testing
scenarios (~1E+06). When analytical formulas are available the main problem is
determining an appropriate threshold for rejecting the null hypotheses and
controlling the family-wide error rate (FWER) or false discovery rate (FDR).
However, often in the case of some test statistics an analytical formula is
unavailable or the assumptions underpinning such a formula are unrealistic (e.g.
small sample sizes or non-Gaussian data). In these cases, a resampling-based
test is desirable. However, this approach can be computationally infeasible: To
reject H0 at a significance threshold of 1E-06, which is common in GWAS studies,
at 1 million SNPs would require 1 trillion resamples. 
# Algorithm
To avoid this pitall, [Wang *et al.*](http://www.biomedcentral.com/1471-2105/10/198/)
proposed a Bayesian scheme  for the differential allocation of resamples. The
algorithm is described in detail in the paper but can succinctly described as
follows- 
> "...we use a Bayesian-inspired approach that assigns resamples to each unit
>based on its individual risk, the chance that the current p-value estimate
>leads to a misclassification of the unit. The goal is to lower the numbers of
>classification errors, since we are giving a higher resolution to the null
>distribution of genes that are more likely to be misclassified in a uniform
>allocation setting. This higher resolution comes at the sacrifice of resamples
>to non-borderline genes that should not need a very resolute null distribution
>for correct inference." 

The intuition for the risk of misclassification can be visually summarized: <img
src="http://www.biomedcentral.com/content/figures/1471-2105-10-198-1.jpg" />
Here, two different densities are visualized; the red density corresponds to a
scenario where the true p-value, p1, is close to p0 and thus there is more area
under the curve; the blue density corresponds to a scenario where the true
p-value, p2, is farther away from p0 and there there is less area under  the
curve. The number of resamples allocated to both p1 and p2 are proportional to
the shaded areas, which visualize the risk of misclassification, respectively.


# Software 
## Description
The code availble in this repository can be used to apply this algorithm to a
matrix with a user-provided function to calculate the test-statistic.
## Usage

