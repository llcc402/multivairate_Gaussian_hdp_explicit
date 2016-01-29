# The infinite mixture of multivariate Gaussians (HDP)

We use the same method as it in the repo of hdp_inf_gaussian_explicit_git. The only different is that we use multivariate Gaussian instead of one variable Gaussians.

## Remark:

The difficulty in the multivariate Gaussian is that when the dimension of the data grows, the pdf of the Gaussians degerate extremely fast when leaving the center. In the sampling of the positions, we generate samples with standard deviation a function of the std of the data set. In the computation of the likelihood, we use Gaussain kernel instead of mvnpdf because it is more steady.

## Functions included:

1. data_generate: Generate multivariate mixture of Gaussians. Each group use some cluster.

2. hdp: The as the one dimensional case.

## Scripts included:

1. main: An example of the repo

## Figures included:

1. whole_set.jpg: The scatter plot of the whole data set.

2. group1.jpg: The scatter plot of the first group.

3. group1_pos.jpg: Use different color to show the posterior of the indicators.


