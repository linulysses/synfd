# synfd
Synthesize dense/sparse functional data/snippets

## Install
devtools::install_github("linulysses/synfd")

## Examples

### generate sparsely observed Gaussian functional data with Matern covariance function
Y <- sparse.fd(mu=1, X=gaussian.process(), n=10, m=5)

### generate samples froma a process defined via K-L representation
Y <- sparse.fd(mu=cos, X=kl.process(eigen.values=1/(2^(1:25)),eigen.functions='FOURIER',distribution='LAPLACE'),n=10, m=5)

### generate samples with mi sample from 2,4,6
Y <- sparse.fd(mu=1, X=gaussian.process(cov=matern), n=10, m=c(2,4,6))
