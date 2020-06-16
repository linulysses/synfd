# synfd
Synthesize dense/sparse functional data/snippets

## Install
devtools::install_github("linulysses/synfd")  or  install.packages('synfd')

## Examples

### generate irregularly observed Gaussian functional data with Matern covariance function
Y <- irreg.fd(mu=1, X=gaussian.process(), n=10, m=5)

### generate samples froma a process defined via K-L representation
Y <- irreg.fd(mu=cos, X=kl.process(eigen.values=1/(2^(1:25)),eigen.functions='FOURIER',distribution='LAPLACE'),n=10, m=5)

### generate regularly observed trajectories
Y <- reg.fd(mu=1, X=gaussian.process(cov=matern), n=10, m=101)
