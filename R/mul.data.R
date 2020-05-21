#' Generate multivariate data with the given mean vector and covariance matrix
#'
#' @param n sample size
#' @param K dimension
#' @param mu mean vector
#' @param Sigma covariance matrix
#' @param distribution distribution
#' @importFrom stats cov2cor rnorm rpois runif
#' @return a \code{n*K} matrix
#' @keywords internal
gen.mul.data <- function(n=100,
                         K=2,
                         mu=rep(0,K),
                         Sigma=diag(K),
                         distribution=c('GAUSSIAN',
                                        'LAPLACE',
                                        'EXPONENTIAL')
)
{
    distribution <- toupper(distribution)
    distribution <- match.arg(distribution)
    
    # fix the bug of lcmix::rmvexp when n=1
    my.rmvexp <- function(n,rate,corr)
    {
        if(n > 1)
            lcmix::rmvexp(n=n,rate=lam,corr=cov2cor(Sigma))
        else
        {
            z <- lcmix::rmvexp(n=2,rate=lam,corr=cov2cor(Sigma))
            z <- z[1,]
        }
    }
    
    if(length(mu) != K){
        stop("Make sure that 'K' and the lenght of 'mu' is equal.")
    }
    
    if(distribution == 'GAUSSIAN')
    {
        X <- MASS::mvrnorm(n=n,mu=mu,Sigma=Sigma)
    }
    else if(distribution == 'LAPLACE')
    {
        X <- L1pack::rmLaplace(n=n,center=mu,Scatter=Sigma)
    }
    else if(distribution == 'EXPONENTIAL')
    {
        lam <- 1/sqrt(diag(Sigma))
        X <- my.rmvexp(n=n,rate=lam,corr=cov2cor(Sigma)) + rep.row(mu,n) - rep.row(lam,n)
    }
    else stop(paste0('Distribution ',distribution,' is not supported.'))
    
    return(X)
}
