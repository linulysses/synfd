#' Generate functional data on a dense grid
#' @param mu a function, scalar or a vector defining the mean function
#' @param X the centered stochastic process, a function that
#'          X(tObs,n) return n*len(tObs) matrix with n trajectories
#' @param n sample size
#' @param m sampling rate, ignored if \code{grid} is specified
#' @param domain the domain
#' @param grid the grid on which the trajectories are observed
#' @param sig0 the std of measurement errors, if NULL, determined by snr
#' @param snr  the signal to noise ratio to determine sig0
#' @return a \code{n*m} matrix with the following attributes
#'      \item{sig0}
#'      \item{snr}{signal-to-nois ratio}
#'      \item{uncontaminated}{n*m matrix of clean observations}
#'      \item{grid}{a grid of points in \code{domain}}
#'      \item{domain}{the domain}
#' @examples
#' Y <- dense.fd(mu=1, X=gaussian.process(), n=10, m=20)
#' Y <- dense.fd(mu=cos, X=kl.process(),n=100, m=20)
#' @export
dense.fd <- function(mu, X, n, m,
                     sig0=NULL, snr=5, domain=c(0,1),
                     grid = seq(domain[1], domain[2], length.out = m))
{
    
    if (m != length(grid))
        warning('m is ignored, since grid is provided.')
    if (min(grid) < domain[1] || max(grid) > domain[2])
        stop('some points in grid are outside of domain.')
    
    
    if (is.function(mu)) mu <- mu(grid)
    else if(length(mu)==1) mu <- rep(mu,m)
    if (length(mu) != m)
        stop('If mu is a vector, it must be of the same length of grid.')
    
    
    y0 <- mcfda::rep.row(mu,n) + X(grid,n)
    
    Z <- scale(y0,center=TRUE,scale=FALSE)
    s <- mean(apply(Z^2,2,mean))
    
    if(is.null(sig0))
    {
        if(!is.infinite(snr))
        {
            sig0 <- sqrt(s/snr)
        }
        else sig0 <- 0
    }
    else
    {
        snr <- s/(sig0^2)
    }
    
    
    y <- y0
    if(sig0 != 0)
    {
        y <- y + sig0 * matrix(rnorm(n*m),nrow=n)
    }
    
    attr(y,'sig0') <- sig0
    attr(y,'snr') <- snr
    attr(y,'uncontaminated') <- y0
    attr(y,'grid') <- grid
    attr(y,'domain') <- domain
    attr(y,'class') <- 'dense.fd'
    return(y)
}

#' plot densely and regularly observed data
#' @param y the dense data object generated from \code{dense.fd}
#' @param ... other parameters passed to \code{matplot}
#' @export
plot.dense.fd <- function(y,...)
{
    matplot(attr(y,'grid'),t(y),...)
}
