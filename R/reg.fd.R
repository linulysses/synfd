#' Sample Regular Functional Data
#' @param mu function, scalar or a vector defining the mean function; default value: \code{0}.
#' @param X centered stochastic process defined by a function of the form 
#'          \code{X(tObs,n)} and returning \code{n*length(tObs)} matrix, where each row represents observations from a trajectory. Default value: \code{wiener.process()}.
#' @param n sample size; default value: \code{100}.
#' @param m sampling rate; ignored if \code{grid} is specified; default value: \code{51}.
#' @param domain the domain; default value: \code{c(0,1)}.
#' @param grid vector of design points; default value: \code{NULL}.
#' @param sig standard deviation of measurement errors; if \code{NULL} then determined by \code{snr}.
#' @param snr  signal to noise ratio to determine \code{sig}; default value: \code{5}.
#' @return a list with the following members 
#'      \describe{
#'          \item{\code{t}}{a vector of design points sorted in increasing order.}
#'          \item{\code{y}}{\code{n*m} matrix; each row represents observations from a trajectory.}
#'      }
#'  and with attributes \code{sig}, \code{snr}, \code{domain}, \code{grid} and 
#'      \describe{
#'          \item{y0}{\code{n*m} matrix of observations without measurement errors.}
#'      }     
#'      
#' @examples
#' Y <- reg.fd()
#' Y <- reg.fd(mu=1, X=gaussian.process(), n=10, m=20)
#' Y <- reg.fd(mu=cos, X=kl.process(),n=100, m=20)
#' Y <- reg.fd(mu=cos, X=kl.process(distribution='EXPONENTIAL'),n=100, m=20)
#' @export
reg.fd <- function(mu=0, X=wiener.process(), n=100, m=51,
                     sig=NULL, snr=5, domain=c(0,1),
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
    
    
    y0 <- rep.row(mu,n) + X(grid,n)
 
    
    Z <- scale(y0,center=TRUE,scale=FALSE)
    s <- mean(apply(Z^2,2,mean))
    
    if(is.null(sig))
    {
        if(!is.infinite(snr))
        {
            sig <- sqrt(s/snr)
        }
        else sig <- 0
    }
    else
    {
        snr <- s/(sig^2)
    }
    
    
    y <- y0
    if(sig != 0)
    {
        y <- y + sig * matrix(rnorm(n*m),nrow=n)
    }
    
    R <- list(t=grid,y=y)
    
    attr(R,'sig') <- sig
    attr(R,'snr') <- snr
    attr(R,'y0') <- y0
    attr(R,'grid') <- grid
    attr(R,'domain') <- domain
    attr(R,'class') <- 'dense.fd'
    return(R)
}

#' Plot Regular Functional Data
#' @param x the data object generated from \code{\link{reg.fd}}.
#' @param ... other parameters passed to \code{matplot}.
#' @importFrom graphics matplot
#' @return a plot of the dataset \code{x}.
#' @export
plot.dense.fd <- function(x,...)
{
    matplot(attr(x,'grid'),t(x),...)
}
