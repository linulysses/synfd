#' Sample Functional Data
#' @param mu function, scalar or a vector defining the mean function; default value: \code{0}.
#' @param X centered stochastic process defined by a function of the form 
#'          \code{X(tObs,n)} and returning \code{n*length(tObs)} matrix, where each row represents observations from a trajectory. Default value: \code{wiener.process()}.
#' @param n sample size; default value: \code{100}.
#' @param m average sampling rate; default value: \code{5}.
#' @param sig standard deviation of measurement errors; if \code{NULL} then determined by \code{snr}.
#' @param snr  signal to noise ratio to determine \code{sig}; default value: \code{5}.
#' @param domain the domain; default value: \code{c(0,1)}.
#' @param delta the proportion of the domain to be observed for each trajectory; only required when \code{type="irregular"}; default value: \code{1}.
#' @param grid vector of design points; only required when \code{type="regular"}. default value: \code{NULL}.
#' @param type the data type, either "regular" or "irregular"; default value: "irregular".
#' @details This is a unified interface for \code{\link{reg.fd}} and \code{\link{irreg.fd}}; see the manual for these functions for details.
#' @return a list with the following members 
#'      \describe{
#'          \item{\code{t}}{design points sorted in increasing order for each trajectory; if \code{type="irregular"} then a list; otherwise a vector.}
#'          \item{\code{y}}{observations for each trajectory; if \code{type="irregular"} then a list; otherwise, a matrix.}
#'      }
#'  and with attributes \code{sig}, \code{snr}, \code{domain}, \code{delta}, \code{grid} and
#'      \describe{
#'          \item{y0}{the measurement-error-free counterpart of \code{y}.}
#'      }    
#'      
#' @references 
#' \insertRef{Lin2020}{synfd}
#' 
#' @examples
#' # irregularly observed Gaussian trajectories with constant mean function 1
#' Y <- rfd(mu=1, X=gaussian.process(), n=10, m=5)
#'
#' # regularly observed trajectories with a K-L representation
#' Y <- rfd(X=kl.process(eigen.functions='FOURIER',distribution='LAPLACE'), type='regular')
#' @export
rfd <- function(mu=0, X=wiener.process(), n=100, m=5,
                sig=NULL, snr=5, domain=c(0,1),delta=1,
                grid = seq(domain[1], domain[2], length.out = m),
                type='irregular')
{
    if(type=='irregular') irreg.fd(mu,X,n,m,sig,snr,domain,delta)
    else reg.fd(mu,X,n,m,sig,snr,domain,grid)
}