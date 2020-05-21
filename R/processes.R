#' create a centered random process objecct
#' @param name the name of the process
#' @param domain the domain where the process is defined
#' @param ... parameters required to define the process
#' \describe{
#'   \item{dispersion}{required by Wiener process}
#'   \item{sig}{required by white noise process}
#'   \item{eigen.values}{required by process defined via K-L representation}
#'   \item{eigen.functions}{required by process defined via K-L representation}
#'   \item{distribution}{the distribution of PC scores, required by process defined via K-L representation}
#' }
#' @return a function hanlde in the form of \code{X(tObs,n)} which generates \code{n} independent trajectories observed at tObs
#' @examples
#' X <- centered.process(name='wiener',dispersion=1)
#' X <- centered.process(name='white.noise',sig=1)
#' X <- centered.process(name='KL',domain=c(0,1),eigen.values=1/(2^(1:25)),eigen.functions='FOURIER',distribution='GAUSSIAN',corr=NULL)
#' X(regular.grid(50),25)
#' @export
centered.process <- function(name=c('WIENER',
                                    'WHITE.NOISE',
                                    'KL',
                                    'KARHUNEN.LOEVE',
                                    'GAUSSIAN'),
                             domain=c(0,1),
                             ...)
{
    name <- toupper(name)
    name <- match.arg(name)
    others <- list(...)
    if (name == 'WIENER')
    {
        X <- wiener.process(dispersion=get.required.param('dispersion',others))
        attr(X,'name') <- name
        attr(X,'dispersion') <- get.required.param('dispersion',others)
        class(X) <- 'random.process'
    }
    else if (name == 'WHITE.NOISE')
    {
        X <- white.noise(sig=get.required.param('sig',others))
        attr(X,'name') <- name
        attr(X,'sig') <- get.required.param('sig',others)
        class(X) <- 'random.process'
    }
    else if (name == 'KL' || name == 'KARHUNEN.LOEVE')
    {
        X <- kl.process(domain=domain,
                        eigen.values=get.required.param('eigen.values',others),
                        eigen.functions=get.required.param('eigen.functions',others),
                        distribution=get.required.param('distribution',others),
                        corr=get.required.param('corr',others))
        attr(X,'name') <- name
        attr(X,'domain') <- domain
        attr(X,'eigen.values') <- get.required.param('eigen.values',others)
        attr(X,'eigen.functions') <- get.required.param('eigen.functions',others)
        attr(X,'distribution') <- get.required.param('distribution',others)
        attr(X,'corr') <- get.required.param('corr',others)
        class(X) <- 'random.process'
    }
    else if (name == 'GAUSSIAN')
    {
        X <- gaussian.process(cov=get.required.param('cov',others))
        attr(X,'name') <- name
        attr(X,'cov') <- get.required.param('cov',others)
        class(X) <- 'random.process'
    }
    else
        stop(paste0('The named process ', name, ' is not supported.'))
    
    return(X)
}



#' create a centered Gaussian process objecct
#' @param cov a function handle that defines covariance function. It shall take two arguments arg1 and arg2, both are vectors, and \code{cov(arg1,arg2)} returns a matrix \code{R} such that \code{R(i,j)} is the value of the covariance function at \code{(arg1[i],arg2[j])}
#' @return a function hanlde in the form of \code{X(tObs,n)} which generates \code{n} independent Gaussian trajectories observed at tObs
#' @examples
#' X <- gaussian.process()
#' X(regular.grid(50),25)
#' @export
gaussian.process <- function(cov=matern)
{
    f <- function(tObs,n)
        gen.mul.data(n,length(tObs),Sigma=cov(tObs))
    return(f)
}


#' create a Wiener process objecct
#' @param dispersion the dispersion parameter of the Wiener process
#' @return a function hanlde in the form of \code{X(tObs,n)} which generates \code{n} independent Wiener trajectories observed at tObs
#' @examples
#' X <- wiener.process()
#' X(regular.grid(50),25)
#' @export
wiener.process <- function(dispersion=1)
{
    f <- function(tObs,n)
    {
        K <- 50 # number of eigenfunctions used (for approximation)
        tObs <- as.matrix(tObs)
        if (dim(tObs)[1] < dim(tObs)[2])
            tObs <- t(tObs)
        basis <- sqrt(2) * sin(tObs %*% matrix(1:K - 1/2, 1, K) * pi)
        samp <- t(basis %*% diag(1/(1:K - 1/2)/pi, K) %*% matrix(rnorm(K * n), K, n))
        return(samp*dispersion)
    }
    return(f)
}


#' create a white noise process objecct
#' @param sig the sigma parameter of the white noise process
#' @return a function hanlde in the form of \code{X(tObs,n)} which generates \code{n} independent Wiener trajectories observed at tObs
#' @examples
#' X <- white.noise()
#' X(regular.grid(50),25)
#' @export
white.noise <- function(sig=1)
{
    f <- function(tObs,n)
        return(matrix(rnorm(n*length(tObs))*sig,n,length(tObs)))
    return(f)
}

#' create a process objecct via K-L representation
#' @param domain the domain
#' @param eigen.values the eigenvalues
#' @param eigen.functions the eigenfunctions
#' @param distribution the distribution of PC scores
#' @param corr a correlation matrix specifying correlation among the random coefficients (default: NULL)
#' @return a function hanlde in the form of \code{X(tObs,n)} which generates \code{n} independent K-L trajectories observed at tObs
#' @examples
#' X <- kl.process()
#' X(regular.grid(50),25)
#' @export
kl.process <- function(domain=c(0,1),
                       eigen.values=1/(2^(1:25)),
                       eigen.functions=c('FOURIER', 'COS', 'SIN', 'LEGENDRE'),
                       distribution=c('GAUSSIAN',
                                      'LAPLACE',
                                      'EXPONENTIAL',
                                      'GAMMA'),
                       corr=NULL)
{
    distribution = match.arg(toupper(distribution),distribution)
    eigen.functions = match.arg(toupper(eigen.functions),eigen.functions)
    
    Sigma <- diag(eigen.values)
    if(!is.null(corr)){
        Sigma <-diag(sqrt(eigen.values))
        Sigma <- Sigma %*% corr %*% Sigma
    }
    
    
    
    f <- function(tObs,n)
        sample.from.kl.process(tObs=tObs,
                               domain=domain,
                               eigen.values=eigen.values,
                               eigen.functions=eigen.functions,
                               distribution=distribution,
                               n=n,
                               Sigma=Sigma)
    
    attr(f,'eigen.values') <- eigen.values
    attr(f,'eigen.functions') <- eigen.functions
    attr(f,'Sigma') <- Sigma
    attr(f, 'corr') <- corr
    
    return(f)
}


sample.from.kl.process <- function(tObs=regular.grid(100),
                                   domain=c(0,1),
                                   eigen.values=1/(2^(1:25)),
                                   eigen.functions='Fourier',
                                   distribution='GAUSSIAN',
                                   n=1,
                                   Sigma=NULL)
{
    
    k <- length(eigen.values)
    m <- length(tObs)
    
    phi <- evaluate.basis(K=k,
                                domain=domain,
                                grid=tObs,
                                type=eigen.functions)
    
    
    
    xi <- gen.mul.data(n=n, K=k,
                       mu=rep(0,k),
                       Sigma=Sigma,
                       distribution=distribution)
    
    y0 <- xi %*% t(phi)
    
    attr(y0,'phi') <- phi
    attr(y0,'xi') <- xi
    class(y0) <- 'kl.process.samples'
    
    return(y0)
}

#' The Matern covariance function
#' @param x a vector
#' @param y If NULL, then y=x.
#' @param nu see https://en.wikipedia.org/wiki/Mat?rn_covariance_function
#' @param rho see https://en.wikipedia.org/wiki/Mat?rn_covariance_function
#' @param sig see https://en.wikipedia.org/wiki/Mat?rn_covariance_function
#' @return A matrix M of len(X)*len(Y) is returned, where M(i,j)=C(X(i),Y(j))
#' @export
matern <- function(x,y=NULL,nu=1,rho=1,sig=1)
{
    if(is.null(y)) y <- x
    
    G <- expand.grid(x,x)
    
    S <- apply(G,1,function(z){
        delta <- abs(z[1]-z[2])/rho
        if(delta == 0)
            0.25^2
        else
            0.25^2 * (sqrt(2*nu)*delta)^nu * besselK(sqrt(2*nu)*delta,nu=nu) / (2^(nu-1)*gamma(nu))
    })
    
    C <- sig^2 * matrix(S,nrow=length(x),ncol=length(y),byrow=F)
    return(C)
}

