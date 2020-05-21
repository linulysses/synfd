get.required.param <- function(name,params)
{
    i <- match(name, names(params)) 
    if (is.na(i)) { 
        stop(paste0(name, ' is not provided.'))
    } 
    params[[i]]
}

#' generate a regular grid of points on a given interval
#' @param M the number of points
#' @param domain the interval
#' @param h margin at boundaries
#' @export
regular.grid <- function(M=100,domain=c(0,1),h=1/(2*M))
{
    seq(domain[1]+h,domain[2]-h,length.out=M)
}



#' replicate a row vector into a matrix
#' @param x the row vector to be replicated
#' @param ... further arguments
#' \describe{
#' \item{n}{The number of times of replication}
#' }
#' @return a matrix of the dimension of \code{n} by \code{length(x)}
#' @keywords internal
rep.row <- function(x,n){
    matrix(rep(x,each=n),nrow=n)
}

#' replicate a column vector into a matrix
#' @param x the column vector to be replicated
#' @param n replicate the column vector n times
#' @return a matrix of the dimension of \code{length(x)} by \code{n}
#' @rdname rep.col
#' @keywords internal
rep.col <- function(x,n){
    matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}


#' Evaluate orthonormal basis functions on a grid
#'
#' @param K A positive integer specifying the number of eigenfunctions to generate.
#' @param domain the domain on which basis functions are defined
#' @param m the number of equispaced points on \code{domain}
#' @param grid A vector specifying the time points to evaluate the basis functions. If \code{grid} is supplied, then \code{m} is ignored
#' @param type A string for the type of orthogonal basis.
#' @return A \code{m} by \code{K} matrix, where rows index basis functions while columns index points in the grid.
#'
#' @examples
#' basis <- evaluate.basis(3, type='fourier')
#' head(basis)
#' 
#' @export evaluate.basis

evaluate.basis <- function(K,
                           m = 51,
                           domain = c(0,1),
                           grid = seq(domain[1], domain[2], length.out = m),
                           type = c('FOURIER', 'COS', 'SIN', 'LEGENDRE'))
{
    type <- toupper(type)
    type <- match.arg(type)
    stopifnot(is.numeric(K) && length(K) == 1 && K > 0)
    
    if (min(grid) < domain[1] || max(grid) > domain[2])
        stop('some points in grid are outside of domain.')
    
    m <- length(grid)
    
    # compute scaling factor accouting for non-canonical domains
    # (x,y) represents canonical domains
    # (a,b) represents the user domain
    a <- domain[1]
    b <- domain[2]
    if (type %in% c('COS','SIN','FOURIER'))
    {
        x <- 0
        y <- 1
    }
    else
    {
        x <- -1
        y <- 1
    }
    
    alpha <- (b*x-a*y)/(x-y)
    beta <- (a-b)/(x-y)
    
    # points in the canonical domain corresponding to points of grid
    pts <- (grid - alpha) / beta
    
    if (type == 'COS') {
        res <- sapply(1:K, function(k) {
            if (k == 1) rep(1, m)
            else sqrt(2) * cos(2 * (k - 1) * pi * pts)
        })
    } else if (type == 'SIN') {
        res <- sapply(1:K, function(k) {
            if (k == 1) rep(1,m)
            else sqrt(2) * sin(2 * k * pi * pts)
        })
    } else if (type == 'FOURIER') {
        res <- sapply(1:K, function(k)
            if (k == 1) rep(1, m)
            else if (k %% 2 == 0)  sqrt(2) * cos(k * pi * pts)
            else sqrt(2) * sin((k - 1) * pi * pts)
        )
    } else if (type == 'LEGENDRE') {
        stop('Legendre basis is not supported yet')
        if (K == 1) {
            res <- matrix(1, length(pts), 1)
        } else if (K > 1) {
            coefMat <- sapply(2:K, function(n) {
                coef <- rep(0, K)
                coef[1] <- 1
                for (k in seq_len(n-1)) {
                    coef[k+1] <- choose(n-1, k) *
                        choose(n-1 + k, k)
                }
                coef * sqrt((2 * n - 1)/2)
            })
            xMat <- cbind(1, stats::poly((pts-1)/2, K - 1, raw = TRUE))
            res <- cbind(rep(sqrt(1/2),length(pts)), xMat %*% coefMat)
        }
        
        if (K >= 25) {
            warning('Numeric unstability may occur. Use K < 25.')
        }
    } else if (type == 'POLY') {
        if (K == 1) {
            res <- matrix(1, length(pts), 1)
        } else if (K > 1) {
            res <- cbind(1, stats::poly(pts, K - 1, raw = TRUE))
        }
        
        if (K >= 25) {
            warning('Numeric unstability may occur. Use K < 25.')
        }
    } else{
        stop('unknown basis type')
    }
    
    res <- res / sqrt(beta) # normalization
    
    res <- matrix(res, ncol = K) # prevent single length pts
    return(res)
}
