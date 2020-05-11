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