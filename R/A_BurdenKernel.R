# Class \code{BurdenKernel}
#
# Class \code{BurdenKernel} holds information when kernel is burden
#   defined as (G W^{1/2} R) (G W^{1/2} R)^T
#
# @name BurdenKernel-class
#
# @keywords internal
#
# @include A_Kernel.R
setClass(Class = "BurdenKernel",
         slots = c(weight = "numeric"),
         contains = c("Kernel"),
         prototype = prototype(weight = numeric()))

# @rdname calcKernel
setMethod(f = ".calcKernel",
          signature = c(kernel = "BurdenKernel"),
          definition = function(kernel, Rmc, snp, ...) {
              tmp <- t(x = snp)*{kernel@weight*Rmc}
              return( as.matrix(x = crossprod(x = tmp)) )
            })

# @rdname point-internal-api
.newBurdenKernel <- function(weight, ...) {
  if (is.null(x = weight)) {
    return( new(Class = "BurdenKernel", weight = 1.0) )
  } else {
    return( new(Class = "BurdenKernel", weight = sqrt(x = weight)) )
  }
}

# Calculate the local kernel matrix and its rank factorization
#
# @param kernel object of class BurdenKernel
# @param Rmc vector subset of variant similarity matrix R for the mth marker 
#   and cvth c value (R[,m,cv])
# @param snp nxM genotype snp matrix
# @param ... ignored
#
# @return {nS x nC} matrix containing rank factorization of kernel matrix
#
# @keywords internal
#
# @rdname calcLocalKernel
setMethod(f = ".calcLocalKernel",
          signature = c(kernel = "BurdenKernel"),
          definition = function(kernel, Rmc, snp, ...) {

              eValues <- {kernel@weight * Rmc}*{kernel@weight * Rmc}
              ieValues <- which(x = eValues > 1e-10)
              if (length(x = ieValues) == 0L) stop("no positive eigen values")
              csnp <- colSums(x = snp)

              eValues <- eValues[ieValues]*csnp[ieValues]
              eVectors <- t(x = snp[,ieValues,drop=FALSE]) / 
                          sqrt(x = csnp[ieValues])

              return( t(x = {sqrt(x = eValues) * eVectors} ) )
  
            })
