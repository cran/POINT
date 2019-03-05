# Class \code{GaussianTrait}
#
# Class \code{GaussianTrait} stores null result for Gaussian trait
#
# @name GaussianTrait-class
# 
# @keywords internal
#
# @include C_BinomialTrait.R
#
setClass(Class = "GaussianTrait",
         slots = c(invA1Hat = "matrix"),
         prototype = prototype(invA1Hat = matrix()))

# @rdname point-internal-api
setMethod(f = ".A2InvA1",
          signature = c(trait = "GaussianTrait"),
          definition = function(trait, X, ZZ1, ...) { 
              a2Hat <- crossprod(x = ZZ1, y = X)
              return( a2Hat %*% trait@invA1Hat )
            })

.newGaussian <- function(fit, X, ...) {

  XX <- crossprod(x = X, y = X)

  A1HatInv <- tryCatch(solve(a = XX),
                       error = function(e){
                                 print(x = e$message)
                                 stop("unable to invert matrix")
                               })

  res <- new(Class = "GaussianTrait",
             invA1Hat = A1HatInv)

  return( res )
}
