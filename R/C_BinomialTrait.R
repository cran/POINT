# Class \code{BinomialTrait}
#
# Class \code{BinomialTrait} stores null result for binomial trait
#
# @name BinomialTrait-class
# 
# @keywords internal
#
setClass(Class = "BinomialTrait",
         slots = c(deriv = "numeric",
                   invA1Hat = "matrix"),
         prototype = prototype(deriv = numeric(), invA1Hat = matrix()))

# @rdname point-internal-api
setGeneric(name = ".A2InvA1",
           def = function(trait, ...) { standardGeneric(f = ".A2InvA1") })

# @rdname point-internal-api
setMethod(f = ".A2InvA1",
          signature = c(trait = "ANY"),
          definition = function(trait, ...) { stop("not allowed") })

# @rdname point-internal-api
setMethod(f = ".A2InvA1",
          signature = c(trait = "BinomialTrait"),
          definition = function(trait, X, ZZ1, ...) { 
              a2Hat <- crossprod(x = ZZ1*trait@deriv, y = X)
              return( a2Hat %*% trait@invA1Hat )
            })

# function to obtain null result and create appropriate subClass object
.newBinomial <- function(fit, X, ...) {

  pred <- predict(object = fit, type = "response")
  deriv <- drop(unname(obj = pred*{1.0-pred}))
  XX <- crossprod(x = X*deriv, y = X)

  A1HatInv <- tryCatch(solve(a = XX),
                       error = function(e){
                                 print(x = e$message)
                                 stop("unable to invert matrix")
                               })

  res <- new(Class = "BinomialTrait",
             deriv = deriv,
             invA1Hat = A1HatInv)

  return( res )
}
