# Class \code{LinearKernel}
#
# Class \code{LinearKernel} holds information when kernel is Linear
#   defined as G W^{1/2} diag(R^2) W^{1/2} G^T
#
# @name LinearKernel-class
#
# @keywords internal
#
# @include A_Kernel.R
setClass(Class = "LinearKernel",
         slots = c(weight = "numeric"),
         contains = c("Kernel"),
         prototype = prototype(weight = numeric()))

# @rdname calcKernel
setMethod(f = ".calcKernel",
          signature = c(kernel = "LinearKernel"),
          definition = function(kernel, Rmc, snp, ...) {
              return( snp %*% { {kernel@weight*Rmc*Rmc} * t(x = snp) } )
            })

# @rdname point-internal-api
.newLinearKernel <- function(weight, ...) {
  if (is.null(x = weight)) {
    return( new(Class = "LinearKernel", weight = 1.0) )
  } else {
    return( new(Class = "LinearKernel", weight = weight) )
  }
}
