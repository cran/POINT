.validity_PolyKernel <- function(object) {

  if (isTRUE(all.equal(current = object@kparam, target = 0.0))) {
    return( "d cannot be zero for polynomial kernel" )
  }

  return(TRUE)

}

# Class \code{PolyKernel}
#
# Class \code{PolyKernel} holds information when kernel is polynomial
#   defined as { 1 + G W^{1/2} diag(R^2) W^{1/2} G^T }^d
#
# @name PolyKernel-class
#
# @slot kparam The order of the polynomial, "d".
# @keywords internal
#
# @include A_Kernel.R
setClass(Class = "PolyKernel",
         slots = c(kparam = "numeric",
                   weight = "numeric"),
         contains = c("Kernel"),
         prototype = prototype("kparam" = numeric(), "weight" = numeric()),
         validity = .validity_PolyKernel)

# @rdname calcKernel
setMethod(f = ".calcKernel",
          signature = c(kernel = "PolyKernel"),
          definition = function(kernel, Rmc, snp, ...) {
              return( { {snp %*% { {kernel@weight*Rmc*Rmc} * t(x = snp) }} +
                        1.0 }^kernel@kparam )
            })

# @rdname point-internal-api
.newPolyKernel <- function(weight, d, ...) {
  if (is.null(x = weight)) {
    return( new(Class = "PolyKernel", kparam = d, weight = 1.0) )
  } else {
    return( new(Class = "PolyKernel", kparam = d, weight = weight) )
  }
}
