.validity_NullResult <- function(object) {

  if (!is(object = object@trait, class2 = "BinomialTrait") &&
      !is(object = object@trait, class2 = "GaussianTrait") &&
      !is.na(x = object@trait)) {
    return("@trait is not of correct type")
  }

  return( TRUE )
}

setClass(Class = "NullResult",
         slots = c(resid = "numeric",
                   trait = "ANY",
                   fit = "ANY"),
         prototype = prototype(resid = numeric(), trait = NA, fit = NA),
         validity = .validity_NullResult)

# Calculate the Influence Function
#
# Given the Null result, X and ZZ1, calculate the influence function
#
# @param nullResult An object of class NullResult
# @param ... X and ZZ1 objects
#
# @name calcPsi
#
# @keywords internal
setGeneric(name = ".calcPsi",
           def = function(nullResult, ...) { standardGeneric(f = ".calcPsi") })

# @rdname calcPsi
setMethod(f = ".calcPsi",
          signature = c(nullResult = "ANY"),
          definition = function(nullResult, ...) { stop("not allowed") })

# @rdname calcPsi
setMethod(f = ".calcPsi",
          signature = c(nullResult = "NullResult"),
          definition = function(nullResult, X, ZZ1) {

              # estimate of psi {Z11 - A2 A1inv X}*eps
              return( (ZZ1 - 
                       tcrossprod(x = X, 
                                  y = .A2InvA1(trait = nullResult@trait, 
                                               X = X,  
                                               ZZ1 = ZZ1)))*
                       nullResult@resid )
            })


# function to obtain null result
#' @importFrom stats glm predict residuals
.newNullResult <- function(trait, X, yy, verbose) {

  fit <- tryCatch(expr = glm(formula = yy ~ X - 1, 
                             data = data.frame(yy,X),
                             family = trait),
                  error = function(e){
                            print(x = e$message)
                            stop("unable to obtain glm fit", call. = FALSE)
                          },
                  warning = function(w){
                            print(x = paste("Warning from glm:", w$message))
                          })

  if (verbose) {
    cat("\nNULL model regression\n")
    print(fit)
  }

  # estimate of A1
  # nC x nC matrix
  if (substr(x = trait, start = 1L, stop = 1L) == "b") {
    traitResult <- .newBinomial(fit = fit, X = X, verbose = verbose)
  } else if (substr(x = trait, start = 1L, stop = 1L) == "g") {
    traitResult <- .newGaussian(fit = fit, X = X, verbose = verbose)
  } else {
    stop("trait type not supported", call. = FALSE)
  }

  res <- new(Class = "NullResult",
             fit = fit,
             resid = residuals(object = fit, type = "working"),
             trait = traitResult)

  return( res )
}
