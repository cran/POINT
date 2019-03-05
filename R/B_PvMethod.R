# Class \code{PvMethod}
#
# Class \code{PvMethod} stores arguments for p-value method
#
# @name PvMethod-class
# 
# @slot args list of arguments to be passed to method
#
# @keywords internal
setClass(Class = "PvMethod",
         slots = c(args = "list"))

# Calculate the p-value
#
# Given influence function and eigenvalues calculate appropriate p-value
#
# @param method An object inheriting from PvMethod
# @param ... psi and ev objects
#
# @name calcPV
#
# @keywords internal
setGeneric(name = ".calcPV",
           def = function(method, ...) { standardGeneric(f = ".calcPV") })

# @rdname calcPV
setMethod(f = ".calcPV",
          signature = c(method = "ANY"),
          definition = function(method, ...) { stop("not allowed") })



# Creates Appropriate PvMethod
#
# Given user specified pv method and arguments, creates appropriate
#   object that inherits from PvMethod.
#
# @name newPvMethodObj
#
# @param pvMethod a character provided by user indicate davies or liu
# @param verbose logical indicating screen printing preferences
# @param ... optional arguments to be provided to p-value method
#
# @return an object of class PvMethod_Davies or PvMethod_Liu
#
# @keywords internal
#
#' @importFrom CompQuadForm liu davies
.newPvMethodObj <- function(pvMethod, verbose, ...) {

  pvMethod <- tolower(x = pvMethod)

  if (substr(x = pvMethod, start = 1L, stop = 1L) == "d") {
    if (verbose) cat("using davies p-value method\n")
    return( new(Class = 'PvMethod_Davies', args = list(...)) ) 
  } else if (substr(x = pvMethod, start = 1L, stop = 1L) == "l") {
    if (verbose) cat("using liu p-value method\n")
    return( new(Class = 'PvMethod_Liu', args = list(...)) ) 
  } else {
    stop( "pvMethod not yet available" )
  }
}

# Calculates Test Statistic
#
# Given influence functions, calculates the test statistic
#
# @name Ts
#
# @param psi an nxd matrix
#
# @return nxn matrix
#
# @keywords internal
.Ts <- function(psi) {
         return( drop(x = crossprod(x = colSums(x = psi)))/nrow(x = psi) )
       }
