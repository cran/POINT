# Hidden methods
#
# @name point-internal-api
# @keywords internal
#' @import methods
tmp <- function(x){}

# Class \code{Kernel}
#
# Class \code{Kernel} holds kernel parameter and weights for kernel calcs
#
# @name Kernel-class
#
# @slot weight The weight vector.
#
# @keywords internal
setClass(Class = "Kernel",
         slots = c(weight = "ANY"),
         prototype = prototype(weight = NA))

# Calculate the kernel matrix
#
# Given a similarity and a genotype snp matrix calculate appropriate kernel
#
# @param kernel An object inheriting from Kernel
# @param ... Rmc and snp objects
#
# @name calcKernel
#
# @keywords internal
setGeneric(name = ".calcKernel",
           def = function(kernel, ...) { standardGeneric(f = ".calcKernel") })

# @rdname calcKernel
setMethod(f = ".calcKernel",
          signature = c(kernel = "ANY"),
          definition = function(kernel, ...) { stop("not allowed") })

# Calculate the local kernel matrix and its rank factorization
#
# @param kernel object of class Kernel
#
# @return {nS x nC} matrix containing rank factorization of kernel matrix
#
# @keywords internal
#
# @name calcLocalKernel
setGeneric(name = ".calcLocalKernel",
           def = function(kernel, ...) { 
               standardGeneric(f = ".calcLocalKernel") 
             })

# Creates Appropriate Kernel Object
#
# Given user specified kernel and weight, creates appropriate
#   object that inherites from Kernel.
#
# @name newKernelObj
#
# @param kernel a character provided by user indicate linear poly or burden
# @param d numeric order for polynomial kernels or NULL for linear and burden
# @param verbose logical indicating screen printing preferences
# @param weight numeric weights or NULL
# @param ... ignored
#
# @return an object of class LinearKernel, PolyKernel, or BurdenKernel
#
# @keywords internal
#
.newKernelObj <- function(kernel, d = NULL, verbose, weight, ...) {

  kernel <- tolower(x = kernel)

  if (substr(x = kernel, start = 1L, stop = 1L) == "l") {
    if (verbose) cat("using linear kernel\n")
    return( .newLinearKernel(weight = weight) ) 
  } else if (substr(x = kernel, start = 1L, stop = 1L) == "p") {
    if (verbose) cat("using polynomial kernel with d=", d, "\n")
    return( .newPolyKernel(weight = weight, d = d) ) 
  } else if (substr(x = kernel, start = 1L, stop = 1L) == "b") {
    if (verbose) cat("using burden kernel\n")
    return( .newBurdenKernel(weight = weight) ) 
  } else {
    stop( "kernel not yet available" )
  }
}

# @param Rmc vector subset of variant similarity matrix R for the mth marker 
#   and cvth c value (R[,m,cv])
# @param snp nxM genotype snp matrix
# @param ... ignored
# @rdname calcLocalKernel
#' @importFrom rARPACK eigs_sym
setMethod(f = ".calcLocalKernel",
          signature = c(kernel = "Kernel"),
          definition = function(kernel, Rmc, snp, ...) {

              localK <- .calcKernel(kernel = kernel, Rmc = Rmc, snp = snp)
              eigSys <- tryCatch(expr = suppressWarnings(expr = eigs_sym(A = localK,
                                                                         k = ncol(x = snp),
                                                                         which = "LM")),
                                 error = function(e){
                                           print(x = e$message)
                                           stop("unable to obtain eigen decomposition")
                                         })

              tmp <- round(x = 1e10*eigSys$values, digits = 0L)/1e10
              pos <- sum(tmp > 0)

              if (pos < 1L) stop("no positive eigen values")

              return( t(x = sqrt(x = eigSys$values[1L:pos]) * 
                      t(x = eigSys$vectors[, 1L:pos, drop=FALSE]) ) )
  
            })
