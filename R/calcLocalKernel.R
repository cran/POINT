# Calculate the local kernel matrix and its rank factorization
#
# @param Rmc vector subset of variant similarity matrix R for the mth marker 
#   and cvth c value (R[,m,cv])
# @param snp nxM genotype snp matrix
# @param kernel object of class Kernel
# @param ... ignored
#
# @return {nS x nC} matrix containing rank factorization of kernel matrix
#
# @keywords internal
#
#' @importFrom rARPACK eigs_sym
calcLocalKernel = function(Rmc, snp, kernel, ...) {

  return( .calcLocalKernel(kernel = kernel, Rmc = Rmc, snp = snp) )

}
