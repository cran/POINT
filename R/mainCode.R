# Calculate the local kernel matrix and its rank factorization
#
# @param Rmc vector subset of variant similarity matrix R for the mth marker 
#   and cvth c value (R[,m,cv])
# @param snp nxM genotype snp matrix
# @param kernel object of class Kernel
# @param nullResult fit results of null model
# @param X matrix of non-genetic covariates
# @param pvMethod object of class pvMethod
# @param eps matrix of perturbations
# @param ... ignored
#
# @return a list containing \cr
#   pvObs : p-value of observed test statistic \cr
#   pvResamp : p-value of perturbed test statistic \cr
#
# @keywords internal
mainCode <- function(Rmc, 
                     snp, 
                     kernel, 
                     nullResult, 
                     X, 
                     pvMethod, 
                     eps, ...) {

  # rank factorization of kernel matrix
  ZZ1 <- calcLocalKernel(Rmc = Rmc,
                         snp = snp,
                         kernel = kernel)

  # influence function
  psi <- .calcPsi(nullResult = nullResult, X = X, ZZ1 = ZZ1)

  # eigen values of test statistic
  ee <- tryCatch(expr = eigen(x = crossprod(x = psi)/nrow(x = snp), 
                              only.values = TRUE),
                 error = function(e){
                           print(x = e$message)
                           stop("unable to obtain eigen decomposition")
                         })

  # non-zero eigen values
  nonzeroev <- ee$value[abs(x = ee$value) > 1e-10]

  # p-value of observed test statistic
  pvObs <- .calcPV(method = pvMethod, 
                   psi = psi, 
                   ev = nonzeroev)

  # perturbed test statistics and p-values
  pvResamp <- pvResampFunc(eps = eps,
                           psi = psi, 
                           ev = nonzeroev, 
                           pvMethod = pvMethod) 

  return( list("pvObs" = pvObs, 
               "pvResamp" = pvResamp) )
       
}
