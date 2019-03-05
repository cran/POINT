# p-value using resampling method
#
# @param eps matrix of perturbations
# @param psi influence function
# @param ev non-zero eigenvalues
# @param pvMethod object of class pvMethod
#
# @return p-value of perturbed test statistic
#
# @keywords internal
pvResampFunc <- function(eps, psi, ev, pvMethod) {

  pvResamp <- apply(X = eps, 
                    MARGIN = 2L, 
                    FUN = function(x){
                            .calcPV(method = pvMethod,
                                    psi = x*psi,
                                    ev = ev)
                          })

  return( pvResamp )
}
