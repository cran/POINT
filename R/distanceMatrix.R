# Calculate the Euclidean distance matrix and its standard deviation
#
# @param coord {nM x 3} matrix of coordinates
#
# @return A list containing: \cr
#  dMatrix : {nM x nM} matrix of Euclidean distances \cr
#  sd_dMatrix : standard deviation of distances \cr
#
# @keywords internal
#
# @name distanceMatrix
#
#' @importFrom stats dist sd
distanceMatrix <- function(coord) {

  Dmat <- as.matrix(x = dist(x = coord,
                             method = "euclidean",
                             diag = TRUE,
                             upper = TRUE))

  sd_dMatrix <- sd(x = Dmat)

  return( list("dMatrix" = Dmat, 
               "sd_dMatrix" = sd_dMatrix) )
}

