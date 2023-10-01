#' Protein Structure Guided Local Test
#'
#' A rare variant association test that utilizes protein tertiary structure to 
#'   increase signal and to identify likely causal variants. Performs 
#'   structure-guided collapsing, which leads to local tests that borrow 
#'   information from neighboring variants on a protein and that provide 
#'   association information on a variant-specific level.
#'
#' @param yy numeric vector; phenotype values.
#' @param X  numeric matrix; non-genetic covariates.
#' @param snp numeric matrix; genotype snp matrix (count of minor alleles).
#'   Matrix cannot contain missing values.
#' @param proteinCoord numeric matrix; columns correspond to 3 dimensional
#'   coordinates (x,y,z) of each variant in the protein tertiary structure.
#' @param ... optional additional arguments for p-value methods 
#'   CompQuadForm::davies and CompQuadForm::liu.
#' @param trait character; type of phenotype data.
#'   Must be one of \{'gaussian','binomial'\}
#'   quantitative or case control data, respectively.
#' @param cValues numeric vector; c values from which to choose the
#'   optimal neighborhood size for borrowing significant information.
#' @param weighted logical; whether or not to weight the local kernel
#'   test using (non-distance based) weights.
#' @param weight numeric vector (optional) If NULL and weighted is TRUE
#'   (1.0-MAF)^24. Ignored if weighted is FALSE.
#' @param kernel character; type of local kernel to use;
#'   Must be one of \{'burden', 'linear', 'polynomial'\}.
#' @param d numeric; If kernel = 'poly', d is the order of the
#'   polynomial kernel.
#' @param pvMethod character; method of calculating the p-value of each
#'   single marker test for fixed c values.
#'   Must be one of \{'davies', 'liu'\}.
#' @param nperturb numeric, number of perturbations/resamples
#'   (perturbed test statistics) to calculate p-value of minP statistic.
#' @param verbose  logical; generate progress screen prints.
#'
#' @return Returns a matrix the rows of which correspond to individual markers.
#'   Columns correspond to: \cr 
#'   (1) minP statistic; \cr 
#'   (2) local kernel test p-value; \cr
#'   (3) optimal scale value from input cValues; \cr
#'   (4) minor allele frequency; and \cr
#'   (5) single variant score test p-value.
#'
#' @examples
#'
#'   # number of subjects
#'   nsubj <- 1000
#'
#'   # number of markers
#'   nm <- 5
#'
#'   # generate coordinates for proteins
#'   protein <- cbind( stats::rnorm(n = nm, mean = 17.6, sd =  6.6),
#'                     stats::rnorm(n = nm, mean =  1.6, sd = 13.6),
#'                     stats::rnorm(n = nm, mean = 22.9, sd = 10.4) )
#'
#'   # generate snp matrix
#'   snp <- matrix(data = rbinom(n = nsubj*nm, size = 1, p = 0.02), 
#'                 nrow = nsubj, ncol = nm)
#'   colnames(snp) = paste0("m",1:nm)
#'
#'   # generate binmoial response
#'   MAF <- colMeans(x = snp)/2
#'   causal <- numeric(nm)
#'   causal[c(2,4)] <- 1.0
#'   betaG <- 0.4*abs(log10(x = MAF))*causal
#'
#'   #no non-genetic covariates
#'   X <- NULL
#'   mu <- -0.05 + snp %*% betaG
#'
#'   pryy <- exp(mu)/(1+exp(mu))
#'   yy <- sapply(X = pryy, FUN = stats::rbinom, n = 1, size = 1)
#'
#'   res <- point(yy = yy, X = X, snp = snp, proteinCoord = protein,
#'                trait = 'binomial', cValues = c(0.1,0.2),
#'                weighted = TRUE, weight = NULL, kernel = 'linear',
#'                pvMethod = 'liu', nperturb = 100,
#'                verbose = FALSE)
#'
#' @importFrom stats rnorm anova
#' @import Matrix
#' @export

point <- function(yy, 
                  X, 
                  snp, 
                  proteinCoord,
                  ..., 
                  trait = 'binomial',
                  cValues = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
                  weighted = TRUE,
                  weight = NULL,
                  kernel = 'linear',
                  d = NULL,
                  pvMethod = 'davies',
                  nperturb = 1000,
                  verbose = TRUE) {

  # ensure that trait is provided as one of {binomial, gaussian}
  trait <- tolower(x = trait)
  if (!{trait %in% c('binomial', 'gaussian')}) {
    stop("trait must be one of {'binomial', 'gaussian'}")
  }

  # create p-value method object
  pvMethod <- .newPvMethodObj(pvMethod = pvMethod, 
                              verbose = verbose, ...)

  # snp must be a matrix
  if (is.data.frame(x = snp)) snp <- data.matrix(frame = snp)
  if (!is.matrix(x = snp)) snp <- matrix(data = snp, ncol = 1L)

  # X must be a matrix
  if (is.data.frame(x = X)) X <- data.matrix(frame = X)
  if (!is.null(x = X) && !is.matrix(x = X)) X <- matrix(data = X, ncol = 1L)

  # verify that dimensions are appropriate
  nsubj <- nrow(x = snp)
  nmarker <- ncol(x = snp)

  if (nsubj < nmarker) {
    stop("number of subjects must be greater than the number of markers")
  }

  if (any(is.na(x = snp))) {
    stop("remove or impute any genotype information ahead of time")
  }
  
  if (nrow(x = snp) != length(x = yy)) {
   stop("phenotype data and genotype matrix must have the same # of individuals")
  }

  if (!is.null(x = X) && {nrow(x = snp) != nrow(x = X)}) {
    stop("genotype and covariate matrices must have the same # of individuals")
  }

  if (ncol(x = snp) != nrow(x = proteinCoord)) {
    stop("# of snp markers and # of marker protein coordinates must be the same")
  }
    
  if (ncol(x = proteinCoord) != 3L) {
    stop("3 dimensional protein coordinate data required")
  }

  # snp is a sparse matrix - convert to sparseMatrix to improve speed
  snp <- Matrix(data = snp, sparse = TRUE)

  if (verbose) cat(nsubj, "subjects\n", nmarker, "markers\n")

  # minor allele frequency
  MAF <- colMeans(x = snp) / 2.0

  if (verbose) {
    cat("\nminor allele frequency\n")
    print(MAF)
  }

  if (!weighted) {
    weight <- NULL
    if (verbose) cat("\nno weighting of markers\n")
  } else {
    if (is.null(x = weight)) {
      wm <- {1.0 - MAF}^24
      weight <- wm / sum(wm)
    } else if (nmarker != length(x = weight)) {
      stop("if provided, dimension of weight must be the # of markers")
    }
    if (verbose) {
      cat("\nmarker weights\n")
      print(weight)
    } 
  }

  # create kernel object
  kernel <- .newKernelObj(kernel = kernel, 
                          d = d,  
                          weight = weight,
                          verbose = verbose)

  # use protein coordinate information to obtain distance matrix and its
  # standard derivative
  dInfo <- distanceMatrix(coord = proteinCoord)

  # ensure intercept is present in non-genetic covariate matrix
  if (is.null(x = X)) {
    X <- rep(x = 1.0, times = nsubj)
  } else {
    ones <- rep(x = 1.0, times = nsubj)
    hasIntercept <- FALSE
    for( i in 1L:ncol(X)) {
      hasIntercept <- isTRUE(all.equal(target = X[,i], current = ones ))
      if (hasIntercept) break
    }
    if (!hasIntercept) {
      if (verbose) cat("\nintercept added to non-genetic covariate matrix\n")
      X <- cbind(1.0, X)
    }
  }

  # fit NULL model
  nullResult <- .newNullResult(trait = trait, 
                               X = X, 
                               yy = yy, 
                               verbose = verbose)

  if (weighted) {
    wgt <- weight
  } else {
    wgt <- rep(x = 1.0, times = nmarker)
  }

  minPResamp <- numeric(length = nperturb)
  minPObs <- numeric(length = nmarker)
  mincv <- numeric(length = nmarker)
  minPpv <- numeric(length = nmarker)
  scoreSingleSNPpv <- numeric(length = nmarker)

  if (verbose) {
    cat("\nscale values considered\n")
    print(cValues)
  }

  # generate perturbation matrix
  eps <- matrix(data = rnorm(n = nsubj*nperturb, 
                             mean = 0.0,  
                             sd = 1.0),
                nrow = nsubj,
                ncol = nperturb)


  for( m in 1L:nmarker) {

    if (verbose) cat("marker", m, "\n")

    InternalGInternal <- snp[,m]*wgt[m]

    fit <- tryCatch(expr = glm(formula = yy~X-1+InternalGInternal, 
                               data = data.frame(X,yy,InternalGInternal), 
                               family = trait),
                    error = function(e){
                              print(x = e$message)
                              stop("unable to fit glm model")
                            },
                    warning = function(w){
                                print(x = paste("Warning from glm:", w$message))
                              })

    if (verbose) print(x = fit)

    scoreSingleSNPpv[m] <- anova(object = nullResult@fit,
                                 fit,
                                 test = "Rao")[2L,6L]

    if (verbose) {
      cat("single variant score test p-value", scoreSingleSNPpv[m], "\n")
    }

    minPObs[m] <- Inf
    minPResamp[] <- Inf

    if (verbose) cat("pvObs ")
    for (cv in 1L:length(x = cValues)) {
      if (cValues[cv] < 1e-6) {
        rVec <- numeric(length = nmarker)
        rVec[m] <- 1.0
      } else {
        h <- cValues[cv]*dInfo$sd_dMatrix
        rVec <- exp(x = -{dInfo$dMatrix[,m]^2}/{2*{h*h}})

        # ensure that NA's or NaN's are not produced
        if (any(is.na(x = rVec)) || any(is.nan(x = rVec))) {
          stop("Variant similarity matrix contains NA or NaN values")
        }
      }

      result <- mainCode(Rmc = rVec, 
                         snp = snp, 
                         kernel = kernel, 
                         nullResult = nullResult, 
                         X = X, 
                         pvMethod = pvMethod, 
                         eps = eps)

      if (result$pvObs < minPObs[m]) {
        minPObs[m] <- result$pvObs
        mincv[m] <- cValues[cv]
      }
      if (verbose) cat(result$pvObs, " ")
      tstMin <- result$pvResamp < minPResamp
      minPResamp[tstMin] <- result$pvResamp[tstMin]

    } #end of loop through c values
    if (verbose) cat("\n")

    minPpv[m] <- sum(minPResamp < minPObs[m]) / nperturb
    if (verbose) cat("min p-value perturbed resampling", minPpv[m], "\n")
    
  } #end of loop through markers
  
  resultsMat <- cbind(minPObs,minPpv,mincv,MAF,scoreSingleSNPpv)
  resultsMat <- cbind(minPObs,minPpv,mincv,MAF,scoreSingleSNPpv)
  prioritizedMarkers <- resultsMat[order(resultsMat[,"minPpv"]),]
  colnames(prioritizedMarkers) <- c("minP","minPpv","mincv","MAF","scoreSingleSNppv")
  
  return( prioritizedMarkers )
  
}
