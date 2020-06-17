#' Clean beta-matrix.
#'
#' \code{clean_beta} reduces the beta-matrix stored in the input
#' SummarizedExperiment object \code{SE} to the right CpGs, imputes missing
#' values if any, and replaces 0 and 1 with min and max values.
#'
#' \code{clean_beta} will transform the the beta-matrix stored in \code{SE} by:
#' 1) reducing it to the 19,401 CpGs used to calibrate DNA methylation profiles
#' to the gold standard
#' 2) checking whether it contains missing values, and impute them with
#' \code{\link[impute]{impute.knn}},
#' 3) check whether it contains 0 and 1 values, and if any, change them to the
#' minimum non-0 and maximum non-1 values in the beta-matrix.
#'
#' @param SE A \code{\link[SummarizedExperiment]{SummarizedExperiment-class}} object.
#' The "assays" component of \code{SE} should contain a beta-matrix of
#' DNA methylation beta-values called "beta", with samples in columns
#' and CpGs in rows.
#' \code{SE} may optionally contain annotation information on the CpGs
#' stored in "rowData" and sample phenotypes stored in "colData".
#'
#' @return A clean version of the input \code{SE} reduced to 19,401 CpGs,
#' with missing values imputed, and without 0 or 1 values.
#' @export
#' @import SummarizedExperiment
#' @importFrom utils data
#' @importFrom stringr str_detect
#' @seealso \code{\link[impute]{impute.knn}} for imputation of missing values,
#' and \code{\link[SummarizedExperiment]{SummarizedExperiment-class}} for more
#' details on how to create and manipulate SummarizedExperiment objects.
#' @seealso
#' @examples
#' # Load matrix of beta-values of two individuals from dataset GSE121961
#' data("GSE121961", envir = environment())
#' # Load phenotypes of the two individuals from dataset GSE121961
#' data("GSE121961_pheno", envir = environment())
#'
#' # Create a SummarizedExperiment object to coordinate phenotypes and
#' # methylation into one object.
#' library(SummarizedExperiment)
#' GSE121961_SE <- SummarizedExperiment(assays=list(beta=GSE121961),
#' colData=GSE121961_pheno)
#'
#' # Run clean_beta() to clean the beta-matrix
#' GSE121961_SE_clean <- clean_beta(SE = GSE121961_SE)
clean_beta <- function(SE=NULL) {

  # Check whether SE is a SummarizedExperiment object
  if (!is(SE, "SummarizedExperiment"))
    stop("Please make sure SE is a SummarizedExperiment object.")

  # Check that the beta-matrix is indeed called "beta"
  if (names(assays(SE))!="beta")
    stop("Please make sure that the beta-matrix stored in the assays component of SE is called beta.")

  # If beta is a matrix, convert to a data frame
  if (is.matrix(assays(SE)$beta))
    assays(SE)$beta <- as.data.frame(assays(SE)$beta)

  # Check beta has row names
  if (!str_detect(rownames(SE)[1], "cg")) {
    stop("Please make sure your the row names of SE correspond to CpGs.")
  }

  # Reduce beta to the 19,401 CpGs used to calibrate the methylation profiles
  message("-------------------Step 1-------------------------------")
  message("Reducing your beta-matrix to the 19,401 CpGs used to calibrate methylation profiles in the muscle clock.")
  gold.mean <- NULL
  data("gold.mean", envir = environment())
  CpGs <- as.character(gold.mean[, "CpGs"])

  # Check that the beta-matrix row names contain the CpGs
  
  L <- length(intersect(CpGs, rownames(SE)))
  message(paste("Your beta-matrix contains", L, "of the 19401 CpGs needed to calibrate methylation profiles."))

  if (L < 19401 * 0.9) {
    message("Your beta-matrix is missing > 10% of the 19,401 CpGs needed to calibrate the methylation profiles.
            Calibration may be off, which may impact the accuracy of epigenetic age estimation.")
  }
  
  SE2 <- SummarizedExperiment(assays=list(beta=assays(SE)$beta[CpGs, ]),
                              colData=colData(SE)) 

  message("-------------------Step 2-------------------------------")
  message("Checking for missing values in the beta-matrix.")
  if (anyNA(assays(SE2)$beta)) {
    L_NA <- length(which(is.na(assays(SE2)$beta)))
    message(paste("Your beta-matrix contains", L_NA, "missing values."))
    assays(SE2)$beta <- impute::impute.knn(as.matrix(assays(SE2)$beta), k = 5)$data
  }

  message("-------------------Step 3-------------------------------")
  message("Checking for the presence of 0 and 1.")
  L0 <- length(which(assays(SE2)$beta == 0))
  L1 <- length(which(assays(SE2)$beta == 1))
  message(paste("Your beta-matrix contains", L0, "0 values and", L1, "1 values."))
  if (L0 > 0) {
    assays(SE2)$beta[assays(SE2)$beta == 0] <- min(assays(SE2)$beta[assays(SE2)$beta != 0])
  }
  if (L1 > 0) {
    assays(SE2)$beta[assays(SE2)$beta == 1] <- max(assays(SE2)$beta[assays(SE2)$beta != 1])
  }
  return(SE2)
}
