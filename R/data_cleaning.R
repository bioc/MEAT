#' Clean beta-matrix.
#'
#' \code{clean_beta} reduces \code{beta} to the right CpGs, imputes missing
#' values if any, and replaces 0 and 1 with min and max values.
#'
#' \code{clean_beta} will
#' 1) reduce \code{beta} to the 19,401 CpGs used
#' to calibrate DNA methylation profiles to the gold standard,
#' 2) check whether \code{beta} contains missing values,
#' and impute them with \code{\link[impute]{impute.knn}},
#' 3) check whether \code{beta} contains 0 and 1 values, and if any,
#' change them to the minimum non-0 and maximum non-1 values in \code{beta}.
#'
#' @param beta A matrix or data frame of DNA methylation beta-values,
#' with samples in columns and CpGs in rows.
#' @return A clean version of the input \code{beta} reduced to 19,401 CpGs,
#' with missing values imputed, and without 0 or 1 values.
#' @export
#' @importFrom utils data
#' @importFrom stringr str_detect
#' @seealso \code{\link[impute]{impute.knn}} for imputation of missing values.
#' @examples
#' data("GSE121961", envir = environment())
#' GSE121961_clean <- clean_beta(beta = GSE121961)
clean_beta <- function(beta) {
  # If beta is a matrix, convert to a data frame
  if (is.matrix(beta))
    beta <- as.data.frame(beta)

  # Check beta has row names
  if (!str_detect(rownames(beta)[1], "cg")) {
    stop("Please make sure your beta-matrix
         has row names corresponding to CpGs.")
  }

  # Reduce beta to the 19,401 CpGs used to calibrate the methylation profiles
  message("-------------------Step 1-------------------------------")
  message("Reducing your beta-matrix to the 19,401 CpGs used to calibrate methylation
    profiles in the muscle clock.")
  gold.mean <- NULL
  data("gold.mean", envir = environment())
  CpGs <- as.character(gold.mean[, "CpGs"])
  # Check that the beta-matrix row names contain the CpGs
  L <- length(intersect(CpGs, rownames(beta)))
  message(paste("Your beta-matrix contains", L, "of the 19,401 CpGs needed to calibrate methylation profiles."))

  if (L < 19401 * 0.9) {
    message("Your beta-matrix is missing > 10% CpGs of the 19,401 CpGs
            needed to calibrate methylation profiles.
            Calibration may be off, which may impact the accuracy of
      age estimation.")
  }

  beta2 <- as.matrix(beta[CpGs, ])
  if (nrow(beta2) == 0) {
    stop("Please make sure the rows of your beta-matrix are named
         after the CpG probes (e.g. cg0012545).")
  }

  message("-------------------Step 2-------------------------------")
  message("Checking for missing values in the beta-matrix.")
  if (anyNA(beta2)) {
    L_NA <- length(which(is.na(beta2)))
    message(paste("Your beta-matrix contains", L_NA, "missing values."))
    beta2 <- impute::impute.knn(beta2, k = 5)$data
  }

  message("-------------------Step 3-------------------------------")
  message("Checking for the presence of 0 and 1.")
  L0 <- length(which(beta2 == 0))
  L1 <- length(which(beta2 == 1))
  message(paste("Your beta-matrix contains", L0, "0 values and", L1, "1 values."))
  if (L0 > 0) {
    beta2[beta2 == 0] <- min(beta2[beta2 != 0])
  }
  if (L1 > 0) {
    beta2[beta2 == 1] <- max(beta2[beta2 != 1])
  }
  return(beta2)
}
