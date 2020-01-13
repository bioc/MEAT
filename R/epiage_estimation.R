#' Estimates age in skeletal muscle from calibrated DNA methylation profiles.
#'
#' \code{epiage_estimation} takes as input a calibrated methylation matrix of
#' beta-values and uses the muscle clock to estimate age.
#'
#' \code{epiage_estimation} estimates age for each sample in the input
#' \code{beta} DNA methylation matrix based on DNA methylation profiles.
#' \code{beta} needs to be a matrix of beta-values that have been calibrated
#' to the gold standard GSE50498 using \code{\link{BMIQcalibration}}.
#' @param beta A matrix with beta values that have been calibrated to
#' the gold standard GSE50498 using \code{\link{BMIQcalibration}},
#' with samples in columns and CpGs in rows.
#' @param pheno An optional data frame, matrix or tibble containing information
#' on the samples such as age, sex, disease status, etc.
#' @param ID_col_name The name of the column in \code{pheno} that contains
#' the sample IDs corresponding to the column names of \code{beta}.
#' @param age_col_name The name of the column in \code{pheno} that contains
#' age (in years).
#'
#' @return A list containing the following items:
#' \enumerate{
#'   \item \code{DNAmage} A numeric vector of predicted age (in years)
#'   for each sample in \code{beta} based on DNA methylation profiles.
#'   \item \code{AAdiff} An optional numeric vector corresponding to the
#'   difference between predicted and actual age (in years) for each sample
#'   in \code{beta}. This is only returned if \code{age_col_name} is not null.
#'   \item \code{AAresid} An optional numeric vector corresponding to
#'   the residuals of linear model (using \code{\link[stats]{lm}}) of
#'   predicted against actual age  for each sample in \code{beta}.
#'   This is only returned if \code{age_col_name} is not null.
#' }
#' @export
#' @import glmnet
#' @importFrom dplyr pull
#' @import tibble
#' @seealso \code{\link[wateRmelon]{BMIQ}} for the original BMIQ algorithm,
#' \url{https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115}
#' for the adapted version of the BMIQ algorithm, and
#' \url{https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115}
#' for the elastic net model of the muscle clock.
#' @examples
#' data("GSE121961", envir = environment())
#' GSE121961_clean <- clean_beta(beta = GSE121961)
#' GSE121961_calibrated <- BMIQcalibration(datM = GSE121961_clean)
#'
#' # Example without specifying phenotypes
#' # (output = predicted age only)
#' GSE121961_epiage <- epiage_estimation(beta = GSE121961_calibrated)
#'
#' # Example with phenotypes
#' # (output = predicted age,
#' # age acceleration difference and age acceleration residual)
#' data("GSE121961_pheno", envir = environment())
#' GSE121961_epiage <- epiage_estimation(
#'   beta = GSE121961_calibrated,
#'   pheno = GSE121961_pheno,
#'   ID_col_name = "ID",
#'   age_col_name = "Age"
#' )
epiage_estimation <- function(beta, pheno = NULL, ID_col_name = NULL, age_col_name = NULL) {
  # If pheno is a matrix or a data frame, convert to tibble
  if (!is_tibble(pheno))
    pheno <- as_tibble(pheno)

  # Load the elastic net model
  elasticnet_model <- NULL
  data("elasticnet_model", envir = environment())
  lambda.glmnet.Training <- 0.025

  # Predict age based on calibrated DNA methylation profile
  DNAmage <- NULL
  DNAmage <- anti.trafo(predict(elasticnet_model, t(beta), type = "response", s = lambda.glmnet.Training))[,
                                                                                                           1]
  if (is.null(age_col_name)) {
    AAdiff <- NULL
    AAresid <- NULL
    return(DNAmage)
  } else {
    Age <- pull(pheno[, age_col_name])
    print(Age)
    AAdiff <- DNAmage - Age
    AAresid <- AAdiff
    AAresid_noNA <- resid(lm(DNAmage ~ Age))
    AAresid[!is.na(AAresid)] <- AAresid_noNA
    return(data.frame(pheno, DNAmage, AAdiff, AAresid))
  }
}

anti.trafo <- function(x, adult.age = 20) {
  ifelse(x < 0, (1 + adult.age) * exp(x) - 1, (1 + adult.age) * x + adult.age)
}
