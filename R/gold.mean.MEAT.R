#' Mean methylation in dataset GSE50498 reduced to the 19,401 CpGs of MEAT
#'
#' Gold standard dataset GSE50498 containing the mean methylation
#' across 24 young and 24 old individuals at the 19,401 CpGs used
#' to calibrate DNA methylation profiles.
#' @format A data frame with 19,401 rows and 2 variables:
#' \describe{
#'   \item{CpGs}{CpG name}
#'   \item{gold.mean}{mean methylation across all samples
#'   at the corresponding CpG (between 0 and 1)}
#'   }
#' @source \url{ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50498/matrix/}
"gold.mean.MEAT"
