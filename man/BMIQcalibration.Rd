% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_calibration.R
\name{BMIQcalibration}
\alias{BMIQcalibration}
\title{Calibrate methylation data to a gold standard.}
\usage{
BMIQcalibration(SE, version = "MEAT2.0")
}
\arguments{
\item{SE}{A \code{\link[SummarizedExperiment]{SummarizedExperiment-class}} object.
The "assays" component of \code{SE} should contain a beta-matrix of
DNA methylation beta-values called "beta" that has been cleaned with
\code{\link{clean_beta}}.
\code{SE} may optionally contain annotation information on the CpGs stored
in "rowData" and sample phenotypes stored in "colData".}

\item{version}{A character specifying which version of the epigenetic clock
you would like to use. Dy default, \code{version} is set to "MEAT2.0" for the
second version of the epigenetic clock. If you would like to use the original
version, set \code{version} to "MEAT".}
}
\value{
A calibrated version of the input \code{SE} calibrated to the gold
standard dataset GSE50498.
}
\description{
\code{BMIQcalibration} uses an adapted version of the BMIQ algorithm to
calibrate the beta-matrix stored in the input SummarizedExperiment object
\code{SE} to the gold standard dataset used in the muscle clock (GSE50498).
}
\details{
\code{BMIQcalibration} was created by Steve Horvath,
largely based on the \code{\link[wateRmelon]{BMIQ}} function from
Teschendorff (2013) to adjust for the type-2 bias in Illumina HM450
and HMEPIC arrays. BMIQ stands for beta mixture quantile normalization.
Horvath fixed minor errors in the v_1.2 version of the BMIQ algorithm
and changed the optimization algorithm to make the code more robust.
He used method = "Nelder-Mead" in \code{\link[stats]{optim}} since
the other optimization method sometimes gets stuck. Toward this end,
the function \code{\link[RPMM]{blc}} was replaced by \code{blc2}.
\code{SE} needs to be a SummarizedExperiment object containing a matrix of
beta-values that has been cleaned using \code{\link{clean_beta}}.
Each sample in \code{SE} is iteratively calibrated to the
gold standard values, so the time it takes to run
\code{BMIQcalibration} is directly proportional to the number
of samples in \code{SE}. This step is essential to estimate
epigenetic age with accuracy.
}
\examples{
# Load matrix of beta-values of two individuals from dataset GSE121961
data("GSE121961", envir = environment())
# Load phenotypes of the two individuals from dataset GSE121961
data("GSE121961_pheno", envir = environment())

# Create a SummarizedExperiment object to coordinate phenotypes and
# methylation into one object.
library(SummarizedExperiment)
GSE121961_SE <- SummarizedExperiment(assays=list(beta=GSE121961),
colData=GSE121961_pheno)

# Run clean_beta() to clean the beta-matrix
GSE121961_SE_clean <- clean_beta(SE = GSE121961_SE, version = "MEAT2.0")

# Run BMIQcalibration() to calibrate the clean beta-matrix
GSE121961_SE_calibrated <- BMIQcalibration(SE = GSE121961_SE_clean, version = "MEAT2.0")

}
\seealso{
\code{\link{clean_beta}} to get the DNA methylation matrix ready
for calibration,
\code{\link[wateRmelon]{BMIQ}} for the original BMIQ algorithm and
\url{https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115}
for the original paper describing Horvath's adapted BMIQ algorithm, and
\code{\link[SummarizedExperiment]{SummarizedExperiment-class}} for more
details on how to create and manipulate SummarizedExperiment objects.
}
