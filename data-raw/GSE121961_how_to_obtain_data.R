# Code to obtain dataset GSE121961
# if (!file.exists("data-raw/GSE121961_series_matrix.txt")) {
#    library(RCurl)
#    library(GEOquery)
#    download.file(url="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121961/matrix/GSE121961_series_matrix.txt.gz",
#                  destfile="data-raw/GSE121961_series_matrix.txt.gz")
#    gunzip(filename='data-raw/GSE121961_series_matrix.txt.gz',
#           destname='data-raw/GSE121961_series_matrix.txt',remove=T)
# }

# GSE121961 <- read.delim('data-raw/GSE121961_series_matrix.txt',
#                        na.strings = c("null","NA"),
#                        comment.char="!")

# Make CpG IDs the row names, and remove the CpG ID as a column
# rownames(GSE121961)=GSE121961$ID_REF
# GSE121961=GSE121961[,-1]

# Select two profiles to use as examples: one profile whose age is known, and one profile whose age is unknown
# GSE121961 <- GSE121961[,c(1,5)]

# Round the methylation values at 2 digits to save memory
# GSE121961 <- round(GSE121961,digits=2)

# usethis::use_data(GSE121961, overwrite = T)
