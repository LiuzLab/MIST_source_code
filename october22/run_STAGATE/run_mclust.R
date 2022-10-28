
.libPaths(c(.libPaths(), "/home/humble_local_25t/alexw/R/x86_64-pc-linux-gnu-library/4.0"))
library(mclust)

setwd("/home/humble_local2_25t/alexw/MIST_round3/MIST_manuscript_revision/run_STAGATE/")

packages <- c("argparser")
# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages], repos = "http://cran.us.r-project.org")
}
# Load packages
invisible(lapply(packages, library, character.only = TRUE))

run_mclust <- function(filepath, G){
  X = read.csv(filepath, row.names=1)
  BIC <- mclustBIC(X)
  mod1 <- Mclust(X, x = BIC, modelNames = 'EEE', G=G)
  res = data.frame(row.names = names(mod1$classification), 
                 'STAGATE' = mod1$classification)
  return(res)
}
# # arguments parser
p <- arg_parser("run mclust")
p <- add_argument(p, "--inpath", help="filepath", type="character")
p <- add_argument(p, "--outpath", help="filepath", type="character")
p <- add_argument(p, "--g", help="number of clusters", type="integer")

argv <- parse_args(p)
input_fn <- argv$inpath
output_fn <- argv$outpath
g <- argv$g
res <- run_mclust(input_fn, G=1:g)
write.csv(res, output_fn)