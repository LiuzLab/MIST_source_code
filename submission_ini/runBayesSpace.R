library(SingleCellExperiment)
library(BayesSpace)
library(ggplot2)

packages <- c("argparser", "BayesSpace")
# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages], repos = "http://cran.us.r-project.org")
}
# Load packages
invisible(lapply(packages, library, character.only = TRUE))

runBayesSpace <- function(fn, q, platform, d, nhvgs=2000){
  data <- read.csv(fn, row.names = 1, stringsAsFactors = F, check.names = T)
  #nhvgs <- as.integer(0.8 * dim(data)[2])
    
  data <- t(data)
  genes <- row.names(data)
  spots <- colnames(data)
  xs <- c()
  ys <- c()
  for (i in 1:length(spots)){
    x <- strsplit(spots[i], split="x")[[1]][1]
    y <- strsplit(spots[i], split="x")[[1]][2]
    xs <- c(xs, as.numeric(x))
    ys <- c(ys, as.numeric(y))
  }
  sce <- SingleCellExperiment(assays=list(counts=as.matrix(data)),
                              colData=DataFrame(spot=spots, row=xs, col=ys),
                              rowData=DataFrame(gene=genes))
  set.seed(2021)
    
  sce<- spatialPreprocess(sce, platform="ST", n.PCs=d, n.HVGs=nhvgs, log.normalize=TRUE)
  sce <- spatialCluster(sce, q=q, platform=platform, d=d,
                        init.method="mclust", model="t", gamma=2,
                        nrep=1000, burn.in=100,
                        save.chain=TRUE)
    
  df <- as.data.frame(sce@colData)
  return(df)
}
# arguments parser
p <- arg_parser("run BayesSpace")
p <- add_argument(p, "--input", help="project directory", type="character")
p <- add_argument(p, "--output", help="data name", type="character")
p <- add_argument(p, "--q", help="number of clusters", type="integer")
p <- add_argument(p, "--platform", help="platform", type="character")
p <- add_argument(p, "--d", help="number of components - PCA", type="integer")
p <- add_argument(p, "--nhvgs", help="number of highly variable genes", type="integer")
argv <- parse_args(p)

input_fn <- argv$input
output_fn <- argv$output
q <- argv$q
platform <- argv$platform
d <- argv$d
nhvgs <- argv$nhvgs
region_df <- runBayesSpace(input_fn, q=q, platform=platform, d=d, nhvgs=nhvgs)
write.csv(region_df, output_fn)