### function to calculate generalized UniFrac using GUniFrac package
library(GUniFrac, quietly = T)

gunifrac <- function(x) {
  y <- GUniFrac(as(t(otu_table(x)), "matrix"), phy_tree(x), alpha = c(0.5))$unifracs
  z <- y[, ,"d_0.5"]
  z.d <- as.dist(z)
  
}
