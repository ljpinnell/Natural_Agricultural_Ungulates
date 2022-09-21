### function to calculate generalized UniFrac using GUniFrac package
library(GUniFrac, quietly = T)

uwunifrac <- function(x) {
  y <- GUniFrac(as(t(otu_table(x)), "matrix"), phy_tree(x), alpha = c(0.5))$unifracs
  z <- y[, ,"d_UW"]
  z.d <- as.dist(z)
  
}