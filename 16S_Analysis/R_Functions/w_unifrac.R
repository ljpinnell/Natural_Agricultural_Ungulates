### function to calculate generalized UniFrac using GUniFrac package
library(GUniFrac, quietly = T)

wunifrac <- function(x) {
  y <- GUniFrac(as(t(otu_table(x)), "matrix"), phy_tree(x), alpha = c(1))$unifracs
  z <- y[, ,"d_1"]
  z.d <- as.dist(z)
  
}
