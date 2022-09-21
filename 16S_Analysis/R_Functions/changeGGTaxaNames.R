changeGGtaxa <- function(x) {
  # remove the D_1__ etc...
  tax.clean <- data.frame(row.names = row.names(x),
                          Kingdom = str_replace(x[,1], "k__",""),
                          Phylum = str_replace(x[,2], "p__",""),
                          Class = str_replace(x[,3], "c__",""),
                          Order = str_replace(x[,4], "o__",""),
                          Family = str_replace(x[,5], "f__",""),
                          Genus = str_replace(x[,6], "g__",""),
                          Species = str_replace(x[,7], "s__",""),
                          stringsAsFactors = FALSE)
}
