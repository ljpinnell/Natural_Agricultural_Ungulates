#load packages
library(phyloseq);library(metagenomeSeq);library(dplyr);library(scales);
library(pairwiseAdonis); library(vegan); library(metagMisc); library(stringr)
library(ggplot2); library(btools); library(randomcoloR); library(cowplot)
library(pairwiseAdonis); library(picante); library(gridExtra); library(grid); library(wrapr)
library(ggalt); library(ggforce); library(concaveman); library(ggdendro)

# setwd
setwd("/Users/ljpinnell/Documents/WTAMU/Project6/16S/phyloseq/")

# source some stuff
source("g_unifrac.R")
source("w_unifrac.R")
source("uw_unifrac.R")
source("changeGGTaxaNames.R")

# import data
qiimedata <- import_biom("table-with-taxonomy.biom", "tree.nwk", "dna-sequences.fasta")
map_file <- import_qiime_sample_data("../../metadata/Project6_metadata.txt") # need to convert the date to date type
map_file$date_collected <- as.Date(map_file$date_collected, format = "%d-%b-%y")
str(map_file)

# combining sample data with the rest
data <- merge_phyloseq(qiimedata,map_file)

## DATA EXPLORATION
data # we have 56 samples, and 36692 ASVs
sum(taxa_sums(data)==0) # no taxa without any counts, good!
sum(sample_sums(data)==0) # no samples without any ASVs, good!

# check the names of our ranks
rank_names(data) # "Rank1" - "Rank7" not ideal, lets change em
colnames(tax_table(data)) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
rank_names(data) # beauty, now they are named properly

# changing the GG style naming (k__Bacteria, etc.)
tax.data <- data.frame(tax_table(data)) # extract the taxonomy table as a data frame
tax.data.names <- changeGGtaxa(tax.data) # this gets rid of the GG format

# now to change the NAs to a better naming scheme
for (i in 1:7){ tax.data.names[,i] <- as.character(tax.data.names[,i])} # converting all columns to characters
tax.data.names[is.na(tax.data.names)] <- "" # replacing the NAs with an empty string

# now filling in the empty slots with the highest assigned taxonomy
for (i in 1:nrow(tax.data.names)){
  if (tax.data.names[i,2] == ""){
    kingdom <- paste("unclassified ", tax.data.names[i,1], sep = "")
    tax.data.names[i, 2:7] <- kingdom
  } else if (tax.data.names[i,3] == ""){
    phylum <- paste("unclassified ", tax.data.names[i,2], sep = "")
    tax.data.names[i, 3:7] <- phylum
  } else if (tax.data.names[i,4] == ""){
    class <- paste("unclassified ", tax.data.names[i,3], sep = "")
    tax.data.names[i, 4:7] <- class
  } else if (tax.data.names[i,5] == ""){
    order <- paste("unclassified ", tax.data.names[i,4], sep = "")
    tax.data.names[i, 5:7] <- order
  } else if (tax.data.names[i,6] == ""){
    family <- paste("unclassified ", tax.data.names[i,5], sep = "")
    tax.data.names[i, 6:7] <- family
  } else if (tax.data.names[i,7] == ""){
    tax.data.names$Species[i] <- paste("unclassified ",tax.data.names$Genus[i], sep = "")
  }
}

head(tax.data.names) # great, no more NAs and no more k__
tax_table(data) <- as.matrix(tax.data.names) # re-insert the taxonomy table into the phyloseq object
tail(tax_table(data), 20) # sweet, lookin good!

## # lets look at the number of reads per sample and the distribution
sample_sum_df <- data.frame(sum = sample_sums(data))
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 1000) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank()) # looks good

# some QC checks
min(sample_sums(data)) #30,035
max(sample_sums(data)) # 70,348
sort(sample_sums(data)) # 30035, 30800, 32583... ...68689,70348; looks great not trimming any

#############################################################################################
##############################         ALPHA DIVERSITY         ##############################
#############################################################################################
#############################################################################################

## data preparation - calculation and adding metadata
alpha_div1 <- estimate_richness(data, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) # calculating most metrics
alpha_div2 <- estimate_pd(data) # calculating Faith's PD

alpha_div <- cbind(alpha_div1, alpha_div2) # combining Faith's and other metrics
alpha_div # looks good, but have some duplicates so lets trim it a bit
alpha_div <- alpha_div[,c(1:5)]
alpha_div
alpha_div.df <- as(sample_data(data), "data.frame") # making into DF for metadata
str(alpha_div.df)
alpha_div_meta <- cbind(alpha_div, alpha_div.df)
alpha_div_meta # great now we've got alpha_div values and metadata

ggplot(alpha_div_meta, aes(x= location_species, y = Observed, fill = location_species, colour = location_species)) + 
  theme_bw() + facet_wrap(~matrix, ncol=2, labeller = labeller(matrix = matrix.labs)) +
  labs(title= "", y= "Richness", x= "") +
  geom_boxplot(alpha = 0.6, size = 1) +
  #geom_point() +
  scale_fill_manual(values = four_colour_palette) +
  scale_colour_manual(values = four_colour_palette) +
  scale_y_continuous(expand = c(0.1,0,0.16,0)) +
  scale_x_discrete(limits = c("Cow; Farm", "Bison; Yellowstone National Park ","Elk; Yellowstone National Park ","Elk; Rocky Mountain National Park")) +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 28, vjust = 1.75),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

ggplot(alpha_div_meta, aes(x= location_species, y = Shannon, fill = location_species, colour = location_species)) + theme_bw() +
  facet_wrap(~matrix, ncol=2, labeller = labeller(matrix = matrix.labs)) +
  labs(title= "", y= "Shannon", x= "") +
  geom_boxplot(alpha = 0.6, size = 1) +
  #geom_point() +
  scale_fill_manual(values = four_colour_palette) +
  scale_colour_manual(values = four_colour_palette) +
  scale_y_continuous(expand = c(0.1,0,0.16,0)) +
  scale_x_discrete(limits = c("Cow; Farm", "Bison; Yellowstone National Park ","Elk; Yellowstone National Park ","Elk; Rocky Mountain National Park")) +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 28, vjust = 1.75),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

ggplot(alpha_div_meta, aes(x= location_species, y = PD, fill = location_species, colour = location_species)) + theme_bw() +
  facet_wrap(~matrix, ncol=1, labeller = labeller(matrix = matrix.labs)) +
  labs(title= "", y= "Faith's PD", x= "") +
  geom_boxplot(alpha = 0.6) +
  #geom_point() +
  scale_fill_manual(values = four_colour_palette) +
  scale_colour_manual(values = four_colour_palette) +
  scale_y_continuous(expand = c(0.1,0,0.16,0)) +
  scale_x_discrete(limits = c("Cow; Farm", "Bison; Yellowstone National Park ","Elk; Yellowstone National Park ","Elk; Rocky Mountain National Park")) +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 28, vjust = 1.75),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
  



richness_stats <- pairwise.wilcox.test(alpha_div_meta$Observed, alpha_div_meta$group, p.adjust.method = "BH")
shannon_stats <- pairwise.wilcox.test(alpha_div_meta$Shannon, alpha_div_meta$group, p.adjust.method = "BH")

#############################################################################################
##############################         BETA-DIVERSITY         ##############################
#############################################################################################
#############################################################################################

data.css <- phyloseq_transform_css(data, log = F) # cumulative sum scaling normalization
data.css.df <- as(sample_data(data.css), "data.frame")

wunifrac_all.dist <- wunifrac(data.css)
gunifrac_all.dist <- gunifrac(data.css)
uwunifrac_all.dist <- uwunifrac(data.css)

wunifrac_all.ord <- ordinate(data.css, method = "NMDS", distance = wunifrac_all.dist)
gunifrac_all.ord <- ordinate(data.css, method = "NMDS", distance = gunifrac_all.dist)
uwunifrac_all.ord <- ordinate(data.css, method = "NMDS", distance = uwunifrac_all.dist) #insufficient data

four_colour_palette <- c("#638652","#4d5660","#8d391e","#ab6116")
matrix.labs <- c(Feces = "FECES", Soil = "SOIL")

plot_ordination(data.css, gunifrac_all.ord, type = "samples", 
                color = "location_species") +
  theme_bw() + facet_wrap(~matrix, ncol=1, labeller = labeller(matrix = matrix.labs)) +
  geom_point(size = 2.5) +
  geom_mark_hull(aes(fill=location_species), concavity = 5, expand = 0, radius = 0, size =0.2, alpha = 0.35) +
  scale_fill_manual(values = four_colour_palette) +
  scale_colour_manual(values = four_colour_palette) +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size =24, colour = "black"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 24, vjust = 1.75),
        axis.title.x = element_text(size = 24, vjust = -1.5))

#stats
gunifrac.adonis <- pairwise.adonis(gunifrac_all.dist, data.css.df$group, perm = 9999, p.adjust.m = "BH")
write.csv(gunifrac.adonis, "16S_gunifrac_adonis.csv")

gunifrac.disper <- betadisper(gunifrac_all.dist, data.css.df$group)
plot(gunifrac.disper)
gunifrac.permdisp <- permutest(gunifrac.disper, permutations = 9999, pairwise = T)
write.csv(gunifrac.permdisp[["pairwise"]][["permuted"]], "16S_gunifrac_permdisp.csv")


## Clustering on gUniFrac
gunifrac.hclust <- hclust(gunifrac_all.dist, method = "ward.D2")
plot(gunifrac.hclust)
gunifrac.dendro <- as.dendrogram(gunifrac.hclust)
gunifrac.dendro.data <- dendro_data(gunifrac.dendro, type = "rectangle")
gunifrac.dendro.data

# want to add a few columns to the dendro.df
## easier to export the data to excel and copy and paste back into df
#write.csv(gunifrac.dendro.data[["labels"]][["label"]], "16S_dendro_sample_order.csv")

# copy and paste from the excel sheet
species_col <- data.frame(species = c("cow","cow","cow","cow","cow","cow","cow","cow",
                                      "elk","elk","elk","elk","elk","elk","elk","elk",
                                      "elk","elk","elk","elk","elk","elk","elk","elk",
                                      "elk","bison","elk","bison","elk","bison","elk","elk",
                                      "bison","bison","elk","elk","elk","bison","elk","elk",
                                      "bison","bison","bison","elk","elk","elk","bison","elk",
                                      "elk","elk","bison","bison","bison","bison","bison","bison"))

matrix_col <- data.frame(matrix = c("feces","feces",	"feces",	"feces",
                                    "feces",	"feces",	"feces",	"feces",
                                    "feces",	"feces",	"feces",	"feces",
                                    "feces",	"feces",	"feces",	"feces",
                                    "feces",	"feces",	"feces",	"feces",
                                    "feces",	"feces",	"feces",	"feces",
                                    "soil",	"feces",	"soil",	"feces",
                                    "soil",	"feces",	"soil",	"soil",	"soil",
                                    "soil",	"soil",	"soil",	"soil",	"soil",
                                    "soil",	"soil",	"feces",	"feces",
                                    "feces",	"soil",	"soil",	"soil",
                                    "feces",	"soil",	"soil",	"soil",
                                    "feces",	"soil",	"soil",	"soil",
                                    "soil",	"soil"))

location_col <- data.frame(location = c("farm","farm","farm","farm","farm","farm","farm","farm",
                                        "RMNP","RMNP","RMNP","YNP","RMNP","RMNP","YNP","YNP",
                                        "YNP","RMNP","YNP","RMNP","YNP","RMNP","YNP","YNP",
                                        "RMNP","YNP","RMNP","YNP","RMNP","YNP","YNP","YNP",
                                        "YNP","YNP","YNP","YNP","YNP","YNP","YNP","RMNP",
                                        "YNP","YNP","YNP","RMNP","RMNP","RMNP","YNP","YNP",
                                        "YNP","RMNP","YNP","YNP","YNP","YNP","YNP","YNP"))

location_species_col <- data.frame(location_species = c("cow_farm","cow_farm","cow_farm","cow_farm",
                                                        "cow_farm","cow_farm",	"cow_farm",	"cow_farm",
                                                        "elk_RMNP",	"elk_RMNP",	"elk_RMNP",	"elk_YNP",
                                                        "elk_RMNP",	"elk_RMNP",	"elk_YNP",	"elk_YNP",
                                                        "elk_YNP",	"elk_RMNP",	"elk_YNP",	"elk_RMNP",
                                                        "elk_YNP",	"elk_RMNP",	"elk_YNP",	"elk_YNP",
                                                        "elk_RMNP",	"bison_YNP",	"elk_RMNP",	"bison_YNP",
                                                        "elk_RMNP",	"bison_YNP",	"elk_YNP",	"elk_YNP",
                                                        "bison_YNP",	"bison_YNP",	"elk_YNP",	"elk_YNP",
                                                        "elk_YNP",	"bison_YNP",	"elk_YNP",	"elk_RMNP",
                                                        "bison_YNP",	"bison_YNP",	"bison_YNP",	"elk_RMNP",
                                                        "elk_RMNP",	"elk_RMNP",	"bison_YNP",	"elk_YNP",
                                                        "elk_YNP",	"elk_RMNP",	"bison_YNP",	"bison_YNP",
                                                        "bison_YNP",	"bison_YNP",	"bison_YNP","bison_YNP"))

gunifrac.dendro.data$labels <- cbind(gunifrac.dendro.data$labels, species_col)
gunifrac.dendro.data$labels <- cbind(gunifrac.dendro.data$labels, matrix_col)
gunifrac.dendro.data$labels <- cbind(gunifrac.dendro.data$labels, location_col)
gunifrac.dendro.data$labels <- cbind(gunifrac.dendro.data$labels, location_species_col)

ggplot(gunifrac.dendro.data$segments) + theme_bw() +
  labs(y= "Ward's Distance") +
  geom_segment(aes(x=x,y=y, xend=xend,yend=yend),
               size = 1, lineend = "round", linejoin = "round") +
  geom_point(data = gunifrac.dendro.data$labels, aes(x,y, shape = location, colour= location_species, fill = location_species), 
             size = c(6,6,6,6,6,6,6,6,
                      8,8,8,8,8,8,8,8,
                      8,8,8,8,8,8,8,8,
                      8,8,8,8,8,8,8,8,
                      8,8,8,8,8,8,8,8,
                      8,8,8,8,8,8,8,8,
                      8,8,8,8,8,8,8,8), position = position_nudge(y= -0.08), stroke = 1.3) +
  scale_y_continuous(expand = c(0.35,0,0.2,0)) +
  scale_colour_manual(values = four_colour_palette) +
  scale_fill_manual(values = alpha(c(four_colour_palette), 0.5)) +
  scale_shape_manual(values = c(23,22,21)) +
  geom_text(data= gunifrac.dendro.data$labels, aes(x,y, label = matrix),
            position = position_nudge(y=-0.2, x=0), size = 8, angle = 90, hjust = 1) +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        panel.border = element_blank(),
        axis.line.y = element_line(size = 1.5, lineend = "square", colour = "black"), 
        axis.ticks.y = element_line(size = 1.5, lineend = "square", colour = "black"),
        axis.ticks.length.y = unit(0.25,"lines"), panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 42, vjust = 3, hjust = 0.59),
        axis.text.y = element_text(size = 28, colour = "black"))

#############################################################################################
#############################         RELATIVE ABUNDANCE        #############################
#############################################################################################
#############################################################################################

# transform to RA
rel_abund <- transform_sample_counts(data.css, function(x) {x/sum(x) * 100})

# glom at taxonomic ranks
rel_abund_kingdom <- tax_glom(rel_abund, taxrank = "Kingdom", NArm = F)
rel_abund_phyla <- tax_glom(rel_abund, taxrank = "Phylum", NArm = F)
rel_abund_class <- tax_glom(rel_abund, taxrank = "Class", NArm = F)
rel_abund_order <- tax_glom(rel_abund, taxrank = "Order", NArm = F)
rel_abund_family <- tax_glom(rel_abund, taxrank = "Family", NArm = F)
rel_abund_genus <- tax_glom(rel_abund, taxrank = "Genus", NArm = F)


### get the proportion of reads assigned to each taxonomic rank
kingdom_otus <- otu_table(rel_abund_kingdom)
write.csv(kingdom_otus, "kingdom_otus.csv")
kingdom_taxa <- tax_table(rel_abund_kingdom)
write.csv(kingdom_taxa, "kingdom_taxa.csv")

phyla_otus <- otu_table(rel_abund_phyla)
write.csv(phyla_otus, "phylum_otus.csv")
phyla_taxa <- tax_table(rel_abund_phyla)
write.csv(phyla_taxa, "phylum_taxa.csv")

class_otus <- otu_table(rel_abund_class)
write.csv(class_otus, "class_otus.csv")
class_taxa <- tax_table(rel_abund_class)
write.csv(class_taxa, "class_taxa.csv")

order_otus <- otu_table(rel_abund_order)
write.csv(order_otus, "order_otus.csv")
order_taxa <- tax_table(rel_abund_order)
write.csv(order_taxa, "order_taxa.csv")

family_otus <- otu_table(rel_abund_family)
write.csv(family_otus, "family_otus.csv")
family_taxa <- tax_table(rel_abund_family)
write.csv(family_taxa, "family_taxa.csv")

genus_otus <- otu_table(rel_abund_genus)
write.csv(genus_otus, "genus_otus.csv")
genus_taxa <- tax_table(rel_abund_genus)
write.csv(genus_taxa, "genus_taxa.csv")

length(unique(rel_abund_genus$Genus)) # 658 genera
genus_palette <- distinctColorPalette(658)

ggplot(rel_abund_genus, aes(x= group, y = Abundance, fill= Genus)) + theme_bw() +
  geom_bar(stat = "summary", colour = "black") +
  scale_fill_manual(values = genus_palette) +
  theme(legend.position = "none")
