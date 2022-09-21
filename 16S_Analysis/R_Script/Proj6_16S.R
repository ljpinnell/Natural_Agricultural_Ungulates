#load packages
library(phyloseq);library(metagenomeSeq);library(dplyr);library(scales);
library(pairwiseAdonis); library(vegan); library(metagMisc); library(stringr)
library(ggplot2); library(btools); library(randomcoloR); library(cowplot)
library(pairwiseAdonis); library(picante); library(gridExtra); library(grid); library(wrapr)
library(ggalt); library(ggforce); library(concaveman); library(ggdendro);
library(microbiome); library(FSA)

# setwd
setwd("/Users/ljpinnell/Documents/VERO/Project6/16S/trimmed_primers/")

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
data # we have 56 samples, and 28804 ASVs
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
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank()) # looks good

# some QC checks
min(sample_sums(data)) #30,035
max(sample_sums(data)) # 70,348
mean(sample_sums(data)) # 139,761
median(sample_sums(data)) # 148,093
sort(sample_sums(data)) # 30035, 30800, 32583... ...68689,70348; skews high but all good

## labels in case I want caps
matrix.labs <- c(Feces = "FECES", Soil = "SOIL")

## colour palette
four_colour_palette <- c("#638652","#4d5660","#8d391e","#ab6116")

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
  labs(title= "", y= "Observed ASVs", x= "") +
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
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
  geom_point(size = 2.5) +
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
richness_stats
write.csv(richness_stats$p.value, "richness_stats.csv")

shannon_stats <- pairwise.wilcox.test(alpha_div_meta$Shannon, alpha_div_meta$group, p.adjust.method = "BH")
write.csv(shannon_stats$p.value, "shannon_stats.csv")
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

plot_ordination(data.css, gunifrac_all.ord, type = "samples", 
                color = "location_species") +
  theme_bw() + facet_wrap(~matrix, ncol=1, labeller = labeller(matrix = matrix.labs)) +
  geom_point(size = 6) +
  #geom_mark_hull(aes(fill=location_species), concavity = 5, expand = 0, radius = 0, size =0.2, alpha = 0.35) +
  stat_ellipse(geom = "polygon", aes(fill=location_species), level = 0.9, lty =4, lwd = 2, alpha= 0.5) +
  scale_fill_manual(values = four_colour_palette) +
  scale_colour_manual(values = four_colour_palette) +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1), "cm"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.border = element_rect(colour = "black", size = 1.5),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        axis.text = element_text(size = 28, colour = "black"),
        axis.title.y = element_text(size = 42, vjust = 1.75),
        axis.title.x = element_text(size = 42, vjust = -1.5))

#stats
gunifrac.adonis <- pairwise.adonis(gunifrac_all.dist, data.css.df$group, perm = 9999, p.adjust.m = "BH")
gunifrac.adonis
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
write.csv(gunifrac.dendro.data[["labels"]][["label"]], "16S_dendro_sample_order.csv")

# copy and paste from the excel sheet
species_col <- data.frame(species = c("cow","cow","cow","cow","cow","cow","cow","cow",
                                      "elk","elk","elk","elk","elk","elk","elk","bison",
                                      "elk","elk","elk","elk","elk","elk","elk","elk",
                                      "elk","elk","elk","elk","elk","elk","elk","bison",
                                      "elk","elk","bison","bison","elk","elk","bison","elk",
                                      "bison","elk","bison","elk","bison","bison","bison","elk",
                                      "bison","bison","elk","elk","bison","bison","bison","bison"))

matrix_col <- data.frame(matrix = c("feces",	"feces",	"feces",	"feces",	"feces",	"feces",	"feces",	"feces",	"feces",	"feces",	"soil",	"soil",	"soil",	"soil",	"soil",	"feces",	"feces",	"feces",	"feces",	"feces",	"feces",	"feces",	"feces",	"feces",	"feces",	"feces",	"feces",	"feces",	"feces",	"feces",	"soil",	"feces",	"soil",	"soil",	"feces",	"feces",	"soil",	"soil",	"soil",	"soil",	"soil",	"soil",	"soil",	"soil",	"soil",	"soil",	"soil",	"soil",	"soil",	"soil",	"soil",	"soil",	"feces",	"feces",	"feces",	"feces"))

location_col <- data.frame(location = c("farm",	"farm",	"farm",	"farm",	"farm",	"farm",	"farm",	"farm",	"RMNP",	"RMNP",	"RMNP",	"RMNP",	"RMNP",	"RMNP",	"RMNP",	"YNP",	"RMNP",	"RMNP",	"RMNP",	"YNP",	"YNP",	"RMNP",	"RMNP",	"YNP",	"RMNP",	"YNP",	"YNP",	"YNP",	"YNP",	"YNP",	"RMNP",	"YNP",	"RMNP",	"RMNP",	"YNP",	"YNP",	"YNP",	"YNP",	"YNP",	"YNP",	"YNP",	"YNP",	"YNP",	"YNP",	"YNP",	"YNP",	"YNP",	"YNP",	"YNP",	"YNP",	"YNP",	"YNP",	"YNP",	"YNP",	"YNP",	"YNP"))

location_species_col <- data.frame(location_species = c("cow_farm",	"cow_farm",	"cow_farm",	"cow_farm",	"cow_farm",	"cow_farm",	"cow_farm",	"cow_farm",	"elk_RMNP",	"elk_RMNP",	"elk_RMNP",	"elk_RMNP",	"elk_RMNP",	"elk_RMNP",	"elk_RMNP",	"bison_YNP",	"elk_RMNP",	"elk_RMNP",	"elk_RMNP",	"elk_YNP",	"elk_YNP",	"elk_RMNP",	"elk_RMNP",	"elk_YNP",	"elk_RMNP",	"elk_YNP",	"elk_YNP",	"elk_YNP",	"elk_YNP",	"elk_YNP",	"elk_RMNP",	"bison_YNP",	"elk_RMNP",	"elk_RMNP",	"bison_YNP",	"bison_YNP",	"elk_YNP",	"elk_YNP",	"bison_YNP",	"elk_YNP",	"bison_YNP",	"elk_YNP",	"bison_YNP",	"elk_YNP",	"bison_YNP",	"bison_YNP",	"bison_YNP",	"elk_YNP",	"bison_YNP",	"bison_YNP",	"elk_YNP",	"elk_YNP",	"bison_YNP",	"bison_YNP",	"bison_YNP",	"bison_YNP"))

dendro_sample_order <- c("BFe16",	"BFe17",	"BFe23",	"BFe3",	"BFe10",	"BFe8",	"BFe28",	"BFe33",	"REF1",	"REF2",	"RES6",	"RES7",	"RES3",	"RES4",	"RES8",	"YBF8",	"REF4",	"REF7",	"REF8",	"YEF7",	"YEF6",	"REF3",	"REF6",	"YEF8",	"REF5",	"YEF4",	"YEF5",	"YEF1",	"YEF2",	"YEF3",	"RES2",	"YBF2",	"RES1",	"RES5",	"YBF1",	"YBF5",	"YES4",	"YES7",	"YBS7",	"YES8",	"YBS8",	"YES3",	"YBS6",	"YES6",	"YBS3",	"YBS1",	"YBS2",	"YES5",	"YBS4",	"YBS5",	"YES1",	"YES2",	"YBF3",	"YBF4",	"YBF6",	"YBF7")

gunifrac.dendro.data$labels <- cbind(gunifrac.dendro.data$labels, matrix_col)
gunifrac.dendro.data$labels <- cbind(gunifrac.dendro.data$labels, location_col)
gunifrac.dendro.data$labels <- cbind(gunifrac.dendro.data$labels, location_species_col)
gunifrac.dendro.data$labels

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
  scale_fill_manual(values = ggplot2::alpha(c(four_colour_palette), 0.5)) +
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

# barplot for under dendro
ggplot(ra_family_filt_melt, aes(x= Sample, y= Abundance, fill= Family)) +
  theme_bw() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "summary", colour = "black") +
  scale_x_discrete(limits = dendro_sample_order) +
  scale_fill_manual(values = ra_family_filt_palette) +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(size = 1.25, lineend = "square", colour = "black"), 
        axis.ticks.y = element_line(size = 1.5, lineend = "square", colour = "black"),
        axis.ticks.length.y = unit(0.25,"lines"),
        axis.title.y = element_text(size = 42),
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

rel_abund_family_filt <- filter_taxa(rel_abund_family, function(x) mean(x) > .15, T)
rel_abund_family_filt

## want to add in the other (> 0.15%) column to OTU and taxonomy
write.csv(otu_table(rel_abund_family_filt),"rel_abund_family_filt_otus.csv")
write.csv(tax_table(rel_abund_family_filt),"rel_abund_family_filt_taxa.csv")

##insert new tables
new_family_filt_otu <- read.table("rel_abund_family_filt_otus.csv", header = T, sep = ",", row.names = 1)
new_family_filt_tax <- read.table("rel_abund_family_filt_taxa.csv", header = T, sep = ",", row.names = 1)
new_OTU = otu_table(as.matrix(new_family_filt_otu), taxa_are_rows = T)
new_TAX = tax_table(as.matrix(new_family_filt_tax))
new_SD = sample_data(rel_abund_family_filt)

new_ra_family_filt <- merge_phyloseq(new_OTU, new_TAX, new_SD)
new_ra_family_filt

ra_phylum_melt <- psmelt(rel_abund_phyla)
ra_order_melt <- psmelt(rel_abund_order)
ra_family_melt <- psmelt(rel_abund_family)
ra_family_filt_melt <- psmelt(new_ra_family_filt)
length(unique(ra_family_filt_melt))

ra_family_filt_palette <- distinctColorPalette(105)
write.csv(ra_family_filt_palette, "family_filt_palette.csv")
new_family_filt_palette <- c("#E8C09E",	"#C4F7EF",	"#5BE797",	"#986C43",	"#B0B65F",	"#AC70A4",	"#B264ED",	"#F6C86C",	"#EDA9EC",	"#C2F0AA",	"#96EF29",	"#49B7C9",	"#4A86E1",	"#A9E3AE",	"#E7F29F",	"#F038DA",	"#B0F2D1",	"#D7CBD0",	"#4AF47B",	"#8EABEE",	"#BC80EE",	"#B132B4",	"#57EFCB",	"#E5F1B9",	"#9CE0E8",	"#CAC7F3",	"#E8CB36",	"#D27ACE",	"#519DC2",	"#4EA5E5",	"#CBBBA9",	"#E4A02C",	"#EFC3F1",	"#F1D9EE",	"#E8EAEE",	"#D2E581",	"#E5D35A",	"#43F5BA",	"#D17B59",	"#C6F643",	"#858D57",	"#79B37C",	"#C6A8E0",	"#46E7E3",	"#5AD093",	"#ECDB7F",	"#E26849",	"#A19583",	"#E793AB",	"#DB44F2",	"#F3D6C8",	"#5AC5F4",	"#F5433F",	"#ABB1BB",	"#E9A36A",	"#7BB25C",	"#CA96E6",	"#C6F266",	"#5AEC37",	"#4F3985",	"#5D4BD1",	"#E4EC33",	"#9FD1ED",	"#C35562",	"#90E78C",	"#8FF273",	"#E5F16D",	"#E3F2DF",	"#62E0F1",	"#682CF1",	"#7086A5",	"#E1E0B7",	"#F29582",	"#44C657",	"#B73E8C",	"#87D348",	"#EC447B",	"#AEE97F",	"#DFAFC5",	"#F181EA",	"#93C8AF",	"#47BB9C",	"#749690",	"#C4D5F0",	"#A7D181",	"#B6D5C9",	"#81E7D9",	"#90F4B9",	"#8FB8E9",	"#9071CB",	"#DCA49C",	"#D46727",	"#A86B6E",	"#3F8E6D",	"#EDCF93",	"#EE7AB0",	"#E8449F",	"#5B71E1",	"#8C2ACD",	"#9695EE",	"#B1C18E",	"#E15BD1",	"#7D7EB6",	"#A9C242",	"#FAF9F6")

ggplot(ra_family_filt_melt, aes(x= location_species, y= Abundance, fill = Family)) +
  facet_wrap(~matrix, labeller = labeller(matrix = matrix.labs)) +
  theme_bw() +
  labs(y= "Relative Abundance(%)") +
  geom_bar(stat = "summary", colour = "black") +
  scale_fill_manual(values = new_family_filt_palette) +
  scale_y_continuous(expand = c(0.0005,0,0.0005,0)) +
  scale_x_discrete(limits = c("Cow; Farm", "Bison; Yellowstone National Park ","Elk; Yellowstone National Park ","Elk; Rocky Mountain National Park"),
                   labels =c("Cow; SW USA","Bison; YNP","Elk; YNP","Elk; RMNP")) +
  theme(#legend.position = "none",
        panel.border = element_rect(colour = "black", size = 1),
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size =24, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size =24, colour = "black",angle = 45, hjust = 0.95, vjust = 0.95),
        axis.title.y = element_text(size = 28, vjust = 1.75),
        axis.title.x = element_blank(),
        plot.title = element_text(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())



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

### procrustes on feces and soil

### we need some kinda core/pan resistome thing for feces
feces.css <- subset_samples(data.css, matrix=="Feces")
soil.css <- subset_samples(data.css, matrix=="Soil")
core_ra <- transform(feces.css,"compositional")
core_ra_family <- aggregate_taxa(core_ra,"Family")
core_ra_family_filt <- core(core_ra_family, detection= 0.001, prevalence =0)

# ordinate for procrustes within soil and feces
feces16s.dist <- gunifrac(feces.css)
soil16S.dist <- gunifrac(soil.css)

feces16S.ord <- ordinate(feces.css, method = "NMDS", distance = feces16s.dist)
soil16S.ord <- ordinate(soil.css, method = "NMDS", distance = soil16S.dist)

# core with compositionals
prevalences <- seq(0.1,1,0.1)
detections <- 10^seq(log10(1e-3), log10(0.1), length = 30)

p1 <- plot_core(core_ra_family,
                plot.type = "heatmap",
                colours = rev(brewer.pal(10, "Spectral")),
                prevalences = prevalences,
                detections = detections, min.prevalence = .90)
p1 + theme_bw() +
  labs(y= "", x= "Detection Threshold \n (Relative Abundance (%))", title = "") +
  theme(legend.key.size = unit(4, "lines"),
        plot.title = element_text(size = 24),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(size = 0.75, colour = "black"),
        axis.text.y = element_text(colour = "black", size = 26),
        axis.text.x = element_text(colour = "black", size = 24),
        axis.title.x = element_text(size = 40, vjust = -0.75))


### bacteroidetes
bacteroidetes_feces <- subset_taxa(rel_abund_feces, Phylum=="Bacteroidetes")

bacteroidetes_feces_phylum <- tax_glom(bacteroidetes_feces, taxrank = "Phylum", NArm = F) %>%
  psmelt()

write.csv(bacteroidetes_feces_phylum,"bacteroidetes_RA.csv")

### bring in bacter + tetra
bacter_tetra <- read.csv("bacteroidetes_tetracyclineRA.csv")
bacter_tetra
plot(bacter_tetra$bacteroidetes_RA, bacter_tetra$tetracyclines_RA)
## hahahaha nothing


###### RA plots for feces only
rel_abund_feces <- subset_samples(rel_abund, matrix=="Feces")
sum(taxa_sums(rel_abund_feces)==0)
rel_abund_feces <- prune_taxa(taxa_sums(rel_abund_feces) > 0, rel_abund_feces)

ra_feces_phyla <- tax_glom(rel_abund_feces, taxrank = "Phylum", NArm = F) %>%
  psmelt() # 58 phyla
feces_phyla_palette <- distinctColorPalette(58)

ggplot(ra_feces_phyla, aes(x= group,y= Abundance, fill = Phylum)) +
  theme_bw() +
  geom_bar(stat = "summary", colour = "black") +
  scale_fill_manual(values = feces_phyla_palette)

# firmicutes, proteobacteria, acidobacteria, actinobacteria

# firmi
feces_firmicutes <- subset_taxa(rel_abund_feces, Phylum=="Firmicutes")

feces_firmicutes_phyla <- tax_glom(feces_firmicutes, taxrank = "Phylum", NArm = F) %>%
  psmelt()
feces_firmicutes_family <- tax_glom(feces_firmicutes, taxrank = "Family", NArm = F) %>%
  psmelt() # 37 families
feces_firm_palette <- c("#FFBF00","#FBCEB1","#F2D2BD","#FFAC1C","#CD7F32",
                        "#DAA06D","#CC5500","#E97451","#E3963E","#eb6a2b",
                        "#D27D2D","#B87333","#FF7F50","#F88379","#8B4000",
                        "#FAD5A5","#E49B0F","#FFC000","#DAA520","#C04000",
                        "#F4BB44","#FFDEAD","#FF5F1F","#CC7722","#FFA500",
                        "#FAC898","#FFE5B4","#EC5800","#F89880","#FF7518",
                        "#E35335","#FF5F15","#F28C28","#FA8072","#F08000",
                        "#E3735E","#FFAA33")

ggplot(feces_firmicutes_family, aes(x= group, y= Abundance)) +
  theme_bw() + coord_flip() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(aes(fill= Family), stat = "summary", colour = "black") +
  geom_errorbar(feces_firmicutes_phyla, mapping= aes(x= group, y= Abundance), stat = "summary", width = 0.55, size = 0.6) +
  scale_fill_manual(values = feces_firm_palette) +
  scale_x_discrete(limits = c("RMNP Elk feces", "YNP Elk feces","YNP Bison feces","Farm Cow feces"),
                   labels = c("Elk; RMNP","Elk; YNP","Bison; YNP","Cow; SW USA")) +
  scale_y_continuous(expand = c(0.002,0,0.1,0)) +
  theme(#legend.position = "none",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    panel.border = element_rect(size = 1.5, colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title.x = element_text(size = 28, vjust = 0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 28, colour = "black"),
    axis.ticks = element_line(size = 0.9, colour = "black"))

##make a csv so i know which to put in legend
feces_firmicutes_family_no_melt <- tax_glom(feces_firmicutes, taxrank = "Family", NArm = F)
write.csv(otu_table(feces_firmicutes_family_no_melt), "firmi_fececs_otus.csv")  
write.csv(tax_table(feces_firmicutes_family_no_melt), "firmi_feces_taxa.csv")

pairwise.wilcox.test(feces_firmicutes_phyla$Abundance, feces_firmicutes_phyla$group, p.adjust.method = "BH")

## bacteroidetes
feces_bacteroidetes <- subset_taxa(rel_abund_feces, Phylum=="Bacteroidetes")

feces_bacteroidetes_phyla <- tax_glom(feces_bacteroidetes, taxrank = "Phylum", NArm = F) %>%
  psmelt()
feces_bacteroidetes_family <- tax_glom(feces_bacteroidetes, taxrank = "Family", NArm = F) %>%
  psmelt() # 35 families
feces_bacter_palette <- c("#5F9EA0","#00A36C","#000080","#89CFF0","#96DED1",
                          "#F0FFFF","#00FFFF","#0818A8","#7393B3","#088F8F",
                          "#0096FF","#0047AB","#6495ED","#00008B","#6F8FAF",
                          "#6082B6","#1434A4","#7DF9FF","#0818A8","#5D3FD3",
                          "#191970","#CCCCFF","#B6D0E2","#4169E1","#ADD8E6",
                          "#1F51FF","#A7C7E7","#0F52BA","#87CEEB","#008080",
                          "#4682B4","#40E0D0","#0437F2","#40B5AD","#0818A8")


ggplot(feces_bacteroidetes_family, aes(x= group, y= Abundance)) +
  theme_bw() + coord_flip() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(aes(fill= Family), stat = "summary", colour = "black") +
  geom_errorbar(feces_bacteroidetes_phyla, mapping= aes(x= group, y= Abundance), stat = "summary", width = 0.55, size = 0.6) +
  scale_fill_manual(values = feces_bacter_palette) +
  scale_x_discrete(limits = c("RMNP Elk feces", "YNP Elk feces","YNP Bison feces","Farm Cow feces"),
                   labels = c("Elk; RMNP","Elk; YNP","Bison; YNP","Cow; SW USA")) +
  scale_y_continuous(expand = c(0.002,0,0.1,0)) +
  theme(legend.position = "bottom",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    panel.border = element_rect(size = 1.5, colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title.x = element_text(size = 28, vjust = 0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 28, colour = "black"),
    axis.ticks = element_line(size = 0.9, colour = "black"))

pairwise.wilcox.test(feces_bacteroidetes_phyla$Abundance, feces_bacteroidetes_phyla$group, p.adjust.method = "BH")

S247 <- subset_taxa(feces_bacteroidetes, Family=="S24-7")
S247_genus <- tax_glom(S247, taxrank = "Genus", NArm = F) %>%
  psmelt()


ggplot(S247_genus, aes(x= group, y= Abundance)) +
         theme_bw() + coord_flip() +
         labs(y= "Relative Abundance (%)") +
  geom_bar(aes(fill= Genus), stat = "summary", colour = "black", alpha = 0.5) + 
  scale_y_continuous(breaks = c(0,0.5,0.75,1,1.5,5,10))

pairwise.wilcox.test(S247_genus$Abundance, S247_genus$group, p.adjust.method = "BH")
  
## proteobacteria
feces_proteobacteria <- subset_taxa(rel_abund_feces, Phylum=="Proteobacteria")

feces_proteobacteria_phyla <- tax_glom(feces_proteobacteria, taxrank = "Phylum", NArm = F) %>%
  psmelt()
feces_proteobacteria_family <- tax_glom(feces_proteobacteria, taxrank = "Family", NArm = F) %>%
  psmelt() #131 families
feces_proteobacteria_order <- tax_glom(feces_proteobacteria, taxrank = "Order", NArm = F) %>%
  psmelt() # 67 orders

feces_firm_palette <- c("#FFBF00","#FBCEB1","#F2D2BD","#FFAC1C","#CD7F32",
                        "#DAA06D","#CC5500","#E97451","#E3963E","#FF4433",
                        "#D27D2D","#B87333","#FF7F50","#F88379","#8B4000",
                        "#FAD5A5","#E49B0F","#FFC000","#DAA520","#C04000",
                        "#F4BB44","#FFDEAD","#FF5F1F","#CC7722","#FFA500",
                        "#FAC898","#FFE5B4","#EC5800","#F89880","#FF7518",
                        "#E35335","#FF5F15","#F28C28","#FA8072","#F08000",
                        "#E3735E","#FFAA33")

ggplot(feces_proteobacteria_order, aes(x= group, y= Abundance)) +
  theme_bw() + coord_flip() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(aes(fill= Order), stat = "summary", colour = "black") +
  geom_errorbar(feces_proteobacteria_phyla, mapping= aes(x= group, y= Abundance), stat = "summary", width = 0.55, size = 0.6) +
  scale_fill_manual(values = distinctColorPalette(131)) +
  scale_x_discrete(limits = c("RMNP Elk feces", "YNP Elk feces","YNP Bison feces","Farm Cow feces"),
                   labels = c("Elk; RMNP","Elk; YNP","Bison; YNP","Cow; SW USA")) +
  scale_y_continuous(expand = c(0.002,0,0.1,0)) +
  theme(legend.position = "none",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    panel.border = element_rect(size = 1.5, colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title.x = element_text(size = 28, vjust = 0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 28, colour = "black"),
    axis.ticks = element_line(size = 0.9, colour = "black"))

### enterobacteriales
feces_enterobacteriales <- subset_taxa(rel_abund_feces, Order=="Enterobacteriales")

feces_entero_order <- tax_glom(feces_enterobacteriales, taxrank = "Order", NArm = F) %>%
  psmelt() # 21 genera
feces_entero_genus <- tax_glom(feces_enterobacteriales, taxrank = "Genus", NArm = F) %>%
  psmelt() # 21 genera

ggplot(feces_entero_genus, aes(x= group, y= Abundance)) +
  theme_bw() + coord_flip() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(aes(fill= Genus), stat = "summary", colour = "black") +
  geom_errorbar(feces_entero_order, mapping= aes(x= group, y= Abundance), stat = "summary", width = 0.55, size = 0.6) +
  scale_fill_manual(values = distinctColorPalette(21)) +
  scale_x_discrete(limits = c("RMNP Elk feces", "YNP Elk feces","YNP Bison feces","Farm Cow feces"),
                   labels = c("Elk; RMNP","Elk; YNP","Bison; YNP","Cow; SW USA")) +
  scale_y_continuous(expand = c(0.002,0,0.1,0)) +
  theme(#legend.position = "none",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    panel.border = element_rect(size = 1.5, colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title.x = element_text(size = 28, vjust = 0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 28, colour = "black"),
    axis.ticks = element_line(size = 0.9, colour = "black"))
