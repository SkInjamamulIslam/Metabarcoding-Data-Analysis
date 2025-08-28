#install.packages("vegan")
#install.packages("tidyverse")
#devtools::install_github("microbiome/microbiome")
#install.packages("phyloseq")
#devtools::install_github("jbisanz/qiime2R")
#install.packages("remotes")
#remotes::install_github("david-barnett/microViz")

library(vegan)
library(tidyverse)
library(microbiome)
library(phyloseq)
library(qiime2R)
library(microViz)
library(ape)
getwd()
list.files()
# Create a new folder
dir.create("results")

# Copy selected files
file.copy(from = c("Figures", "file2.csv"), to = "results/")

# Change working directory to the new folder
setwd("results")

# Confirm
getwd()  # Shows the current working directory


getwd()

setwd("D:/Biomac_training/Raw_reads/R_analysis") #set directory where you have data
list.files()



# create phyloseq object with qiime2 artifacts
physeq_16S<-qiime2R::qza_to_phyloseq(features="table.qza", #OTU table
                                     tree="rooted-tree.qza", # rooted tree
                                     taxonomy="taxonomy.qza",# taxonomy
                                     metadata = "sample-metadata.tsv") # Metadata

# data in phyloseq object can be examined with different accessor functions
otu_table(physeq_16S)
tax_table(physeq_16S)
sample_data(physeq_16S)


### Defining a sample order (sample_order) as a vector of sample names
sample_order <- c("L1S8", "L1S57", "L1S76", "L1S105", "L2S155", "L2S175", "L2S204",
                  "L2S222", "L3S242", "L3S294", "L3S313", "L3S341", "L3S360",
                  "L5S104", "L5S155", "L5S174", "L5S203", "L5S222", "L1S140", "L1S208", "L1S257", "L2S204",
                  "L1S281", "L2S240", "L2S309", "L2S357", "L2S382", "L3S378",
                  "L4S63", "L4S112", "L4S137", "L5S240", "L6S20", "L6S68", "L6S93")

physeq_16S <- microViz::ps_reorder(physeq_16S, sample_order)
sample_data(physeq_16S)


# Inspect results ----
print(physeq_16S)
head(otu_table(physeq_16S))
head(tax_table(physeq_16S))
head(sample_data(physeq_16S))


# rarefaction curve
tab <- otu_table(physeq_16S)
class(tab)<-"matrix"
tab<-t(tab)
rarecurve(tab,step=50,cex=0.9)
rarecurve(tab,step=50,cex=0.9, label=FALSE)

#{raefaction curve optional
# Define colors for each sample
cols <- rainbow(nrow(tab))  # Generate unique colors for each sample

# Define layout: 1st panel for rarefaction plot, 2nd panel for legend
layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1))  # Adjust widths to separate figure areas

# Plot the rarefaction curve in the first panel
par(mar = c(5, 4, 4, 2))  # Normal margins for main plot
rarecurve(tab, step = 50, cex = 0.9, col = cols, label = FALSE)

# Move to the second panel and create an empty plot for the legend
par(mar = c(0, 0, 0, 0))  # Remove margins for clean legend
plot.new()  # Create a blank plot for the legend
legend("center", legend = rownames(tab), col = cols, lty = 1, cex = 0.8, 
       title = "Samples", bty = "n")
# Reset layout after plotting
layout(1) #}


# row_sums <- rowSums(tab)
colSums(otu_table(physeq_16S))

# remove samples with lowest number of reads from phyloseq object
cols_to_remove<-c("L3S341")
# Create a logical vector indicating which samples to keep
samples_to_keep <- !sample_names(physeq_16S) %in% cols_to_remove
physeq_16S.pruned = prune_samples(samples_to_keep, physeq_16S)

# rarefaction
physeq_16S.pruned.rarefied <- rarefy_even_depth(physeq_16S.pruned, rngseed = 5, replace = FALSE)
colSums(otu_table(physeq_16S.pruned.rarefied))


# alpha rarefaction curve after subsampling
rarecurve(t(data.frame(otu_table(physeq_16S.pruned.rarefied))), step=50, cex=0.9)

#######################
# tax glom function. taxonomic aggregation (grouping)
physeq.species  <- tax_glom(physeq_16S.pruned.rarefied, NArm = FALSE, taxrank = rank_names(physeq_16S.pruned.rarefied)[7])
tax_table(physeq.species)
physeq.genus  <- tax_glom(physeq_16S.pruned.rarefied, taxrank = rank_names(physeq_16S.pruned.rarefied)[6],
                          NArm = FALSE)
tax_table(physeq.genus) ## View species-level taxonomy

otu_table(physeq.species)## View species-level OTU table (counts summed by species)
otu_table(physeq.genus)  

# fill in NAs in tax_table
# this is optional
physeq.genus <- microViz::tax_fix(physeq.genus, sep = "|", unknowns = NA)
tax_table(physeq.genus)
physeq.species <- microViz::tax_fix(physeq.species, sep = "|", unknowns = NA)
tax_table(physeq.species)

# relative abundance
psi_rela.ASV  = phyloseq::transform_sample_counts(physeq_16S.pruned.rarefied, function(x) x / sum(x))
#psi_rela.species  = phyloseq::transform_sample_counts(physeq_16S.pruned.rarefied, function(x) x / sum(x)) 
#psi_rela.genus  = phyloseq::transform_sample_counts(physeq.genus, function(x) x / sum(x)) 


##check the dataset before run the below command
head(sample_data(physeq_16S.pruned.rarefied))

##check the category eg. Treatment if needed change with sample location

# create new column in sample-metadata to compare subject and different bodysite
sample_data(physeq_16S.pruned.rarefied)$Species.Region <- paste(
  sample_data(physeq_16S.pruned.rarefied)$subject,
  sample_data(physeq_16S.pruned.rarefied)$body.site, sep = "|"
)
sample_data (physeq_16S.pruned.rarefied)

#Export the OTU/ASV Table, taxonomy files, and sample metadata
otu_df <- as.data.frame(otu_table(physeq_16S.pruned.rarefied))
write.csv(otu_df, "otu_table.csv", quote = FALSE)

tax_df <- as.data.frame(tax_table(physeq_16S.pruned.rarefied))
write.csv(tax_df, "taxonomy_table.csv", quote = FALSE)

meta_df <- as.data.frame(sample_data(physeq_16S.pruned.rarefied))
write.csv(meta_df, "sample_metadata.csv", quote = FALSE)


####Alpha diversity analysis

##Plot Chao1 richness estimator and Shannon diversity estimator
plot_richness(physeq_16S.pruned.rarefied, measures=c("Chao1", "Shannon"))


##Regroup together samples from the same fraction

plot_richness(physeq_16S.pruned.rarefied, measures=c("Chao1", "Shannon"), x="body.site", color="subject")

##alpha diversity

alpha_meas = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
(p <- plot_richness(physeq_16S.pruned.rarefied, "subject", "body.site", measures=alpha_meas))

##Add a ggplot2 box plot layer to the previous plot

p + geom_boxplot(data=p$data, aes(x=subject, y=value, color=NULL), alpha=0.1)

# Different options are available (vegan, phyloseq, etc.)
# Shannon index using vegan package
count = as.data.frame(otu_table(physeq_16S.pruned.rarefied))
x = t(count)
head(x)

Shannon = vegan::diversity(x)
Shannon
# species number, for eveness calculation
Richness <- vegan::specnumber(x)
Pielou_evenness <- Shannon/log(Richness)
report = cbind(Shannon, Pielou_evenness, Richness) 
head(report)
report
mapping = phyloseq::sample_data(physeq_16S.pruned.rarefied)
index = merge(mapping,report , by="row.names",all=F)
head(index)

# alpha richness indices using phyloseq function
richness <- phyloseq::estimate_richness(physeq_16S.pruned.rarefied)

#############################

data <- index
data_long <- data %>%
  pivot_longer(cols = c(Pielou_evenness, Richness, Shannon), names_to = "Index", values_to = "Value")

#check your data before running
str(data_long)
unique(data_long$body.site)  # Check unique treatment groups


# Create boxplots
#very important to identify actual experimental group

ggplot(data_long, aes(x = body.site, y = Value, fill = body.site)) +  
  geom_boxplot() +  
  facet_wrap(~ Index, scales = "free") +  
  labs(title = "Comparison of Alpha Diversity Across Treatments",  
       x = "Body Sites", y = "Alpha Diversity Value") +  
  theme_minimal() +  
  theme(legend.position = "bottom") # Moves legend below the figure


ggplot(data_long, aes(x = body.site, y = Value, fill = body.site)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.6) +
  facet_wrap(~ Index, scales = "free") +
  labs(title = "Comparison of Alpha Diversity Across body-site", x = "body.site", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##check available variable
colnames(sample_data(psi_rela.ASV))
# beta diversity
# psi_rela.ASV
ord <- ordinate(psi_rela.ASV, "MDS", distance = "unifrac")
plot_ordination(psi_rela.ASV, ord, color = "subject",shape = "reported.antibiotic.usage", label = "body.site", title = "Unifrac distance") +
  geom_point(size = 4) 

ord <- ordinate(psi_rela.ASV, "PCoA", distance = "bray")
plot_ordination(psi_rela.ASV, ord, color = "subject",shape = "reported.antibiotic.usage", label = "body.site", title = "Bray distance") +
  geom_point(size = 4) 


##Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(physeq_16S.pruned.rarefied))
standf = function(x, t=total) round(t * (x / sum(x)))
Microbiome_clean = transform_sample_counts(physeq_16S.pruned.rarefied, standf)

#### fill in NAs in tax_table
# this is optional

# Subset taxa where Species is not NA and not empty
physeq_species_clean <- subset_taxa(
  physeq_16S.pruned.rarefied,
  !is.na(Species) & Species != ""
)

# Create bar plot using only taxa with known species
plot_bar(physeq_species_clean, fill = "Phylum") +
  theme_minimal() +
  labs(title = "Barplot of Abundance by Species", x = "sample-id", y = "Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Check your available taxonomy ranks
rank_names(physeq_16S.pruned.rarefied)

##Basic bar graph based on taxa
plot_bar(physeq_16S.pruned.rarefied, fill = "Species")

##Basic bar graph based on Division
plot_bar(physeq_16S.pruned.rarefied, fill = "Phylum")

#Make the bargraph nicer by removing OTUs boundaries. This is done by adding ggplot2 modifier
plot_bar(physeq_16S.pruned.rarefied, fill = "Family") + 
  geom_bar(aes(color=Phylum, fill=Kingdom), stat="identity", position="stack")


#Regroup together control vs treatment samples (regroup sample category wise)
Microbiome_group <- merge_samples(physeq_16S.pruned.rarefied, "body.site")
plot_bar(Microbiome_group, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")


#Heatmaps##continue
##A basic heatmap using the default parameters

plot_heatmap(physeq_16S.pruned.rarefied, method = "NMDS", distance = "bray")

##It is very very cluttered. It is better to only consider the most abundant OTUs for heatmaps. 
#For example one can only take OTUs that represent at least 20% of reads in at least one sample. Remember we normalized all the sampples to median number of reads (total)
Micro_abund <- filter_taxa(physeq_16S.pruned.rarefied, function(x) sum(x > total*0.2) > 0, TRUE)
Micro_abund
otu_table(Micro_abund)[1:8, 1:5]

plot_heatmap(Micro_abund, method = "NMDS", distance = "bray")

##It is possible to use different distances and different multivaraite methods. 
##For example Jaccard distance and MDS and label OTUs with Class, order by Class. We can also change the Palette (the default palette is a bit ugly…).

plot_heatmap(Micro_abund, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
             taxa.label = "Class", taxa.order = "Phylum", 
             trans=NULL, low="beige", high="red", na.value="beige") ##change the distance 


##Another strategy is to do a heatmap for a specific taxonomy group

plot_heatmap(Micro_abund, method = "NMDS", distance = "bray", 
             taxa.label = "Phylum", taxa.order = "Phylum", 
             low="beige", high="red", na.value="beige")



##Exploratory tree plots
Micro_Proto <- subset_taxa(physeq_16S.pruned.rarefied, Genus=="Bacillus")
plot_tree(Micro_Proto, color="body.site", shape="subject", label.tips="Genus", size="Abundance")
sample_variables(Micro_Proto)
str(Micro_Proto)
View(Micro_Proto)
Micro_Proto

##Ordination
#multivariate analysis based on Bray-Curtis distance and NMDS ordination

Micro.ord <- ordinate(physeq_16S.pruned.rarefied, "NMDS", "bray")

##Plot OTUs
plot_ordination(physeq_16S.pruned.rarefied, Micro.ord, type="taxa", color="Order", shape= "Phylum", 
                title="OTUs")

##A bit confusing, 
#so make it more easy to visualize by breaking according to taxonomic division
plot_ordination(physeq_16S.pruned.rarefied, Micro.ord, type="taxa", color="Class", 
                title="OTUs", label="Phylum") + 
  facet_wrap(~Phylum, 3)


##Now display samples and enlarge the points to make it more easy to read
sample_variables(physeq_16S.pruned.rarefied)

plot_ordination(physeq_16S.pruned.rarefied, Micro.ord, type="samples", color="subject", 
                shape="body.site", title="Samples") + geom_point(size=3)


##Display both samples and OTUs but in 2 different panels
plot_ordination(physeq_16S.pruned.rarefied, Micro.ord, type="split", color="Phylum", 
                shape="body.site", title="biplot", label = "subject") +  
  geom_point(size=3)


##### Principal Coordinates Analysis (PCoA)+unweighted-UniFrac distance

GPUF <- UniFrac(physeq_16S.pruned.rarefied)
GloPa.pcoa = ordinate(physeq_16S.pruned.rarefied, method="PCoA", distance=GPUF) #Calculate the PCoA on this distance matrix, GPUF.

plot_scree(GloPa.pcoa, "Scree plot for Microbiome, UniFrac/PCoA")

(p12 <- plot_ordination(physeq_16S.pruned.rarefied, GloPa.pcoa, "samples", color="body.site") + 
    geom_point(size=5) + geom_path() + scale_colour_hue(guide = FALSE) )

(p13 <- plot_ordination(physeq_16S.pruned.rarefied, GloPa.pcoa, "samples", axes=c(1, 3),
                        color="subject") + geom_line() + geom_point(size=5) )

##non-metric Multi-Dimensional Scaling (NMDS)
# perform NMDS, set to 2 axes
GP.NMDS <- ordinate(physeq_16S.pruned.rarefied, "NMDS", GPUF)
(p <- plot_ordination(physeq_16S.pruned.rarefied, GP.NMDS, "samples", color="body.site") +
    geom_line() + geom_point(size=5) )




##This part is for assignment
## some extra (https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121)
#Reading the table of sequence variants (SVs)

install.packages('gplots')
library("gplots")


SVs<-read_qza("table.qza")

names(SVs)
SVs$data[1:5,1:5]
#We can also look at the unique identifier for this object
SVs$uuid
#We can see the type of artifact:
SVs$type
#We can also get a complete list of the files within the artifact and their sizes
SVs$contents

#We can also print the providence
print_provenance(SVs)

#Reading Metadata
metadata<-read_q2metadata("sample-metadata.tsv")
head(metadata)
#Reading Taxonomy
taxonomy<-read_qza("taxonomy.qza")
head(taxonomy$data)


#Creating a Phyloseq Object

physeq<-qza_to_phyloseq(
  features="table.qza",
  tree="rooted-tree.qza",
  "taxonomy.qza",
  metadata = "sample-metadata.tsv"
)
physeq

#Alpha Diversity Over Time

library(tidyverse)
library(qiime2R)

metadata<-read_q2metadata("sample-metadata.tsv")
shannon<-read_qza("shannon_vector.qza")

# this moves the sample names to a new column that matches the metadata and allows them to be merged

shannon<-shannon$data %>% rownames_to_column("SampleID")

gplots::venn(list(metadata=metadata$SampleID, shannon=shannon$SampleID))


metadata<-
  metadata %>% 
  left_join(shannon)
head(metadata)
metadata <- metadata %>% rename(shannon = shannon_entropy)
metadata %>%
  filter(!is.na(shannon)) %>%
  ggplot(aes(x=`days-since-experiment-start`, y=shannon, color=`body-site`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se) +
  xlab("Days") +
  ylab("Shannon Diversity") +
  theme_q2r() + # try other themes like theme_bw() or theme_classic()
  scale_color_viridis_d(name="Body Site") # use different color scale which is color blind friendly
ggsave("Shannon_by_time.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches

##if there was an effect of antibiotics on diversity.

metadata %>%
  filter(!is.na(shannon)) %>%
  ggplot(aes(x=`reported-antibiotic-usage`, y=shannon, fill=`reported-antibiotic-usage`)) +
  stat_summary(geom="bar", fun.data=mean_se, color="black") + #here black is the outline for the bars
  geom_jitter(shape=21, width=0.2, height=0) +
  coord_cartesian(ylim=c(2,7)) + # adjust y-axis
  facet_grid(~`body-site`) + # create a panel for each body site
  xlab("Antibiotic Usage") +
  ylab("Shannon Diversity") +
  theme_q2r() +
  scale_fill_manual(values=c("cornflowerblue","indianred")) + #specify custom colors
  theme(legend.position="none") #remove the legend as it isn't needed
ggsave("../../../images/Shannon_by_abx.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches

metadata %>%
  filter(!is.na(shannon)) %>%
  ggplot(aes(x=subject, y=shannon, fill=`subject`)) +
  stat_summary(geom="bar", fun.data=mean_se, color="black") + #here black is the outline for the bars
  geom_jitter(shape=21, width=0.2, height=0) +
  coord_cartesian(ylim=c(2,7)) + # adjust y-axis
  facet_grid(~`body-site`) + # create a panel for each body site
  xlab("Antibiotic Usage") +
  ylab("Shannon Diversity") +
  theme_q2r() +
  scale_fill_manual(values=c("cornflowerblue","indianred")) + #specify custom colors
  theme(legend.position="none") #remove the legend as it isn't needed
ggsave("Shannon_by_person.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches


##Plotting PCoA

library(ggplot2)
library(gplots)

# Read in the metadata
metadata <- read_q2metadata("sample-metadata.tsv")

# Read in the PCoA results
uwunifrac <- read_qza("unweighted_unifrac_pcoa_results.qza")

# Read and clean up Shannon diversity vector
shannon <- read_qza("shannon_vector.qza")$data %>%
  rownames_to_column("SampleID") %>%
  rename(shannon = shannon_entropy)

# Step 1: Get intersecting SampleIDs from all sources
common_ids <- Reduce(intersect, list(
  uwunifrac$data$Vectors$SampleID,
  metadata$SampleID,
  shannon$SampleID
))

# Step 2: Filter all datasets to the common sample set
vectors_filtered <- uwunifrac$data$Vectors %>%
  filter(SampleID %in% common_ids) %>%
  select(SampleID, PC1, PC2)

metadata_filtered <- metadata %>%
  filter(SampleID %in% common_ids)

shannon_filtered <- shannon %>%
  filter(SampleID %in% common_ids)

# Step 3: Join filtered data
plot_data <- vectors_filtered %>%
  left_join(metadata_filtered, by = "SampleID") %>%
  left_join(shannon_filtered, by = "SampleID")

# Step 4: Plot
p <- ggplot(plot_data, aes(
  x = PC1,
  y = PC2,
  color = `body-site`,
  shape = `reported-antibiotic-usage`,
  size = shannon
)) +
  geom_point(alpha = 0.5) +
  theme_q2r() +
  scale_shape_manual(values = c(16, 1), name = "Antibiotic Usage") +
  scale_size_continuous(name = "Shannon Diversity") +
  scale_color_discrete(name = "Body Site")

# Step 5: Display and save
print(p)
ggsave("PCoA.pdf", p, height = 4, width = 5, device = "pdf")


##Plotting a Heatmap

library(tidyverse)
library(qiime2R)

metadata<-read_q2metadata("sample-metadata.tsv")
SVs<-read_qza("table.qza")$data
taxonomy<-read_qza("taxonomy.qza")$data

SVs<-apply(SVs, 2, function(x) x/sum(x)*100) #convert to percent

SVsToPlot<-  
  data.frame(MeanAbundance=rowMeans(SVs)) %>% #find the average abundance of a SV
  rownames_to_column("Feature.ID") %>%
  arrange(desc(MeanAbundance)) %>%
  top_n(30, MeanAbundance) %>%
  pull(Feature.ID) #extract only the names from the table

SVs %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  gather(-Feature.ID, key="SampleID", value="Abundance") %>%
  mutate(Feature.ID=if_else(Feature.ID %in% SVsToPlot,  Feature.ID, "Remainder")) %>% #flag features to be collapsed
  group_by(SampleID, Feature.ID) %>%
  summarize(Abundance=sum(Abundance)) %>%
  left_join(metadata) %>%
  mutate(NormAbundance=log10(Abundance+0.01)) %>% # do a log10 transformation after adding a 0.01% pseudocount. Could also add 1 read before transformation to percent
  left_join(taxonomy) %>%
  mutate(Feature=paste(Feature.ID, Taxon)) %>%
  mutate(Feature=gsub("[kpcofgs]__", "", Feature)) %>% # trim out leading text from taxonomy string
  ggplot(aes(x=SampleID, y=Feature, fill=NormAbundance)) +
  geom_tile() +
  facet_grid(~`body-site`, scales="free_x") +
  theme_q2r() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_viridis_c(name="log10(% Abundance)")
ggsave("heatmap.pdf", height=4, width=11, device="pdf") # save a PDF 3 inches by 4 inches

#Plotting Per-Feature Abundances

clr<-apply(log2(SVs+0.5), 2, function(x) x-mean(x))
clr %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  gather(-Feature.ID, key=SampleID, value=CLR) %>%
  filter(Feature.ID=="4b5eeb300368260019c1fbc7a3c718fc") %>%
  left_join(metadata) %>%
  filter(`body-site`=="gut") %>%
  ggplot(aes(x=subject, y=CLR, fill=subject)) +
  stat_summary(geom="bar", color="black") +
  geom_jitter(width=0.2, height=0, shape=21) +
  theme_q2r() +
  theme(legend.position="none")
ggsave("aldexbar.pdf", height=2, width=1.5, device="pdf")




