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
