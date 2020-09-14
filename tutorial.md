---
title: "Microbiome Tutorial"
authors: Vithia Gunalan, Carmen Saenz, Mani Arumugam
date: "15/09/2020"
output:
  html_document:
    keep_md: true
---


## Welcome to the Microbiome Tutorial!
Today we will analyse some microbiome data and look for a microbiome signature that distinguishes a cohort of patients between those with colorectal cancer (CRC) and healthy controls.

To run through this tutorial, follow the text and instructions below.
To run the code below, copy and paste it to the **Console** window in R studio and hit *Enter*.
You can use the *up* and *down* arrows on your keyboard to cycle through previous commands.
To view plots in a larger format, use the *Zoom* button in the lower right panel of Rstudio, on the tab marked **Plots**.
The *left* and *right* arrows in that same tab help you to cycle through plots you have already made.

### Load Packages

First step is to load the relevant R packages we will require for this exercise:

```r
library(phyloseq)
library(dplyr)
library(ggplot2)
library(DESeq2)
```

### Load Data

Now that we have our packages loaded, we can load the data files for this exercise.
Data files are prepared for you as *phyloseq* objects.
We load the phyloseq object required for this tutorial:

```r
thomasB2 <- readRDS("thomasB2.rds")
thomasBtree <- read_tree("ThomasAM_2018b.tree")
```

you can see from the code that we also load a tree file (_ThomasAM_2018b.tree_) with the phyloseq object.
This is a phylogenetic tree which shows the relationship between all of the taxa in our data.

The first step will be to merge this phylogenetic tree with the loaded phyloseq object.
This will be needed later when we work with distance measures (in particular Weighted UniFrac):

```r
thomasB2phylo = phyloseq::merge_phyloseq(thomasB2, thomasBtree)
```

### Set Colour Theme

The plots we will generate require the colourspace to be defined:

```r
theme_set(theme_bw())
pal = "Set1"
scale_colour_discrete <-  function(palname=pal, ...){
scale_colour_brewer(palette=palname, ...)
}
scale_fill_discrete <-  function(palname=pal, ...){
scale_fill_brewer(palette=palname, ...)
}
```

### Taxa Pruning

Prune all taxa across samples with 0 counts:

```r
thomasB2prune0 <- phyloseq::prune_taxa(taxa_sums(thomasB2phylo) > 0, thomasB2phylo)
```

### Alpha-Diversity

let's have a look at the alpha diversity between CRC and healthy samples.
Alpha diversity shows us the within sample diversity, i.e. how many taxa there are within each sample.
We can compare this for all CRC samples against all samples in healthy individuals:

```r
phyloseq::plot_richness(thomasB2prune0, measures=c("Observed", "Shannon"), x="disease", color="disease")
```

Are there any differences between the alpha diversity between the CRC and healthy groups?
What does this tell us about the difference between healthy individuals and individuals with CRC?

How many samples do we have in the dataset for each condition?

```r
phyloseq::nsamples(thomasB2prune0)
```

Now that we know how many samples to expect, does the plot accurately reflect the number of samples in the dataset?
Even bioinformaticians sometimes just count the dots.


```r
p = phyloseq::plot_richness(thomasB2prune0, measures=c("Observed", "Shannon"), x="disease", color="disease")
p$layers <- p$layers[-1]
p + ggplot2::geom_jitter(width=0.2)
```

Check the number of samples again now.
Adding "jitter" to the plot helps us examine the distribution of datapoints.
We can also use a violin plot (that shows probablility density) as a nice visualisation of the distribution of datapoints:

```r
p + aes(fill = disease) + geom_violin(trim=FALSE, colour = "black", weight = 1) + geom_jitter(width = 0.2, colour = "black")
```

### Beta-Diversity

Beta-diversity is the diversity in taxa composition for each group where alpha was about examining within-sample diversity
One useful measure of beta-diversity is the Weighted UniFrac distance:

```r
ord = ordinate(thomasB2prune0, "PCoA", "unifrac", weighted=TRUE)
ordplot <- plot_ordination(thomasB2prune0, ord, color="disease", title="PCoA Weighted Unifrac Distance CRC vs Healthy")
ordplot +
stat_ellipse(type = "norm") +
theme_bw()
```

Another useful measure of beta-diversity is Bray-Curtis Dissimilarity:

```r
ordbray = ordinate(thomasB2prune0, "PCoA", "bray")
ordbrayplot <- plot_ordination(thomasB2prune0, ordbray, color="disease", title="PCoA Bray-Curtis Dissimilarity CRC vs Healthy")
ordbrayplot +
stat_ellipse(type = "norm") +
theme_bw()
```
You can use the left and right arrows in the plot window to cycle between this and the Weighted Unifrac plot you made before.
What do these plots show us in terms of the taxonomic diversity between the CRC vs the healthy samples?
How different are the two plots? Look up these dissimilarity measures to learn more.
Is there a consistent difference between the CRC and healthy groups?
What would you expect to see if the 2 groups were completely different? Discuss.

### Differential Abundance

Differential abundance analysis lets us know if there are any taxa which might be significantly different in abundance between the CRC and healthy groups.
In order to perform this analysis, we will use the DESeq2 package.



```r
thomasB2ds = phyloseq_to_deseq2(thomasB2prune0, ~ disease)
```

One of the first things to do here is to set the baseline against which we calculate the log2 fold change of abundances.
This is called 'reveling'.
In this instance, we relevel the groups (CRC and healthy) in order to calculate log2 fold changes of taxa abundance for CRC against the healthy group:

```r
thomasB2ds$disease <- relevel( thomasB2ds$disease, "healthy")
```

Then calculate the log2 fold change:

```r
thomasB2ds = DESeq(thomasB2ds, test="Wald", fitType="parametric")
thomasB2res = results(thomasB2ds, cooksCutoff = FALSE)
alpha = 0.005
sigtab = thomasB2res[which(thomasB2res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(thomasB2prune0)[rownames(sigtab), ], "matrix"))
View(sigtab)
```
Examine the table that comes up in Rstudio.
the log2FoldChange column shows the difference in abundance as log2(CRC/healthy) for all taxa in the dataset.
You can use the arrows for each column header to sort table in different ways.
For example, sort the table by log2FoldChange in descending order. What does this sort show?
There are more species in the dataset that we see here. Why are only a few shown in the table?

Run the code below to get a plot of differentially abundant genera:

```r
library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set2", ...) {
+ scale_fill_brewer(palette = palname, ...)
}
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum, )) + geom_point(size=6) +
theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust=0.5)) + ggtitle("Log2 Fold Change Genus Abundance (CRC vs Healthy)")
```
Data points are sorted by log2FoldChange in descending order and coloured by phylum. Which phyla are highly represented amongst the colorectal cancer biomarkers?

Viewing at the Genus level is one way to identify CRC biomarkers.
To see this at species level:
```{r]}
species_order = sigtab$Species[order(sigtab$log2FoldChange, decreasing = TRUE)]
sigtab$Species = factor(sigtab$Species, levels = species_order, ordered = T)
ggplot(sigtab, aes(x=Species, y=log2FoldChange, color=Phylum)) + geom_point(size=6) +
theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+ ggtitle("Log2 Fold Change Species Abundance (CRC vs Healthy)")
```
