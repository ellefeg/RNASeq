# Differential Cell Abundance

5th October, Laura Grice

# Method 1: EdgeR

Adapted from OSCA book (`https://bioconductor.org/books/release/OSCA/multi-sample-comparisons.html`)

This method uses EdgeR to perform differential cell abdundance

## Step 0: Setup

```{r}
library(Seurat)
library(edgeR)
library(ggplot2)

outdir <- "/path/to/outdir"
df <- readRDS("/path/to/SeuratObj.RDS")
```

## Step 1: Set up our EdgeR object

First, we generate our count matrix of cell types per sample. The easiest way to do that is to extract two columns of metadata from the Seurat object and use `table` to create an indicence matrix. Here, `sample_ident` should be whatever you are calling your samples while `annotations` should be the cell annotations.

```{r}
counts <- as.data.frame.matrix(table(df[[c("annotations", "sample_ident")]]))
```

Now we will extract the metadata which will later be used to create a design matrix. It is very important to ensure that the metadata samples are sorted in the same order as in the counts matrix, or your design matrix will not be correct.

```{r}
meta <- unique(df[[c("sample_ident", "orig_ident", "patient")]])
rownames(meta) <- meta$sample_ident
meta$sample_ident <- NULL
meta <- meta[colnames(counts),] #very important!
meta$disease <- factor(meta$disease, levels = c("Normal", "Disease")) #Put whatever the "reference" condition is first
```

Now we can make our EdgeR object

```{r}
y.ab <- DGEList(specific_counts, samples = meta)
y.ab
```

## Step 2: Filter low abundance labels

This is less important than when we are performing DGE analysis, but if we have very rare cell types we may want to remove them here. By default, `filterByExpr` removes cell types if there are less than 15 cells total or every sample has fewer than 10 cells for that cell type.

```{r}
keep <- filterByExpr(y.ab, group=y.ab$samples$disease)
y.ab <- y.ab[keep,]
summary(keep)
```

OSCA book: _"Unlike DE analyses, we do not perform an additional normalization step with calcNormFactors(). This means that we are only normalizing based on the “library size”, i.e., the total number of cells in each sample. Any changes we detect between conditions will subsequently represent differences in the proportion of cells in each cluster. The motivation behind this decision is discussed in more detail in Section 14.5.3 [of the OSCA book]"._

EdgeR manual: _"In this design, Patient is included in the design matrix to correct for baseline differences between the Patients, but we will not be testing for differential expression between the Patients. The filtering should therefore be based soley Treatment rather than on Patient."_

```{r}
design <- model.matrix(~patient + disease, y.ab$samples)
```

## Step 3: Dispersion

We use the estimateDisp() function to estimate the NB dispersion for each cluster (Figure 14.11). We turn off the trend as we do not have enough points for its stable estimation.

```{r}
y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion) #all values will be the same
plotBCV(y.ab, cex=1)
```

```{r}
fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior) #all values will be the same
summary(fit.ab$df.prior)  #all values will be the same
plotQLDisp(fit.ab, cex=1)
```
 
## Step 4: Run DA analysis

```{r}
res <- glmQLFTest(fit.ab, coef=ncol(design))
summary(decideTests(res))
```

## Step 5: Look at results and save
```{r}
topTags(res, n = 20)
tagtable <- topTags(res, n = 20)
tagtable <- tagtable$table
tagtable$stats <- gtools::stars.pval(tagtable$FDR)
write.table(tagtable, file = paste0(outdir, "DA_statTable.txt"), sep = "\t", quote = FALSE, col.names = NA)
```

# Method 2: Milo

Github: `https://github.com/MarioniLab/miloR`

Bioconductor: `http://www.bioconductor.org/packages/release/bioc/vignettes/miloR/inst/doc/milo_demo.html`

Paper: Dann, E., Henderson, N.C., Teichmann, S.A. et al. Differential abundance testing on single-cell data using k-nearest neighbor graphs. Nat Biotechnol (2021). https://doi.org/10.1038/s41587-021-01033-z

## Step 0: Setup

```{r}
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)

outdir <- "/path/to/outdir"
df <- readRDS("/path/to/SeuratObj.RDS")
```

Convert to SCE object and confirm the UMAP is correctly transferred

```{r}
mySCE <- Seurat::as.SingleCellExperiment(df)
plotUMAP(mySCE)
```

Now we have to convert the object again into a Milo format object (essentially an SCE object with slots for neighbourhoods/KNN graphs)

```{r}
myMilo <- Milo(mySCE)
reducedDim(myMilo, "UMAP") <- reducedDim(mySCE, "UMAP")
myMilo
```

## Step 1: Construct a KNN graph

Here we construct a KNN graph, which gets stored in igraph format in the `graph` slot of our Milo object. We just use the suggested parameters here. Note: The authors state that they are "perfecting the functionality to add a precomputed KNN graph (for example constructed with Seurat or scanpy) to the graph slot using the adjacency matrix". I would have to look deeper into the object format to see if it would be possible/easy to do this manually...

```{r}
myMilo <- buildGraph(myMilo, k = 10, d = 30)
```

## Step 2: Define representative neighbourhoods

_"We define the neighbourhood of a cell, the index, as the group of cells connected by an edge in the KNN graph to the index cell. For efficiency, we don’t test for DA in the neighbourhood of every cell, but we sample as indices a subset of representative cells, using a KNN sampling algorithm used by Gut et al. 2015._

_For sampling you need to define a few parameters:_

* prop: the proportion of cells to randomly sample to start with (usually 0.1 - 0.2 is sufficient)
* k: the k to use for KNN refinement (we recommend using the same k used for KNN graph building)
* d: the number of reduced dimensions to use for KNN refinement (we recommend using the same d used for KNN graph building)
* refined: indicated whether you want to use the sampling refinement algorith, or just pick cells at random. The default and recommended way to go is to use refinement. The only situation in which you might consider using random instead, is if you have batch corrected your data with a graph based correction algorithm, such as BBKNN, but the results of DA testing will be suboptimal

Here we just use the suggested values:

```{r}
myMilo <- makeNhoods(myMilo, prop = 0.1, k = 10, d=30, refined = TRUE)
```

_Once we have defined neighbourhoods, it’s good to take a look at how big the neighbourhoods are (i.e. how many cells form each neighbourhood). This affects the power of DA testing. We can check this out using the plotNhoodSizeHist function. Empirically, we found it’s best to have a distribution peaking between 50 and 100. Otherwise you might consider rerunning makeNhoods increasing k and/or prop_

```{r}
plotNhoodSizeHist(myMilo)
```

This is what I get when I run this:
![image](https://user-images.githubusercontent.com/22607689/135974366-69a8b04c-5d9f-4c32-96d9-90446c0a3390.png)

## Step 3: Count cells in neighbourhoods

_Now we have to count how many cells from each sample are in each neighbourhood. We need to use the cell metadata and specify which column contains the sample information_

```{r}
myMilo <- countCells(myMilo, meta.data = data.frame(colData(myMilo)), samples="sample_ident")
```

_This adds to the Milo object a n \times m matrix, where n is the number of neighbourhoods and m is the number of experimental samples. Values indicate the number of cells from each sample counted in a neighbourhood. This count matrix will be used for DA testing._

```{r}
head(nhoodCounts(myMilo))
```

## Step 4: Differential abundance testing

Now we are all set to test for differential abundance in neighbourhoods. We implement this hypothesis testing in a generalized linear model (GLM) framework, specifically using the Negative Binomial GLM implementation in edgeR.

We first need to think about our experimental design. The design matrix should match samples to a condition of interest. In this case the Condition is the covariate we are going to test for.

Here I just use the same method I used above to make my design, this is simpler than the method in the vignette.

```{r}
meta <- unique(df[[c("sample_ident", "orig_ident")]])
meta <- data.frame(samples = meta$sample_ident,
                   disease = meta$orig_ident,
                   patient = meta$sample_ident)
meta$patient <- gsub("_.*", "", meta$patient)
rownames(meta) <- meta$samples
meta$disease <- factor(meta$disease, levels = c("Normal", "Cancer"))
traj_design <- meta
```

_Milo uses an adaptation of the Spatial FDR correction introduced by cydar, which accounts for the overlap between neighbourhoods. Specifically, each hypothesis test P-value is weighted by the reciprocal of the kth nearest neighbour distance. To use this statistic we first need to store the distances between nearest neighbors in the Milo object._

```{r}
myMilo <- calcNhoodDistance(myMilo, d=30)
```

_Now we can do the test, explicitly defining our experimental design_. As before, I include Patient to correct for "baseline differences between the patients"

```{r}
da_results <- testNhoods(myMilo, design = ~patient + disease, design.df = traj_design)

da_results %>%
  arrange(- SpatialFDR) %>%
  head() 
```

## Step 5:Visualise neighbourhoods displaying DA

_"To visualize DA results relating them to the embedding of single cells, we can build an abstracted graph of neighbourhoods that we can superimpose on the single-cell embedding. Here each node represents a neighbourhood, while edges indicate how many cells two neighbourhoods have in common. Here the layout of nodes is determined by the position of the index cell in the UMAP embedding of all single-cells. The neighbourhoods displaying singificant DA are colored by their log-Fold Change."_

```{r}
myMilo <- buildNhoodGraph(myMilo)

plotUMAP(myMilo) + plotNhoodGraphDA(myMilo, da_results, alpha=0.05) +
  plot_layout(guides="collect")
```


## Step 6: Additional visualisations

See vignette: https://rawcdn.githack.com/MarioniLab/miloR/7c7f906b94a73e62e36e095ddb3e3567b414144e/vignettes/milo_gastrulation.html#5_Finding_markers_of_DA_populations
