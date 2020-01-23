Laura Grice
_adapted from 20191213_HowToScrnaseq.Rmd_

Seurat is a very powerful R package for performing quality control, analysis and exploration of scRNASeq data. You can perform unsupervised clustering and cell type/state discovery, spatial reconstruction of single cell data, and perform data "anchoring/integration" to integrate sequencing data from different sources. This latter function is quite amazing - you can take samples from not only different conditions, but also different labs, different sequencing technologies (scRNASeq, ATACSeq, RNASeq, CELSeq...) or even different species. You can also perform label transfer, where the annotations from one dataset can be projected onto a new unknown dataset of interest. 

According to two clustering benchmarking papers (Freytag et al. 2018 and Duo et al. 2018), Seurat is one of the best clustering algorithms currently available for scRNASeq data. However, Quan is concerned that all of the Seurat output figures look the same, so while it is a good "starter" tool, additional tools and more creative graphics are required to push our publications to the next level.

The Seurat website is very good, it contains many [vignettes](https://satijalab.org/seurat/vignettes.html) to guide you through the Seurat functionalities, and also a useful ["cheat sheet" of Seurat commands](https://satijalab.org/seurat/essential_commands.html).

There are other alternative tools to Seurat; Monocle, by Cole Trapnell's group, is amongst the most popular alternatives. It is not currently (January 2020) possible to (automatically) convert between Seurat3 and Monocle3 objects.

# Seurat objects

A `Seurat` object is essentially the R object that holds all of your experimental input data and the results of all your downstream analyses in Seurat. If you've run Cellranger, you can make a Seurat object out of the `filtered_feature_bc_matrix` folder output by Cellranger, which contains our Cellranger count data but with the non-cellular redundancy removed from it. To load in a dataset, just point Seurat  to the `filtered_feature_bc_matrix` folder produced by cellranger and Seurat will deal with the files inside; you do not need to unzip the individual files.

Because the vast majority of read counts in the data matrix are 0s, the `Read10X` command stores the data as a "sparse matrix", a way of representing all of the necessary information in a more compact form. The `CreateSeuratObject` command makes a class `Seurat` object (see **Seurat 101** above). We also do some very basic pre-filtering, to eliminate any very rarely expressed genes (<4 cells) and any cells with very low expression (<201 genes).

```
df.data <- Read10X(data.dir = "/path/to/filtered_feature_bc_matrix")
df <- CreateSeuratObject(counts = df.data, project = "someProjectName", min.cells = 3, min.features = 200)
```

For each dataset, we see the number of genes (features), the number of cells (samples) and the number of datasets (assays):
```
df
# output:
An object of class Seurat 
14363 features across 12782 samples within 1 assay 
Active assay: RNA (14363 features)
```

Seurat objects are S4 objects in R, so it is a little more difficult to interact with than your standard S3 R objects. You can learn more about `Seurat` objects at the [Seurat Wiki page](https://github.com/satijalab/seurat/wiki)

A `Seurat` object has 10 data slots, each of which contains some information such as "assays", "metadata", etc. You can access this information with the `@` - for instance, `pbmc@metadata`. Some of these slots have their own sub-slots (for instance, `assays` has 6 sub-slots). You can access these sub-slots with the `$`. So, if you want to see the RNA assays you have run, you'd type `pbmc@assays$RNA`.  To make things more complicated, some sub-slots have THEIR OWN sub-sub-slots. These can be accessed with the `@` again. So if you want to see a small selection of counts of your RNA assays, you can use `pbmc@assays$RNA@counts[1:6,1:6]`.

You can also use `[[]]` symbols as a "hack" to access the sub-slots. There are also some special shortcuts for key datatypes. This complicates the system further:

* `pbmc@reductions$pca` = `pbmc[["pca"]] #[[]] workaround`
* `pbmc@assays$RNA` = `GetAssay(pbmc) #shortcut`
* `pbmc@reductions$pca@cell.embeddings %%>% as.tibble()` = `Embeddings(pbmc@reductions$pca %>% as.tibble() #shortcut`
* `pbmc@assays$RNA@counts[1:6,1:6]` = `GetAssayData(pbmc, slot = "counts")[1:6, 1:6] #shortcut`

You can learn more about these shortcuts (and other Seurat commands) at the [Seurat Command List](https://satijalab.org/seurat/essential_commands.html)

But what sorts of information can you find in the Seurat object anyway? 

* The `assays` slot has six sub-slots: `counts` (unnormalised counts), `data` (normalised data matrix), `scale.data` (a scaled matrix), `key` (so you can look up features for a given assay), `var.features` (a vector of variable genes) and `meta.features` (gene-level metadata). 
* The `reductions` slot holds `DimReduc` objects - i.e. dimension reduction results from tSNE, UMAP, PCA. Each of these DimReduc objects has 8 slots: cell embeddings, feature loadings, projected, assay.used, standard deviation, key, jackstraw and misc.

You can add in your own metadata too, such as information about the timepoint or the treatment (this will come in handy if you plan to integrate multiple samples). To do so, you simply run: `df$myNewColumn <- 'someMetadata'` (this would add "someMetadata" to every single cell, but if you had a list of metadata that was the same length as the list of cells, it'd work just the same - for instance, see "Measure 1: mitochondria" below).

Once you integrate 2+ samples, you'll need to tell Seurat whether you want to work on the integrated or unintegrated (raw) data - some analyses are incompatible with one or the other. To do this, run `DefaultAssay(myObject) <- "integrated" # or "RNA"`.

# Pre-processing and QC visualisation

## Measure 1: mitochondria

A common step in QC is to look at is the percentage of mitochondrial genes per cell. A high % mitochondrial genes may indicate technical problems or the biological state of your sample, for a few reasons:

* A high number of lysed cells. As cells burst, their cytoplasmic transcripts spill out and are lost, but the mitochondrial transcripts remain inside the mitochondria and are captured by scRNASeq. Thus, we see a higher percentage of mitochondrial transcripts than normal
* A high number of apoptotic cells, which express mitochondrial genes and export the transcripts to the cytoplasm
* If you see a cell cluster upregulating mitochondrial genes, this cluster may represent dead/dying cells
* A stressful biological state due to your experimental/medical conditions. For instance, diseases with elevated metabolic activity or necrosis may have increased mitochondrial gene expression

If your model species is well-annotated, hopefully the mitochondrial genes have some common naming scheme that allows you to identify them. For instance, all mitochondrial gene IDs in mouse begin with the letters "mt-". We can then calculate our percentage of mitochondrial genes per cell like this:

```
df[["percent.mt"]] <- PercentageFeatureSet(df, pattern = "^mt-")
```
Lysed cells may lose a higher proportion of nuclear genes compared to mitochondrial genes. We might then expect that cells with a high percentage of mitochondrial genes may have a low number of total expressed genes. We can visualise this by plotting percent.mt vs nFeature_RNA. This may help you determine your percentage threshold for filtering cells with an overabundance of mitochondrial genes.


## Measure 2: Number of expressed genes or number of mRNA sequences

A common step in scRNASeq pipelines is to identify cells that have extreme (low or high) numbers of expressed genes. Low gene counts may indicate low-quality cells (e.g. cell lysing) or empty GEM droplets that did not take up a cell, and high gene counts may indicate that multiple cells were taken up by a single GEM droplet ("cell doublets" or "multiplets"). 

(Note that we already did some filtering for low gene counts in our `CreateSeuratObject` step - we removed genes expressed in fewer than 3 cells, and cells expressing fewer than 200 genes.)

A common heuristic is to filter cells with a high `nFeature_RNA` or `nCount_RNA` value. However, it can be challenging to distinguish between doublets (high nCount_RNA value due to the presence of two cells) vs large or mRNA-rich cell types. Thus, it may be of interest to run tools which detect doublets or empty cells more directly. 

Doublet detection tools include:

* demuxlet (shell)
* DoubletFinder (R)
* DoubletDecon (R)
* DoubletDetection (R, Python)
* Scrublet (Python)
* scds
* scran's DoubletCluster/DoubletCell (R)
* Solo

Empty droplet detection tools include:

* soupX

Note that depending on the results of this detection, you may choose to take a combinatorial approach to filtering. For instance, you might remove cells with a postive doublet call AND an nFeature_RNA value over X.

# Seurat Normalisation

First we normalise the raw read counts to acount for differences in sequencing depth per cell **for each sample** (hence, we perform this step prior to anchoring). There are two main normalisation strategies that you can use in Seurat, `NormalizeData` or `sctransform`. SCTransform is a newer normalisation strategy but according to the [Seurat FAQ #4](https://satijalab.org/seurat/faq), "DE and integration are not currently supported with sctransform but will be soon". For this reason, I will use the `NormalizeData` procedure here. However, see _**Appendix 3: scTransform as a normalisation strategy**_ for more information about `sctransform`, including a "beta" method to perform integration with sctransform.

`NormalizeData` is the strategy used in the base vignette. By default it uses the "global scaling normalisation" strategy, `LogNormalize`. It (1) normalises gene expression measurements for each cell by total expression, (2) multiplies it by a scale factor (10k by default), and (3) log-transforms the result. Normalized values are stored in `pbmc[["RNA"]]@data`. There are other normalisation strategies available in `NormalizeData` (CLR = centered log ratio transformation, RC = relative counts). 

`RC` is basically the same method as `LogNormalize` but the values aren't normalised:

* logNormalize: ln(1 + 10000 x [cell count / total count])
* RC: 10000 x [cell count / total count] (or you can change 10k to 1 million for CPM)

# Find variable features

Here we find gene subsets that exhibit high variability across cells (highly expressed in some cells, lowly expressed in others). These genes are heterogenous features that are of interest in downstream analysis. To account for the "mean-variance relationship inherent in single-cell data" (See _*Appendix 4: The mean-variance relationship in RNASeq and scRNASeq*_), the data is first transformed; the authors use `vst` (variance-stabilising transformation) but alternatively, you can use `mvp` (mean.var.plot) or `disp` (dispersion - genes with highest dispersions).

`vst` fits a curve to the relationship of log(variance) and log(mean). Features are standardised with observed mean and expected variance (given fitted line). These standardised values are used to compute variance if these values, across all cells. Then the genes are ranked based on standardised variance. By default, the top 2000 genes are selected for downstream PCA/clustering.

# Scale the data

Now we linear transform (scale) our data. Scaling can be slow, but we only need to scale the genes used for PCA input (i.e. the 2000 variable features). `ScaleData` alters the mean expression and variance so that all genes have equal weight and highly-expressed genes don't dominate:

* Shifts expression of each gene so mean expression across cells is 0
* Scales expression of each gene, so variance across cells is 1

To learn more about the general stats of scaling the mean and the variance, see Ming Tang's explanation in "[Scaling the data](https://crazyhottommy.github.io/scRNA-seq-workshop-Fall-2019/scRNAseq_workshop_1.html)".

If you incorporate the `vars.to.regress` flag, you can use ScaleData to regress out unwanted or uninteresting sources of variation (e.g. technical noise, batch effect, cell cycle stage, % mitochondrial expression, etc...). _"As suggested in Buettner et al, NBT, 2015, regressing these signals out of the analysis can improve downstream dimensionality reduction and clustering"_. _"Seurat constructs linear models to predict gene expression based on user-defined variables. The scaled z-scored residuals of these models are stored in the scale.data slot, and are used for dimensionality reduction and clustering."_ However, the authors [recommend](https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html) that if you want to do this regression, you could instead use SCTransform.  

# PCA


Principal component analysis (PCA) is a method of linear dimension reduction. If you want to learn more about how PCA works, you can watch a [short](https://www.youtube.com/watch?v=HMOI_lkzW08) or [long](https://www.youtube.com/watch?v=_UVHneBUBW0) video from StatQuest, or read [Ming Tang's blog post](https://divingintogeneticsandgenomics.rbind.io/post/pca-in-action/) about PCA.

"Also, Seurat uses irlba::irlba() to compute PCA (links to paper and package vignette). Typical output difference between prcomp and irlba with my data is that irlba results in PCs with lower cell embedding values (and a flip in the sign for some PCs) compared to prcomp." [ref](https://github.com/satijalab/seurat/issues/751)

# Cluster the cells

NB: This section adapted from [Ming Tang](https://crazyhottommy.github.io/scRNA-seq-workshop-Fall-2019/scRNAseq_workshop_1.html#cluster_the_cells)

Seurat implements graph-based clustering. Clustering embeds cells in a graph structure (e..g a k-nearest neighbour - KNN - graph); edges are drawn between cells iwth similar feature expression patterns. Then the graph is partitioned into highly inter-connected "quasi-cliques" or "communities". The KNN graph is built based on euclidian distance in PCA space (hence you specify `reduction = "pca"` in `FindNeighbors`). Edge weights between cells are refined based on shared overlap in their local neighbourhood (Jaccard similarty. Then we cluster cells. Cells are iteratively grouped together. When we run `FindClusters`, we need to set the resolution value. This is the granularity of the clusters - the higher the value, the more clusters you get. The bigger the dataset, the more resolution you probably want to use. The Seurat authors recommend a range of 0.4 - 1.2 for ~3000 cell datasets. (For reference, I have 1223 cells in `treat`, 1221 cells in `cntl` and 2444 cells in `dpo.1.combined`)

# What clustering resolution to use?

[Biogitte asks](https://github.com/satijalab/seurat/issues/1565), _"Are there functions in Seurat 3 where it is possible to compare the different cluster resolutions? Or a pipeline?"_ User mojaveazure suggests a tool called `clustree` which can work with Seurat 3 objects, and tilofreiwalk suggests `IKAP`.

Let's have a look at clustree. From the [vignette](https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html): 

>_"Statistics designed to help you (choose the clustering parameter - here, resolution) typically either compare two clusterings or score a single clustering. A clustering tree is different in that it visualises the relationships between at a range of resolutions."_

>_"To build a clustering tree we need to look at how cells move as the clustering resolution is increased. Each cluster forms a node in the tree and edges are constructed by considering the cells in a cluster at a lower resolution (say k=2) that end up in a cluster at the next highest resolution (say k=3). By connecting clusters in this way we can see how clusters are related to each other, which are clearly distinct and which are unstable. Extra information about the cells in each node can also be overlaid in order to help make the decision about which resolution to use. For more information about clustering trees please refer to our associated publication (Zappia and Oshlack 2018)."_

## Identify conserved cell markers

Let's look for canonical cell type marker genes that are conserved across conditions. The `FindConservedMarkers` function performs DE testing for each dataset/group and combines the p-values using the MetaDE R package. The purpose is to find genes that are markers in all conditions. This amounts to splitting the dataset based on condition, running separate DE tests, then combining the results. Thus, in the output you get [conditionX]_p_val and [conditionY]_p_val and so on. You also get `max_pval` (the maximum p-value across all groups) and `minimump_p_val` is the combined p-value, which prepresents a meta-analysis of significance values (combining p-values across different tests) ([Andrew Butler](https://github.com/satijalab/seurat/issues/1164))

`FindMarkers` is used with the RNA assay and `FindConservedMarkers` is used with the integrated set (it is like performing `FindMarkers` for each dataset separately and calculating a combined p-value). `FindConservedMarkers` is run on a single cluster at a time but see [this link](https://hbctraining.github.io/scRNA-seq/lessons/sc_exercises_integ_marker_identification.html) for how to loop through it.

## Cell annotation

I wanted to run monocle3/garnett with my existing seurat object, but while Monocle3 remains in Beta, [Seurat isn't supporting](https://github.com/satijalab/seurat/issues/1658#issuecomment-502763190) converting SEurat objects to monocle objects. If I can find another tool that can directly be loaded with my Seurat object, that'd be good. Here's a command to load across data but I think it's to an early Monocle verison: https://github.com/cole-trapnell-lab/monocle-release/issues/311

I'm getting problems becuase you can run Monocle2 (and thus convert Seurat --> Monocle2 object called a `CellDataSet`) but Garnett currently loads Monocle3 as a dependency, and Monocle3/Garnett can't read the `CellDataSet` object. It wants a `cell_data_set` object.

## de analysis
https://satijalab.org/seurat/v3.1/de_vignette.html


not compatible with integration

See appendix 3 for a note about why DE is of limited interest in scRNASeq.


# UMAP

PCA is a linear dimension reduction technique. UMAP and tSNE are two non-linear dimension reduction techniques. I don't focus on tSNE here, but see _**Appendix 5: tSNE**_ for some more information (in particular about some of the problems with tSNE and _why you shouldn't cluster on tSNE-separated data!_). UMAP is underlaid by complicated maths, but to learn more you can [watch this video by the UMAP creator](https://www.youtube.com/watch?v=nq6iPZVUxZU) or look at Laura's (hard copy) notebook #1, page 183.


# Data integration

## Estimating the dimensionality of the individual datasets in preparation for merging

In the anchoring step (Part Five), you need to tell Seurat how many CCA (canonical correlation analysis) dimensions to use when specifying the neighbour search space (more details in Part Five). The Seurat team [recommends](https://github.com/satijalab/seurat/issues/1248) _"using whatever dimensionality you would use in a standard clustering analysis (i.e. from PCA, for example by exploring an elbow plot, or a statistical resampling analysis)"_. 

How precisely to choose this dimensionality value is not entirely clear. In [Stuart et al. 2019](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8), they state that _"a robust fully unsupervised procedure to identify this value remains a fundamental challenge in the analysis of high-dimensional data"_ but that _"larger datasets will typically have increased dimensionality, particularly if they represent increasingly heterogeneous populations"_. They don't worry too much about the problem and instead _"neglect to finely tune this parameter for each dataset"_, because they _"still observe robust performance over diverse use cases"_. They choose particular values that work, mostly based on increases in dataset size and heterogeneity:

* Neuronal, bipolar and pancreatic analyses: dimensionality = 30
* scATAC-seq in mouse cortex: dimensionality = 20 
* Human bone marrow: dimensionality = 50
* Integration of mouse cell atlas: dimensionality = 100

We will run PCA later on our integrated data, but for now let's have a quick look at the dimensionality of each sample using PCA, to help us choose our parameters for anchoring. We will use the 2000 variable genes we just calculated in the `FindVariableFeatures()` step. We need to work with scaled data for our PCA, but we haven't run this yet because we need to do it after anchoring. So we'll make a temporary object we can scale. Also, below we will explore how to run our own PCA analysis on the data using base R. However, for now we will just use Seurat's `RunPCA()` command.

## Merge the data

Seurat has a method called "anchoring" or "integration", to allow the analysis of multiple scRNASeq samples. _"These methods aim to identify shared cell states that are present across different datasets, even if they were collected from different individuals, experimental conditions, technologies, or even species."_ [ref](https://satijalab.org/seurat/v3.1/integration.html). The method determines pairwise correspondances between individual pairs of cells in the database, which are thought to originate from the same biological state (i.e. represent the same cell type/state). Alternatively, you can use a well-annotated reference transcriptome to transfer labels (metadata, e.g. cell type annotations) to your unknown dataset of interest.

The anchoring step requires you to provide an estimation of the "dimensionality" of the dataset (default = 30) - _"Which dimensions to use from the CCA to specify the neighbor search space"_ (`?FindIntegrationAnchors()`). This is why we did some preliminary PCA analyses in Part Four.

After we run `IntegrateData`, the Seurat object contains a new `assay` which has the integrated/batch-corrected expression matrix. You can switch the default assay back and forth between `integrated` and `RNA` (the raw/uncorrected expression values).

## Appendix 3: scTransform as a normalisation strategy

_My question 1: with LogNormalize, you run NormalizeData and FindVariableFeatures before anchoring and ScaleData after anchoring. What order do you do SCTransform and anchoring?_

_My question 2: SCTransform can "remove confounding sources of variation, for example, mitochondrial mapping percentage" - does this mean you don't filter on %mt before SCTransform?_

Seurat's new normalisation strategy is to run `SCTransform`, which is described in the [sctransform preprint](https://www.biorxiv.org/content/10.1101/576827v1), [the sctransform Github](https://github.com/ChristophH/sctransform), and [a Seurat vignette](https://satijalab.org/seurat/v3.0/sctransform_vignette.html). You can also use sctransform to remove confounding sources of variation, such as mitochondrial mapping percentage. If you use `SCTransform`, you do NOT need to run `NormalizeData`, `ScaleData` or `FindVariableFeatures`. Instead of using log normalised expression values, they "correct/harmonise Pearson residuals that are output from SCTransform". **SCTransform and DE analysis/integration are not currently supported in the main version of Seurat.**

**Data integration (anchoring) and SCTransform**

Although SCTransform and sample integration are not technically compatible in the current version of Seurat, there are ways around it. For instance, Github user [KoichiHashikawa](https://github.com/ChristophH/sctransform/issues/4) proposes to run SCTransform, integrate the data, and then run ScaleData on the data. According to the satijalab admin, "This is a reasonable approach, and uses the Pearson residuals to find anchors across datasets. So your exact code is valid, and for users that are anxious to try this out right away, this is how we would suggest to proceed". This is their code:

```{r eval=FALSE, include=TRUE}
# not run

# load the data, run SCTransform, make a project list
cntl<- CreateSeuratObject(counts = cntl.data, min.cells = 3, min.features = 200, project = "10X")
stim<- CreateSeuratObject(counts =stim.data, min.cells = 3, min.features = 200, project = "10X")
cntl@meta.data$stim <- "cntl"
stim@meta.data$stim <- "stim"
cntl<- SCTransform(object = cntl , verbose = FALSE)
stim <- SCTransform(object = stim , verbose = FALSE)
integrate.list<-objects()
integrate.list$cntl<-cntl
integrate.list$cntl<-stim

# integrate the list as normal
reference.list <- integrate.list[c("cntl","stim")]
integrate.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:40)
integrated.data <- IntegrateData(anchorset = integrate.anchors, dims = 1:40)

DefaultAssay(object = integrated.data) <- "integrated"

# scale data
integrated.data <- ScaleData(object = integrated.data, verbose = FALSE) # do I need this?
integrated.data <- RunPCA(object = integrated.data, npcs = 40, verbose = FALSE)
integrated.data <- RunUMAP(object = integrated.data, reduction = "pca", dims = 1:40)

integrated.data <- FindNeighbors(object = integrated.data)
integrated.data <- FindClusters(integrated.data, resolution = 0.6, dims.use = 1:40)
```

But the seuratlab admin adds that: "the absolute best thing to do in our view would be to correct the Pearson residuals themselves (which are stored in the @scale.data slot of the SCT assay). In this case, you would not run ScaleData after integration, as the corrected residuals would already be placed in the @scale.data slot of the 'corrected' assay. This is not currently implemented in the public version of v3, but will be soon. The results are quite similar both ways, but we prefer the latter approach for a few technical reasons (and in some cases, it does improve the biological results slightly as well)."

Since this question was posed on Github, the lab have released an [SCTransform + integration vignette](https://satijalab.org/seurat/v3.0/integration.html) which requires you to install the Seurat developers version. This allows you to SCTransform normalisation and integration automagically:
1. Each data set is run separately through sctransform
2. Pre-integration is performed using `PrepSCTIntegration`
3. Integration is run as normal, using `normalization.method = "SCT"`
4. Normal downstream analysis (visualisation, clustering...). *DO NOT run `ScaleData`*

Apparently, you can switch back and forth between development mode and normal CRAN packages with the following code:
```{r eval=FALSE, include=TRUE}
# https://stackoverflow.com/questions/9656016/how-to-install-development-version-of-r-packages-github-repository
# not run

install.packages("devtools")
library(devtools)
dev_mode(on=T)
install_github("hadley/ggplot2")
# use dev ggplot2 now
# when finished do:
dev_mode(on=F)  #and you are back to having stable ggplot2
```

See also:

A [question thread](https://github.com/ChristophH/sctransform/issues/4) from March about how to do integration and SCTransform.

Another [question thread](https://github.com/satijalab/seurat/issues/1330) from April.

The [Seurat vignette](https://satijalab.org/seurat/v3.0/integration.html) about integration and label transfer with SCTransform (note: needs the Seurat development version).

**Differential gene expression and SCTransform**

There is not currently a workaround for running differential expression on SCTransform data (or indeed on integrated data). Currently, the Seurat team [officially recommend](https://satijalab.org/seurat/faq) running differential expression tests on unintegrated data, which by default is stored in the RNA assay. This is because "integration inherently produces dependencies between data points, violating the assumotions of the stats tests used for DE"

**SCTransform and downstream analysis**
If you decided to use SCTransform for normalisation, you can use more PCs in downstream steps. This is because sctransform performs more effective normalsation, strongly removing technical effects from the data
(see "why can we choose more PCs when using sctransform" in https://satijalab.org/seurat/v3.0/sctransform_vignette.html)

## Appendix 4: The mean-variance relationship in RNASeq and scRNASeq

See Valentine Svensson's blog about [variance stabilising scRNASeq counts](http://www.nxn.se/valent/2017/10/15/variance-stabilizing-scrna-seq-counts).

In normal RNASeq data, "the process of counting however implies that variation will propagate as the number of events increase. The effect of this is that there will be an inherent relation between mean (expected value) and variance of counts." ([Svensson](http://www.nxn.se/valent/2017/10/15/variance-stabilizing-scrna-seq-counts)). 

![_In normal RNASeq, the variance increases as mean expression level increases_](https://images.squarespace-cdn.com/content/v1/5797927459cc680fdc992469/1508097104065-0AX21JYRBJQ8AIOSM04O/ke17ZwdGBToddI8pDm48kDs9ZAg36B0pAVjkbJIv4u9Zw-zPPgdn4jUwVcJE1ZvWEtT5uBSRWt4vQZAgTJucoTqqXjS3CfNDSuuf31e0tVEKE0062xkgjbYH_h_hFoxVucLFHtLqT6iNIYfJ0jiYlltO8nJtk629tZGIWiyY3XQ/1.png?format=500w)

In scRNASeq, low expressing genes have high variance. We need to stabilise the variance before selecting genes based on their variance.

Clustering is very important in the scRNASeq analysis pipeline (see Appendix 3). However, most existing clustering algorithms expect normally distributed data. We can transform our scRNASeq data to make it more similar to normal data, often with logs. However, since we have so many zero counts in scRNASeq data, it's common to add a count of 1 to all values to prevent errors associated with log-transforming 0s - `log(count + 1)`. Look what this does to the mean-variance relationship (black plots):

![_In log-transformed scRNASeq, the mean-variance relationship changes_](https://images.squarespace-cdn.com/content/v1/5797927459cc680fdc992469/1508097382788-A015UUMM6VWFY36W1JYG/ke17ZwdGBToddI8pDm48kDHKp_lUVa_H233WCatbnitZw-zPPgdn4jUwVcJE1ZvWEtT5uBSRWt4vQZAgTJucoTqqXjS3CfNDSuuf31e0tVFvx5GmvJfD1FgYrih0Tie1Aefz9qGT1nKbdUtdPppoSFtO8nJtk629tZGIWiyY3XQ/2.png?format=500w)

In the black plots, we can see that the linear mean-variance relationship persists with low-count genes, but is lost with higher-mean genes.

Our goal is to transform our data and remove the mean-variance relationship effectively, which we can do with variance stabilising transformation (vst). 

TO BE CONTINUED.

## Appendix 5: tSNE

Dimension reduction is the process of representing high dimensional data (e.g. data from many many thousands of cells) into a low dimensional graph, and yet still retain a lot of the same information. tSNE should only be used for visualisation - and you should NEVER cluster on your tSNE results. This is because while the distances/relationships between points in the same cluster is meaningful, the distances between the different clusters are not really meaningful or easy to interpret.

tSNE is a non-linear dimension reduction technique (as is UMAP; PCA is a linear dimension reduction method). You can learn about how tSNE works in this [StatQuest video](https://www.youtube.com/watch?v=NEaUSP4YerM); see also Laura's (hardcopy) notebook #1, page 163 amd 189. Here are some problems with tSNE, [according to Nikolay Oskolkov](https://www.youtube.com/watch?v=NEaUSP4YerM):

* It is slow and doesn't work well with large datasets. The algorithm fitSNE attempts to speed it up but it uses heaps of memory - you need to run it on a computer cluster.
* It doesn't preserve the global data structure - within-cluster distances are meaningful but between-cluster distances are not. So don't cluster on tSNE!
* It is a visualisation-only tool, because you can only look at 2-3 dimensions (think ~"principal components") - you can't use tSNE to look at 20-50 PCs like you can in PCA
* tSNE performs "non-parametric mapping" from high to low dimensions. It doesn't use PCA loadings (features) that drive the observed clustering
* You need to run PCA (or autoencoder) first anyway to perform pre-dimensionality reduction. tSNE can't work with high dimensional data directly
* It uses too much memory

UMAP is faster than tSNE if you have a large number of data points, >2-3 embedding dimensions, or large number of ambient dimensions in the dataset

Some benefits of UMAP over tSNE, according to [Dan Allison](https://medium.com/@dan.allison/dimensionality-reduction-with-umap-b081837354dd):

* UMAP is faster than tSNE
* UMAP captures global structure better than tSNE
* tSNE doesn't have much use outside data visualisation, but UNAP is a general-purpose dimension reduction technique
* UMAP has solid theoretical/mathematical background as a manifold approximation technique. tSNE is mostly a visualisation heuristic

## Appendix 3: Why is DE a lower priority in scRNASeq than in RNASeq?

A blogger called Valentine Svensson wrote an interesting post about [variance stabilising scRNASeq counts](http://www.nxn.se/valent/2017/10/15/variance-stabilizing-scrna-seq-counts). They write:

>_Single cell RNA-seq data is different. Not necessarily because the data wouldn't be suitible for these tools, but rather because **differential expression is a minor question of limited interest in single cell studies**. By far the most popular use of scRNA-sequencing is to **identify groups of cells which are similar to each other** and might correspond to functionally distinct cell types. In addition to such clustering analysis, **inference of developmental trajectories are popular**, as well as quantifying the **degree of variation between conditions**.
So unlike bulk RNA-sequencing, the key analysis modality is in terms of **multivariate analysis such as clustering or "dimensionality reduction" like PCA**. In the coming years I believe figuring out effective ways to think about these issues for count data will be important, especially for sparse counts from low depth!_

# Seurat resources

See also [this readme](https://github.com/ellefeg/RNASeq/blob/master/README.md)

* Seurat has many [tutorials](https://satijalab.org/seurat/vignettes.html); of particular interest is the [basic clustering vignette](https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html), the [label transfer vignette](https://satijalab.org/seurat/v3.1/integration.html) and the [control vs. stimulated PBMCs vignette](https://satijalab.org/seurat/v3.1/immune_alignment.html). Also potentially of interest are the data downloads linked with the various vignettes, e.g. the [Microwell Mouse Cell Atlas data](https://satijalab.org/seurat/v3.1/mca.html). See also this useful ["cheat sheet" of Seurat commands](https://satijalab.org/seurat/essential_commands.html).

* Ming Tang is a bioinformatician at Harvard. He often discusses scRNASeq on [twitter](https://twitter.com/tangming2005). He also has a [blog](https://divingintogeneticsandgenomics.rbind.io/) and a very useful [scRNAseq tutorial](https://crazyhottommy.github.io/scRNA-seq-workshop-Fall-2019/index.html). This tute runs over the same basic content as the Seurat vignettes, but also includes useful links, explanations and "DIY versions" of different plots. He also has his own long list of scRNASeq resources on [github](https://github.com/crazyhottommy/scRNAseq-analysis-notes), and a [separate repo of genomics tools and resources](https://github.com/crazyhottommy/getting-started-with-genomics-tools-and-resources).

* The Broad Institute has published a long set of notes for an scRNASeq workshop. Of particular interest is [section 9.4]((https://broadinstitute.github.io/2019_scWorkshop/data-wrangling-scrnaseq.html#beginning-with-seurat-httpsatijalab.orgseurat) which discuss Seurat.

* Another [Seurat tutorial](https://www.fimm.fi/sites/default/files/Seurat-guideline-10x.pdf) from FIMM (Institute for Molecular Medicine Finland); this tutorial is not very detailed.

* Yet [another tutorial](https://nbisweden.github.io/excelerate-scRNAseq/session-integration/Data_Integration.html); this one involves data integration.

* [And another](https://hbctraining.github.io/scRNA-seq/lessons/sc_exercises_clustering_analysis.html)

* [Hemberg lab scRNASeq course notes](https://github.com/hemberg-lab/scRNA.seq.course) or in [book format](https://scrnaseq-course.cog.sanger.ac.uk/website/index.html) - Chapter 9 is about Seurat

**Other software**

* [schex](https://github.com/SaskiaFreytag/schex), a tool to plot `FeaturePlot()` as hexagons, avoiding the problem of overlapping cells.

* Benchmarking papers looking at (1) different clustering tools, by [Freytag et al. 2018](https://f1000research.com/articles/7-1297) and [Duo et al. 2018](https://f1000research.com/articles/7-1141) or (2) cell type assignment tools, by [Abdelaal et al. 2018](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1795-z). See also this [Medium post](https://medium.com/@HeleneOMICtools/your-top-3-single-cell-rna-sequencing-analysis-tools-221b65fbc57e) ranking various scRNASeq analysis tools.

* [Scanpy](https://scanpy.readthedocs.io/en/stable/) is another scRNASeq analysis toolkit; I haven't used it but some of the scRNASeq data I've been sent has already been analysed with this software.

* [Monocle3](https://cole-trapnell-lab.github.io/monocle3/) is tool which is similar to Seurat. It is possible to convert Seurat data to an object compatible with Monocle2, but not yet for Monocle3.

* SingleCellExperiment is a generic R object type used for scRNASeq experiments. The Seurat object is a different object type that is based off SCE types, but is explored in a different way. [This blog post](http://lazappi.id.au/2018/06/exploring-the-sce-verse/) is an explainer for the SCE object

* [Cerebro](https://github.com/romanhaa/Cerebro/blob/master/README.md) "allows users to interactively visualize various parts of single cell transcriptomics data without requiring bioinformatic expertise". You can export a Seurat object in the right format using `cerebroApp` R package, and then do downstream visualisation. See also [the paper](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz877/5640500).

* [schex](https://github.com/SaskiaFreytag/schex), a tool to plot `FeaturePlot()` as hexagons, avoiding the problem of overlapping cells.
