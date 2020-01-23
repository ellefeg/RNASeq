_adapted from 20191213_HowToScrnaseq.Rmd_

Seurat is a very powerful R package for performing quality control, analysis and exploration of scRNASeq data. You can perform unsupervised clustering and cell type/state discovery, spatial reconstruction of single cell data, and perform data "anchoring/integration" to integrate sequencing data from different sources. This latter function is quite amazing - you can take samples from not only different conditions, but also different labs, different sequencing technologies (scRNASeq, ATACSeq, RNASeq, CELSeq...) or even different species. You can also perform label transfer, where the annotations from one dataset can be projected onto a new unknown dataset of interest. 

According to two clustering benchmarking papers (Freytag et al. 2018 and Duo et al. 2018), Seurat is one of the best clustering algorithms currently available for scRNASeq data. However, Quan is concerned that all of the Seurat output figures look the same, so while it is a good "starter" tool, additional tools and more creative graphics are required to push our publications to the next level.

The Seurat website is very good, it contains many [vignettes](https://satijalab.org/seurat/vignettes.html) to guide you through the Seurat functionalities, and also a useful ["cheat sheet" of Seurat commands](https://satijalab.org/seurat/essential_commands.html).

There are other alternative tools to Seurat; Monocle, by Cole Trapnell's group, is amongst the most popular alternatives. It is not currently (January 2020) possible to (automatically) convert between Seurat3 and Monocle3 objects.

**A note about Seurat objects**

A `Seurat` object is essentially the R object that holds all of your experimental input data and the results of all your downstream analyses in Seurat. It is an S4 object in R, so it is a little more difficult to interact with than your standard S3 R objects. You can learn more about `Seurat` objects at the [Seurat Wiki page](https://github.com/satijalab/seurat/wiki)

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

Once you integrate 2+ samples, you'll need to tell Seurat whether you want to work on the integrated or unintegrated (raw) data - some analyses are incompatible with one or the other. To do this, run `DefaultAssay(myObject) <- "integrated" # or "RNA"`.

# Seurat resources

See [RNASeq repo README](https://github.com/ellefeg/RNASeq/blob/master/README.md)
