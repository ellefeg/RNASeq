Laura Grice
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

* Seurat has many [tutorials](https://satijalab.org/seurat/vignettes.html); of particular interest is the [basic clustering vignette](https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html), the [label transfer vignette](https://satijalab.org/seurat/v3.1/integration.html) and the [control vs. stimulated PBMCs vignette](https://satijalab.org/seurat/v3.1/immune_alignment.html). Also potentially of interest are the data downloads linked with the various vignettes, e.g. the [Microwell Mouse Cell Atlas data](https://satijalab.org/seurat/v3.1/mca.html). See also this useful ["cheat sheet" of Seurat commands](https://satijalab.org/seurat/essential_commands.html).

* Ming Tang is a bioinformatician at Harvard. He often discusses scRNASeq on [twitter](https://twitter.com/tangming2005). He also has a [blog](https://divingintogeneticsandgenomics.rbind.io/) and a very useful [scRNAseq tutorial](https://crazyhottommy.github.io/scRNA-seq-workshop-Fall-2019/index.html). This tute runs over the same basic content as the Seurat vignettes, but also includes useful links, explanations and "DIY versions" of different plots. He also has his own long list of scRNASeq resources on [github](https://github.com/crazyhottommy/scRNAseq-analysis-notes), and a [separate repo of genomics tools and resources](https://github.com/crazyhottommy/getting-started-with-genomics-tools-and-resources).

* The Broad Institute has published a long set of notes for an scRNASeq workshop. Of particular interest is [section 9.4]((https://broadinstitute.github.io/2019_scWorkshop/data-wrangling-scrnaseq.html#beginning-with-seurat-httpsatijalab.orgseurat) which discuss Seurat.

* Another [Seurat tutorial](https://www.fimm.fi/sites/default/files/Seurat-guideline-10x.pdf) from FIMM (Institute for Molecular Medicine Finland); this tutorial is not very detailed.

* Yet [another tutorial](https://nbisweden.github.io/excelerate-scRNAseq/session-integration/Data_Integration.html); this one involves data integration.

* [And another](https://hbctraining.github.io/scRNA-seq/lessons/sc_exercises_clustering_analysis.html)

* [Hemberg lab scRNASeq course notes](https://github.com/hemberg-lab/scRNA.seq.course) or in [book format](https://scrnaseq-course.cog.sanger.ac.uk/website/index.html)

**Other software**

* [schex](https://github.com/SaskiaFreytag/schex), a tool to plot `FeaturePlot()` as hexagons, avoiding the problem of overlapping cells.

* Benchmarking papers looking at (1) different clustering tools, by [Freytag et al. 2018](https://f1000research.com/articles/7-1297) and [Duo et al. 2018](https://f1000research.com/articles/7-1141) or (2) cell type assignment tools, by [Abdelaal et al. 2018](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1795-z). See also this [Medium post](https://medium.com/@HeleneOMICtools/your-top-3-single-cell-rna-sequencing-analysis-tools-221b65fbc57e) ranking various scRNASeq analysis tools.

* [Scanpy](https://scanpy.readthedocs.io/en/stable/) is another scRNASeq analysis toolkit; I haven't used it but some of the scRNASeq data I've been sent has already been analysed with this software.

* [Monocle3](https://cole-trapnell-lab.github.io/monocle3/) is tool which is similar to Seurat. It is possible to convert Seurat data to an object compatible with Monocle2, but not yet for Monocle3.

* SingleCellExperiment is a generic R object type used for scRNASeq experiments. The Seurat object is a different object type that is based off SCE types, but is explored in a different way. [This blog post](http://lazappi.id.au/2018/06/exploring-the-sce-verse/) is an explainer for the SCE object

* [Cerebro](https://github.com/romanhaa/Cerebro/blob/master/README.md) "allows users to interactively visualize various parts of single cell transcriptomics data without requiring bioinformatic expertise". You can export a Seurat object in the right format using `cerebroApp` R package, and then do downstream visualisation. See also [the paper](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz877/5640500).
