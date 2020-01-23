Laura Grice
**A general list of (sc)RNASeq resources**

# scRNASeq 101

* [Awesome Single Cell](https://github.com/seandavi/awesome-single-cell) is a more comprehensive version of this list, listing software packages, web portals/databases, key journal articles, researchers...

* An extensive [Introduction to single-cell RNA-seq analysis](http://barc.wi.mit.edu/education/hot_topics/scRNAseq_March2019/SingleCellRNAseq.pdf) from MIT (powerpoint PDF) which introduces the theory (with a bit of code) through scRNASeq technology, experimental design, analysis pipeline, pseudotime...

* [SimpleSingleCell](https://bioconductor.org/packages/release/workflows/html/simpleSingleCell.html)
 - a step-by-step workflow for basic scRNASeq analysis

* A blog post by Nikolay Oskolkov about general methods [to normalise scRNASeq data](https://towardsdatascience.com/how-to-normalize-single-cell-a438281ea654)

* Ming Tang is a bioinformatician at Harvard. He often discusses scRNASeq on [twitter](https://twitter.com/tangming2005). He also has a [blog](https://divingintogeneticsandgenomics.rbind.io/) and a very useful [scRNAseq tutorial](https://crazyhottommy.github.io/scRNA-seq-workshop-Fall-2019/index.html). This tute runs over the same basic content as the Seurat vignettes, but also includes useful links, explanations and "DIY versions" of different plots. He also has his own long list of scRNASeq resources on [github](https://github.com/crazyhottommy/scRNAseq-analysis-notes), and a [separate repo of genomics tools and resources](https://github.com/crazyhottommy/getting-started-with-genomics-tools-and-resources).

* The Broad Institute has published a long set of notes for an [scRNASeq workshop](https://broadinstitute.github.io/2019_scWorkshop/). Of particular interest is [section 9.4]((https://broadinstitute.github.io/2019_scWorkshop/data-wrangling-scrnaseq.html#beginning-with-seurat-httpsatijalab.orgseurat) which discuss Seurat.

* [Hemberg lab scRNASeq course notes](https://github.com/hemberg-lab/scRNA.seq.course) or in [book format](https://scrnaseq-course.cog.sanger.ac.uk/website/index.html)

* [A book about scRNASeq in Bioconductor](https://osca.bioconductor.org/), which mostly focusses on using SingleCellExperiment objects. It contains tutorials/code but also the thinking behind the different analyses. The book was presented in a [Nature Methods publication](https://www.nature.com/articles/s41592-019-0654-x)

# Benchmarking

* Benchmarking papers looking at (1) different clustering tools, by [Freytag et al. 2018](https://f1000research.com/articles/7-1297) and [Duo et al. 2018](https://f1000research.com/articles/7-1141) or (2) cell type assignment tools, by [Abdelaal et al. 2018](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1795-z). See also this [Medium post](https://medium.com/@HeleneOMICtools/your-top-3-single-cell-rna-sequencing-analysis-tools-221b65fbc57e) ranking various scRNASeq analysis tools.

# Tools
* [schex](https://github.com/SaskiaFreytag/schex), a tool to plot `FeaturePlot()` as hexagons, avoiding the problem of overlapping cells.

* [Scanpy](https://scanpy.readthedocs.io/en/stable/) is another scRNASeq analysis toolkit; I haven't used it but some of the scRNASeq data I've been sent has already been analysed with this software.

* [Monocle3](https://cole-trapnell-lab.github.io/monocle3/) is tool which is similar to Seurat. It is possible to convert Seurat data to an object compatible with Monocle2, but not yet for Monocle3.

* SingleCellExperiment is a generic R object type used for scRNASeq experiments. The Seurat object is a different object type that is based off SCE types, but is explored in a different way. [This blog post](http://lazappi.id.au/2018/06/exploring-the-sce-verse/) is an explainer for the SCE object

* [Cerebro](https://github.com/romanhaa/Cerebro/blob/master/README.md) "allows users to interactively visualize various parts of single cell transcriptomics data without requiring bioinformatic expertise". You can export a Seurat object in the right format using `cerebroApp` R package, and then do downstream visualisation. See also [the paper](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz877/5640500).


# Datasets and databases

See "List of single cell data and databases" on Lab Archives

# Statistics

* The [StatQuest Youtube page](https://www.youtube.com/channel/UCtYLUTtgS3k1Fg4y5tAhLbw) has a lot of statistical explanation videos, and in particular a playlist of videos about [high throughput sequencing](https://www.youtube.com/playlist?list=PLblh5JKOoLUJo2Q6xK4tZElbIvAACEykp). For instance, see his short videos about [PCA](https://www.youtube.com/watch?v=HMOI_lkzW08) or [tSNE](https://www.youtube.com/watch?v=NEaUSP4YerM).

* UMAP algorithm resources: [How Exactly UMAP Works](https://towardsdatascience.com/how-exactly-umap-works-13e3040e1668 ) blog post, [Understanding UMAP](https://pair-code.github.io/understanding-umap/) overview, introduction to the UMAP algorithm at [SciPy 2018](https://www.youtube.com/watch?v=nq6iPZVUxZU) by the original author.

* Dave Tang is a computational biologist in Perth at UWA. He has written about scRNASeq [on his computational biology and genomics blog](https://davetang.org/muse/category/single-cell-2/). However, these posts (especially the ones about [Seurat](https://davetang.org/muse/2017/08/01/getting-started-seurat/) and [Monocle](https://davetang.org/muse/2017/10/01/getting-started-monocle/) are quite old (2017-2018).

* Sajita Lab's [Annual Single Cell Genomics Day Workshops](https://satijalab.org/scgd/). See also the workshops from [2019](http://www.satijalab.org/scgd19/) and [2018](http://www.satijalab.org/scgd18/). Slides and videos are available for previous years.

* Principal component analysis (PCA) is a method of linear dimension reduction. If you want to learn more about how PCA works, you can watch a [short](https://www.youtube.com/watch?v=HMOI_lkzW08) or [long](https://www.youtube.com/watch?v=_UVHneBUBW0) video from StatQuest, or read [Ming Tang's blog post](https://divingintogeneticsandgenomics.rbind.io/post/pca-in-action/) about PCA or [permutation tests for PCA components](https://divingintogeneticsandgenomics.rbind.io/post/permute-test-for-pca-components/). 

* How to perform [PCA in R](http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/). More [PCA in R and Python](https://www.analyticsvidhya.com/blog/2016/03/practical-guide-principal-component-analysis-python/), with suggested plots.


# Computing

* The [difference](https://blogs.nvidia.com/blog/2018/08/02/supervised-unsupervised-learning/) between supervised and unsupervised learning (relevant to many aspects of computational biology, but especially to cell type annotation)

* How to [structure R projects](https://chrisvoncsefalvay.com/2018/08/09/structuring-r-projects/) so you don't go crazy. [Another view](https://nicercode.github.io/blog/2013-04-05-projects/) and [another by Software Carpentry](https://swcarpentry.github.io/r-novice-inflammation/06-best-practices-R/)

* A [long thread](https://community.rstudio.com/t/data-science-project-template-for-r/3230/13) of suggestions for R coding project folder structure, and a [suggested folder structure](https://github.com/pavopax/new-project-template) that I like. Another long list of Twitter-sourced ideas about [R workflows](https://maraaverick.rbind.io/2017/09/r-workflow-fun/)

* For visualisation, the [R Graph Gallery](https://www.r-graph-gallery.com/) has different R/ggplot2 charts and the code to make them. The [ggplot flipbook](https://evamaerey.github.io/ggplot_flipbook/) shows how to build plots from the ground up, with step-by-step code. Or there's the [R colour picker](https://deanattali.com/blog/colourpicker-ggmarginal-gadgets/) which gives you the R HEX/name codes for colours.

* [A Scientist's Guide to R](https://craighutton.netlify.com/post/2019-05-17-asgr-basic-workflow/) is a series of blog posts "on how to use R as an analytical and productivity tool in the process of conducting scientific research"

* Introductory book, [A Primer for Computational Biology](https://open.oregonstate.education/computationalbiology/). Part 1 - intro to Unix, part 2 - programming in Python, part 3 - programming in R.
