**A general list of (sc)RNASeq resources**

Laura Grice

# scRNASeq 101

* [Awesome Single Cell](https://github.com/seandavi/awesome-single-cell) is a much more comprehensive version of this list, listing software packages, web portals/databases, key journal articles, researchers...

* An extensive [Introduction to single-cell RNA-seq analysis](http://barc.wi.mit.edu/education/hot_topics/scRNAseq_March2019/SingleCellRNAseq.pdf) from MIT (powerpoint PDF) which introduces the theory (with a bit of code) through scRNASeq technology, experimental design, analysis pipeline, pseudotime...

* [SimpleSingleCell](https://bioconductor.org/packages/release/workflows/html/simpleSingleCell.html)
 - a step-by-step workflow for basic scRNASeq analysis

* A blog post by Nikolay Oskolkov about general methods [to normalise scRNASeq data](https://towardsdatascience.com/how-to-normalize-single-cell-a438281ea654)

* Ming Tang is a bioinformatician at Harvard. He often discusses scRNASeq on [twitter](https://twitter.com/tangming2005). He also has a [blog](https://divingintogeneticsandgenomics.rbind.io/) and a very useful [scRNAseq tutorial](https://crazyhottommy.github.io/scRNA-seq-workshop-Fall-2019/index.html). This tute runs over the same basic content as the Seurat vignettes, but also includes useful links, explanations and "DIY versions" of different plots. He also has his own long list of scRNASeq resources on [github](https://github.com/crazyhottommy/scRNAseq-analysis-notes), and a [separate repo of genomics tools and resources](https://github.com/crazyhottommy/getting-started-with-genomics-tools-and-resources).

* The Broad Institute has published a long set of notes for an [scRNASeq workshop](https://broadinstitute.github.io/2019_scWorkshop/). Of particular interest is [section 9.4]((https://broadinstitute.github.io/2019_scWorkshop/data-wrangling-scrnaseq.html#beginning-with-seurat-httpsatijalab.orgseurat) which discuss Seurat.

* [Hemberg lab scRNASeq course notes](https://github.com/hemberg-lab/scRNA.seq.course) or in [book format](https://scrnaseq-course.cog.sanger.ac.uk/website/index.html)

* [A book about scRNASeq in Bioconductor](https://osca.bioconductor.org/), which mostly focusses on using SingleCellExperiment objects. It contains tutorials/code but also the thinking behind the different analyses. The book was presented in a [Nature Methods publication](https://www.nature.com/articles/s41592-019-0654-x)

# Benchmarking

* Benchmarking papers looking at (1) different clustering tools, by [Freytag et al. 2018](https://f1000research.com/articles/7-1297), [Duo et al. 2018](https://f1000research.com/articles/7-1141) and [Menon 2017](https://academic.oup.com/bfg/article/17/4/240/4728639) or (2) cell type assignment tools, by [Abdelaal et al. 2018](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1795-z). See also this [Medium post](https://medium.com/@HeleneOMICtools/your-top-3-single-cell-rna-sequencing-analysis-tools-221b65fbc57e) ranking various scRNASeq analysis tools.

* Genome Biology [Benchmarking issue](https://www.biomedcentral.com/collections/benchmarkingstudies)

* [Essential guidelines for computational method benchmarking](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1738-8)

* A quest for benchmarking papers

<blockquote class="twitter-tweet"><p lang="en" dir="ltr">Dear <a href="https://twitter.com/hashtag/academictwitter?src=hash&amp;ref_src=twsrc%5Etfw">#academictwitter</a>, I&#39;m looking for a list of benchmarks for single cell RNA-seq data analysis. Generally, so normalization, DE, clustering, cell assignment, etc. I&#39;m OFC aware of the &quot;Methods comparisons&quot; section in <a href="https://twitter.com/seandavis12?ref_src=twsrc%5Etfw">@seandavis12</a>&#39;s awesome list (<a href="https://t.co/lT7DivJdAZ">https://t.co/lT7DivJdAZ</a>) ..</p>&mdash; Mark Robinson (@markrobinsonca) <a href="https://twitter.com/markrobinsonca/status/1220048175320440832?ref_src=twsrc%5Etfw">January 22, 2020</a></blockquote> <script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script> 

# Clustering

* Thank you to Ming Tang's resources for [general bioinformatics](https://github.com/crazyhottommy/getting-started-with-genomics-tools-and-resources) and [RNASeq](https://github.com/crazyhottommy/RNA-seq-analysis) - several of of these links were  taken from these lists

* Clustering benchmark by [Freytag et al. 2018](https://f1000research.com/articles/7-1297) or [Duo et al. 2018](https://f1000research.com/articles/7-1141) or [Menon 2017](https://academic.oup.com/bfg/article/17/4/240/4728639)

* Seurat [clustering tutorial](https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html) and [paper](www.cell.com/abstract/S0092-8674(15)00549-8). Seurat implements graph-based clustering. Clustering embeds cells in a graph structure (e.g a k-nearest neighbour - KNN - graph); edges are drawn between cells with similar feature expression patterns. Then the graph is partitioned into highly inter-connected "quasi-cliques" or "communities". The KNN graph is built based on euclidian distance in PCA space (hence you specify `reduction = "pca"` in `FindNeighbors`). Edge weights between cells are refined based on shared overlap in their local neighbourhood (Jaccard similarty. Then we cluster cells. Cells are iteratively grouped together. When we run `FindClusters`, we need to set the resolution value. This is the granularity of the clusters - the higher the value, the more clusters you get. The bigger the dataset, the more resolution you probably want to use. The Seurat authors recommend a range of 0.4 - 1.2 for ~3000 cell datasets. (For reference, I have 1223 cells in `treat`)

* SC3 [package](http://bioconductor.org/packages/release/bioc/html/SC3.html) and [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5410170/) - consensus clustering of single-cell RNA-Seq data. A non-graph-based clustering approach. 
SC3 achieves high accuracy and robustness by consistently integrating different clustering solutions through a consensus approach. Determines optimal K value, then iteratively clusters and generates consensus set. Can run as a shiny object to visualise data. Not particularly fast. 

* IKAP [paper](https://academic.oup.com/gigascience/article/8/10/giz121/5579995) and [github](https://github.com/NHLBI-BCB/IKAP). Runs based on Seurat object. For a range of PCs and K-values, performs clustering and determines the optimal combination.



* [sincell](http://bioconductor.org/packages/devel/bioc/html/sincell.html): R package for the statistical assessment of cell state hierarchies from single-cell RNA-seq data

* [optCluster](https://cran.r-project.org/web/packages/optCluster/index.html): An R Package for Determining the Optimal Clustering Algorithm. Not sure if it's designed for scRNASeq

* [ConsensusClusterPlus algorithm](https://bioconductor.org/packages/release/bioc/html/ConsensusClusterPlus.html) for determining cluster count and membership by stability evidence in unsupervised analysis. Not sure if it's designed for scRNASeq

* [Interactive visualisation and fast computation of the solution path: convex bi-clustering](https://www.youtube.com/watch?v=2g-akN6q8aI) by Genevera Allen cvxbiclustr and the clustRviz package coming.

* Paper: [Fast and accurate single-cell RNA-Seq analysis by clustering of transcript-compatibility counts](http://biorxiv.org/content/early/2016/01/15/036863)

* [CIDR](https://github.com/VCCRI/CIDR): Ultrafast and accurate clustering through imputation for single-cell RNA-seq data

* [pcaReduce](http://biorxiv.org/content/early/2015/09/08/026385): Hierarchical Clustering of Single Cell Transcriptional Profiles.

* [CountClust](https://www.bioconductor.org/packages/3.3/bioc/html/CountClust.html): Clustering and Visualizing RNA-Seq Expression Data using Grade of Membership Models. Fits grade of membership models (GoM, also known as admixture models) to cluster RNA-seq gene expression count data, identifies characteristic genes driving cluster memberships, and provides a visual summary of the cluster memberships

* [FastProject](http://biorxiv.org/content/early/2016/03/12/043463): A Tool for Low-Dimensional Analysis of Single-Cell RNA-Seq Data

* [Compare clusterings for single-cell sequencing bioconductor package](http://bioconductor.org/packages/devel/bioc/html/clusterExperiment.html).The goal of this package is to encourage the user to try many different clustering algorithms in one package structure. We give tools for running many different clusterings and choices of parameters. We also provide visualization to compare many different clusterings and algorithm tools to find common shared clustering patterns.



**theory**

* [PCA, MDS, k-means, Hierarchical clustering and heatmaps](https://rpubs.com/crazyhottommy/PCA_MDS) by Ming Tang

* [Cluster Analysis in R](https://www.datanovia.com/en/blog/types-of-clustering-methods-overview-and-quick-start-r-code/#at_pco=smlre-1.0&at_si=58765a95fcb21379&at_ab=per-2&at_pos=3&at_tot=4) - Unsupervised machine learning very practical intro on STHDA website.

* Clustering analysis for high-dimentional biological data: [Avoiding common pitfalls when clustering biological data](http://stke.sciencemag.org/content/9/432/re6)

* [How does gene expression clustering work?](http://www.nature.com/nbt/journal/v23/n12/full/nbt1205-1499.html)

* [Understanding UMAP](https://pair-code.github.io/understanding-umap/)

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
