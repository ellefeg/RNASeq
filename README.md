**A general list of (sc)RNASeq resources**

Laura Grice

# scRNASeq 101

* [Awesome Single Cell](https://github.com/seandavi/awesome-single-cell) is a much more comprehensive version of this list, listing software packages, web portals/databases, key journal articles, researchers...

* Ming Tang's lists of [scRNASeq resources](https://github.com/crazyhottommy/scRNAseq-analysis-notes) and [genomics tools and resources](https://github.com/crazyhottommy/getting-started-with-genomics-tools-and-resources)

* [scRNA-tool](https://www.scrna-tools.org/) database - "a catalogue of tools analysing single-cell RNA sequencing data"

* Ming Tang is a bioinformatician at Harvard. He often discusses scRNASeq on [twitter](https://twitter.com/tangming2005). He also has a [blog](https://divingintogeneticsandgenomics.rbind.io/)

* An extensive [Introduction to single-cell RNA-seq analysis](http://barc.wi.mit.edu/education/hot_topics/scRNAseq_March2019/SingleCellRNAseq.pdf) from MIT (powerpoint PDF) which introduces the theory (with a bit of code) through scRNASeq technology, experimental design, analysis pipeline, pseudotime...

* The [DataCamp interactive course](https://www.datacamp.com/courses/single-cell-rna-seq-workflows-in-r?tap_a=5644-dce66f&tap_s=411670-1f1ebc) is mostly paid but Chapter 1 (What is Single-Cell RNASeq) is free

* Ming Tang's very useful [scRNAseq tutorial](https://crazyhottommy.github.io/scRNA-seq-workshop-Fall-2019/index.html). Most of the tutorial covers the basic Seurat vignette, but it also includes useful links, explanations and "DIY versions" of different plots.

* [SimpleSingleCell](https://bioconductor.org/packages/release/workflows/html/simpleSingleCell.html)
 - a step-by-step workflow for basic scRNASeq analysis

* The Broad Institute has published a long set of notes for an [scRNASeq workshop](https://broadinstitute.github.io/2019_scWorkshop/). Of particular interest is [section 9.4]((https://broadinstitute.github.io/2019_scWorkshop/data-wrangling-scrnaseq.html#beginning-with-seurat-httpsatijalab.orgseurat) which discuss Seurat.

* [Hemberg lab scRNASeq course notes](https://github.com/hemberg-lab/scRNA.seq.course) or in [book format](https://scrnaseq-course.cog.sanger.ac.uk/website/index.html)

* eBook: [scRNASeq in Bioconductor](https://osca.bioconductor.org/). Mostly focusses on using SingleCellExperiment objects. It contains tutorials/code but also the thinking behind the different analyses. The book was presented in a [Nature Methods publication](https://www.nature.com/articles/s41592-019-0654-x).

* [Stuart & Satika 2019](https://www.nature.com/articles/s41576-019-0093-7) - Integrative single-cell analysis. *"In this Review, we discuss the recent advances in the collection and integration of different data types at single-cell resolution with a focus on the integration of gene expression data with other types of single-cell measurement"*

* [Holland et al. 2019](https://www.biorxiv.org/content/10.1101/753319v1) - Robustness and applicability of functional genomics tools on scRNA-seq data. *"Our analyses suggest that bulk functional genomics tools can be applied to scRNA-seq data, outperforming dedicated single cell tools. Furthermore we provide a benchmark for further methods development by the community."*

# Benchmarking

* Genome Biology [Benchmarking issue](https://www.biomedcentral.com/collections/benchmarkingstudies)

* Weber et al. 2019: [Essential guidelines for computational method benchmarking](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1738-8)

* **Dimensionality reduction:** [Sun et al. 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1898-6)

* **Clustering:** [Freytag et al. 2018](https://f1000research.com/articles/7-1297)

* **Clustering:** [Duo et al. 2018](https://f1000research.com/articles/7-1141)

* **Clustering:** [Menon 2017](https://academic.oup.com/bfg/article/17/4/240/4728639)

* **Clustering:** [Krzak et al. 2019](https://www.frontiersin.org/articles/10.3389/fgene.2019.01253/full?&utm_source=Email_to_authors_&utm_medium=Email&utm_content=T1_11.5e1_author&utm_campaign=Email_publication&field=&journalName=Frontiers_in_Genetics&id=486077)

* **Clustering:** [Qi et al. 2019](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbz062/5528236)

* **Cell type assignment:** [Abdelaal et al. 2018](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1795-z)

* **False discovery control:** [Korthauer et al. 2019](https://www.ncbi.nlm.nih.gov/pubmed/?term=A+practical+guide+to+methods+controlling+false+discoveries+in+computational+biology)

* [Tian et al. 2019](https://www.nature.com/articles/s41592-019-0425-8) - Benchmarking single cell RNA-sequencing analysis pipelines using mixture control experiments - *"Here, we generated a realistic benchmark experiment that included single cells and admixtures of cells or RNA to create ‘pseudo cells’ from up to five distinct cancer cell lines. In total, 14 datasets were generated using both droplet and plate-based scRNA-seq protocols. We compared 3,913 combinations of data analysis methods for tasks ranging from normalization and imputation to clustering, trajectory analysis and data integration. Evaluation revealed pipelines suited to different types of data for different tasks."*

* A [Medium post](https://medium.com/@HeleneOMICtools/your-top-3-single-cell-rna-sequencing-analysis-tools-221b65fbc57e) ranking various scRNASeq analysis tools.

* A Twitter quest for benchmarking papers:

<blockquote class="twitter-tweet"><p lang="en" dir="ltr">Dear <a href="https://twitter.com/hashtag/academictwitter?src=hash&amp;ref_src=twsrc%5Etfw">#academictwitter</a>, I&#39;m looking for a list of benchmarks for single cell RNA-seq data analysis. Generally, so normalization, DE, clustering, cell assignment, etc. I&#39;m OFC aware of the &quot;Methods comparisons&quot; section in <a href="https://twitter.com/seandavis12?ref_src=twsrc%5Etfw">@seandavis12</a>&#39;s awesome list (<a href="https://t.co/lT7DivJdAZ">https://t.co/lT7DivJdAZ</a>) ..</p>&mdash; Mark Robinson (@markrobinsonca) <a href="https://twitter.com/markrobinsonca/status/1220048175320440832?ref_src=twsrc%5Etfw">January 22, 2020</a></blockquote> <script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script> 

* And a [Google Docs spreadsheet](https://docs.google.com/spreadsheets/d/1Gqn0eZ8oiNh8-9ovyh4D_ZsoWcKYoKTJ0C1AfrJrdRA/edit#gid=0) summarising the results:

<blockquote class="twitter-tweet"><p lang="en" dir="ltr">OK, my current list is 31 scRNA-seq (analysis method) benchmarks:<a href="https://t.co/gfan8Kupcp">https://t.co/gfan8Kupcp</a><br><br>I consolidated <a href="https://twitter.com/seandavis12?ref_src=twsrc%5Etfw">@seandavis12</a>&#39;s <a href="https://twitter.com/_lazappi_?ref_src=twsrc%5Etfw">@_lazappi_</a>&#39;s <a href="https://twitter.com/JMA_Data?ref_src=twsrc%5Etfw">@JMA_Data</a>&#39;s and went through the <a href="https://twitter.com/GenomeBiology?ref_src=twsrc%5Etfw">@GenomeBiology</a> benchmark issue.<br><br>Any others I missed? <a href="https://t.co/AEpHCWcYsS">https://t.co/AEpHCWcYsS</a></p>&mdash; Mark Robinson (@markrobinsonca) <a href="https://twitter.com/markrobinsonca/status/1221903274678333441?ref_src=twsrc%5Etfw">January 27, 2020</a></blockquote> <script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script> 

* Benchmarking: [Machine Learning and Statistical Methods for Clustering Single Cell RNA-sequencing Data](https://github.com/kuanglab/single-cell-review)


# Clustering

* Thank you to Ming Tang's resources for [general bioinformatics](https://github.com/crazyhottommy/getting-started-with-genomics-tools-and-resources) and [RNASeq](https://github.com/crazyhottommy/RNA-seq-analysis) and also to [awesome single cell](https://github.com/seandavi/awesome-single-cell) - several of of these links were  taken from these lists

* See "Benchmarking" section above for Benchmarking papers

* [Kiselev et al. 2019](https://www.nature.com/articles/s41576-018-0088-9) - Challenges in unsupervised clustering of single-cell RNA-seq data

* [Lee & Hemberg 2019](https://www.nature.com/articles/s41592-019-0534-4) - Supervised clustering for single-cell analysis

* [Petegrosso et al. 2019](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbz063/5519426) - Machine learning and statistical methods for clustering single-cell RNA-sequencing data 

* [Kim et al. 2019](https://www.ncbi.nlm.nih.gov/pubmed/30137247) - Impact of similarity metrics on single-cell RNA-seq data clustering.

* PanoView - [Hu et al. 2019](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007040) - PanoView: An iterative clustering method for single-cell RNA sequencing data

* SOUP - [Zhu et al. 2019](https://www.pnas.org/content/116/2/466) - Semisoft clustering of single-cell data

* SAME-Clustering - [Huh et al. 2019](https://academic.oup.com/nar/article/48/1/86/5644992) - SAME-clustering: Single-cell Aggregated Clustering via Mixture Model Ensemble 

* VPAC - [Chen et al. 2019](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2742-4) - VPAC: Variational projection for accurate clustering of single-cell transcriptomic data

* PARC - [Stassen et al. 2020](https://www.ncbi.nlm.nih.gov/pubmed/31971583?dopt=Abstract&utm_source=dlvr.it&utm_medium=twitter) - PARC: ultrafast and accurate clustering of phenotypic data of millions of single cells.

* Seurat [clustering tutorial](https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html) and [paper](www.cell.com/abstract/S0092-8674(15)00549-8). Seurat implements graph-based clustering. Clustering embeds cells in a graph structure (e.g a k-nearest neighbour - KNN - graph); edges are drawn between cells with similar feature expression patterns. Then the graph is partitioned into highly inter-connected "quasi-cliques" or "communities". The KNN graph is built based on euclidian distance in PCA space (hence you specify `reduction = "pca"` in `FindNeighbors`). Edge weights between cells are refined based on shared overlap in their local neighbourhood (Jaccard similarty. Then we cluster cells. Cells are iteratively grouped together. When we run `FindClusters`, we need to set the resolution value. This is the granularity of the clusters - the higher the value, the more clusters you get. The bigger the dataset, the more resolution you probably want to use. The Seurat authors recommend a range of 0.4 - 1.2 for ~3000 cell datasets. (For reference, I have 1223 cells in `treat`)

* SC3 [package](http://bioconductor.org/packages/release/bioc/html/SC3.html) and [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5410170/) - consensus clustering of single-cell RNA-Seq data. A non-graph-based clustering approach. 
SC3 achieves high accuracy and robustness by consistently integrating different clustering solutions through a consensus approach. Determines optimal K value, then iteratively clusters and generates consensus set. Can run as a shiny object to visualise data. Not particularly fast. 

* IKAP [paper](https://academic.oup.com/gigascience/article/8/10/giz121/5579995) and [github](https://github.com/NHLBI-BCB/IKAP). Runs based on Seurat object. For a range of PCs and K-values, performs clustering and determines the optimal combination.

* DC3 - [Zeng et al. 2019](https://www.ncbi.nlm.nih.gov/pubmed/31601804?dopt=Abstract&utm_source=dlvr.it&utm_medium=twitter) - DC3 is a method for deconvolution and coupled clustering from bulk and single-cell genomics data.

* [GOAE and GONN](https://www.biorxiv.org/content/10.1101/437020v1) - Combining Gene Ontology with Deep Neural Networks to Enhance the Clustering of Single Cell RNA-Seq Data

* [celaref](http://bioconductor.org/packages/release/bioc/html/celaref.html) - Single-cell RNAseq cell cluster labelling by reference

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

* Galapagos - [Wagner 2019](https://www.biorxiv.org/content/10.1101/770388v1) - Straightforward clustering of single-cell RNA-Seq data with t-SNE and DBSCAN. *"Here, I propose Galapagos, a simple and effective clustering workflow based on t-SNE and DBSCAN that does not require a gene selection step."*

* [clusterExperiment](https://github.com/epurdom/clusterExperiment) - [R] - Functions for running and comparing many different clusterings of single-cell sequencing data. Meant to work with SCONE and slingshot.

* Some other (more general-purpose?) packages that include clustering are: [ascend](https://github.com/IMB-Computational-Genomics-Lab/ascend) (R), [CALISTA](https://github.com/CABSEL/CALISTA) (R), [celda](https://github.com/campbio/celda) (R), [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) (10x Genomics - linux)

* [CountClust](https://github.com/kkdey/CountClust) - [R] - Functions for fitting Grade-of-Membership models, also known as "Topic models", to RNA-seq counts. These models generalize clustering methods to allow that each cell may belong to more than one cluster/topic.

* [dropClust](https://github.com/debsin/dropClust) - [R/Python] - Efficient clustering of ultra-large scRNA-seq data.

* [GiniClust](https://github.com/lanjiangboston/GiniClust) - [Python/R] - GiniClust is a clustering method implemented in Python and R for detecting rare cell-types from large-scale single-cell gene expression data. GiniClust can be applied to datasets originating from different platforms, such as multiplex qPCR data, traditional single-cell RNAseq or newly emerging UMI-based single-cell RNAseq, e.g. inDrops and Drop-seq.

* [netSmooth](https://github.com/BIMSBbioinfo/netSmooth) - [R] - netSmooth is a network-diffusion based method that uses priors for the covariance structure of gene expression profiles on scRNA-seq experiments in order to smooth expression values. We demonstrate that netSmooth improves clustering results of scRNA-seq experiments from distinct cell populations, time-course experiments, and cancer genomics.

* [pcaReduce](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0984-y) - [R] - hierarchical clustering of single cell transcriptional profiles.

* [SAKE](https://github.com/naikai/sake) - [R] - Single-cell RNA-Seq Analysis and Clustering Evaluation.

* [SC3](https://github.com/hemberg-lab/sc3) - [R] - SC3 is a tool for the unsupervised clustering of cells from single cell RNA-Seq experiments.

* [SCell](https://github.com/diazlab/SCell) - [matlab] - SCell is an integrated software tool for quality filtering, normalization, feature selection, iterative dimensionality reduction, clustering and the estimation of gene-expression gradients from large ensembles of single-cell RNA-seq datasets. SCell is open source, and implemented with an intuitive graphical interface.

* [SCENIC](https://github.com/aertslab/SCENIC) - [R] - SCENIC is an R package to infer Gene Regulatory Networks and cell types from single-cell RNA-seq data. [SCENIC: single-cell regulatory network inference and clustering](https://www.nature.com/articles/nmeth.4463)

* SCMarker - [R] - SCMarker is a method performing ab initial marker gene set selection from scRNA-seq data to achieve improved clustering/cell-typing results. SCMarker: ab initio marker selection for single cell transcriptome profiling.

* SCUBA - [matlab/R] - SCUBA stands for "Single-cell Clustering Using Bifurcation Analysis." SCUBA is a novel computational method for extracting lineage relationships from single-cell gene expression data, and modeling the dynamic changes associated with cell differentiation.

SIMLR - [R, matlab] - SIMLR (Single-cell Interpretation via Multi-kernel LeaRning) learns an appropriate distance metric from the data for dimension reduction, clustering and visualization. SIMLR is capable of separating known subpopulations more accurately in single-cell data sets than do existing dimension reduction methods.

sincera - [R] - R-based pipeline for single-cell analysis including clustering and visualization.

scClustViz - An interactive R Shiny tool for visualizing single-cell RNAseq clustering results from common analysis pipelines (SingleCellExperiment or Seurat, currently). Its main goal is two-fold: A: to help select a biologically appropriate resolution or K from clustering results by assessing differential expression between the resulting clusters; and B: help annotate cell types and identify marker genes. See the demo app here! scClustViz can also be used to generate R data packages for sharing published data - see the website for details and a list of published datasets.

SeuratWizard - a web-based (wizard style) interactive R Shiny application to perform guided single-cell RNA-seq data analysis and clustering. demo

SeuratV3Wizard - a web-based (wizard style) interactive R Shiny application to perform guided single-cell RNA-seq data analysis and clustering based on Seurat v3. demo

**clustering theory**

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

* A blog post by Nikolay Oskolkov about general methods [to normalise scRNASeq data](https://towardsdatascience.com/how-to-normalize-single-cell-a438281ea654)

# Computing

* The [difference](https://blogs.nvidia.com/blog/2018/08/02/supervised-unsupervised-learning/) between supervised and unsupervised learning (relevant to many aspects of computational biology, but especially to cell type annotation)

* How to [structure R projects](https://chrisvoncsefalvay.com/2018/08/09/structuring-r-projects/) so you don't go crazy. [Another view](https://nicercode.github.io/blog/2013-04-05-projects/) and [another by Software Carpentry](https://swcarpentry.github.io/r-novice-inflammation/06-best-practices-R/)

* A [long thread](https://community.rstudio.com/t/data-science-project-template-for-r/3230/13) of suggestions for R coding project folder structure, and a [suggested folder structure](https://github.com/pavopax/new-project-template) that I like. Another long list of Twitter-sourced ideas about [R workflows](https://maraaverick.rbind.io/2017/09/r-workflow-fun/)

* For visualisation, the [R Graph Gallery](https://www.r-graph-gallery.com/) has different R/ggplot2 charts and the code to make them. The [ggplot flipbook](https://evamaerey.github.io/ggplot_flipbook/) shows how to build plots from the ground up, with step-by-step code. Or there's the [R colour picker](https://deanattali.com/blog/colourpicker-ggmarginal-gadgets/) which gives you the R HEX/name codes for colours.

* [A Scientist's Guide to R](https://craighutton.netlify.com/post/2019-05-17-asgr-basic-workflow/) is a series of blog posts "on how to use R as an analytical and productivity tool in the process of conducting scientific research"

* Introductory book, [A Primer for Computational Biology](https://open.oregonstate.education/computationalbiology/). Part 1 - intro to Unix, part 2 - programming in Python, part 3 - programming in R.
