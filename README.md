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

Thank you to [Mark Robinson](https://twitter.com/markrobinsonca/status/1221903274678333441?ref_src=twsrc%5Etfw) for part of this list

* Genome Biology [Benchmarking issue](https://www.biomedcentral.com/collections/benchmarkingstudies)

* Weber et al. 2019: [Essential guidelines for computational method benchmarking](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1738-8)

* [Zhang & Zhang 2018](https://www.ncbi.nlm.nih.gov/pubmed/29994128) - **Comparison of computational methods for imputing single-cell RNA-sequencing data** - Studies the dropout effect. *"To this end, we compared eight imputation methods, evaluated their power in recovering original real data, and performed broad analyses to explore their effects on clustering cell types, detecting differentially expressed genes, and reconstructing lineage trajectories in the context of both simulated and real data. Simulated datasets and case studies highlight that there are no one method performs the best in all the situations. Some defects of these methods such as scalability, robustness and unavailability in some situations need to be addressed in future studies."*

* **General pipelines:** [Vieth et al. 2019](https://www.nature.com/articles/s41467-019-12266-7) -  A systematic evaluation of single cell RNA-seq analysis pipelines - *"We find that choices of normalisation and library preparation protocols have the biggest impact on scRNA-seq analyses. Specifically, we find that library preparation determines the ability to detect symmetric expression differences, while normalisation dominates pipeline performance in asymmetric DE-setups"*

* **Bulk RNASeq to scRNASeq:** [Holland et al. 2019](https://www.biorxiv.org/content/10.1101/753319v1) - Robustness and applicability of functional genomics tools on scRNA-seq data - *"Our benchmarks on both the simulated and real data revealed comparable performance to the original bulk data. Additionally, we showed that the TF and pathway activities preserve cell-type specific variability by analysing a mixture sample sequenced with 13 scRNA-seq different protocols. Our analyses suggest that bulk functional genomics tools can be applied to scRNA-seq data, outperforming dedicated single cell tools"*

* **General pipelines:** [Tian et al. 2019](https://www.nature.com/articles/s41592-019-0425-8) - Benchmarking single cell RNA-sequencing analysis pipelines using mixture control experiments - *"Evaluation revealed pipelines suited to different types of data for different tasks. Our data and analysis provide a comprehensive framework for benchmarking most common scRNA-seq analysis steps."*

* **Normalisation:** [Ding et al. 2017](https://www.ncbi.nlm.nih.gov/pubmed/28468817) - Assessment of Single Cell RNA-Seq Normalization Methods - *"Our analyses showed that methods considering spike-in External RNA Control Consortium (ERCC) RNA molecules significantly outperformed those not considering ERCCs"*

* **Normalisation:** [Cole et al. 2019](https://www.sciencedirect.com/science/article/abs/pii/S2405471219300808) - Performance Assessment and Selection of Normalization Procedures for Single-Cell RNA-Seq - *"We have developed “scone”— a flexible framework for assessing performance based on a comprehensive panel of data-driven metrics. Through graphical summaries and quantitative reports, scone summarizes trade-offs and ranks large numbers of normalization methods by panel performance ... top-performing normalization methods lead to better agreement with independent validation data for a collection of scRNA-seq datasets"*

* **Dimensionality reduction:** [Sun et al. 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1898-6) - Accuracy, robustness and scalability of dimensionality reduction methods for single-cell RNA-seq analysis - *"we provide important guidelines for choosing dimensionality reduction methods for scRNA-seq data analysis"*

* **Dimensionality reduction:** [Heiser and Lau 2019](https://www.biorxiv.org/content/10.1101/684340v1) - A quantitative framework for evaluating single-cell data structure preservation by dimensionality reduction techniques - *"Using discrete and continuous scRNA-seq datasets, we find that input cell distribution and method parameters are largely determinant of global, local, and organizational data structure preservation by eleven published dimensionality reduction methods"*

* **Dimensionality reduction (PCA):** [Tsuyuzaki et al. 2020](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1900-3) - Benchmarking principal component analysis for large-scale single-cell RNA-sequencing - *"PCA algorithms based on Krylov subspace and randomized singular value decomposition are fast, memory-efficient, and more accurate than the other algorithms."*

* **Dimensionality reduction:** [Becht et al. 2018](https://www.nature.com/articles/nbt.4314) - Dimensionality reduction for visualizing single-cell data using UMAP - *"UMAP provides the fastest run times, highest reproducibility and the most meaningful organization of cell clusters."*

* **Similarity metrics for clustering:** [Kim et al. 2019](https://www.ncbi.nlm.nih.gov/pubmed/30137247) - Impact of similarity metrics on single-cell RNA-seq data clustering. - *"These findings demonstrate the importance of similarity metrics in clustering scRNA-seq data and highlight Pearson's correlation as a favourable choice. Further comparison on different scRNA-seq library preparation protocols suggests that they may also affect clustering performance"*

* **Clustering:** [Freytag et al. 2018](https://f1000research.com/articles/7-1297) -  Comparison of clustering tools in R for medium-sized 10x Genomics single-cell RNA-sequencing data - *"Seurat outperformed other methods, although performance seems to be dependent on many factors, including the complexity of the studied system. Furthermore, we found that solutions produced by different methods have little in common with each other."*

* **Clustering:** [Duo et al. 2018](https://f1000research.com/articles/7-1141) -  A systematic performance evaluation of clustering methods for single-cell RNA-seq data - *"We found substantial differences in the performance, run time and stability between the methods, with SC3 and Seurat showing the most favorable results. Additionally, we found that consensus clustering typically did not improve the performance compared to the best of the combined methods, but that several of the top-performing methods already perform some type of consensus clustering"*

* **Clustering:** [Menon 2017](https://academic.oup.com/bfg/article/17/4/240/4728639)

* **Clustering:** [Krzak et al. 2019](https://www.frontiersin.org/articles/10.3389/fgene.2019.01253/full?&utm_source=Email_to_authors_&utm_medium=Email&utm_content=T1_11.5e1_author&utm_campaign=Email_publication&field=&journalName=Frontiers_in_Genetics&id=486077) - Benchmark and Parameter Sensitivity Analysis of Single-Cell RNA Sequencing Clustering Methods - *"great variability in the performance of the models is strongly attributed to the choice of the user-specific parameter settings. We describe several tendencies in the performance attributed to their modes of usage and different types of datasets, and identify which methods are strongly affected by data dimensionality in terms of computational time"*

* **Clustering:** [Qi et al. 2019](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbz062/5528236)

* **Cell type assignment:** [Kohler et al. 2019](https://www.biorxiv.org/content/10.1101/653907v1) - Deep learning does not outperform classical machine learning for cell-type annotation - *"deep learning does not outperform classical machine-learning methods in the task. Thus, cell-type prediction based on gene-signature derived cell-type labels is potentially too simplistic a task for complex non-linear methods, which demands better labels of functional single-cell readouts."*

* ?? [Andrews et al. 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6415334/) - False signals induced by single-cell imputation - *"Imputation of single-cell RNA-seq data introduces circularity that can generate false-positive results. Thus, statistical tests applied to imputed data should be treated with care. Additional filtering by effect size can reduce but not fully eliminate these effects. Of the methods we considered, SAVER was the least likely to generate false or irreproducible results, thus should be favoured over alternatives if imputation is necessary."*

* **Cell type assignment:** [Abdelaal et al. 2018](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1795-z) - A comparison of automatic cell identification methods for single-cell RNA sequencing data - *"We find that most classifiers perform well on a variety of datasets with decreased accuracy for complex datasets with overlapping classes or deep annotations. The general-purpose support vector machine classifier has overall the best performance across the different experiments."*

* **Cell type assignment:** [Diaz-Mejia et al. 2019](https://f1000research.com/articles/8-296) -  Evaluation of methods to assign cell type labels to cell clusters from single-cell RNA-sequencing data - *"GSVA was the overall top performer and was more robust in cell type signature subsampling simulations, although different methods performed well using different datasets. METANEIGHBOR and GSVA were the fastest methods. CIBERSORT and METANEIGHBOR were more influenced than the other methods by analyses including only expected cell types"*

* **Cell type assignment:** [Zhao et al. 2019](https://www.biorxiv.org/content/10.1101/827139v1) -  Evaluation of Cell Type Deconvolution R Packages on Single Cell RNA-seq Data - *"Overall, methods such as Seurat, SingleR, CP, RPC and SingleCellNet performed well, with Seurat being the best at annotating major cell types. Also, Seurat, SingleR and CP are more robust against down-sampling. However, Seurat does have a major drawback at predicting rare cell populations, and it is suboptimal at differentiating cell types that are highly similar to each other, while SingleR and CP are much better in these aspects."*

* **Cell type assignment:** [Huang et al. 2019](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbz096/5593804) -  Evaluation of single-cell classifiers for single-cell RNA sequencing data sets - *"Results showed that Seurat based on random forest, SingleR based on correlation analysis and CaSTLe based on XGBoost performed better than others. A simple ensemble voting of all tools can improve the predictive accuracy. Under nonideal situations, such as small-sized and class-imbalanced reference data sets, tools based on cluster-level similarities have superior performance. However, even with the function of assigning ‘unassigned’ labels, it is still challenging to catch novel cell types by solely using any of the single-cell classifiers"*

* **Cell type proportions:** [Cobos et al. 2020](https://www.biorxiv.org/content/10.1101/2020.01.10.897116v1) - Comprehensive benchmarking of computational deconvolution of transcriptomics data - *"overall, single-cell methods have comparable performance to the best performing bulk methods and bulk methods based on semi-supervised approaches showed higher error and lower correlation values between the computed and the expected proportions. Moreover, failure to include cell types in the reference that are present in a mixture always led to substantially worse results, regardless of any of the previous choices. "*

* **Marker selection:** [Gilbert and Vargo 2019](https://www.biorxiv.org/content/10.1101/679761v1) - Comparison of marker selection methods for high throughput scRNA-seq data - *"most existing marker selection methods show similar performance on experimental scRNA-seq data; thus, the speed of the algorithm is the most important consid-eration for large data sets. With this in mind, we introduce RANKCORR, a fast marker selection method with strong mathematical underpinnings that takes a step towards sensible multi-class marker selection."*

* **False discovery control:** [Korthauer et al. 2019](https://www.ncbi.nlm.nih.gov/pubmed/?term=A+practical+guide+to+methods+controlling+false+discoveries+in+computational+biology)

* **Batch effect correction:** [Tran et al. 2020](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1850-9) - A benchmark of batch-effect correction methods for single-cell RNA sequencing data - *"Harmony, LIGER, and Seurat 3 are the recommended methods for batch integration. Due to its significantly shorter runtime, Harmony is recommended as the first method to try, with the other methods as viable alternatives."* (NB also suggested scMerge for DE)

* **Differential expression:** [Jaakkola et al. 2017](https://www.ncbi.nlm.nih.gov/pubmed/27373736) - *"We compared five statistical methods to detect differentially expressed genes between two distinct single-cell populations. Currently, it remains unclear whether differential expression methods developed originally for conventional bulk RNA-seq data can also be applied to single-cell RNA-seq data analysis. Our results in three diverse comparison settings showed marked differences between the different methods... They, however, did not reveal systematic benefits of the currently available single-cell-specific methods."*

* **Differential expression:** [Soneson & Robinson 2018](https://www.ncbi.nlm.nih.gov/pubmed/29481549) - *"found considerable differences in the number and characteristics of the genes that are called differentially expressed. Prefiltering of lowly expressed genes has important effects, particularly for some of the methods developed for bulk RNA-seq data analysis. However, we found that bulk RNA-seq analysis methods do not generally perform worse than those developed specifically for scRNA-seq"*

* **Differential expression:** [Dal Monil et al. 2017](https://www.ncbi.nlm.nih.gov/pubmed/28588607) - *"We observed marked differences between the selected methods in terms of precision and recall, the number of detected differentially expressed genes and the overall performance. Globally, the results obtained in our study suggest that is difficult to identify a best performing tool and that efforts are needed to improve the methodologies for single-cell RNA-sequencing data analysis and gain better accuracy of results."*

* **Differential expression:** [Wang et al. 2019](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2599-6) - Comparative analysis of differential gene expression analysis tools for single-cell RNA sequencing data - *"current methods designed for scRNAseq data do not tend to show better performance compared to methods designed for bulk RNAseq data. Data multimodality and abundance of zero read counts are the main characteristics of scRNAseq data, which play important roles in the performance of differential gene expression analysis methods and need to be considered in terms of the development of new methods."*

* **Trajectory analysis:** [Saelens et al. 2019](https://www.nature.com/articles/s41587-019-0071-9) - *"Our results highlight the complementarity of existing tools, and that the choice of method should depend mostly on the dataset dimensions and trajectory topology."*

* **Isoform quantification:** [Westoby et al. 2018](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1571-5) - Simulation-based benchmarking of isoform quantification in single-cell RNA-seq - *"Performance is generally good for simulated data based on SMARTer and SMART-seq2 data. The reduction in performance compared with bulk RNA-seq is small"*

* **SNVs:** [Liu et al. 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1863-4) - Systematic comparative analysis of single-nucleotide variant detection methods from single-cell RNA sequencing data - *"We recommend SAMtools, Strelka2, FreeBayes, or CTAT, depending on the specific conditions of usage"*

* **Gene regulatory networks:** [Pratapa et al. 2020](https://www.nature.com/articles/s41592-019-0690-6) - Benchmarking algorithms for gene regulatory network inference from single-cell transcriptomic data - *"We find that the area under the precision-recall curve and early precision of the algorithms are moderate. The methods are better in recovering interactions in synthetic networks than Boolean models. The algorithms with the best early precision values for Boolean models also perform well on experimental datasets. Techniques that do not require pseudotime-ordered cells are generally more accurate"*

* **Gene regulatory networks:** [Chen and Mar 2018](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2217-z) - Evaluating methods of inferring gene regulatory networks highlights their lack of performance for single cell gene expression data - *"most of these assessed network methods are not able to predict network structures from single cell expression data accurately, even if they are specifically developed for single cell methods. Also, single cell methods, which usually depend on more elaborative algorithms, in general have less similarity to each other in the sets of edges detected"*

* **scATAC-seq:** [Chen et al. 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1854-5) - Assessment of computational methods for the analysis of single-cell ATAC-seq data - *"Despite variation across methods and datasets, SnapATAC, Cusanovich2018, and cisTopic outperform other methods in separating cell populations of different coverages and noise levels in both synthetic and real datasets. Notably, SnapATAC is the only method able to analyze a large dataset (> 80,000 cells)."*

* **single cell potency:** [Shi et al. 2018](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bby093/5115275) -  Quantifying Waddington’s epigenetic landscape: a comparison of single-cell potency measures - *"integration of RNA-Seq data with a protein interaction network dramatically improves the robustness and reliability of single-cell potency estimates. Thus, this study provides novel systems-biological insight into cellular potency and may provide a foundation for improved models of differentiation potency with far-reaching implications for the discovery of novel stem cell or progenitor cell phenotypes."*

* **Visualisation:** [Cakir et al. 2020](https://www.biorxiv.org/content/10.1101/2020.01.24.918342v1) - Comparison of visualisation tools for single-cell RNAseq data - *"We review and compare several of the currently available analysis and visualisation tools and benchmark those that allow to visualize the scRNAseq data on the web and share it with others. To address the problem of format compatibility for most visualisation tools, we have also developed a user-friendly R package, sceasy, which allows users to convert their own scRNAseq datasets into a specific data format for visualisation."*

* **Differential state analyses:** [Crowell et al. 2019](https://www.biorxiv.org/content/10.1101/713412v1) - On the discovery of population-specific state transitions from multi-sample multi-condition single-cell RNA sequencing data - *"In this work, we surveyed the methods available to perform cross-condition differential state analyses, including cell-level mixed models and methods based on aggregated “pseudobulk” data"*

* **Association (gene-gene or cell-cell relationships)** - [Skinnider et al. 2019](https://www.nature.com/articles/s41592-019-0372-4) - Evaluating measures of association for single-cell transcriptomics - *"Our analysis provides data-driven guidance for gene and cell network analysis in single-cell transcriptomics."*

* **Highly variable genes:** - [Yip et al. 2019](https://www.ncbi.nlm.nih.gov/pubmed/29481632) - Evaluation of tools for highly variable gene discovery from single-cell RNA-seq data. - *"we compare seven HVG methods from six software packages, including BASiCS, Brennecke, scLVM, scran, scVEGs and Seurat. Our results demonstrate that reproducibility in HVG analysis requires a larger sample size than DEG analysis"*

* [Tian et al. 2019](https://www.nature.com/articles/s41592-019-0425-8) - Benchmarking single cell RNA-sequencing analysis pipelines using mixture control experiments - *"Here, we generated a realistic benchmark experiment that included single cells and admixtures of cells or RNA to create ‘pseudo cells’ from up to five distinct cancer cell lines. In total, 14 datasets were generated using both droplet and plate-based scRNA-seq protocols. We compared 3,913 combinations of data analysis methods for tasks ranging from normalization and imputation to clustering, trajectory analysis and data integration. Evaluation revealed pipelines suited to different types of data for different tasks."*

* A [Medium post](https://medium.com/@HeleneOMICtools/your-top-3-single-cell-rna-sequencing-analysis-tools-221b65fbc57e) ranking various scRNASeq analysis tools.

* A Twitter quest for benchmarking papers:

<blockquote class="twitter-tweet"><p lang="en" dir="ltr">Dear <a href="https://twitter.com/hashtag/academictwitter?src=hash&amp;ref_src=twsrc%5Etfw">#academictwitter</a>, I&#39;m looking for a list of benchmarks for single cell RNA-seq data analysis. Generally, so normalization, DE, clustering, cell assignment, etc. I&#39;m OFC aware of the &quot;Methods comparisons&quot; section in <a href="https://twitter.com/seandavis12?ref_src=twsrc%5Etfw">@seandavis12</a>&#39;s awesome list (<a href="https://t.co/lT7DivJdAZ">https://t.co/lT7DivJdAZ</a>) ..</p>&mdash; Mark Robinson (@markrobinsonca) <a href="https://twitter.com/markrobinsonca/status/1220048175320440832?ref_src=twsrc%5Etfw">January 22, 2020</a></blockquote> <script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script> 

* And a [Google Docs spreadsheet](https://docs.google.com/spreadsheets/d/1Gqn0eZ8oiNh8-9ovyh4D_ZsoWcKYoKTJ0C1AfrJrdRA/edit#gid=0) summarising the results:

<blockquote class="twitter-tweet"><p lang="en" dir="ltr">OK, my current list is 31 scRNA-seq (analysis method) benchmarks:<a href="https://t.co/gfan8Kupcp">https://t.co/gfan8Kupcp</a><br><br>I consolidated <a href="https://twitter.com/seandavis12?ref_src=twsrc%5Etfw">@seandavis12</a>&#39;s <a href="https://twitter.com/_lazappi_?ref_src=twsrc%5Etfw">@_lazappi_</a>&#39;s <a href="https://twitter.com/JMA_Data?ref_src=twsrc%5Etfw">@JMA_Data</a>&#39;s and went through the <a href="https://twitter.com/GenomeBiology?ref_src=twsrc%5Etfw">@GenomeBiology</a> benchmark issue.<br><br>Any others I missed? <a href="https://t.co/AEpHCWcYsS">https://t.co/AEpHCWcYsS</a></p>&mdash; Mark Robinson (@markrobinsonca) <a href="https://twitter.com/markrobinsonca/status/1221903274678333441?ref_src=twsrc%5Etfw">January 27, 2020</a></blockquote> <script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script> 

* Benchmarking: [Machine Learning and Statistical Methods for Clustering Single Cell RNA-sequencing Data](https://github.com/kuanglab/single-cell-review)


* **Conquer:** [Soneson & Robinson 2018](https://www.ncbi.nlm.nih.gov/pubmed/29481549) - in their DE benchmarking paper they also *"present conquer, a repository of consistently processed, analysis-ready public scRNA-seq data sets that is aimed at simplifying method evaluation and reanalysis of published results"*


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
