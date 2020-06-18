# Imputation

scRNASeq data is "sparse" - that is, many genes in the count matrix are given a zero count. While some of these zeroes may be biologically meaningful, they may also be technical artefacts of the sequencing protocol. Data sparsity is a much greater problem for scRNASeq than it is for bulk RNASeq, due to the small amount of starting material (single cells); this problem is also known as the "dropout effect".

Imputation is a process used to "rescue" these zero-count genes; there are a variety of different techniques to do this that fall under the imputation umbrella.

## Which tools to use?

Hou et al. 2020: _"Of the methods considered, **MAGIC, kNN-smoothing and SAVER** were found to outperform the other methods most consistently (Figure 6). However, the performance of methods varied across evaluation criteria, experimental protocols, datasets, and downstream analyses. For example, scVI was one of the top performers in terms of the similarity between the imputed single-cell and bulk expression profiles (Figure 2), but it did not perform among the top in clustering and trajectory analysis (Figures 4,5). SAVER-X performed consistently well using UMI-based method, but less well with non-UMI based methods (Figure 6). While MAGIC was one of the best performer overall (Figure 6), it performed worse than many other methods when identifying differentially expressed genes in hypothesis testing type settings that take into account cell variability (Figure 3)."_

Hou et al. 2020: _"In terms of computation, **MAGIC, DCA and DeepImpute are among the most efficient** methods. **kNN-smoothing, ALRA, bayNorm, scImpute, SAUCIE, scScope exhibit high scalability**. SAVER-X, SAVER and SAVER-X are intermediary while the remaining methods do not scale well for large datasets (Figure S9)."_

[Table 2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1926-6) in Lahnemann et al. 2020 lists many imputation tools, split by their general method of action: 
* A: model-based imputation
* B: data smoothing
* C: data reconstruction, matrix factorization
* T: using external information

## To impute or not impute?

Hou et al. 2020: _"Broadly, we found significant variability in the performance of the methods across evaluation settings. While most scRNA-seq imputation methods recover biological expression observed in bulk RNA-seq data, **the majority of the methods do not improve performance in downstream analyses** compared to no imputation, in particular for clustering and trajectory analysis, and thus should be used with caution. Furthermore, we find that the performance of scRNA-seq imputation methods depends on many factors including the experimental protocol, the sparsity of the data, the number of cells in the dataset, and the magnitude of the effect sizes."_

Hou et al. 2020: _"One important observation is that, while the majority of imputation methods outperformed no imputation in recovering bulk expression (16/18 methods) and log fold changes of individual genes between cell types without considering cell variability within each cell type (13/18 methods) (Figure 2), much fewer methods performed better than no imputation for identifying differentially expressed genes after considering cell variability (1-10/18 depending on the test scenario), clustering cells (5-6/21 methods) or inferring pseudotemporal trajectories (4-11/21 methods). Thus, the current imputation methods as a whole seem to be **most effective for providing a point estimate of the activity of individual genes, and they become less effective when coupled with various downstream analysis tasks**. ... **Thus, imputation could be more helpful for analyzing individual genes rather than cell-to-cell relationship.**"_ i.e. if you want to detect key marker genes / genes of interest, imputation may be helpful, but it probably will not be beneficial to run at the start of every pipeline before all other analyses.

Zhang et al. 2018: _"Simulated datasets and case studies highlight that there are no one method performs the best in all the situations. Some defects of these methods such as scalability, robustness and unavailability in some situations need to be addressed in future studies."_

## [Incomplete] Bibliography

### 2020

A Systematic Evaluation of Single-cell RNA-sequencing Imputation Methods. Wenpin Hou, Zhicheng Ji, Hongkai Ji, Stephanie C. Hicks. bioRxiv **2020**.01.29.925974; doi: https://doi.org/10.1101/2020.01.29.925974 

Lähnemann, D., Köster, J., Szczurek, E. et al. Eleven grand challenges in single-cell data science. Genome Biol 21, 31 (2020). https://doi.org/10.1186/s13059-020-1926-6 (Section: Challenge I: Handling sparsity in single-cell RNA sequencing)

Hyundoo Jeong, Zhandong Liu, PRIME: a probabilistic imputation method to reduce dropout effects in single-cell RNA sequencing, Bioinformatics, **2020**, btaa278, https://doi.org/10.1093/bioinformatics/btaa278

G2S3: a gene graph-based imputation method for single-cell RNA sequencing data. Weimiao Wu, Qile Dai, Yunqing Liu, Xiting Yan, Zuoheng Wang. bioRxiv **2020**.04.01.020586; doi: https://doi.org/10.1101/2020.04.01.020586 

### 2019

Peng, T., Zhu, Q., Yin, P. et al. SCRABBLE: single-cell RNA-seq imputation constrained by bulk RNA-seq data. Genome Biol 20, 88 (2019). https://doi.org/10.1186/s13059-019-1681-8

Mongia, Sengupta, Majumdar. McImpute: Matrix Completion Based Imputation for Single Cell RNA-seq Data. Frontiers in Genetics 2019 Vol 10. DOI: 10.3389/fgene.2019.00009.

### 2018

Andrews TS, Hemberg M. False signals induced by single-cell imputation. F1000Res. 2018;7:1740. Published 2018 Nov 2. doi:10.12688/f1000research.16613.2

Gong, W., Kwak, I., Pota, P. et al. DrImpute: imputing dropout events in single cell RNA sequencing data. BMC Bioinformatics 19, 220 (2018). https://doi.org/10.1186/s12859-018-2226-y

[MAGIC](https://github.com/KrishnaswamyLab/MAGIC) (**2018**)

Li, W.V., Li, J.J. An accurate and robust imputation method scImpute for single-cell RNA-seq data. Nat Commun 9, 997 (**2018**). https://doi.org/10.1038/s41467-018-03405-7 [See also, Github](https://github.com/Vivianstats/scImpute)

L. Zhang and S. Zhang, "Comparison of Computational Methods for Imputing Single-Cell RNA-Sequencing Data," in IEEE/ACM Transactions on Computational Biology and Bioinformatics, vol. 17, no. 2, pp. 376-389, **2018**, doi: 10.1109/TCBB.2018.2848633. [or freely available at: Comparison of computational methods for imputing single-cell RNA-sequencing data Lihua Zhang, Shihua Zhang bioRxiv 241190; doi: https://doi.org/10.1101/241190]

### Unknown

[NetDECODE](https://github.com/shmohammadi86/NetDECODE)

[VIPER](https://github.com/ChenMengjie/VIPER)
