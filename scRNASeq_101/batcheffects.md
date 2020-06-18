# Batch effects

Batch effect correction in scRNA-seq refers to processes undertaken to correct artefacts in scRNA-Seq data when the same cells are processed differently (e.g. major differences such as different protocols, chemistries, sequencing technologies, etc, or seemingly minor changes such as run dates, lab personnel, well/plates, etc.). A major step in scRNASeq analysis is grouping cells into cell types, but if the batch effect is stronger than the cell type / biological effect, it will overwhelm your biological signal and your "cell types" may actually be the result of something irrelevant. 

The batch effect might affect your **whole** dataset (i.e. a genome-wide effect) or it might affect only certain **genes**.

There are two steps in the process:
* Step 1: **detect** the batches
    * **supervised** batch detection is where we know what the batches could be _a priori_, but we don't know if they affect anything
    * **unsupervised** batch detection is where we can find batches without knowing what they are, similar to how SVA is used with RNASeq data. However, [Nikolay Oskolkov does not recommend these methods for scRNASeq](https://towardsdatascience.com/how-to-batch-correct-single-cell-7bad210c7ae1)
* Step 2: **correct** the batches



"Removing (technical) variability means changing the data for individual samples. For expression data, this means that **you should never take individual samples from a batch corrected set for a separate analysis** such as differential expression."

## Detecting batch effects

Nikolay Oskolkov has written a [blog post](https://towardsdatascience.com/how-to-batch-correct-single-cell-7bad210c7ae1) about batch correction. One of his useful visualisations is to "display how much of variation in each Principal Component (PC) is explained by the Donor, Plate and Cell Type". The code to do this analysis can be found [here](https://github.com/NikolayOskolkov/HowToBatchCorrectSingleCell/blob/master/batch_analysis.md).

He also shows how you can "rank all genes influenced by the batch by correlating their expression with the Donor variable." He doesn't remove these genes, but does note that you should keep an eye on them down the track.

## Correcting batch effects

The following three tools all worked well in [Nikolay Oskolkov's trials](https://towardsdatascience.com/how-to-batch-correct-single-cell-7bad210c7ae1):

**Combat**

The simplest of the three tools tested is **Combat**. Combat uses a linear algorithm, and was originally developed to work with microarray data. _"Combat performs batch correction through the Bayesian linear regression and therefore offers a considerable improvement compared to the simple aka Frequentist linear batch removal implemented e.g. in Limma. Besides better performance for small sample sizes, Combat uses the Bayesian shrinkage towards the mean, this implies that correcting each gene we use a “mean” information from other genes while Limma does the correction gene-by-gene independently. Nevertheless, Combat applies a linear batch-effect correction."_ [link](https://towardsdatascience.com/how-to-batch-correct-single-cell-7bad210c7ae1)

As an aside, Combat is the tool that Maika found worked the best with her data, too.

**MNN**

MNN uses a non-linear algorithm. _"Mutual Nearest Neighbors (MNN) algorithm uses a non-linear batch-effect correction inspired by the idea from K-Nearest Neighbors (KNN). The method tries to find most similar cells (mutual neighbors) between the batches, those cells are assumed to belong to the same type and the systematic differences between the MNN cells quantify the strength of the batch, this information is used to scale the rest of the cells in the batches. My general problem with this method is that **I always failed to understand the assumption that the MNN cells should belong to the same type**. My intuition is that a batch can completely mess up the concept of distance between the cells in the batches, so **two cells that are closest neighbors between the batches can actually belong to different cell types but seem to be very similar because of the batch itself**. Nevertheless, the **MNN batch-effect correction seems to work well in practice** and we will see later that it does a good job on our benchmark."_ [link](https://towardsdatascience.com/how-to-batch-correct-single-cell-7bad210c7ae1)

**CCA**

CCA uses a non-linear algorithm. _"Another popular scRNAseq specific batch correction method which **sometimes is seen as an across samples integration is the Seurat Canonical Correlation Analysis (CCA) technique**. ... The data from the batches are projected into a low-dimensional space by maximizing correlation (or covariance) between the data sets from different batches. This implies that the projections of the data sets are correlated but not necessarily aligned, i.e. do not overlap quite well in the low-dimensional space. The latter is solved by applying the Dynamic Time Warping (DTW) technique that locally stretches and squeezes the CCA data projections in order to further align them._"

Tran et al. benchmarked 14 batch effect tools in their 2020 paper. They _"found that each batch-effect removal method has its advantages and limitations, with no clearly superior method. "_ They selected three tools that performed well overall, and gave recommendations for which tools are better in certain scenarios. Their top-three picks were:

**Harmony**

_"Due to its significantly shorter runtime, Harmony is recommended as the first method to try"_. _"Harmony first employs PCA for dimensionality reduction. In the PCA space, Harmony iteratively removes batch effects present"_

**LIGER**

LIGER's point of difference is that it can overcome the problematic assumption that _"differences between datasets are entirely due to technical variations and not of biological origins, thus aiming to remove all of them"_

**Seurat 3**

As noted abouve, Seurat's method uses CCA _"to reduce data dimensionality and capture the most correlated data features to align the data batches"_. CCA is primarily used in Seurat to perform data integration by finding _"correlations across datasets"_. Seurat's method also incorporates MNNs (mutual nearest neighbours) to form anchors between the data.

# References

Tran, H.T.N., Ang, K.S., Chevrier, M. et al. A benchmark of batch-effect correction methods for single-cell RNA sequencing data. Genome Biol 21, 12 (2020). https://doi.org/10.1186/s13059-019-1850-9

BLOG: Nikolay Oskolkov. How to Batch Correct Single Cell. [URL](https://towardsdatascience.com/how-to-batch-correct-single-cell-7bad210c7ae1)

[Systems Biology Group: Batch effect analysis](https://sysbiowiki.soe.ucsc.edu/node/323)
