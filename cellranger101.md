Laura Grice

_Adapted from 20191213_HowToScrnaseq.Rmd_

## Appendix 2: Cellranger 101

Cellranger is software released by 10x Genomics for dealing with scRNAseq output. It contains tools for various stages of the pre-processing and analysis pathway, and overlaps somewhat with Seurat. Some of the Cellranger tools include (but are not limited to):

* `cellranger mkfastq` demultiplexes the raw base call files produced by Illumina, and produces fastq files. This tool is actually just a wrapper of Illumina's `bcl2fastq` tool.
* `cellranger mkref` takes a reference genome fasta file and gtf file and produces a genome reference for read mapping
* `cellranger count` performs alignment, filtering, barcode counting, UMI counting. It can also work out clusters and perform gene expression analysis, but here we will use Seurat to do this.

10x Genomics provides a [download page](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest) to download prebuilt human or mouse reference datasets. The "build steps" for each prebuilt set describe the process: GTF and fasta files were downloaded from Ensembl, some attributes are removed (e.g. VDJ genes and pseudogenes), and the `mkref` command is run. If you want to modify a dataset (for example, add in the zs-green sequence), you need to edit the original fasta and gtf files to add in the fasta and gtf information for your sequences. Here, each fasta sequence should represent a SCAFFOLD (e.g. an entire chromosome) and the GTF file is responsible for providing information about the various genomic coordinates of your scaffold. The downstream `cellranger count` command will only count features that are annotated as `exon` type, though. To add zs-green into the mouse genome, I just added a fasta sequence of the single-exon zs-green sequence, and a GTF line where the only feature was one exon, running from position 1 to the end of the fasta sequence.

For more about how I built the new mouse + zsgreen reference genome, see the `SCI_7dpi_notes` page on LabArchives.

**Output**

The output of the `cellranger count` analysis will be contained in a folder called, surprisingly, `counts`, which in turn contains a subfolder called `outs`. It contains many different output files, including:

* `web_summary.html` - this is a useful file which opens in a web browser, showing summary statistics for the mapping
* `metrics_summary.csv` - a tabular format of the above summary stats
* `possorted_genome_bam.bam` and `possorted_genome_bam.bam.bai` - a BAM file (and index file) of reads aligned to genome/transcription, with barcode info
* `filtered_feature_bc_matrix` and `...matrix_h5.h5` - filtered feature-barcode matrices for ONLY cellular barcodes in MEX format (and *.h5 is the same thing in HDF5 format). `filtered_feature_bc_matrix` is the file input into Seurat. This contains only detected CELLULAR barcodes.
* `raw_feature_bc_matrices` and `...matrix_h5.h5` - filtered feature-barcode matrices for ALL BARCODES. MEX and HDF5 formats. This contains EVERY barcode from a fixed list of known barcode sequences, including background and non-cellular barcodes.
* `analysis` - secondary analysis data - dimensionality reduction, cell clustering, differential expression. This is run for you in cellranger, but you can also run it yourself in Seurat. It would be useful to compare the clusterings. The `clustering` folder contains csv files listing the cluster IDs of each cell, for different k-means values. The `diffexp` folder contains csv files listing various DE data, like p-values and log fold changes. The `PCA` folder contains various components information and the `tSNE` folder contains tSNE coordinates.
* `molecule_info.h5` - this is useful if you want to run `cellranger aggr`, which aggregates these `cellranger count` output files, normalises them, and then re-generates an aggregate count matrix. This is useful if you want to combine data from multiple samples into an experiment-wide feature-barcode matrix and analysis.
* `cloupe.cloupe` - this is useful if you want to visualise your data in 10x Genomics' Loupe cell browser

Seurat works with the `filtered_feature_bc_matrix` folder. There are three .gz files in the folder (you do not need to unzip them). The `barcodes.tsv.gz` file contains all the cell barcodes that passed cellranger's filter. The `features.tsv.gz` file contains the ENSEMBL accession numbers, the corresponding gene IDs (aka "symbols" in Seurat), and the type of observation (in this case, every line is `Gene Expression`. The third file, `matrix.mtx.gz`, is the a sparse `gene x cell` matrix showing only non-zero counts per cell per gene, which contains a header and then information that looks like this:

`33538 11769 24825783`

`33509 1 1`

`33506 1 4`

`...`

The first line summarises the dimensionality of the matrix: There are `33538` genes and `11769` cells in the matrix, and `24825783` non-zero entries. In rows 2 and 3 (and so on), the first number is the row (gene) index, the second number is the column (cell) index, and the third number is the count number. Thus, gene 33506 for cell 1 has a count of 4.

Note that most gene-cell combinations have a count value of 0; these have been omitted from this matrix to save disk space (hence, "sparse matrix")
