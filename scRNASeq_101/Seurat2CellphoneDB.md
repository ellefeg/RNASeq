# Seurat2CellphoneDB

This code will show you how to extract data from a mouse Seurat object for use in CellphoneDB. Note that if you are running mouse data with gene symbols you will need to convert your genes to human Ensembl IDs. 

For code to run CellphoneDB once your input files are ready, see the very end of this document.

**Pitfalls**

1. If CellphoneDB tells you the cells in your metadata don'tt match the cells in your count data, check the symbols in your names. I had problems where the underscores in my metadata got converted to full stops. For some reason it would only work if I converted all underscores in the count data to full stops (not the other way around!)
2. In theory CellphoneDB can accept (human) gene symbols, you just have to specify this in your command. I couldn't get it to work though, so I always convert to ENSEMBL IDs. Regardless, you will get CellphoneDB output with gene symbols not ENSEMBL IDs

**Load the libraries**

```{r}
#library(EWCE)
library(tibble)
library(biomaRt)
library(tidyr)
library(dplyr)
library(Seurat)
library("org.Hs.eg.db") # remember to install it if you don't have it already
```

**Load the data**

How to access `biomartMouse2HumanOrthos.txt` annotations:
1. Go to EnsemblBiomart
2. Dataset: Mouse genes (GRCm38.p6)
3. Filters: Orthologous human genes only
4. Attributes: gene stable ID, human gene stable ID, human gene name, gene name
This gives a table with 26708 rows

```{r}
outdir <- "/path/to/outdir/"
alldata <- readRDS("/path/to/indir/SeuratObject.RDS")
mouse2human <- read.delim("/Users/uqlgrice/Documents/IMB/Research/Data/orthologues/mouse2human/biomartMouse2HumanOrthos.txt")
```

**In your mouse-to-human orthologue conversion, remove many-to-\* or \*-to-many hits**

```{r}
mouse2human_uni <- unique(mouse2human[,3:4])
mouse2human_uni <- mouse2human_uni[mouse2human_uni$Gene.name %in% rownames(alldata@assays$RNA@counts),] #only keep genes in seurat obj
human_keep <- names(which(table(mouse2human_uni$Human.gene.name) == 1)) #only want human genes which appear once (i.e. not many-to-*)
mouse_keep <- names(which(table(mouse2human_uni$Gene.name) == 1)) #only want mouse genes which appear once (i.e. not many-to-*)
mouse2human_uni <- mouse2human_uni[mouse2human_uni$Human.gene.name %in% human_keep,] #remove dups
mouse2human_uni <- mouse2human_uni[mouse2human_uni$Gene.name %in% mouse_keep,] #remove dups
colnames(mouse2human_uni) <- c("HumanID", "MouseID")
rownames(mouse2human_uni) <- mouse2human_uni$MouseID
mouse2human_uni$MouseID <- NULL
```

**If you have multiple timepoints, subset the data**

It does not make sense to include cells that don't temporally co-exist.

```{r}
d1_cells <- names(which(alldata$time == "d1"))
d1_data <- subset(alldata, cells = d1_cells)
# repeat for other timepoints
```

**Define a function to convert mouse gene IDs to human**
```{r}
# convertMouseGeneList <- function(x){
# require("biomaRt")
# human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
# humanx <- unique(genesV2[, 2])
# # Print the first 6 genes found to the screen
# #print(head(humanx))
# return(humanx)
# }
```

**For each timepoint, convert the counts**
```{r}
counts_d1 <- d1_data@assays$RNA@data
counts_d1 <- counts_d1[rowSums(counts_d1) > 0,] #remove genes with 0 counts
m2h_d1 <- mouse2human_uni[rownames(mouse2human_uni) %in% rownames(counts_d1), , drop = FALSE] #make sure mouse conversion and counts have same genes
counts_d1 <- counts_d1[rownames(counts_d1) %in% rownames(m2h_d1), , drop = FALSE] #make sure mouse conversion and counts have same genes
m2h_d1 <- m2h_d1[rownames(counts_d1), , drop = FALSE]  #make sure mouse conversion and counts are in the same order
counts_d1 <- cbind(m2h_d1, counts_d1)
rownames(counts_d1) <- counts_d1$HumanID
counts_d1$HumanID <- NULL
colnames(counts_d1) <- gsub("_", ".", colnames(counts_d1)) #Important! For some reason my metadata got converted to having . not _ and your count data needs to match. For me, t didn't work the other way around i.e. converting metadata from . to _
#counts_d1_2 = Matrix::Matrix(as.matrix(counts_d1),sparse=TRUE)
write.table(counts_d1, file = paste0(outdir, "d1", "_counts.txt"), sep = "\t", quote = FALSE, col.names = NA)

# Convert to Ensembl IDs (CPDB needs this, but it will output the results as gene symbols
# Technically it is possible to use symbols and tell CPDB that you are doing so, but it failed for me
library("org.Hs.eg.db") # remember to install it if you don't have it already
symbols <- mapIds(org.Hs.eg.db, keys = rownames(counts_d1), keytype = "SYMBOL", column="ENSEMBL")
counts_d1 <- as.data.frame(counts_d1)
counts_d1$Gene <- unname(symbols)
counts_d1_ens <- counts_d1
counts_d1_ens <- counts_d1_ens[!(is.na(counts_d1_ens$Gene)),]
# counts_d1_ens %>%
#   group_by(Gene) %>%
#   summarise_all(sum) %>%
#   data.frame() -> counts_d1_ens_consolidate #merge duplicates
counts_d1_ens_consolidate <- counts_d1_ens
rownames(counts_d1_ens_consolidate) <- counts_d1_ens_consolidate$Gene
counts_d1_ens_consolidate$Gene <- NULL

meta_d1 <- d1_data[[c("cellclust", "labels")]]
# remove NAs
meta_d1 <- meta_d1[!(is.na(meta_d1$cellclust) | is.na(meta_d1$labels)),]
# extract specific labels (cellclust)
meta_spec_d1 <- meta_d1[,"cellclust", drop = FALSE]
rownames(meta_spec_d1) <- gsub("_", ".", rownames(meta_spec_d1))
# extract broad labels (labels)
meta_broad_d1 <- meta_d1[,"labels", drop = FALSE]
rownames(meta_broad_d1) <- gsub("_", ".", rownames(meta_broad_d1))

# we removed some NAs so make sure to filter the ensembl counts too
rownames(meta_d1) <- gsub("_", ".", rownames(meta_d1))
counts_d1_ens_consolidate_v2 <- counts_d1_ens_consolidate[,(colnames(counts_d1_ens_consolidate) %in% rownames(meta_d1))]

write.table(meta_spec_d1, file = paste0(outdir, "d1", "_specific_meta.txt"), sep = "\t", quote = FALSE, col.names = NA)
write.table(meta_broad_d1, file = paste0(outdir, "d1", "_broad_meta.txt"), sep = "\t", quote = FALSE, col.names = NA)
write.table(counts_d1_ens_consolidate_v2, file = paste0(outdir, "d1", "_counts_ensembl.txt"), sep = "\t", quote = FALSE, col.names = NA)
```

**repeat for other timepoints**

**Now run CellphoneDB in the terminal**

```
#cellphonedb method statistical_analysis "$data_path"/d7_broad_meta.txt "$data_path"/d7_counts_ensembl.txt --output-path outdir_d7_broad --threshold 0.2 --subsampling --subsampling-log FALSE --pvalue 0.01

```
