---
title:  "How to Use Alevin with Seurat"
date:   2018-05-22 15:04:23
categories: [tutorial]
tags: [alevin]
---
## Alevin-Seurat Connection

```R
library(Seurat)
library(dplyr)
```

```R
# Parts of the function is taken from Seurat's Read10x parsing function
ReadAlevin <- function( base.path = NULL ){
    if (! dir.exists(base.path )){
      stop("Directory provided does not exist")
    }

    barcode.loc <- paste0( base.path, "alevin/quants_mat_rows.txt" )
    gene.loc <- paste0( base.path, "alevin/quants_mat_cols.txt" )
    matrix.loc <- paste0( base.path, "alevin/quants_mat.csv" )
    if (!file.exists( barcode.loc )){
      stop("Barcode file missing")
    }
    if (! file.exists(gene.loc) ){
      stop("Gene name file missing")
    }
    if (! file.exists(matrix.loc )){
      stop("Expression matrix file missing")
    }
    matrix <- t(as.matrix(read.csv( matrix.loc, header=FALSE))[,-1])

    cell.names <- readLines( barcode.loc )
    gene.names <- readLines( gene.loc )

    colnames(matrix) <- cell.names
    rownames(matrix) <- gene.names
    matrix[is.na(matrix)] <- 0
    return(matrix)
}
```


```R
# path to the given output directory when ran Alevin
base.path <- "/mnt/scratch6/avi/data/cgat/tom_pipe/alevin/t_3k-wo_whitelist-20-0/"
alv.data <- ReadAlevin(base.path)
```


```R
"and we are good to go !! Cells after this has been taken from Seurat tutrial:
https://satijalab.org/seurat/pbmc3k_tutorial.html
Below lines are for example purposes only and could be suboptimal. We recommend
checking out Seurat tool for more detailed tutorial of the downstream analysis."
```

```R

```


```R
pbmc <- CreateSeuratObject(raw.data = alv.data, min.cells = 3, min.genes = 200, project = "10X_PBMC")
```


```R
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
```


```R
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
```


![png](../../images/output_7_0.png)



```R
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
```

    [1] "Scaling data matrix"
      |======================================================================| 100%
    [1] "PC1"
    [1] "ENSG00000070081.15" "ENSG00000256128.5"  "ENSG00000166866.12"
    [4] "ENSG00000199691.1"  "ENSG00000229508.2" 
    [1] ""
    [1] "ENSG00000134684.10" "ENSG00000166887.15" "ENSG00000165806.19"
    [4] "ENSG00000172731.13" "ENSG00000251586.1" 
    [1] ""
    [1] ""
    [1] "PC2"
    [1] "ENSG00000106460.18" "ENSG00000184162.14" "ENSG00000199691.1" 
    [4] "ENSG00000163634.11" "ENSG00000211843.1" 
    [1] ""
    [1] "ENSG00000126768.12" "ENSG00000233430.3"  "ENSG00000201622.1" 
    [4] "ENSG00000230076.1"  "ENSG00000253322.1" 
    [1] ""
    [1] ""
    [1] "PC3"
    [1] "ENSG00000165806.19" "ENSG00000266984.1"  "ENSG00000256389.1" 
    [4] "ENSG00000267221.2"  "ENSG00000185085.2" 
    [1] ""
    [1] "ENSG00000256128.5" "ENSG00000242028.6" "ENSG00000240499.7"
    [4] "ENSG00000217514.1" "ENSG00000136490.8"
    [1] ""
    [1] ""
    [1] "PC4"
    [1] "ENSG00000237149.5" "ENSG00000223335.1" "ENSG00000225327.3"
    [4] "ENSG00000263593.1" "ENSG00000250342.1"
    [1] ""
    [1] "ENSG00000240499.7"  "ENSG00000126768.12" "ENSG00000172731.13"
    [4] "ENSG00000242028.6"  "ENSG00000225022.1" 
    [1] ""
    [1] ""
    [1] "PC5"
    [1] "ENSG00000235303.1"  "ENSG00000144043.11" "ENSG00000136490.8" 
    [4] "ENSG00000242371.1"  "ENSG00000283633.1" 
    [1] ""
    [1] "ENSG00000166664.13" "ENSG00000225479.1"  "ENSG00000260594.1" 
    [4] "ENSG00000244159.1"  "ENSG00000255422.1" 
    [1] ""
    [1] ""



```R
pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, 
                     save.SNN = TRUE)
```


```R
pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)
```


```R
TSNEPlot(object = pbmc)
```


![png](../../images/output_11_0.png)