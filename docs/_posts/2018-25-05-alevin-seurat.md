---
title:  "How to Use alevin with Seurat"
date:   2018-05-25 15:04:23
categories: [tutorial]
tags: [alevin]
---
## Alevin-Seurat Connection

```R
# Seurat >3.0 and tximport >1.13.0
library(Seurat)
library(tximport)
```

```R
# path to the output directory of Alevin run
files <- file.path("alevin_quants/alevin/quants_mat.gz")
file.exists(files)
```
TRUE


```R
# Reading in the alevin quants quants
txi <- tximport(files, type="alevin")
```


```R
"and we are good to go !! Cells after this has been taken from Seurat tutorial:
https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html
Below lines are for example purposes only and could be suboptimal. We recommend
checking out Seurat tool for more detailed tutorial of the downstream analysis."
```

```R

```


```R
pbmc <- CreateSeuratObject(counts = txi$counts , min.cells = 3, min.features = 200, project = "10X_PBMC")
```


```R
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
```


```R
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
```


```R
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
```


```R
pbmc <- RunUMAP(pbmc, dims = 1:10)
```


```R
DimPlot(pbmc, reduction = "umap")
```


![png](../../images/output_11_0.png)
