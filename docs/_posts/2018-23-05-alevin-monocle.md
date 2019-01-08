---
title:  "How to Use Alevin with Monocle"
date:   2018-05-23 15:04:23
categories: [tutorial]
tags: [alevin]
---
## Alevin-Monocle Connection





```R
library(monocle)
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
    matrix <- as.matrix(read.csv( matrix.loc, header=FALSE))
    matrix <- t(matrix[,1:ncol(matrix)-1])

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
base.path <- "/mnt/scratch6/avi/data/cgat/tom_pipe/alevin/neurons_900-wo_whitelist-20-0/"
alv.data <- ReadAlevin(base.path)
```


```R
# and we are good to go !! Cells after this has been taken from Monocle tutrial:
# http://cole-trapnell-lab.github.io/monocle-release/docs/#constructing-single-cell-trajectories
# Below lines are for example purposes only and could be suboptimal. We recommend
# checking out Monocle tool for more detailed tutorial of the downstream analysis.
```


```R

```


```R
alv.monocle.data <- newCellDataSet(alv.data,
                                   lowerDetectionLimit = 0.5,
                                   expressionFamily = negbinomial.size())
```

    Warning message in newCellDataSet(alv.data, lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size()):
    “Warning: featureData must contain a column verbatim named 'gene_short_name' for certain functions”Warning message in newCellDataSet(alv.data, lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size()):
    “Warning: featureData must contain a column verbatim named 'gene_short_name' for certain functions”Warning message in newCellDataSet(alv.data, lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size()):
    “Warning: featureData must contain a column verbatim named 'gene_short_name' for certain functions”


```R
alv.monocle.data
```


    CellDataSet (storageMode: environment)
    assayData: 52325 features, 931 samples 
      element names: exprs 
    protocolData: none
    phenoData
      sampleNames: TACTCGCAGATAGTCA AGGGAGTTCCTCAATT ... GACTGCGAGACACGAC
        (931 total)
      varLabels: Size_Factor
      varMetadata: labelDescription
    featureData: none
    experimentData: use 'experimentData(object)'
    Annotation:  



```R
print("estimating Size/Dispersion Input -> ")
alv.monocle.data <- estimateSizeFactors(alv.monocle.data)
alv.monocle.data <- estimateDispersions(alv.monocle.data)

#adding genes present and num-cells expressed features
print("Detect genes -> ")
alv.monocle.data <- detectGenes(alv.monocle.data, min_expr = 0.1)

#making a list of expressed genes
print("getting a list of expressed genes in at least 10 cells -> ")
# expressed_genes.tap73 <- row.names(subset(fData(alv.data), num_cells_expressed >= 2))

#count number of mRna in each cell
print("counting mRNA in each cell ")
pData(alv.monocle.data)$Total_mRNAs <- Matrix::colSums(exprs(alv.monocle.data))
```

    [1] "estimating Size/Dispersion Input -> "


    Warning message in log(ifelse(y == 0, 1, y/mu)):
    “NaNs produced”Warning message:
    “step size truncated due to divergence”Warning message in log(ifelse(y == 0, 1, y/mu)):
    “NaNs produced”Removing 195 outliers


    [1] "Detect genes -> "
    [1] "getting a list of expressed genes in at least 10 cells -> "
    [1] "counting mRNA in each cell "



```R
alv.monocle.data <- reduceDimension(alv.monocle.data, max_components = 2, method = 'DDRTree', cores=5)
```


```R
alv.monocle.data <- orderCells(alv.monocle.data)
```


```R
plot_cell_trajectory(alv.monocle.data)
```




![png](../../images/output_10_1.png)
