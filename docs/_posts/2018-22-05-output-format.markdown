---
title:  "Understanding Output"
date:   2018-05-22 15:04:23
categories: [tutorial]
tags: [alevin]
---
## Output Format:

Traditionally single-cell tools dumps the Cell-v-Gene count matrix in various formats. Typical 10x experiment can range form hundreds to tens of thousand of cells -- resulting in huge size of count-matrix. Although, this itself is an open area of research but by default alevin dumps a per-cell level gene-count matrix in a binary-compressed format with the row and column indexes in a separate file.

A typical run of alevin will generate 3 files:

* *quants_mat.gz* -- Compressed count matrix.
* *quants\_mat\_cols.txt* -- Column Header (Gene-ids) of the matrix.
* *quants\_mat\_rows.txt* -- Row Index (CB-ids) of the matrix.

### Pandas Compatibility:

Alevin can also dump the count-matrix in a human readable -- comma-separated-value (_CSV_) format, if given flag `--dumpCsvCounts`. The new output file `quants_mat.csv` can be easily read-in using the following python-pandas based function:

```python
def read_quant_alevin(base):
    alv = pd.read_table(base+"quants_mat.csv", sep=",", header=None)
    names = pd.read_table(base+"quants_mat_rows.txt", header=None)
    genes = pd.read_table(base+"quants_mat_cols.txt", header=None)
    alv.drop([len(alv.columns)-1], axis=1, inplace=True)
    alv.columns = genes[0].values
    alv.index = names[0].values
    alv = alv.loc[:, (alv != 0).any(axis=0)]
    return alv
```

The above function takes the location of the `alevin` directory inside the specified -- output directory ( while running the alevin-quantification tool ). The function returns a python `dataframe` for the count matrix with Cellular-Barcodes as the index and Gene-id as the header which can be used for the downstream analysis.

### Ipython Notebook
Prefer to read ipython notebook ?
Check out the gist [here](https://gist.github.com/k3yavi/c501705ed2d29b12b0d10cf78b3ed001) for the python-code to parse the binary/csv count-matrix fomrat generated through alevin pipeline.