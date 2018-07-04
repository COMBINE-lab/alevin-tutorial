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

## Python Import
Alevin's per-cell level count matrix can be imported directly into a python dataframe. The following python3 dependency is needed and can be installed using pip as follows:

```python
pip3 install vpolo
```

### Reading Binary format
Alevin's `quants_mat.gz` file can be easily imported to generate Cell v Gene dataframe using the following piece of python3 code:

``` python
from vpolo.alevin import parser
alevin_df = parser.read_quants_bin("<PATH TO ALEVIN output folder>")
```

### Pandas Compatibility:

Alevin can also dump the count-matrix in a human readable -- comma-separated-value (_CSV_) format, if given flag `--dumpCsvCounts`. The new output file `quants_mat.csv` can be easily read-in using the following python-pandas based function:

```python
from vpolo.alevin import parser
alevin_df = parser.read_quants_csv("<PATH TO ALEVIN output folder>")
```

The above function takes the path of the directory specified in --output directory ( while running the alevin-quantification tool ). Note that if alevin was finished successfully then the above specified directory will contain a directory inside of it with the name `alevin` . The function returns a python `dataframe` for the count matrix with Cellular-Barcodes as the index and Gene-id as the header which can be used for the downstream analysis.

### Ipython Notebook
Prefer to read ipython notebook ?
Check out the gist [here](https://gist.github.com/k3yavi/c501705ed2d29b12b0d10cf78b3ed001) for the python-code to parse the binary/csv count-matrix fomrat generated through alevin pipeline.
