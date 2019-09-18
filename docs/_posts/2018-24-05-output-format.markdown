---
title:  "Understanding Output"
date:   2018-05-24 15:04:23
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

Alternatively, vpolo can be installed using conda as follows:

```python
conda install -c bioconda vpolo
```

If you want the python code to understand the schema of the binary output, it can be found [here](https://github.com/k3yavi/vpolo/blob/master/vpolo/alevin/parser.py)

### Reading Binary format
Alevin's `quants_mat.gz` file can be easily imported to generate Cell v Gene dataframe using the following piece of python3 code:

``` python
from vpolo.alevin import parser
alevin_df = parser.read_quants_bin("<PATH TO ALEVIN output folder>")
```

### Pandas Compatibility:

Alevin can also dump the count-matrix in a human readable -- matrix market exchange (_MTX_) format, if given flag `--dumpMtx`. The new output file `quants_mat.mtx.gz` can be easily read-in using the following python-pandas based function:

```python
from scipy.io import mmread
import gzip
with gzip.open(mat_file) as f:
  alevin_df = mmread("<PATH TO ALEVIN quants_mat.mtx.gz file>").toarray()
```

The above function takes the path of the directory specified in --output directory ( while running the alevin-quantification tool ). Note that if alevin was finished successfully then the above specified directory will contain a directory inside of it with the name `alevin` . The function returns a python `dataframe` for the count matrix with Cellular-Barcodes as the index and Gene-id as the header which can be used for the downstream analysis.

### Reading Alevin UMI graphs

When run with the command line flag `--dumpUmiGraph` alevin generates the per cell level Parsimonious Umi Graph (PUGs) into a compressed binary file `cell_umi_graphs.gz`. The file can be prased and generate a per cell level `dot` graph file using following python script:

```python
from vpolo.alevin import parser
parser.read_umi_graph("<PATH to alevin output folder>", "<output folder>")
```

### Reading Alevin's bfh (big freaking hash) file

When run with the command line flag `--dumpBfh` alevin generates the big hash file used by alevin for performing the deduplication, along with the mapped equivalence classes. The file can be parsed by the following python script:

```python
from vpolo.alevin import parser
parser.read_bfh("<PATH to alevin output folder>")
```
