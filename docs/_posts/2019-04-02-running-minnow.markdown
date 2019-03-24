---
title:  "Running Minnow"
date:   2019-02-04 15:04:23
categories: [tutorial]
tags: [minnow]
---
## Installation 


Minnow is written in C++14 and tested in a ubuntu server, please let us know if you have difficulty compiling it in your own machine.

```bash
git clone https://github.com/COMBINE-lab/minnow.git
cd minnow
mkdir build
cd build
cmake ..
make
```
Given the above steps run without error, the binary will be stored within the `build/src`.

## How to run ?

### Simple runs

Minnow has many options some are basic and some are advanced, the most basic one would require to have a 
transcript to gene map (a `tsv` file) and the reference file in addition to that it needs three different  
files to generate data from: `quants_mat.csv`, a gene to cell count file (comma separated), `quants_mat_rows.txt`, 
and `quants_mat_cols.txt`. To list the required files, a complete file tree would look like following, 

```
+-- reference
|   +-- transcript.fasta
|   +-- t2g.tsv
+-- experiment
|   +-- quants_mat.csv
|   +-- quants_mat_rows.txt
|   +-- quants_mat_cols.txt 
```
The t2g.tsv file contains the transcript to gene information. The csv file is the main matrix file and assumed to be 
in gene to cell format and comma separated of dimension `m * n` where there are `m` genes and `n` cells. In the given setting 
one must use a `quants_mat_rows.txt` with `m` lines and `quants_mat_cols.txt` with `n` lines. Not having this convention
might lead to failure and segmentation faults.

The matrix can be generated in many ways, like using splatter or any generative model user might want to implement. Given 
the above format minnow will generate the reads keeping the matrix intact. If not specified otherwise minnow generates 
`100,000` molecules per cell after PCR, but this can be customized by passing the required flag with the program.  

Below we will go over the options that can be provided to minnow in order to get different functionality. Minnow 
primarily works in two different modes. One where you can provide minnow with a simple matrix and it will with gene 
and cell information according to above format. Additional files can be obtained from https://github.com/COMBINE-lab/minnow

We have also linked different files in the main website. 

## Basic mode
In basic mode minnow takes a set of genes to assign them to whatever gene names that are provided to it. Please note that 
if you are using real gene names that are present in the reference then please don't use this mode. Since in that case 
this mode would reassign those real gene names by it's own algorithm and the purpose is defeated. 

One the genes are assigned, minnow uses the `weibull` distribution to assign reads to the candidate transcripts. A typical 
command in normal mode might look as follows, 

``` console
src/minnow simulate --normal-mode --g2t /mnt/scratch6/avi/data/cgat/references/metadata/hg_t2g.tsv -m example_data/normal_data_100_Cells_50K_Genes/ --PCR 7 -r /mnt/scratch6/avi/data/cgat/references/txome/hg_transcriptome.fasta -e 0.01 -p 25 -o /mnt/scratch1/avi/minnow/benchmark/ismb_19_runs/splatter_test/100_50K_Genes_pbmc/defaultnorm --useWeibull

```
Here the folder `example_data/normal_data_100_Cells_50K_Genes/` contains the matrix and the corresponding gene names 
and cell names in the format shown above. 

## Splatter-mode
The `splatter-mode` is very similar to the `normal-mode` except, the provided matrix comes from `splatter`. 

In either mode user can make use of the `de-Bruijn` graph provided in the github repo in order to generate realistic reads 
in the following way 

```console
src/minnow simulate --splatter-mode --g2t /mnt/scratch7/hirak/mm/mm_t2g.tsv -m /mnt/scratch7/hirak/descend_matrices/parametric_sim/ --PCR 7 -r /mnt/scratch7/hirak/mm/mm_transcriptome_fixed.dedup.fasta  -e 0.01 -p 25 -o /mnt/scratch7/hirak/descend_matrices/parametric_sim/dbg2  --dbg --gfa /mnt/scratch7/hirak/mm/mm.gfa --bfh /mnt/scratch7/hirak/mm/bfh.txt --uniq /mnt/scratch7/hirak/mm/mm_101_stranded_gene_uniqueness.txt

```


Alternatively a `gfa` file can be provided to represent `de-Bruijn` file.

## Alevin mode 
Without illumina model

```console
src/minnow simulate --splatter-mode --g2t /mnt/scratch7/hirak/mm/mm_t2g.tsv -m /mnt/scratch7/hirak/descend_matrices/parametric_sim/ --PCR 7 -r /mnt/scratch7/hirak/mm/mm_transcriptome_fixed.dedup.fasta  -e 0.01 -p 25 -o /mnt/scratch7/hirak/descend_matrices/parametric_sim/dbg2  --dbg --gfa /mnt/scratch7/hirak/mm/mm.gfa --bfh /mnt/scratch7/hirak/mm/bfh.txt --uniq /mnt/scratch7/hirak/mm/mm_101_stranded_gene_uniqueness.txt

```

With illumina model

```console
src/minnow simulate --splatter-mode --g2t /mnt/scratch7/hirak/mm/mm_t2g.tsv -m /mnt/scratch7/hirak/descend_matrices/parametric_sim/ --PCR 7 -r /mnt/scratch7/hirak/mm/mm_transcriptome_fixed.dedup.fasta  -e 0.01 -p 25 -o /mnt/scratch7/hirak/descend_matrices/parametric_sim/dbg2  --dbg --gfa /mnt/scratch7/hirak/mm/mm.gfa --bfh /mnt/scratch7/hirak/mm/bfh.txt --uniq /mnt/scratch7/hirak/mm/mm_101_stranded_gene_uniqueness.txt --illum ../data/illumina_error_models/ill100v5_mate2

```

