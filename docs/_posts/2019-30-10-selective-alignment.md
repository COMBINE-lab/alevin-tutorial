---
title:  "Selective Alignment"
date:   2019-10-30 15:04:23
categories: [tutorial]
tags: [selective-alignment]
---

## Fast is Good but Fast and accurate is better !

The accuracy of transcript quantification using RNA-seq data depends on many factors, such as the choice of alignment or mapping method and the quantification model being adopted. After investigating the influence of mapping and alignment on the accuracy of transcript quantification in both simulated and experimental data, as well as the effect on subsequent differential expression analysis, we designed selective alignment method. Selective Alignment overcomes the shortcomings of lightweight approaches without incurring the computational cost of traditional alignment. Here we give a short tutorial on how to index your genome and transcriptome to get the most accurate quantification estimates.

### Downloading Reference

We are first going to download the reference transcriptome and genome for salmon index. As an example we are downloading the gencode mouse reference

```python
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.transcripts.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz
```

### Installing Salmon

Although there are mutliple ways to download salmon (ex: binary from github, docker image), we are going to install it throug conda. Assuming a conda environment is already set up, we can install salmon through following command:

```python
conda install --channel bioconda salmon
```

Make sure you have the latest version of salmon (v1.0 as on November 1st, 2019) by using `salmon --version`

### Preparing metadata

Salmon indexing requires the names of the genome targets, which is extractable by using the `grep` command:

```python
grep "^>" <(zcat GRCm38.primary_assembly.genome.fa.gz) > decoys.txt
sed -i -e 's/>//g' decoys.txt
```

Along with the list of decoys salmon also needs the concatenated transcriptome and genome reference file for index.
NOTE: the genome targets (decoys) should come after the transcriptome targets in the reference

```python
cat gencode.vM23.transcripts.fa.gz GRCm38.primary_assembly.genome.fa.gz > gentrome.fa.gz
```

### Salmon Indexin

We have all the ingredients ready for the salmon recipe. We can run salmon indexing step as follows:

```python
salmon index -t gentrome.fa.gz -d decoys.txt -p 12 -i salmon_index
```

### Ipython Notebook
Prefer to read ipython notebook ?
Check out the gist [here](https://gist.github.com/k3yavi/a486647c35158a8296cec543ed9b526f).
