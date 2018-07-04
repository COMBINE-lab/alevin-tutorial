---
title:  "Setting Up Resources"
date:   2018-05-22 15:04:23
categories: [tutorial]
tags: [alevin]
---
## Example on 10x PBMC data

## Basic Requirements:
* Salmon v0.10 Binary.
* Reference transcriptome.
* Transcript to Gene Mapping.
* Raw _fastq_ files.

### Salmon Binary
There are multiple ways to install alevin:

* Using Conda (Recommended)

```python
conda install salmon
```

* Using prebuild binary (only for linux and Mac-osx):

```python
wget https://github.com/COMBINE-lab/salmon/releases/download/v0.10.0/Salmon-0.10.0_linux_x86_64.tar.gz
tar -xvzf Salmon-0.10.0_linux_x86_64.tar.gz
```

* Compiling from Source:

```python
git clone https://github.com/COMBINE-lab/salmon.git
cd salmon; mkdir build; cd build
cmake ..
make install
```

*Note:* if you find problem installing salmon binary, please raise an issue (following the issue-template) on github [here](https://github.com/COMBINE-lab/salmon/issues).

### Reference Transcriptome

Alevin uses the same framework as Salmon to make index of the reference transcriptome and hence requires no extra flags or customization for the indexing stage.

In our tutorial we are working on `PBMC` data and will download human transcriptome. For example, we are downloading trancriptome with transcripts from protein-coding genes only as follows:

```python
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.pc_transcripts.fa.gz
echo "Done Downloading"
```

Once we have the reference-transcriptome, alevin index can be created using the following command:

```python
./bin/salmon index -i index -k 31 --gencode -p 4 -t gencode.v28.pc_transcripts.fa.gz
```
This command will build the alevin index inside the folder `index` in your current working directory.

### Transcript to Gene Mapping
Alevin works on transcript level equivalence classes to resolves potential UMI collision, while it also benefits from transcript to gene relation by sharing the information among the equivalene classes form one gene -- hence the need for a map from transcript id to gene-ids. Alevin requires the user to input a *tab* separated (one transcript-gene pair per line) file. For our-example we can extract the exact file using the following command for the `GTF` file downloaded from [here](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.basic.annotation.gtf.gz).

```python
bioawk -c gff '$feature=="transcript" {print $group}' human/gtf/gencode.v26.primary_assembly.annotation.gtf | awk -F ' ' '{print substr($4,2,length($4)-3) "\t" substr($2,2,length($2)-3)}' - > txp2gene.tsv
```

The above script, all it does is subsample the `transcript` feature from the *GTF* and dumps the corresponding txp-gene-ids pair in a _tab separated file (tsv)_.

### Dowloading Raw-Fastq

Raw Fastq file can be downloaded from the 10x-Genomics Support website from [here](https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k). You might have to provide your email-id and register to download the raw-fastqs.

### Ipython Notebook
Prefer to read ipython notebook ?
Check out the gist [here](https://gist.github.com/k3yavi/c501705ed2d29b12b0d10cf78b3ed001).
