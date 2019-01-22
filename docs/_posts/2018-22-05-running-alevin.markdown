---
title:  "How to Run Alevin"
date:   2018-05-22 15:04:23
categories: [tutorial]
tags: [alevin]
---
## Running Alevin:

Alevin has multiple features and modes which based on the requirements and use-case can toggled on and off.

### Default Mode
Once all the above resource are available alevin can be run using the following command:

```bash
./bin/salmon alevin -lISR -1 fastqs/pbmc4k_S1_L001_R1_001.fastq.gz -2 fastqs/pbmc4k_S1_L001_R2_001.fastq.gz --chromium -i index -p 8 -o alevin_output --tgMap txp2gene.tsv
```

### External Whitelist 
Many times (specifically in case of 10x data), the whitelist is already available and based on the use-case, one don't wan't alevin to perform whitelisting of it's own. This mode can be activated by giving `--whitelist` flag to alevin and providing a list of known whitelisted CB (Cellular Barcodes). Given, this flag alevin will sequence correct CB towards the given CB and generated the gene-count matrix accordingly.

### 10x v3 Data
The data format is almost similar to v2, the only chnage we have to do is to use `--chromiumV3` instead of `--chromium`.

### Drop-seq Data
Alevin can be easily ported to quantify the data from *drop-seq* protocol as well. Since the data format is almost similar to `10x`, the only chnage we have to do is use `--dropseq` instead of `--chromium`.

### 10x v1 Data
Alevin is designed to primarily work with the file-format having CB-UMI in file while the corresponding read-sequence in the other. However, the 10x's `v1` chemistry does not follow the same convention and primarily have UMI and read-sequence in the same file, breaking the parsing format for alevin. To support the working of alevin we have written a wrapper script which takes in the `v1` chemistry data and feed it to alevin in the required format. Since the wrapper script is not optimized for performance one can observe a time-hit compared to analysis done on `v2` chemistry data.

To run alevin in `v1` mode the following three changes are required:

1. Using the flag `--gemcode` instead of `--chromium`.
2. NO `-1`, `-2`; flags should be given instead `-b` flag which specifies the path to the parent directory containing the reads (\*I1\* and \*RA\*).
3. The wrapper script has to be downloaded and compile using the following command:

```python
git clone git@github.com:COMBINE-lab/salmon.git
cd salmon/scripts/v1_10x;
g++ -std=c++11 -O3 -I ../../include -o wrapper wrapper.cpp -lz
```

Once compiled, the above command will generate a binary with the name `wrapper` and the path to this binary should be updated in the `run.sh` file. The command to quantify the `v1` data would be as follows (from inside the `v1_10x` folder):

```bash
./run.sh ./salmon alevin -lISR -b pbmc3k_fastqs/ --gemcode -i index -p 8 -o alevin_output --tgMap txp2gene.tsv
```

### Ipython Notebook
Prefer to read ipython notebook ?
Check out the gist [here](https://gist.github.com/k3yavi/c501705ed2d29b12b0d10cf78b3ed001).
