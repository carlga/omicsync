# [![Language](https://img.shields.io/badge/Language-Python-yellow.svg?style=flat-square)](http://dx.doi.org/10.1093/bioinformatics/btw354) &emsp; <img src="docs/images/omicsync_logo_whitebg.png#gh-light-mode-only" alt="omicsync" width="240"> <img src="docs/images/omicsync_logo_darkbg.png#gh-dark-mode-only" alt="omicsync" width="240">

---

`omicsync` is an integrated framework that aims to streamline the processing and analysis
of publicly available or newly generated high-throughput sequencing datasets. Its
modular architecture makes it applicable to a variety of omics technologies and by
using the snakemake workflow management system it allows to efficiently process files.
`omicsync` provides a powerful and reproducible execution environment which makes it 
particularly suitable for multi-omics integration and large-scale meta-analysis projects.

Relevant information for execution and reproducibility is stored in the configuration 
file. To create and clear the configuration file use `omicsync makeconfig` and 
`omicsync clearconfig`, respectively. `omicsync init` helps generate required index files
and databases. Retrieval of reference genomes and annotations is facilitated with
`omicsync getref`, and publicly available datasets with `omicsync getfq`. Quality control
steps with extended functionality can be accessed via `omicsync runqc`. Tailored
bioinformatics workflows can be executed with `omicsync run_rnaseq` and `omicsync run_atacseq`.

For more details and use case examples check the [documentation]().


## Installation



## Usage



## Citation



## Development & Contributions

