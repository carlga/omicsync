# ![osync](docs/images/osync_logo_whitebg.png#gh-light-mode-only) ![osync](docs/images/osync_logo_darkbg.png#gh-dark-mode-only)

[![lang](https://img.shields.io/badge/Language-Python-yellow.svg?style=flat-square)](http://dx.doi.org/10.1093/bioinformatics/btw354)

---

`osync` is an integrated framework that aims to streamline the processing and analysis
of publicly available or newly generated high-throughput sequencing datasets. Its
modular architecture makes it applicable to a variety of omics technologies and by
using the snakemake workflow management system it allows to efficiently process files.
`osync` provides a powerful and reproducible execution environment which makes it 
particularly suitable for multi-omics integration and large-scale meta-analysis projects.

Relevant information for execution and reproducibility is stored in the configuration 
file. To create and clear the configuration file use `osync makeconfig` and 
`osync clearconfig`, respectively. `osync init` helps generate required index files
and databases. Retrieval of reference genomes and annotations is facilitated with
`osync getref`, and publicly available datasets with `osync getfq`. Quality control
steps with extended functionality can be accessed via `osync runqc`. Tailored
bioinformatics workflows can be executed with `osync run_rnaseq`, `osync run_atacseq`,
`osync run_scrnaseq` and `osync run_scatacseq`.

For more details and use case examples check the [documentation]().


## Installation



## Usage



## Citation



## Development & Contributions

