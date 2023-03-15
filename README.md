# 3'aQTL-pipe

[![Github Release](https://img.shields.io/badge/release-v1.1-brightgreen)](https://github.com/3UTR/3aQTL-pipe)
[![python Release](https://img.shields.io/badge/python-3.8-brightgreen)](https://www.python.org/downloads/)
[![R Release](https://img.shields.io/badge/R-3.6.3-brightgreen)](https://cran.r-project.org/)
[![DOI](https://zenodo.org/badge/480019097.svg)](https://zenodo.org/badge/latestdoi/480019097)

**Abbreviation** 
* APA: alternative polyadenylation
* 3'aQTL: 3′UTR alternative polyadenylation quantitative trait loci

This pipeline describes the step-by-step methods for analyzing dynamics alternative polyadenylation events across population-scale samples and performing association analysis between common genetic variants and APA usages to obtain a map of genetic regulation of APA. 

The scripts in this repository have been tested on 89 samples from Geuvadis RNA-seq Project and GTEx Project.
For conditions to reuse of these scripts please refer to LICENSE file.

## Using this pipeline
Details on how to prepare environment and use the script can be found on [GitHub wiki](https://github.com/3UTR/3aQTL-pipe/wiki) pages for this repository.

This pipeline relies on [Dapars2](https://github.com/3UTR/DaPars2) for APA quantification, [Matrix-eQTL](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/) for association mapping, and [SuSieR](https://github.com/stephenslab/susieR![image](https://user-images.githubusercontent.com/10413520/160762171-0a0e0d3c-f3ee-43a5-8b12-0920eba2dfac.png)
) for fine-mapping.

## Authors

Xudong Zou, Ruofan Ding, Wenyan Chen, Gao Wang, Shumin Cheng, Wei Li, Lei Li

Institute of Systems and Physical Biology, Shenzhen Bay Laboratory, Shenzhen 518055, China

## Citation
* Code and Execution:

**Using population-scale transcriptomic and genomic data to map 3' UTR alternative polyadenylation quantitative trait loci**
Xudong Zou, Ruofan Ding, Wenyan Chen, Gao Wang, Shumin Cheng, Qin Wang, Wei Li, Lei Li. ***STAR Protocols***,3(3):101566 **(2022)**.
DOI: https://doi.org/10.1016/j.xpro.2022.101566
https://www.sciencedirect.com/science/article/pii/S2666166722004464?via%3Dihub

* The first 3'aQTL atlas of human tissues:

**An atlas of alternative polyadenylation quantitative trait loci contributing to complex trait and disease heritability**

Lei Li, Kai-Lieh Huang, Yipeng Gao, Ya Cui, Gao Wang, Nathan D. Elrod, Yumei Li, Yiling Elaine Chen, Ping Ji, Fanglue Peng, William K. Russell, Eric J. Wagner & Wei Li. ***Nature Genetics***,53,994-1005 **(2021)**. DOI:https://doi.org/10.1038/s41588-021-00864-5

https://www.nature.com/articles/s41588-021-00864-5

## Contact
For any issues, please create a GitHub Issue.

## Funding
This work was supported by National Natural Science Foundation of China (no. 32100533) and startup funds from Shenzhen Bay Laboratory to L.L.
