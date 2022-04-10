# Mapping 3′UTR alternative polyadenylation quantitative trait loci through population-scale transcriptomic and genomic data

**Abbreviation** 
* APA: alternative polyadenylation
* 3'aQTL: 3′UTR alternative polyadenylation quantitative trait loci

This pipeline describes the detailed steps for analyzing dynamics alternative polyadenylation events across lymphoblastoid cell lines (LCL) samples from 445 unrelated individuals and performing association analysis between common genetic variants and APA events to obtain a map of genetic regulation of APA. The whole pipeline includes APA quantitative analysis across samples, association test between common genetic variants and APA usage (mapping 3'aQTL), fine-mapping 3'aQTLs, and other steps for preparing phenotype, genotype data and processing output of above analyses. The final outputs of this pipeline including the matrix of APA usage profile across samples, the table of association between common genetic variants and APA usage (3'aQTLs), the table of fine-mapped 3'aQTLs.

The scripts in this repository were tested on Geuvadis RNA-seq Project dataset which contains RNA sequencing data of 445 unrelated individuals (belong to five sub-populations, ~90 samples in each) and corresponding genotype data from 1000 Genome Project. We recommend analyzing each of the five sub-populations separately.

For conditions to reuse of these scripts please refer to LICENSE file.

## Using this pipeline
Details on how to prepare environment and use the script can be found on [GitHub wiki](https://github.com/Xu-Dong/3aQTL_pipe/wiki) pages for this repository.

This pipeline relies on [Dapars2](https://github.com/3UTR/DaPars2) for APA quantification, [Matrix-eQTL](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/) for association mapping, and [SuSieR](https://github.com/stephenslab/susieR![image](https://user-images.githubusercontent.com/10413520/160762171-0a0e0d3c-f3ee-43a5-8b12-0920eba2dfac.png)
) for fine-mapping.

## Authors

Xudong Zou, Ruofan Ding, Wenyan Chen, Gao Wang, Shumin Cheng, Wei Li, Lei Li

Institute of Systems and Physical Biology, Shenzhen Bay Laboratory, Shenzhen 518055, China

## Citation
* Code and Execution:

[Ref TBD]

* The first 3'aQTL atlas of human tissues:

**An atlas of alternative polyadenylation quantitative trait loci contributing to complex trait and disease heritability**

Lei Li, Kai-Lieh Huang, Yipeng Gao, Ya Cui, Gao Wang, Nathan D. Elrod, Yumei Li, Yiling Elaine Chen, Ping Ji, Fanglue Peng, William K. Russell, Eric J. Wagner & Wei Li. ***Nature Genetics***,53,994-1005 **(2021)**. DOI:https://doi.org/10.1038/s41588-021-00864-5

https://www.nature.com/articles/s41588-021-00864-5

## Contact
For any issues, please create a GitHub Issue.

## Funding
This work was supported by National Natural Science Foundation of China (no. 32100533) and startup funds from Shenzhen Bay Laboratory to L.L.
