[![DOI](https://zenodo.org/badge/621298541.svg)](https://zenodo.org/doi/10.5281/zenodo.7924348)

# Malzl_Peycheva_et_al_2023
Repository for code of analyses presented in the paper

This repository contains code for all non-trivial analysis and plots presented in the paper. Most files used can be found in the GEO repository GSE228880.

## Software requirements and installation guide
Required software tools and packages used to run the code in this repository, including their versioning, can be found in the [`environment.yml`]() file which is also used to install them in an easy fashion. A quick rundown on how to do this from the commandline is provided as follows:

1. install a version of the anaconda package manager (we recommend using [miniconda](https://docs.conda.io/en/latest/miniconda.html))
2. clone the code from this repository using `git clone https://docs.conda.io/en/latest/miniconda.html`
3. change the cloned directory and run the following command to set up the environment `conda env create -f environment.yml` (this step usually takes a couple of minutes)
4. activate the environment using `conda activate analysis` and prepare it for use in jupyter with `python3 -m ipykernel install --name analysis --user`
5. for executing the code in the `*.ipynb` notebooks run `jupyter lab`

Versions of other software and how they were run can be found in the manuscript. Most of the code used to analyse the Hi-C data can be found at [pavrilab/hicer-nf](https://github.com/PavriLab/hicer-nf) and [pavrilab/hic_analysis](https://github.com/PavriLab/hic_analysis)

## Hardware requirements
All code found in the jupyter notebooks should be able to run on standard machines (i.e. laptops or PCs). For more demanding tasks like the MEF Repli-seq analysis, which requires running bowtie2 for alignment of reads and parsing the results we recommend running it on work stations or clusters with at least 16 cores and 80GB of RAM (the analysis in this work was run on the [CLIP cluster](https://www.clip.science/))

## HMM partitioning of RT profiles
This analysis was done using [@gspracklin's](https://github.com/gspracklin) [hmm_bigwigs](https://github.com/gspracklin/hmm_bigwigs) tool which was really helpful. However, there was a small caveat in that the default of this tool was to ignore chromosomes X, Y and MT. Since in our case X had quite a substantial coverage in Repli-Seq data we modified the the code such before running such that chromosome X was retained. In particular this concerned the `get_chroms` and the `create_df` functions where for `get_chroms` we simply deleted the chromosome filtering part and for `create_df` we added a safeguard for fetching values from the bigwig files which would otherwise fail when attempting to fetch a chromosome that is not contained in the file. This modification can be found below
```python
chromsizes = bbi.chromsizes(inputfile)
for item in chroms:
    if not item in chromsizes:
        continue
```

## MEF replication timing
By the time of writing the only MEF WT Repli-Seq data known to us was the [Rivera-Mulia et al. dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114747) which was generated from a 129/sv CAST/Ei hybrid mouse. Thus, it was necessary to extract the 129/sv reads only using a Python implementation of their [HARP algorithm](https://github.com/dvera/harp) (see also [Rivera-Mulia et al.](www.doi.org/10.1101/gr.232561.117)). In brief, 129/sv and CAST/Ei reference genomes were downloaded from Ensembl and reads were aligned to both using bowtie2. Only reads that aligned well to 129/sv were kept. Please see [mef_repliseq](https://github.com/dmalzl/Malzl_Peycheva_et_al_2023/tree/main/mef_repliseq) for details.
