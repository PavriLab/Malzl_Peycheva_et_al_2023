# Malzl_Peycheva_et_al_2023
Repository for code of analyses presented in the paper

This repository contains code for all non-trivial analysis and plots presented in the paper. Most files used can be found in the GEO repository.

## HMM partitioning of RT profiles
This analysis was done using [@gspracklin's](https://github.com/gspracklin) [hmm_bigwigs](https://github.com/gspracklin/hmm_bigwigs) tool which was really helpful. However, there was a small caveat in that the default of this tool was to ignore chromosomes X, Y and MT. Since in our case X had quite a substantial coverage in Repli-Seq data we modified the the code such before running such that chromosome X was retained. In particular this concerned the `get_chroms` and the `create_df` functions where for `get_chroms` we simply deleted the chromosome filtering part and for `create_df` we added a safeguard for fetching values from the bigwig files which would otherwise fail when attempting to fetch a chromosome that is not contained in the file. This modification can be found below
```python
chromsizes = bbi.chromsizes(inputfile)
for item in chroms:
    if not item in chromsizes:
        continue
```

## MEF replication timing
By the time of writing the only MEF WT Repli-Seq data known to us was the [Rivera-Mulia et al. dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114747) which was generated from a 129/sv CAST/Ei hybrid mouse. Thus, it was necessary to extract the 129/sv reads only using a Python implementation of their [HARP algorithm](https://github.com/dvera/harp) (see also [Rivera-Mulia et al.](www.doi.org/10.1101/gr.232561.117)). In brief, 129/sv and CAST/Ei reference genomes were downloaded and reads were aligned to both using bowtie2. Only reads that aligned well to 129/sv were kept. Please see [mef_repliseq](https://github.com/dmalzl/Malzl_Peycheva_et_al_2023/tree/main/mef_repliseq) for details.
