## Overview

Hydra is a non-parametric clustering pipeline to identify differentially expressed pathways. 
The hydra pipeline includes routines for identifying multimodally distributed genes, scanning for 
differentially expressed gene sets, and identifying enriched gene sets from multimodally expressed 
genes. Hydra is distributed as a docker container and includes sever

## Identify Multimodal Genes
Use the filter tool to identify multimodally expressed genes

```
docker 
```

## Post-Analysis
Each gene set will get its own directory in the output directory. Start a jupyter notebook session to investigate expression clusters. Open the *.ipynb files to begin analyzing gene set clustering.  




## Installation (optional)
Download hydra
```
git clone https://github.com/jpfeil/hydra.git
cd hydra
```

Hydra requires python 2.7. We recommend creating a hydra specific conda environment.
```
conda create -n hydra python=2.7
source activate hydra
conda install -y --file requirements.tsv
```

Install bnpy clustering tools:
```
git clone https://github.com/bnpy/bnpy.git
cd bnpy/
pip install -e .
```

Install fgsea (optional)
Open an R session
```
library(devtools)
install_github("ctlab/fgsea")
```


Gene Set Sources:

The hallmarks of cancer gene sets were curated by Ted Goldstein:

https://github.com/tedgoldstein/hallmarks

The MSigDB gene sets were downloaded from the MSigDB portal:

http://software.broadinstitute.org/gsea/msigdb/index.jsp
