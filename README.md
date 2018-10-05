## Overview

Hydra is a semi-parametric clustering pipeline to identify differentially expressed pathways. The hydra pipeline iterates over pathway gene sets, removes unimodally expressed genes, and then performs a multivariate clustering analysis to identify robust gene expression signatures. Hydra incorporates domain knowledge to improve clustering of disease populations. The trained model can then be used for N-of-1 tumor classification and online learning applications.

## Installation
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

## Run
Hydra comes pre-packaged with useful gene sets including the MSigDB HALLMARK gene sets. User generated gene sets must be placed in the gene-sets/misc/ directory to be included in the analysis.

```
cd hydra/hydra/
./run --expression test/test-exp.tsv --CPU 4 --all-genes
```

## Post-Analysis
Each gene set will get its own directory in the output directory. Start a jupyter notebook session to investigate expression clusters. Open the *.ipynb files to begin analyzing gene set clustering.  


Gene Set Sources:

The hallmarks of cancer gene sets were curated by Ted Goldstein:

https://github.com/tedgoldstein/hallmarks

The MSigDB gene sets were downloaded from the MSigDB portal:

http://software.broadinstitute.org/gsea/msigdb/index.jsp
