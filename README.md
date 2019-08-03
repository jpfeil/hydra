## Overview

Hydra is a non-parametric clustering pipeline to identify 
differentially expressed pathways. The hydra pipeline includes 
routines for identifying multimodally distributed genes, scanning for 
differentially expressed gene sets, and identifying enriched gene 
sets from multimodally expressed genes. Hydra is available as a 
docker container for easy deployment.

## Gene Set Analysis
Use the sweep command to search for differential expression across a database 
of gene set annotations. Make sure the GMT file gene annotations matches the 
annotations in the expression matrix.
```
docker run -it -v $PWD:/data jpfeil/hydra:0.2.4 sweep \
                                                -e test/test-exp.tsv \
                                                -o test-sweep \
                                                --min-mean-filter 1.0 \
                                                --gmt /opt/hydra/gene-sets/h.all.v6.2.symbols.gmt \
                                                --gmt-regex IL2_STAT5
                                                --CPU 15 
```

## Unsupervised Enrichment Analysis
Use the filter tool to identify multimodally expressed genes

```
docker run -it -v $PWD:/data jpfeil/hydra:0.2.4 filter \
                                                -e test/test-exp.tsv \ 
                                                -o test-filter \
                                                --CPU 15 
```
This will generate a MultiModalGenes directory. The next step in the pipeline 
is to cluster enriched gene sets gene set genes.  

Perform GO enrichment clustering across multimodally expressed genes
```
docker run -it -v $PWD:/data jpfeil/hydra:0.2.4 enrich \
                                                -e <PATH to expression tsv file> \
                                                -m <PATH to MultiModalGenes dir> \
                                                --min-prob-filter 0.1 \
                                                --go-enrichment \
                                                -o <output directory> 
```

Or you can perform enrichment analysis using a user-specified gene set with the --gmt flag. The enrichment analysis uses the clusterProfiler tool, which requires the gene set database use entrez ids. The gene expression input matrix should still use gene symbols.
```
docker run -it -v $PWD:/data jpfeil/hydra:0.2.4 enrich \
                                                -e test/test-exp.tsv \
                                                -m <PATH to MultiModalGenes dir> \
                                                --min-prob-filter 0.1 \
                                                --gmt /opt/hydra/gene-sets/h.all.v6.2.entrez.gmt \
                                                -o test-enrich
```

## Jupyter Notebook
Interactive environment for investigating expression data. This comes with all of the 
hydra code and dependencies pre-installed.

```
docker run -it -v $(pwd):/data/ -p 8889:8888 jpfeil/hydra:0.2.4 notebook -e None
```

Perform interactive analysis through a browser
```
http://127.0.0.1:8889
```
The token for accessing the Jupyter notebook is printed to stdout.

Multimodal gene analysis
```
import sys
sys.path.append('/opt/hydra/')
import library.analysis as hydra


mm = hydra.EnrichmentAnalysis(mm_path,   # PATH to MultiModelGenes dir
                              exp_path,  # PATH to gene expression input
                              min_prob_filter=0.2,
                              gmt_path='GO')
                              
mm.get_enriched_terms()

genes = mm.get_enriched_term_genes()

clus = hydra.MultivariateMixtureModel(data=exp.reindex(genes),
                                      center=True,
                                      gamma=5.0,
                                      variance=2.0,
                                      K=5)
                                   
assignments = clus.get_assignments(exp.reindex(genes))
```

The minimum component probability filter can be used to tune 
the resolution of the clustering analysis with respect to the 
number of samples available. We provide a method `ScanEnrichmentAnalysis`
to explore how the minimum probability thresholds influences 
gene set enrichment and the number of clusters.
```
scan = hydra.ScanEnrichmentAnalysis(mm_path, 
                                    exp_path, 
                                    'GO', 
                                    CPU=7
                                    cluster=True).scan()
```
