
Run:
```
docker run -it -v $(pwd):/data jpfeil/hydra:0.1.0 --hallmark --expression exp.tsv --CPU 30
```

Output:
The output includes the bnpy models and a jupyter notebook for exploring clusters. The jupyter 
notebook depends on several python 2.7 libraries:
python 2.7
bnpy
seaborn
numpy
matplotlib
pandas 
scipy

Gene Set Sources:

The hallmarks of cancer gene sets were curated by Ted Goldstein:

https://github.com/tedgoldstein/hallmarks

The MSigDB gene sets were downloaded from the MSigDB portal:

http://software.broadinstitute.org/gsea/msigdb/index.jsp
