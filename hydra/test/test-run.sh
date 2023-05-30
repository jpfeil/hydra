docker run -it -v $PWD:/data jpfeil/hydra:0.3.4 filter -e test/test-exp.tsv.gz --output-dir hydra-test-filter --overwrite

docker run -it -v $PWD:/data jpfeil/hydra:0.3.4 sweep --output-dir hydra-test-sweep -e test/test-exp.tsv.gz --gmt /opt/hydra/gene-sets/h.all.v6.2.symbols.gmt --gmt-regex IL2_STAT5 --overwrite

docker run -it -v $PWD:/data jpfeil/hydra:0.3.4 enrich --output-dir hydra-test-enrich -e test/test-exp.tsv.gz --gmt /opt/hydra/gene-sets/h.all.v6.2.entrez.gmt --overwrite 
