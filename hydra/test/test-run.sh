docker run -it -v $PWD:/data jpfeil/hydra:0.3.0 filter -e test/test-exp.tsv --output-dir hydra-test-filter --overwrite

docker run -it -v $PWD:/data jpfeil/hydra:0.3.0 sweep --output-dir hydra-test-sweep -e test/test-exp.tsv --gmt /opt/hydra/gene-sets/h.all.v6.2.symbols.gmt --gmt-regex IL2_STAT5 --overwrite

docker run -it -v $PWD:/data jpfeil/hydra:0.3.0 enrich --output-dir hydra-test-enrich -e test/test-exp.tsv --gmt /opt/hydra/gene-sets/h.all.v6.2.entrez.gmt --overwrite 
