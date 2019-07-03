docker run -it -v $PWD:/data jpfeil/hydra:0.1.0 filter -e test/test-exp.tsv --output-dir hydra-test-filter

docker run -it -v $PWD:/data jpfeil/hydra:0.1.0 sweep --output-dir hydra-test-sweep -e test/test-exp.tsv --gmt /opt/hydra/gene-sets/h.all.v6.2.symbols.gmt --gmt-regex IL2_STAT5

docker run -it -v $PWD:/data jpfeil/hydra:0.1.0 enrich --output-dir hydra-test-enrich -e test/test-exp.tsv --gmt /opt/hydra/gene-sets/h.all.v6.2.entrez.gmt --enrichment

