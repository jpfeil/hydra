#!/usr/bin/env/Rscript

library('base')
library('stats')

args = commandArgs(trailingOnly = TRUE)


m <- as.matrix(read.table(args[1],
                          header=FALSE,
                          sep='\t'))

res <- stats::fisher.test(m)

write(res$p.value, "")
