#!/usr/bin/env/Rscript

library(data.table)
library(fgsea)
library(ggplot2)

args = commandArgs(trailingOnly = TRUE)

gmt.pth <- args[1]

rnk.pth <- args[2]

ranks <- read.table(rnk.pth,
                    header=FALSE,
                    col.names =c('Gene', 'Z'),
                    sep='\t')

ranks <- setNames(ranks$Z, ranks$Gene)

pathways <- gmtPathways(gmt.pth)

fgseaRes <- fgsea(pathways, 
                  ranks, 
                  maxSize=500, 
                  nperm=15000,
                  gseaParam = 0.5)

fwrite(fgseaRes, 
       file=args[3])
