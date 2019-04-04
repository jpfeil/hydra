require(org.Hs.eg.db)
require(clusterProfiler)

args <- commandArgs(trailingOnly=T)

mm <- read.table(args[1], header=F)

bb <- read.table(args[2], header=F)

gene_map <- bitr(mm$V1,
                 fromType='SYMBOL',
                 toType='ENTREZID',
                 OrgDb='org.Hs.eg.db')

gene_background <- bitr(bb$V1,
                        fromType='SYMBOL',
                        toType='ENTREZID',
                        OrgDb='org.Hs.eg.db')

cu <- read.gmt(args[4])

res <- enricher(gene=gene_map$ENTREZID,
                 universe=gene_background$ENTREZID,
                 TERM2GENE = cu,
                 pAdjustMethod = 'BH',
                 qvalueCutoff = 0.01)

res <- setReadable(hall, 
                   keytype='ENTREZID',
                   OrgDb = org.Hs.eg.db)

df <- data.frame(res)

write.csv(df, file=args[3])
