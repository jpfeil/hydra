require(org.Hs.eg.db)
require(clusterProfiler)

args <- commandArgs(trailingOnly=T)

mm <- read.table(args[1], header=F)
head(mm)

bb <- read.table(args[2], header=F)
head(bb)

gene_map <- bitr(mm$V1,
                 fromType='SYMBOL',
                 toType='ENTREZID',
                 OrgDb='org.Hs.eg.db')

gene_background <- bitr(bb$V1,
                        fromType='SYMBOL',
                        toType='ENTREZID',
                        OrgDb='org.Hs.eg.db')

go <- enrichGO(gene          = gene_map$ENTREZID,
               universe      = gene_background$ENTREZID,
               OrgDb         = org.Hs.eg.db,
               ont           = "BP",
               pAdjustMethod = "BH",
               qvalueCutoff  = 0.01,
               readable      = TRUE)

simple <- clusterProfiler::simplify(go,
                                    cutoff=0.25,
                                    by='p.adjust',
                                    select_fun=min)

df <- data.frame(simple)

write.csv(df, file=args[3])