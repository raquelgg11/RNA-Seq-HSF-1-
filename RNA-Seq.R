library(NOISeq)
library(GO.db)
library(AnnotationHub)
library(org.Hs.eg.db)
library(org.Ce.eg.db)
library(clusterProfiler)
library(ggplot2)
library(gplots)
library(plotly)
library(topGO)
library(Rgraphviz)
library(limma)
library(gProfileR)
library(biomaRt)


# Reading data and creating a dataframe.
datos_a_filtrar<- read.csv('datosomicos/NOISeq/FPKM.tabular', header = FALSE, sep="\t",stringsAsFactors = FALSE)

colnames(datos_a_filtrar)<-(c("gene","WT_0","WT_1", "WT_2", "Mutant_0", "Mutant_1", "Mutant_2"))
rownames(datos_a_filtrar)<- datos_a_filtrar$gene

#Rows in which the difference in expression between all WTs and all mutants is less than 30% are removed..
dif <- abs(datos_a_filtrar[,"WT_0"] - datos_a_filtrar[,"Mutant_0"])
dif1 <- abs(datos_a_filtrar[,"WT_0"] - datos_a_filtrar[,"Mutant_1"])
dif2<- abs(datos_a_filtrar[,"WT_0"] - datos_a_filtrar[,"Mutant_2"])
dif3 <- abs(datos_a_filtrar[,"WT_1"] - datos_a_filtrar[,"Mutant_0"])
dif4 <- abs(datos_a_filtrar[,"WT_1"] - datos_a_filtrar[,"Mutant_1"])
dif5 <- abs(datos_a_filtrar[,"WT_1"] - datos_a_filtrar[,"Mutant_2"])
dif6 <- abs(datos_a_filtrar[,"WT_2"] - datos_a_filtrar[,"Mutant_0"])
dif7 <- abs(datos_a_filtrar[,"WT_2"] - datos_a_filtrar[,"Mutant_1"])
dif8 <- abs(datos_a_filtrar[,"WT_2"] - datos_a_filtrar[,"Mutant_2"])
filas_a_eliminar <- which(dif < datos_a_filtrar[,"WT_0"]*0.30 | dif < datos_a_filtrar[,"Mutant_0"]*0.30 & 
                            dif1 < datos_a_filtrar[,"WT_0"]*0.30 | dif1 < datos_a_filtrar[,"Mutant_1"]*0.3 & 
                            dif2 < datos_a_filtrar[,"WT_0"]*0.30 | dif2 < datos_a_filtrar[,"Mutant_2"]*0.3 &
                            dif3 < datos_a_filtrar[,"WT_1"]*0.30 | dif3 < datos_a_filtrar[,"Mutant_0"]*0.3 & 
                            dif4 < datos_a_filtrar[,"WT_1"]*0.30 | dif4 < datos_a_filtrar[,"Mutant_0"]*0.3 & 
                            dif5 < datos_a_filtrar[,"WT_1"]*0.30 | dif5 < datos_a_filtrar[,"Mutant_1"]*0.3 & 
                            dif6 < datos_a_filtrar[,"WT_2"]*0.30 | dif6 < datos_a_filtrar[,"Mutant_0"]*0.3 & 
                            dif7 < datos_a_filtrar[,"WT_2"]*0.30 | dif7 < datos_a_filtrar[,"Mutant_1"]*0.3 & 
                            dif8 < datos_a_filtrar[,"WT_2"]*0.30 | dif8 < datos_a_filtrar[,"Mutant_2"]*0.3  )
filtrados <- datos_a_filtrar[-filas_a_eliminar,]

#Now, those rows in which the expression values lower than 10 are removed.
menores_diez <- which(filtrados$Mutant_1 < 10 & filtrados$WT_0 < 10 &filtrados$WT_1 < 10 &filtrados$WT_2 < 10 &filtrados$Mutant_0 < 10 &filtrados$Mutant_2 < 10 )
filtrados_2 <- filtrados[-menores_diez,]

filtrados_reducidos <- as.data.frame(filtrados_2[2:7])


#factores_new = data.frame(Muestras = c("WT_0","WT_1", "WT_2", "Mutant_0", "Mutant_1", "Mutant_2"))

replicates = data.frame(Replicates = c(rep("WT", 3), rep("Mut", 3)))
mydata <- NOISeq::readData(data = filtrados_reducidos,  factors = replicates)
head(assayData(mydata)$exprs)



#Now, we do a PCA per sample.

png("PCA-newdata.png",width = 600, height = 600)
myPCA = dat(mydata, type = "PCA")
explo.plot(myPCA, factor = "Replicates")

dev.off()

#COMPARACIÓN WT,Mut
comparation=c("WT,Mut")
mynoiseq1 = noiseqbio(mydata,k = 0.5, norm = "n", nclust = 15, plot = FALSE,
                      factor="Replicates", conditions = comparation, lc = 1, r = 50, adj = 1.5,
                      a0per = 0.9, random.seed = 12345, filter = 3, depth = NULL,
                      cv.cutoff = NULL, cpm = 1)


# Differentially expresed genes 
mynoiseq1.deg = degenes(mynoiseq1, q = 0.9, M = NULL)
mynoiseq1.deg
# Upregulated genes
mynoiseq1.deg1 = degenes(mynoiseq1, q = 0.9, M = "up")

# Downregulated genes
mynoiseq1.deg2 = degenes(mynoiseq1, q = 0.9, M = "down")

#Expression plot
png("expresionplot-newdata.png",width = 600, height = 600)
DE.plot(mynoiseq1, q = 0.9, graphic = "expr", log.scale = TRUE)
dev.off()

#MD plot.
png("MD-newdata.png",width = 600, height = 600)
DE.plot(mynoiseq1, q = 0.9, graphic = "MD")
dev.off()


geneList<-rownames(mynoiseq1.deg)


# Removing the .1 from the end otherwise it does not recognize the identifier.
geneList <- gsub("\\.[0-9]+", "", geneList)
# Removing "rna-" from the beginning of each gene name in the list
geneList <- sub("^rna-", "", geneList)

# Converting identifiers
gene.df <- bitr(geneList, fromType = "REFSEQ", 
                toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
                OrgDb = org.Ce.eg.db )

gene.df <- bitr(geneList, fromType = "REFSEQ", 
                toType = c("ENSEMBL", "SYMBOL", "ENTREZID"), OrgDb = org.Ce.eg.db)



# Grouping the GO terms to the desired level of depth, the larger it is, the more specific terms will appear (we are not enriching).
ggo <- groupGO(gene     =  geneList,
               OrgDb    = org.Ce.eg.db,
               ont      = "MF",
               keyType = "REFSEQ",
               level    = 4,
               readable = TRUE)


png("barplot-ggo-newdata.png",width = 1200, height = 600)
barplot(ggo, drop=TRUE, showCategory=20)
dev.off()



# GO Enrichment
ego2 <- enrichGO(gene         = geneList,
                 OrgDb         = org.Ce.eg.db,
                 keyType       = 'REFSEQ',
                 ont           = "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05)

png("barplot-ego2-newdata.png",width = 600, height = 600)
barplot(ego2, drop=TRUE, showCategory=20)
dev.off()
# Mapping ids to SYMBOL
ego2 <- setReadable(ego2, OrgDb = org.Ce.eg.db)


#Enrichment plots
png("barplot-ego2-newdata-2.png",width =1200 , height = 600)
barplot(ego2)
dev.off()

png("dotplot-ego2-newdata.png", width =1200 , height = 600)
dotplot(ego2)
dev.off()


gsea_genes=data.frame(rownames(mynoiseq1.deg),mynoiseq1.deg[3:3])
gsea_genes[,1] <-sub("^rna-", "", gsea_genes[,1])
gsea_genes[,1] <-sub("\\.[0-9]+", "", gsea_genes[,1])
geneList = gsea_genes[,2]
## feature 2: named vector
names(geneList) = as.character(gsea_genes[,1])
## feature 3: decreasing orde
geneList = sort(geneList, decreasing = TRUE)

ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Ce.eg.db,
              keyType       = 'REFSEQ',
              ont          = "MF",
              #      nPerm        = 1000,
              minGSSize    = 10,
              maxGSSize    = 500,
              pvalueCutoff = 0.1,
              verbose      = FALSE, eps = 0)

ego3 <- setReadable(ego3, OrgDb = org.Ce.eg.db)
png("dotplot-ego3-newdata.png", width =1200 , height = 600)

dotplot(ego3)
dev.off()


# Retrieving Ensemble GO annotations

ensembl=useMart("ENSEMBL_MART_ENSEMBL") 
ensembl=useMart("ensembl")
ensembl = useDataset("celegans_gene_ensembl",mart=ensembl)



datos_genes <-getBM(attributes = c("ensembl_gene_id","namespace_1003","go_id","name_1006"), 
                    filters="ensembl_gene_id", 
                    values=gene.df$ENSEMBL, 
                    mart=ensembl)

#Heatmap
x<-as.matrix(filtrados_2[2:7])
datos_heatmap<-apply(x, 2, as.double)
rownames(datos_heatmap)<-rownames(filtrados_2)
group <- as.factor(c("WT_0","WT_1", "WT_2", "Mutant_0", "Mutant_1", "Mutant_2"))
# creates a own color palette from red to green
mycol <- colorpanel(9,"red","green")
png("heatmap-newdata.png",width = 1024, height = 1024)
#rowside_colors=c(rep("black", 29), rep("blue", 38))
col_breaks = c(seq(-1,0,length=5),  # for red
               seq(0.01,1,length=5) )

heatmap.2(datos_heatmap,
          col = mycol,
          xlab = "Experiment",
          ylab = "Genes",
          main = "Heatmap by expression",
          labRow = rownames(datos_heatmap),
          labCol = group,
          dendrogram = c("both"),
          margins=c(12,8),
          trace="column",
          tracecol="black",
          #  breaks=col_breaks,
          density.info="none",
          scale = "row",
          notecol="black",
          reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
          distfun=function(x) as.dist(1-cor(t(x))),
          #  distfun = dist(datos_fusionados_FC, method = "manhattan"),
          hclustfun=function(x) hclust(x, method="complete")
          #     RowSideColors = rowside_colors ,   # line width
          
)

dev.off()