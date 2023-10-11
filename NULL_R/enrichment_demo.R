##https://ycl6.github.io/GO-Enrichment-Analysis-Demo/2_topGO.html#Start_R

library("topGO")
BiocManager::install("org.Mm.eg.db")
library("org.Mm.eg.db")
library("ggplot2")

#Load data
data <- data.table::fread("https://raw.githubusercontent.com/ycl6/GO-Enrichment-Analysis-Demo/master/DESeq2_DEG.txt")
data$GeneID <- substr(data$GeneID, 1, 18)
data

#Define significance threshold
up.idx <- which(data$padj < 0.05 & data$log2fc > 0) # FDR < 0.05 and logFC > 0
dn.idx <- which(data$padj < 0.05 & data$log2fc < 0) # FDR < 0.05 and logFC < 0
dim(data)

length(up.idx)
length(dn.idx)

#Define significant genes
all.genes <- data$GeneID
up.genes <- data[up.idx,]$GeneID
dn.genes <- data[dn.idx,]$GeneID

#Decide the sub-ontology to test
ontology <- "BP"

#Decide test algorithm
algorithm <- "weight01"

#Definde the statistical test used
statistic <- "fisher"

#set outfile prefix
outTitle <- paste0("topGO_GO-", ontology, "_ORA_", algorithm,"_", statistic)
outTitle

#Prepare input data
upList <- factor(as.integer(all.genes %in% up.genes))
names(upList) <- all.genes

head(upList, 30)
table(upList)

all.genes

#Create topGOdata object
upGOdata <- new("topGOdata",
                ontology = ontology,
                allGenes = upList,
                geneSel = function(x)(x == 1), 
                nodeSize = 10,
                annot = annFUN.org,
                mapping = "org.Mm.eg.db",
                ID = "ensembl")
org.Mm.eg()

#Test for enrichment
upRes <- runTest(upGOdata,
                 algorithm = algorithm,
                 statistic = statistic)
upRes
