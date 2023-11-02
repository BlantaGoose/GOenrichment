##TopGO_2.Rの続き
library(tidyverse)
library(topGO)

path = "231016"

OG2GO <- topGO::readMappings(file = paste("..", path, "data/Processed_data/GOterm.tsv", sep = "/"))
Int <- read_tsv(paste("..", path, "data/Processed_data/Interested.tsv", sep = "/"))

background <- names(OG2GO)
Interested <- as.character(Int$FamilyID)

input <- factor(as.integer(background %in% Interested))
names(input) <- background
##確認
table(input)

##topGOdata making
library(topGO)
topgodata <- new("topGOdata", 
                 ontology = "BP", 
                 allGenes = input,
                 geneSel = function(x) {x == 1},  #Interested
                 nodeSize = 10,
                 annot = annFUN.gene2GO,
                 gene2GO = OG2GO)  #background

##topGOdataの中の遺伝子を調べる
sg <- sigGenes(topgodata)
str(sg)
numSigGenes(topgodata)

##Enrichment analysis
##"weight01"はGO hierarchyを考えてくれるが、"classic"は考えてくれない
##algorithm ... classic, elim, weight, weight01
##statistic ... ks, fisher
alg <- "classic"
sta <- "fisher"
res_en <- runTest(topgodata, algorithm = alg, statistic = sta)
res_en

Res <- GenTable(topgodata, res_en, orderBy = 'res_en', topNodes = length(res_en@score))
Res

##6. Create Result Table
resTable <- GenTable(topgodata, res_en, topNodes = 20)
resTab <- resTable %>% dplyr::rename(pvalue = "result1")

##多重検定補正
fdr <- p.adjust(resTab[,6], method = "BH")

resTable_fdr <- resTab %>% 
  bind_cols(fdr) %>%
  dplyr::rename(FDR = ...7)
resTable_fdr

##図示する
showSigOfNodes(topgodata,
               score(resFisher),
               firstSigNodes = 10,
               useInfo = "all")

