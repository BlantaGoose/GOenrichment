##Orthogroup内の遺伝子全てに紐づいたGO termをまとめて、OG単位で解析
library(tidyverse)
library(topGO)
library(org.Gg.eg.db)

##Please enter the directory containing interested dataframe
path <- "5000_Paleognate3_28"

##全体のデータが入っているdataframe読み込み
##OrthogroupごとにbiomaRtでGO termを紐付けたバージョン
file <- "distinct_Gagal_norapid.tsv"
OG2GO <- readMappings(paste("input", path, file, sep = "/"))

##Please enter your interested gene dataframe
int <- "intersectGF.txt"
Interested <- read_tsv(paste("../cafe/summary", path, int, sep = "/"),
                       col_names = "Orthogroup")

##background。OG2GOでもいいし、commonでもいい
background <- names(OG2GO)
Interested <- as.character(Interested$Orthogroup)

##backgroundのどこにInterested
input <- factor(as.integer(background %in% Interested))
names(input) <- background
##確認
table(input)

##topGOdata making
topgodata <- new("topGOdata", 
                ontology = "BP", 
                allGenes = input,
                geneSel = function(x) {x == 1},  #Interested
                nodeSize = 10,
                annot = annFUN.gene2GO,
                gene2GO = OG2GO)  #background

##topGOdataにはどんな遺伝子が入っているのかな
sg <- sigGenes(topgodata)
str(sg)
numSigGenes(topgodata)

##Enrichment analysis
##"weight01"はGO hierarchyを考えてくれるが、"classic"は考えてくれない
##algorithm ... classic, elim, weight, weight01
##statistic ... ks, fisher
alg <- "weight01"
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

##GO termで優位なのがきたら、大事そうなやつを絞れる

