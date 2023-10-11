library(tidyverse)
library(topGO)
library(org.Gg.eg.db)

##Please enter the directory containing interested dataframe
path <- "Paleognate3_28"

##全体のデータが入っているdataframe読み込み
##OrthogroupごとにbiomaRtでGO termを紐付けたバージョン
file <- "distinct_Gagal_overall.tsv"
OG2GO <- readMappings(paste("input", path, file, sep = "/"))

##Please enter your interested gene dataframe
int <- "specificGF.txt"
Interested <- read_tsv(paste("../cafe/summary", path, int, sep = "/"),
                       col_names = "Orthogroup")

##background。OG2GOでもいいし、commonでもいい
background <- names(OG2GO)
Interested <- as.character(Interested$Orthogroup)

##backgroundのどこにInterestedがアルカナ
input <- factor(as.integer(background %in% Interested))
names(input) <- background

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
alg <- "weight01"
resFisher <- runTest(topgodata, algorithm = alg, statistic = "fisher")
resFisher

Res <- GenTable(topgodata, resFisher, orderBy = 'resFisher', topNodes = length(resFisher@score))


##6. Create Result Table
resTable <- GenTable(topgodata, resFisher, topNodes = 20)
resTable

##図示する
showSigOfNodes(topgodata,
               score(resFisher),
               firstSigNodes = 10,
               useInfo = "all")

##GO termで優位なのがきたら、大事そうなやつを絞れる
