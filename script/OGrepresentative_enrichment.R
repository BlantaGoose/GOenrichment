##Orthogroup内の適当な1遺伝子を代表にして解析するとき
library(tidyverse)
library(topGO)
library(org.Gg.eg.db)

##Please enter the directory containing interested dataframe
directory <- "Paleognate3_28"

##全体のデータが入っているdataframe読み込み
df <- read_csv(paste("input", directory, "background.csv", sep = "/")) %>% 
  dplyr::rename(Entrez = V2)


##Please enter your interested gene dataframe
Interested <- read_tsv(paste
                       ("../cafe/summary/", directory, "/commonGF.txt", sep = ""),
                       col_names = "Orthogroup")

newInterested <- Interested %>% 
  left_join(df)

##Preparation for topGO data
##OrthogroupごとにGO termを取ってきているので
##df$Entrez %in% newInterested$Entrez
input <- as.integer(df$Entrez %in% newInterested$Entrez)
names(input) <- df$Entrez
table(input)

##WORK ON TOPGO, Entrez
topgodata <- new("topGOdata",
                 ontology = "MF",
                 allGenes = input,
                 geneSel = function(x) {x == 1},
                 nodeSize = 10,
                 annot = annFUN.org,
                 mapping = "org.Gg.eg.db",
                 ID = "ENTREZ")

##
whichAlgorithms()
whichTests()
#resKS <- runTest(topgodata, algorithm = "elim", statistic = "ks")
resFisher <- runTest(topgodata, algorithm = "classic", statistic = "Fisher")

GenTable(topgodata, resFisher, orderBy='KS', topNodes=length(resFisher@score))


##6. Create Result Table
resTable <- GenTable(topgodata, resFisher, topNodes = 20)
resTable

showSigOfNodes(topgodata,
               score(resFisher),
               firstSigNodes = 10,
               useInfo = "all")

##GO termで優位なのがきたら、大事そうなやつを絞れる