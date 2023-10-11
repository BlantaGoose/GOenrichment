##7CorvidaeのtopGO解析を行う。カラスノードにおけるGFの機能は、全７種の鳥類を使用したときとは異なる偏りをしているのでは。
library(topGO)
library(tidyverse)

#BiocManager::install("org.Gg.eg.db") #ニワトリゲノムへのアノテーションデータ
#library(org.Gg.eg.db)

#Let's assign Gene ID into interested GF!
#まず、バックグラウンドを読み込む。
#entrez_ids = mappedkeys(org.Gg.egGO) #entrezはgene identifierの一種。他にもGenBank, Aliasなど
#scores = runif(length(entrez_ids), 0, 1)  # p-value-like

named_scores = setNames(scores, entrez_ids)

#creating topGOdata
#・gene identifier
#・mapping between identifiers and GO term (Bioconductorのmicroarray specific annotation packageで利用可能なことが多い)
#・GOの階層構造（GO.db package）
tg_data = new("topGOdata", #クラスはtopGOdata
              ontology = "BP", #Biology Process, Cellular Component, Molecular Functionの中から気になるのを選ぼう
              allGenes = named_scores, #backgroundの選択。
              geneSelectionFun = function(x) {x < 0.01}, #allGenesがnumericな場合のみ可能。allGenesの数値を引数として、興味ある遺伝子にTRUEとなる
              nodeSize = 10, #nodeSizeよりも少ない数の遺伝子にしかアノテーションできなかったGOtermの排除。1以上のinteger
              annotationFun = annFUN.org, #遺伝子IDとGO termをつなぐ。annFun.orgの他にも色々
              mapping = "org.Gg.eg.db", #BioConductor AnnotationData Packagesにあるような種の遺伝子IDなら、mappingを使う
              ID = "entrez") #allGenesすなわちnamed_scoresに与えた名前。


#conducting test
resElimFisher = runTest(tg_data, algorithm = "elim", statistic = "fisher")
str(resElimFisher) #str()はオブジェクトの表示をする。
geneData(resElimFisher) #結果が知りたい→geneData()で要約統計量入りのリストを表示。
score(resElimFisher) %>% head()

#Summarying results






##geneID2GOが完成するまでの物語。あほくさい
##出力時に1列目OGs、2列目対応するGO termsになっているようにしたい。
##ここの作業はPythonでやろう。あほくさい
ANN1300 <- ANN[1,] %>% left_join(ANN[2,], by = "Orthogroup") %>%
  left_join(ANN[3,], by = "Orthogroup") %>%
  left_join(ANN[4,], by = "Orthogroup") %>%
  left_join(ANN[5,], by = "Orthogroup") %>%
  left_join(ANN[6,], by = "Orthogroup") %>%
  left_join(ANN[7,], by = "Orthogroup") %>%
  left_join(ANN[8,], by = "Orthogroup") %>%
  left_join(ANN[9,], by = "Orthogroup") %>%
  left_join(ANN[10,], by = "Orthogroup") %>%
  left_join(ANN[11,], by = "Orthogroup") %>%
  left_join(ANN[12,], by = "Orthogroup") %>%
  left_join(ANN[13,], by = "Orthogroup") %>%
  left_join(ANN[14,], by = "Orthogroup") %>%
  left_join(ANN[15,], by = "Orthogroup") %>%
  left_join(ANN[16,], by = "Orthogroup") %>%
  left_join(ANN[17,], by = "Orthogroup") %>%
  left_join(ANN[18,], by = "Orthogroup") %>%
  left_join(ANN[19,], by = "Orthogroup")
ANN1300 <- ANN1300 %>% unite(go_id.x, go_id.y, go_id.x.x, go_id.y.y, go_id.x.x.x, go_id.x.x.x.x, go_id.x.x.x.x.x, go_id.y.y.y, go_id.x.x.x.x.x.x, go_id.x.x.x.x.x.x.x, go_id.x.x.x.x.x.x.x.x, go_id.x.x.x.x.x.x.x.x.x, go_id.y.y.y.y, go_id.y.y.y.y.y, go_id.y.y.y.y.y.y, go_id.y.y.y.y.y.y.y, go_id.y.y.y.y.y.y.y.y, go_id.y.y.y.y.y.y.y.y.y, go_id, sep = " ")
ANN1300
ANN1373 <- ANN[20,]
ANN1923 <- ANN[21,]
names(ANN1373)[2] <- "go_id.x"
names(ANN1923)[2] <- "go_id.x"


ANN0819 <- ANN[22,] %>% left_join(ANN[23,], by = "Orthogroup") %>%
  left_join(ANN[24,], by = "Orthogroup") %>%
  left_join(ANN[25,], by = "Orthogroup") %>%
  left_join(ANN[26,], by = "Orthogroup") %>%
  left_join(ANN[27,], by = "Orthogroup") %>%
  left_join(ANN[28,], by = "Orthogroup") %>%
  left_join(ANN[29,], by = "Orthogroup") %>%
  left_join(ANN[30,], by = "Orthogroup") %>%
  left_join(ANN[31,], by = "Orthogroup") %>%
  left_join(ANN[32,], by = "Orthogroup") %>%
  left_join(ANN[33,], by = "Orthogroup") %>%
  left_join(ANN[34,], by = "Orthogroup") %>%
  left_join(ANN[35,], by = "Orthogroup") %>%
  left_join(ANN[36,], by = "Orthogroup") %>%
  left_join(ANN[37,], by = "Orthogroup") %>%
  left_join(ANN[38,], by = "Orthogroup") %>%
  left_join(ANN[39,], by = "Orthogroup") %>%
  left_join(ANN[40,], by = "Orthogroup") %>%
  left_join(ANN[41,], by = "Orthogroup") %>%
  left_join(ANN[42,], by = "Orthogroup") %>%
  left_join(ANN[43,], by = "Orthogroup") %>%
  left_join(ANN[44,], by = "Orthogroup") %>%
  left_join(ANN[45,], by = "Orthogroup") %>%
  left_join(ANN[46,], by = "Orthogroup") %>%
  left_join(ANN[47,], by = "Orthogroup") %>%
  left_join(ANN[48,], by = "Orthogroup") %>%
  left_join(ANN[49,], by = "Orthogroup") %>%
  left_join(ANN[50,], by = "Orthogroup") %>%
  left_join(ANN[51,], by = "Orthogroup") %>%
  left_join(ANN[52,], by = "Orthogroup") %>%
  left_join(ANN[53,], by = "Orthogroup") %>%
  left_join(ANN[54,], by = "Orthogroup") %>%
  left_join(ANN[55,], by = "Orthogroup") %>%
  left_join(ANN[56,], by = "Orthogroup") %>%
  left_join(ANN[57,], by = "Orthogroup") %>%
  left_join(ANN[58,], by = "Orthogroup") %>%
  left_join(ANN[59,], by = "Orthogroup") %>%
  left_join(ANN[60,], by = "Orthogroup") %>%
  left_join(ANN[61,], by = "Orthogroup") %>%
  left_join(ANN[62,], by = "Orthogroup") %>%
  left_join(ANN[63,], by = "Orthogroup") %>%
  left_join(ANN[64,], by = "Orthogroup") %>%
  left_join(ANN[65,], by = "Orthogroup") %>%
  left_join(ANN[66,], by = "Orthogroup") %>%
  left_join(ANN[67,], by = "Orthogroup")
ANN0819 <- ANN0819 %>% unite(go_id.x, go_id.y,
                  go_id.x.x, go_id.y.y, 
                  go_id.x.x.x, go_id.y.y.y,
                  go_id.x.x.x.x, go_id.y.y.y.y,
                  go_id.x.x.x.x.x, go_id.y.y.y.y.y,
                  go_id.x.x.x.x.x.x, go_id.y.y.y.y.y.y,
                  go_id.x.x.x.x.x.x.x, go_id.y.y.y.y.y.y.y,
                  go_id.x.x.x.x.x.x.x.x, go_id.y.y.y.y.y.y.y.y,
                  go_id.x.x.x.x.x.x.x.x.x, go_id.y.y.y.y.y.y.y.y.y,
                  go_id.x.x.x.x.x.x.x.x.x.x, go_id.y.y.y.y.y.y.y.y.y.y,
                  go_id.x.x.x.x.x.x.x.x.x.x.x, go_id.y.y.y.y.y.y.y.y.y.y.y,
                  go_id.x.x.x.x.x.x.x.x.x.x.x.x, go_id.y.y.y.y.y.y.y.y.y.y.y.y,
                  go_id.x.x.x.x.x.x.x.x.x.x.x.x.x, go_id.y.y.y.y.y.y.y.y.y.y.y.y.y,
                  go_id.x.x.x.x.x.x.x.x.x.x.x.x.x.x, go_id.y.y.y.y.y.y.y.y.y.y.y.y.y.y,
                  go_id.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x, go_id.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y, 
                  go_id.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x,go_id.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y,
                  go_id.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x,go_id.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y,
                  go_id.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x, go_id.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y,
                  go_id.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x, go_id.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y,
                  go_id.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x,
                  go_id.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y,
                  go_id.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x,
                  go_id.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y,
                  go_id.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x,
                  go_id.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y,
                  go_id.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x.x,
                  go_id.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y.y)


ANN0059 <- ANN[68,] %>% left_join(ANN[69,], by = "Orthogroup") %>%
  left_join(ANN[70,], by = "Orthogroup") %>%
  left_join(ANN[71,], by = "Orthogroup") %>%
  left_join(ANN[72,], by = "Orthogroup")
ANN0059 <- ANN0059 %>% unite(go_id.x, go_id.y, go_id.x.x, go_id.y.y, go_id)

ANN0436 <- ANN[73,] %>% left_join(ANN[74,], by = "Orthogroup") %>%
  left_join(ANN[75,], by = "Orthogroup") %>%
  left_join(ANN[76,], by = "Orthogroup") %>%
  left_join(ANN[77,], by = "Orthogroup")
ANN0436 <- ANN0436 %>% unite(go_id.x, go_id.y, go_id.x.x, go_id.y.y, go_id)


ANN0798 <- ANN[78,] %>% left_join(ANN[79,], by = "Orthogroup") %>%
  left_join(ANN[80,], by = "Orthogroup") %>%
  left_join(ANN[81,], by = "Orthogroup")
ANN0798 <- ANN0798 %>% unite(go_id.x, go_id.y, go_id.x.x, go_id.y.y)

ANN <- merge(ANN0059, ANN0436, all = T)
ANN
ANN <- merge(ANN, ANN0798, all = T)
ANN
ANN <- merge(ANN, ANN0819, all = T)
ANN <- merge(ANN, ANN1300, all = T)
ANN <- merge(ANN, ANN1373, all = T)
ANN <- merge(ANN, ANN1923, all = T)
ANN
names(ANN)[2] <- "go_id"

ANN %>% write.table("/Users/araisota/research/enrichment/topGO/ANN.txt",
                    quote = F, col.names = T)

# http://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html
#1.まずOg to GOされたbackgroundを用意（うまくいかない）
library(tidyverse)
library(topGO)

#Loading your data containing GO-to-gene。file=your data name
#一列目は遺伝子名、二列目はGO terms
geneID2GO <- readMappings(file="ANN.txt", sep = "") #タブ区切りになっていないのならしよう。
geneID2GO #なんかリストになってるけど、仕方なし

names(geneID2GO)[1] <- "Orthogroup"
names(geneID2GO)[2] <- "go_id"

#geneListの用意。比較した全てのOGs to GO termsを


myGOdata <- new("topGOdata", 
                ontology = "BP",
                allGenes = geneList, #geneListは、興味ある系統群と比較したい全てのアノテーション情報が入っている。
                annot = annFUN.gene2GO,
                gene2GO = geneID2GO)
#





#興味ある遺伝子リストと、gene universe内の比較したい遺伝子リストを指摘。gene universeにはyour fileの遺伝子のGO annotationが入っている。
geneUniverse <- names(geneID2GO)
genesOfInterest <- read.table("ANN.txt", header=TRUE)
genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse #geneListは、gene universeの中のどれが興味ある遺伝子なのか指摘している


geneUniverse <- names(geneID2GO) #興味ある遺伝子群と比較する用のものをgeneUniverseから取り出す
genesOfInterest <- read.table("ANN.txt", sep = "\t", quote = "", comment.char = "",) #引数をこうするとなんかうまくいく
genesOfInterest <- as.character(genesOfInterest)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

#topGOdataの準備
myGOdata <- new("topGOdata",
                ontology="BP", 
                allGenes=geneList,
                annot = annFUN.gene2GO,
                gene2GO = geneID2GO)

