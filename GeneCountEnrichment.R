library(tidyverse)
library(topGO)
library(org.Gg.eg.db)

##Please enter the directory containing interested dataframe
directory <- "Paleognate3_28"

##dataframe読み込み
df <- read_csv(paste("input", directory, "background.csv", sep = "/"))
df <- df %>% dplyr::rename(Entrez = V2)

##otamesi
orthologues <- read_tsv(paste("input", directory, "Orthogroups.tsv", sep = "/")) %>% 
  relocate(c("Orthogroup", "Ga.gal"))

##Interestedがspecificの時
Interested <- read_tsv(paste("../cafe/summary/", directory, "/SpecificGF.txt", sep = ""), col_names = "Orthogroup") %>%
  separate(Orthogroup, into=c("number", "Orthogroup")) %>% 
  dplyr::select(Orthogroup) %>%
  filter(Orthogroup != "object")

##Interestedがcommonの時
Interested <- read_tsv(paste("../cafe/summary/", directory, "/commonGF.txt", sep = ""), col_names = "Orthogroup")

newInterested <- Interested %>% left_join(df)

##Preparation for topGO data
##df$Entrez %in% newInterested$Entrez
input <- as.integer(df$Entrez %in% newInterested$Entrez)
names(input) <- df$Entrez
table(input)

##WORK ON TOPGO
topgodata <- new("topGOdata",
                 ontology = "BP",
                 allGenes = input,
                 geneSel = function(x) {x == 1},
                 nodeSize = 10,
                 annot = annFUN.org,
                 mapping = "org.Gg.eg.db",
                 ID = "ENTREZ")

##
whichAlgorithms()
whichTests()
resKS <- runTest(topgodata, algorithm = "elim", statistic = "ks")
resFisher <- runTest(topgodata, algorithm = "classic", statistic = "Fisher")
GenTable(topgodata, resFisher, orderBy='Fisher', topNodes=length(resFisher@score))


##6. Create Result Table
resTable <- GenTable(topgodata, resFisher, topNodes = 20)
resTable

showSigOfNodes(topgodata,
               score(resFisher),
               firstSigNodes = 10,
               useInfo = "all")

##For PANTHER
##the input file includes a column of numerical values
otam <- inputtoEntrez %>% as_data_frame() %>% 
  mutate(GO = nouvelledf3$Entrez) %>% relocate(GO)
write_tsv(otam, "test.tsv")