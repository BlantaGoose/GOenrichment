##After introduction4topgo.R
##もしinputデータの遺伝子IDがEntrez_idなら、biomaRt不要なので、OG_representative_enrichment.Rに移る
library(tidyverse)
library(biomaRt)

##Please enter the directory containing interested dataframe
directory <- "Paleognate3_28"

##data 読み込み
df <- read_csv(paste("input", directory, "df.csv", sep = "/"))

##listMarts() ##使えるデータ確認
db <- useMart("ENSEMBL_MART_ENSEMBL")
Gg <- useDataset("ggallus_gene_ensembl", mart = db)
filters <- listFilters(Gg)
attributes <- listAttributes(Gg)

##欲しいbiomart attributesを検索
filters[grep("id", filters[,1]),]
fil <- "refseq_peptide_predicted"
attributes[grep("id", filters[,1]),]
att <- c("entrezgene_id")

rep <- df$V1

##valuesには紐付けしたいinput dataを用意。
##resNP,XPは元データがRefSeqの時
resXPtoEntrez <- getBM(attributes = att,
                       filters = fil,
                       values = rep,
                       mart = Gg)

fil <- "refseq_peptide"
resNPtoEntrez <- getBM(attributes = att,
                       filters = fil,
                       values = rep,
                       mart = Gg)

##dfのVI列に注目して結合。nouvelledfはVIのRefSeqをEntrezに変えたdf
##以下はRefSeqの時
newdf <- df %>% 
  left_join(resXPtoEntrez, 
            by = c(V1 = "refseq_peptide_predicted")) %>%
  left_join(resNPtoEntrez, 
            by = c(V1 = "refseq_peptide")) %>%
  relocate(Orthogroup, V1, entrezgene_id.x, entrezgene_id.y) %>% 
  dplyr::select(c(Orthogroup, V1, entrezgene_id.x, entrezgene_id.y)) %>%
  unite(Entrez2, entrezgene_id.x, entrezgene_id.y, sep = "") %>%
  mutate(Entrez = str_replace(Entrez2, pattern = "NA", replacement = "")) %>%
  filter(Entrez != "NA") %>% dplyr::select(!Entrez2)


##topgoのbackgroundデータ完成
write_csv(newdf, paste("input", directory, "representative_Gagal.csv", sep = "/"))