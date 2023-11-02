##TopGO_1.Rの続き
##Prepare all GO term data per each OG and eliminate the duplication of GO terms.
BiocManager::install("biomaRt")
BiocManager::install("AnnotationDbi")
library(org.Gg.eg.db)
library(tidyverse)

path = "231016"

##Fams_total is the file which is processed in CAFE5 output (REGfamsignificance.R, REGfamextraction.R)
##ただし、df.csvはEntrez_idになっている必要がある。
Ggal_base <- read_csv(paste("..", path, "data/Processed_data/df_dotedited.csv", sep = "/"), quote = "")
back <- read_csv(paste("..", path, "data/Processed_data/Fams_total.csv", sep = "/")) %>% 
  dplyr::select(1) %>% distinct()
Ggal_before <- back %>% left_join(Ggal_base, by = c(FamilyID = "Orthogroup"))

##Please watch Gga via environment pane. Where is the column which have any elements(that is, all element are NA)?
##The number is the very input below!! (suji <- )
suji <- 47
Ggal <- Ggal_before[1:suji]
Ggal[[2]][1]

db <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL")
data <- "ggallus_gene_ensembl"

Gg <- biomaRt::useDataset(data, mart = db)
filters <- biomaRt::listFilters(Gg)
attributes <- biomaRt::listAttributes(Gg)

##Search for biomaRt attributes which you want
filters[grep("id", filters[,1]),]
fil <- "entrezgene_id"
attributes[grep("id", filters[,1]),]
att <- c("entrezgene_id", "go_id")

##for文で回す
neodf = data_frame()
i = 2
while (i <= ncol(Ggal)) {
  rep <- Ggal[,i]
  resEntreztoGO <- biomaRt::getBM(attributes = att,
                                  filters = fil,
                                  values = rep[[1]],
                                  mart = Gg)
  while (class(resEntreztoGO[[2]]) != "character") {
    i <- i + 1
    rep <- Ggal[,i]
    resEntreztoGO <- biomaRt::getBM(attributes = att,
                                    filters = fil,
                                    values = rep[[1]],
                                    mart = Gg)
  }
  num <- assign(paste("V", i, sep = ""), colnames(Ggal[,i]))
  newdf <- Ggal %>% 
    dplyr::rename(entrezgene_id = num) %>%
    left_join(resEntreztoGO, by = "entrezgene_id") %>%
    dplyr::select(c(FamilyID, go_id)) %>%
    filter(go_id != "NA") %>%
    filter(go_id != "")
  neodf <- neodf %>% bind_rows(newdf)
  print(i)
  print(neodf)
  i <- i + 1
}
print(neodf)
##neodfはOGごとにGO termを出した結果

##OG GO Go GOみたいなdfにしたい。pivot_wider用の連番列を作成
neodf2 <- neodf %>% 
  group_by(FamilyID) %>% 
  mutate(number = row_number())
##pivot_widerして、冗長性排除
nouvelledf <- neodf2 %>% 
  pivot_wider(names_from = number, values_from = go_id)

ot <- nouvelledf %>% 
  t() %>%
  as_data_frame()

datf = data_frame()
for (m in 1:ncol(ot)) {
  result <- ot %>% 
    distinct(ot[m]) %>%
    na.omit() %>% 
    t() %>% 
    as_data_frame()
  datf <- datf %>% bind_rows(result)
  print(m)
}
datf2 <- datf %>% ungroup()

write_tsv(datf2, paste("..", path, "data/Processed_data/GOterm.tsv", sep = "/"))
