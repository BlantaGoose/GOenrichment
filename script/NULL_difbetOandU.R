##Difference between Overall and Union
##Overallとunionで何が違うのか調べる
library(tidyverse)
library(biomaRt)
library(org.Gg.eg.db)

path = "Paleognate3_28"

##introduction4topgo.Rで処理したニワトリデータを読み込み
Ggal_base <- read_csv(paste("input", path, "df.csv", sep = "/"))

##backgroundによって、Ggal_baseに色々つける
file <- "overallGF.txt"

back <- read_tsv(paste("../cafe/summary", path, file, sep = "/"),
                 col_names = "Orthogroup")
back2 <- read_tsv("../cafe/summary/Paleognate3_28/unionGF.txt", col_names = "Orthogroup")

##back...overallのなかに、back2...unionが持たないやつがどのくらいある？
nouni <- back[[1]][!back[[1]] %in% back2[[1]]]
df_nouni <- nouni %>% as_data_frame() %>%
  dplyr::rename(Orthogroup = value)

##biomaRtで紐付け。
db <- useMart("ENSEMBL_MART_ENSEMBL")
Gg <- useDataset("ggallus_gene_ensembl", mart = db)
filters <- listFilters(Gg)
attributes <- listAttributes(Gg)

##欲しいbiomart attributesを検索
filters[grep("id", filters[,1]),]
fil <- "entrezgene_id"
attributes[grep("id", filters[,1]),]
att <- c("go_id", "entrezgene_id")

#関数定義
g <- function(x){
  res <- getBM(att, fil, x, Gg)
  res %>% return()
}

##unionに存在しないOverallのOrthogroupの中身を結合
Gga <- df_nouni %>% left_join(Ggal_base)
##すべてNAなOrthogroupの行を消去


suji <- 75
Ggal <- Gga[1:suji]
Ggal[suji]


naggal <- is.na(Ggal)
##naggalで、2列目以降すべての列がTRUEになっている行を削除
naggal_df <- naggal %>% as_data_frame()
naggal_df %>% 

neodf = data_frame()

##1列目はOrthogroupの名前ベクトルなので、2列目から
i = 2
while (i <= ncol(Ggal_num5)) {

  ##exはGgal_baseの列ベクトル
  ex <- Ggal_num5[,i]

  ##関数の結果をaに入れる。
  a <- try(g(ex), silent = FALSE)
  
  ##もしgetBMで紐付けができなかったら、以下の処理を行う
  if (class(a) == "try-error") {
    go_id = "0"
    entrezgene_id = 0
    a <- data_frame(go_id, entrezgene_id)
  }

  while (class(a[[1]]) != "character") {
    i <- i + 1
    ex <- Ggal_num5[,i]
    a <- g(ex)
  }
    
  ##numはGgal_baseの各列名
  num <- assign(paste("V", i, sep = ""), colnames(Ggal_num5[,i]))
  
  newdf <- Ggal_num5 %>% 
    dplyr::rename(entrezgene_id = num) %>%
    left_join(a, by = "entrezgene_id") %>%
    dplyr::select(c(Orthogroup, go_id)) %>%
    filter(go_id != "NA") %>%
    filter(go_id != "")
  neodf <- neodf %>% bind_rows(newdf)
  
  print(i)
  i <- i + 1
}

neodf2 <- neodf %>% distinct(Orthogroup, go_id)
##紐づいている...

neodf2 %>% distinct()

##GOannotationされないときは、getBMででたGOtermの結果は全部numericになる
##characterじゃないときは無視して、iを大きくし、再度getBMでannotationする