##introduction4topgo.Rのつづき
##OGはそれぞれ含まれる遺伝子数が異なる。
##各OGに含まれるニワトリオーソログ全てのGOtermをまとめ、冗長性をなくしたデータを用意する
library(tidyverse)
library(biomaRt)
library(org.Gg.eg.db)
library(org.Hs.eg.db)

path = "5000_Paleognate3_28"
file <- "intersectGF.txt"

##introduction4topgo.Rで処理したニワトリデータを読み込み
Ggal_base <- read_csv(paste("input", path, "df.csv", sep = "/"))

##backgroundによって、Ggal_baseに色々つける
back <- read_tsv(paste("../cafe/summary", path, file, sep = "/"),
                 col_names = "Orthogroup")
Gga <- back %>% left_join(Ggal_base)

#EnvironmentからGgaを確認し、何列目までアノテーションが続いているか確認
#Overall ... 75列目までだったから、75にする
#Union ... 58
#Intersect ... 18
suji <- 18
Ggal <- Gga[1:suji]
Ggal[suji]

##BIOMARTでひもづけの準備
db <- useMart("ENSEMBL_MART_ENSEMBL")
##data <- "hsapiens_gene_ensembl"
data <- "ggallus_gene_ensembl"

Gg <- useDataset(data, mart = db)
filters <- listFilters(Gg)
attributes <- listAttributes(Gg)

##欲しいbiomart attributesを検索
filters[grep("id", filters[,1]),]
fil <- "entrezgene_id"
attributes[grep("id", filters[,1]),]
att <- c("go_id", "entrezgene_id")

##for文で、tidyGgal_baseの各列にgo_idを紐付けしていく


##関数定義
g <- function(x){
  res <- getBM(att, fil, x, Gg)
  res %>% return()
}

##whileで各列繰り返し
##まずは、繰り返しで蓄積させる用のneodfを定義
neodf = data_frame()

##1列目はOrthogroupの名前ベクトルなので、2列目から
i = 2
while (i <= ncol(Ggal)) {
  ##exはGgal_baseの列ベクトル
  ex <- Ggal[,i]
  ##関数の結果をaに入れる。
  a <- g(ex)
  while (class(a[[1]]) != "character") {
    i <- i + 1
    ex <- Ggal[,i]
    a <- g(ex)
  }
  ##numはGgal_baseの各列名
  num <- assign(paste("V", i, sep = ""), colnames(Ggal[,i]))
  
  newdf <- Ggal %>% 
    dplyr::rename(entrezgene_id = num) %>%
    left_join(a, by = "entrezgene_id") %>%
    dplyr::select(c(Orthogroup, go_id)) %>%
    filter(go_id != "NA") %>%
    filter(go_id != "")
  neodf <- neodf %>% bind_rows(newdf)
  
  print(i)
  i <- i + 1
}

##OrthogroupごとにGO termがたくさん取れるので、group_byでまとめ、連番を振る
neodf2 <- neodf %>% 
  group_by(Orthogroup) %>% 
  mutate(number = row_number())

##nouvelledfは連番ごとにpivot_widerして横に長くした、各OGごとのGO
nouvelledf <- neodf2 %>% 
  pivot_wider(names_from = number, values_from = go_id)


##冗長性排除
##完成したものを入れるdfを作る。列数は一応nouveledfに準拠
datf <- nouvelledf %>%
  filter(between(Orthogroup, 1, nrow(nouvelledf))) %>% 
  dplyr::select(Orthogroup)

ot <- nouvelledf %>% 
  t() %>%
  as_data_frame()

for (m in 1:ncol(ot)) {
  result <- ot %>% 
    distinct(ot[m]) %>%
    na.omit() %>% 
    t() %>% 
    as_data_frame()
  datf <- datf %>% bind_rows(result)
  print(m)
}
datf2 <- datf %>% ungroup() %>% dplyr::select(!"Orthogroup")

#outname <- "distinct_Gagal_union.tsv"
#outname <- "distinct_Gagal_union_11.tsv"
#outname <- "distinct_Gagal_intersect.tsv"
#outname <- "distinct_Gagal_overall.tsv"
#outname <- "distinct_Gagal_norapid.tsv"
outname <- "5000distinct_Gagal_intersect.tsv"
#outname <- "5000distinct"

datf2 %>% write_tsv(paste("input", path, outname, sep = "/"))

##distinct_Gagalから、Interested:都市鳥拡大、Background:全OGというふうにenrichmentしてみる
##OGdistinct_enrichmentに続く