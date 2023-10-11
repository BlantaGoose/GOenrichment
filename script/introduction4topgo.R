##目的はenrichment analysis。ニワトリオーソログを紐付けなければならない。
##このため、fams_lineage_specific.pyではなく、fams_lineage_common.py出なければならない。

library(tidyverse)
library(biomaRt)

##Please enter the directory containing interested dataframe
path <- "Paleognate3_28"
file <- "Orthogroups.tsv"

##data読み込み
##OrthoFinder outputのOrthogroups.tsvを移動させておこう。
orthologues <- read_tsv(paste("input", path, file, sep = "/")) %>% 
  relocate(c("Orthogroup", "Ga.gal"))
Ggalgene <- orthologues %>% 
  dplyr::select(c("Orthogroup", "Ga.gal")) %>% 
  na.omit()

##dataframe展開
tidyGgal = data.frame()
for (i in 1:nrow(Ggalgene)) {
  henkan =  as.character(Ggalgene[i, 2]) %>% 
    strsplit("\\D+\\.\\D+&") %>%
#    str_remove_all("\\.\\d") %>% ##遺伝子の後ろにある.1を消すやつ。
#    strsplit(", ") %>%
    as.data.frame() %>%
    t() %>%
    as_data_frame()
  tidyGgal = tidyGgal %>% bind_rows(henkan)
}

tidyGgal2 <- tidyGgal %>% 
  cbind(Ggalgene[1]) %>%
  relocate("Orthogroup") %>% 
  dplyr::select(!V1)
##遺伝子名が"種名&遺伝子ID"になっている場合


##ファイル書き込み
tidyGgal2 %>% write_csv(paste("input", path, "df.csv", sep = "/"))

##Orthogroups.tsv全体で20159個のファミリーが抽出でき、そのうちニワトリが持っていたのは15252だった
##OGdistinct（おすすめ）やOGrepresentativeに続く