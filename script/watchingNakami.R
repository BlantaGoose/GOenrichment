library(tidyverse)
library(topGO)

##introduction4topgo.Rで作成したdf.csvを読み込む
lib <- read_csv("input/5000_Paleognate3_28/df.csv")

##抽出したforegroundを読み込む
fore <- read_table("../cafe/summary/5000_Paleognate3_28/intersectGF_passerines.txt", col_names = FALSE)

foreGgal <- fore %>% left_join(lib, by = c(X1 = "Orthogroup"))

##foreGgalを確認して、何列目からNAが現れているか調べる
na <- 5
foreGgal <- foreGgal[1:na]

file <- "intersectOrthologs_passerines.csv"
foreGgal %>% write_csv(paste("../cafe/summary/sotsuken", file, sep = "/"))

##Inteesectをbackgroundにして、どんなGO termがあるのか調べる
ot <- read_tsv("input/5000_Paleognate3_28/5000distinct_Gagal_intersect.tsv")

##ot行はfamily、列はそれぞれ入っているGOterm。
##GOtermそれぞれが何このファミリーにくっついているのか確認
ott <- ot %>% pivot_longer(2:ncol(ot), names_to = "retsu", values_to = "GOterm")
ottt <- ott %>% 
  group_by(GOterm) %>% 
  count() %>%
  na.omit()
