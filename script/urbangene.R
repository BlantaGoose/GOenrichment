library(tidyverse)

ot <- read_tsv("input/Paleognate3_19/Orthogroups.tsv")
ot2 <- ot %>% dplyr::select(c("Orthogroup" , "Ac.gen", "At.cun", "Ca.ann", "Ca.ust", "Pa.mon", "Py.ruf", "Zo.alb"))

Interested <- read_tsv("../cafe/summary/1214/SpecificGF.txt", col_names = "Orthogroup") %>%
  separate(Orthogroup, into=c("number", "Orthogroup")) %>% 
  dplyr::select(Orthogroup) %>%
  filter(Orthogroup != "object")

neo <- Interested %>% left_join(ot2, by = "Orthogroup")
namae <- colnames(neo)[2:length(colnames(neo))]

##cafeで抽出できた特異的ファミリーに含まれる遺伝子を可視化
##遺伝子名の手前についてる種名を消したい。
for (m in 1:length(namae)) {
  n <- neo[m+1]
  n2 <- Interested
  henkan <- data_frame(matrix(nrow = nrow(Interested)), (ncol = ncol(Interested)))
  for (i in 1:nrow(n)) {
    henkan[i] <- as.character(n[i, 1]) %>%
      strsplit("\\D+\\.\\D+&")
    print(henkan)
    break
    
  }
}

