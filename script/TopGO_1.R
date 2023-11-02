####topGO_1
##Use Gallus gallus orthologs from OrthoFinder output
##This script modify chicken's ortholog file.
library(tidyverse)

path = "231016"
orthologues <- read_tsv(paste("..", path, "data/Original_data/Orthogroups.tsv", sep = "/")) %>%
  relocate(c("Orthogroup", "Ga.gal"))
Ggalgene <- orthologues %>% 
  dplyr::select(c("Orthogroup", "Ga.gal")) %>% 
  na.omit()

tidyGgal = data.frame()
for (i in 1:nrow(Ggalgene)) {
  henkan =  as.character(Ggalgene[i, 2]) %>% 
    strsplit("\\D+\\.\\D+&") %>%
#    str_remove_all("\\.\\d") %>% ##遺伝子の後ろにある.1を消すやつ。
#    strsplit(", ") %>%
    as.data.frame() %>% 
    t() %>%
    as.data.frame()
  tidyGgal = tidyGgal %>% bind_rows(henkan)
}

tidyGgal2 <- tidyGgal %>% 
  cbind(Ggalgene[1]) %>%
  relocate("Orthogroup") %>% 
  dplyr::select(!V1)

tidyGgal2 %>% write_csv(paste("..", path, "data/Processed_data/df_dotedited.csv", sep = "/"))
