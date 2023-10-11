library(tidyverse)
library(topGO)
library(org.Gg.eg.db)
library(biomaRt)

##Entrez IDでtopGO
orthologues <- read_tsv("input/crowsniwatori/Orthogroups.tsv") %>% 
  relocate(c("Orthogroup", "Gallus_gallus"))
Ggalgene <- orthologues %>% dplyr::select(c("Orthogroup", "Gallus_gallus")) %>% 
  na.omit()
Ggalgene_orthogroup <- Ggalgene %>% dplyr::select("Orthogroup")

df = data.frame()
for (i in 1:nrow(Ggalgene)) {
  henkan =  as.character(Ggalgene[i, 2]) %>% 
    str_remove_all("\\.\\d") %>%
    strsplit(", ") %>%
    as.data.frame() %>%
    t() %>%
    as_data_frame()
  df = bind_rows(df, henkan)
}

df <- df %>% cbind(Ggalgene_orthogroup) %>%
  relocate("Orthogroup")
##dfはニワトリオーソログとオーソグループ
dfshort <- df %>% dplyr::select("Orthogroup", "V1")

##スデにわかってる場合
Interested <- read_tsv("../cafe/summary/0720summary/CommonGF.txt", col_names = "Orthogroup") %>%
  separate(Orthogroup, into=c("number", "Orthogroup")) %>% 
  dplyr::select(Orthogroup) %>%
  filter(Orthogroup != "object")

Interested
#Interestedにニワトリおーソロ具を配分
Interested2 <- Interested %>% left_join(dfshort, by = "Orthogroup") %>% dplyr::select("Orthogroup", "V1")
Interested2$V1

##何このOrthogroupに、ニワトリおーソロ具が紐づいたかメモ

geneListniwatori <- factor(as.integer(dfshort$V1 %in% Interested2$V1))
names(geneListniwatori) <- dfshort$V1
geneListniwatori

##topGO
topgodata <- new("topGOdata",
                 ontology = "BP",
                 allGenes = geneListniwatori,
                 geneSel = function(x)(x == 1),
                 nodeSize = 5,
                 annot = annFUN.org,
                 mapping = "org.Gg.eg.db",
                 ID = "genbank")



##custom annotation
##biomartでgene-to-GOs & GO-to-genes
resXPGO <- getBM(attributes = c("go_id", "refseq_peptide_predicted"),
                 filters = "refseq_peptide_predicted",
                 values = replace,
                 mart = Gg)
resNPGO <- getBM(attributes = c("go_id", "refseq_peptide"),
                 filters = "refseq_peptide",
                 values = replace,
                 mart = Gg)
resGO <- bind_rows(resXPGO, resNPGO)
resGO[is.na(resGO)] <- ""
resGO <- resGO %>% 
  unite(refseq, refseq_peptide_predicted, refseq_peptide, sep = "")

sukoshi <- resGO %>% distinct(refseq, .keep_all = TRUE) %>% na_if("") %>% na.omit()
sukoshi

nouvelledf <- df %>% 
  left_join(resGO, by = c(V1 = "refseq")) %>%
  dplyr::select(c("Orthogroup", "V1", "go_id"))




atarasi <- data_frame()
for (i in 1:ncol(genetogo)) {
#  print(genetogo[i-1, i])
  list = (c(colnames(genetogo[i-1,i]), genetogo[i-1,i]))
  atarasi <- atarasi %>% bind_rows(list)
}
