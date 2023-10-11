library(tidyverse)
library(topGO)
library(org.Gg.eg.db)
library(biomaRt)

##1. REFSEQ IDでtopGOするパターン
##OG内の遺伝子抽出
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
##とりあえず、OGの中に含まれるGeneのピリオドを消し、tidyにしたdf

##Prepare input data
#dfの中からbackgroundとinterestedを抽出
allgen <- df$V1
allgen


#head(Interested_orthogroup)をしても何も表示されない→backgroundにInterestedのOGがない！訳ではない
#→backgroundとInterestedのOGのニワトリオーソログを紐付けする。17/39に紐付けした
#遺伝子を選ぶ基準がわからないので、とりあえず一列目のものを使おう。

#次に、interestedを作成。今回は、カラス特異的で、カラス全種に共通するGF。
Interested <- read_tsv("../cafe/summary/0720summary/CommonGF.txt", col_names = "Orthogroup") %>%
  separate(Orthogroup, into=c("number", "Orthogroup")) %>% 
  dplyr::select(Orthogroup) %>%
  filter(Orthogroup != "object")
Himoduke_int <- Interested %>% 
  left_join(df, by = "Orthogroup")
intgen <- Himoduke_int$V1
intgen

#Input data←多分ここが変なのよ
#all.genes、up.genesみたいにcharacter型にしよう
input <- factor(as.integer(allgen %in% intgen))
names(input) <- allgen

head(input)
table(input)

#topGO for refseq ID
topgo <- new("topGOdata",
             ontology = "BP",
             allGenes = input,
             geneSel = function(x)(x == 1),
             nodeSize = 10,
             annot = annFUN.org,
             mapping = "org.Gg.eg.db",
             ID = "entrez")

##IDに入れるやつってどうやって調べればいいの？


#2. topGO for ensembl IDのパターン
#before topGO, we replace refseq ID for ensembl ID
replace <- df[2] %>% as.list()

listMarts()
db <- useMart("ENSEMBL_MART_ENSEMBL")
Gg <- useDataset("ggallus_gene_ensembl", mart = db)

filters <- listFilters(Gg)
filters[grep("refseq", filters[,1]),]
attributes <- listAttributes(Gg)
attributes[grep("refseq", filters[,1]),]
attributes[grep("ensembl", filters[,1]),]

##OrthoFinderによって取り出したニワトリOrthologuesには、NP_xxxxxx、XP_xxxxxxの両方が取り出されていた。
##とりまどっちも作成し、元々のallgen(=allGenes)にleft_join
resXP <- getBM(attributes = c("ensembl_peptide_id", "refseq_peptide_predicted"),
               filters = "refseq_peptide_predicted",
               values = replace,
               mart = Gg)

resNP <- getBM(attributes = c("ensembl_peptide_id", "refseq_peptide"),
               filters = "refseq_peptide",
               values = replace,
               mart = Gg)

nouvelledf <- df %>% 
  left_join(resXP, by = c(V1 = "refseq_peptide_predicted")) %>%
  left_join(resNP, by = c(V1 = "refseq_peptide")) %>%
  relocate(Orthogroup, V1, ensembl_peptide_id.x, ensembl_peptide_id.y)

df2 <- nouvelledf %>% dplyr::select(c(Orthogroup, V1, ensembl_peptide_id.x, ensembl_peptide_id.y))
df2[is.na(df2)] <- ""
df2 <- df2 %>% unite(Ensembl, ensembl_peptide_id.x, ensembl_peptide_id.y, sep = "") %>%
  filter(Ensembl != "")
##紐づけ完了

##input dataの準備
Interested2 <- Interested %>% left_join(df2)
Interested2$Ensembl

input2 <- factor(as.integer(df2$Ensembl %in% Interested2$Ensembl))
names(input2) <- df2$Ensembl

head(input2)
table(input2)

##create topGOdata object

topgo <- new("topGOdata",
             ontology = "BP",
             allGenes = input2,
             geneSel = function(x)(x == 1),
             nodeSize = 5,
             annot = annFUN.org,
             mapping = "org.Gg.eg.db",
             ID = "ENSEMBL")
inputtogene
input2 %>% class()
#おてあげ
#→protein idでもenrichment analysisできるの？あとraw data変えたらできないのでは





##input2のクラスをfactorジャなくてvectorにしてnumericにする
List <- InterestedtoEntrez$Entrez[!is.na(InterestedtoEntrez$Entrez)] #興味あるEntrezのやつ

li = c()
for (i in 1:length(List)) {
  if (List[i] %in% dftoEntrez$Entrez) {
    li <- c(li, which(dftoEntrez$Entrez == List[i]))
  }
}
print(li) #liに入っている番号は、dftoEntrez$Entrezの、興味ある遺伝子がある番号
##ここに1をふり、残りに0を降ったものを作りたい@9/22
li
for (i in 1:length(li)) {
  names(dftoEntrez$Entrez[i]) <- 1
}
copy <- dftoEntrez$Entrez
head(copy) ##copyに名前をつけよう


##まず、copyの要素数と同じリストを作成
analogue <- vector("list", length = length(copy))
for (i in 1:length(li)) {
  analogue[li[i]] <- copy[li[i]]
}

analogue2 <- str_replace_all(analogue, pattern = "NULL", "0") %>% as.numeric()
analogue2 %>% print()
##ついに、お名前シールができた。

named_scores_Gg <- setNames(analogue2, copy)
named_scores_Gg %>% class()

##topGodata
tg_data_Gg2 = new("topGOdata",
              ontology="BP",
              allGenes=named_scores_Gg,
              geneSelectionFun=function(x) {x > 0},
              nodeSize=10,
              annotationFun=annFUN.org,
              mapping="org.Gg.eg.db",
              ID="entrez")
##できた！！
##############################################################

##org.Gg.eg.dbは正常か？
scores_Gg = runif(length(mapped_genesGg), 0, 1)  # p-value-like
named_scores_Gg = setNames(scores_Gg, mapped_genesGg)

tg_data_Gg = new("topGOdata",
              ontology="BP",
              allGenes=named_scores_Gg,
              geneSelectionFun=function(x) {x < 0.01},
              nodeSize=10,
              annotationFun=annFUN.org,
              mapping="org.Gg.eg.db",
              ID="entrez")

##できてる...

##試しに、org.Gg.egdbからのEntrezをbiomartでEnsemblに変えてみる
attributes[grep("ensembl", filters[,1]),]


resEnsem <- getBM(attributes = "",
               filters = "entrezgene_id",
               values = replace,
               mart = Gg)
resEnsem
##これはできない。Ggmartのなかでentrezgeneが紐付けされてない？

"""
##.1をkesanai
dfkesanai = data.frame()
for (i in 1:nrow(Ggalgene)) {
  henkan = as.character(Ggalgene[i, 2]) %>% 
  #  str_remove_all("\\.\\d") %>%
    strsplit(", ") %>%
    as.data.frame() %>%
    t() %>%
    as_data_frame()
  dfkesanai = bind_rows(dfkesanai, henkan)
}

dfkesanai <- dfkesanai %>% 
  cbind(Ggalgene_orthogroup) %>% 
  relocate("Orthogroup")

allgenkesanai <- dfkesanai$V1

Himoduke_intkesanai <- Interested %>% 
  left_join(dfkesanai, by = "Orthogroup")
intgenkesanai <- Himoduke_intkesanai$V1

inputkesanai <- factor(as.integer(allgenkesanai %in% intgenkesanai))
names(inputkesanai) <- allgenkesanai

head(inputkesanai)
table(inputkesanai)

topgo <- new("topGOdata",
             ontology = "BP",
             allGenes = inputkesanai,
             geneSel = function(x)(x == 1),
             nodeSize = 10,
             annot = annFUN.org,
             mapping = "org.Gg.eg.db",
             ID = "ensembl")

##protein_idではなく、Gene idにしてみよう
attributes[grep("ensembl_gene_id", filters[,1]),]

resXPtogene <- getBM(attributes = c("ensembl_gene_id", "refseq_peptide_predicted"),
               filters = "refseq_peptide_predicted",
               values = replace,
               mart = Gg)

resNPtogene <- getBM(attributes = c("ensembl_gene_id", "refseq_peptide"),
               filters = "refseq_peptide",
               values = replace,
               mart = Gg)

nouvelledftogene <- df %>% 
  left_join(resXPtogene, by = c(V1 = "refseq_peptide_predicted")) %>%
  left_join(resNPtogene, by = c(V1 = "refseq_peptide")) %>%
  relocate(Orthogroup, V1, ensembl_gene_id.x, ensembl_gene_id.y)

dftogene <- nouvelledftogene %>% dplyr::select(c(Orthogroup, V1, ensembl_gene_id.x, ensembl_gene_id.y))
dftogene[is.na(dftogene)] <- ""
dftogene <- dftogene %>% unite(Ensembl, ensembl_gene_id.x, ensembl_gene_id.y, sep = "") %>%
  filter(Ensembl != "")

Interestedtogene <- Interested %>% left_join(dftogene)
Interestedtogene$Ensembl

inputtogene <- factor(as.integer(dftogene$Ensembl %in% Interestedtogene$Ensembl))
names(inputtogene) <- dftogene$Ensembl

head(inputtogene)
table(inputtogene)

topgo <- new("topGOdata",
             ontology = "BP",
             allGenes = inputtogene,
             geneSel = function(x)(x == 1),
             nodeSize = 5,
             annot = annFUN.org,
             mapping = "org.Gg.eg.db",
             ID = "ENSEMBL")
"""
##EntrezIDだとできる？
attributes[grep("entrez", filters[,1]),]

resXPtoEntrez <- getBM(attributes = c("entrezgene_id", "entrezgene_accession", "refseq_peptide_predicted"),
                     filters = "refseq_peptide_predicted",
                     values = replace,
                     mart = Gg)

resNPtoEntrez <- getBM(attributes = c("entrezgene_id", "entrezgene_accession", "refseq_peptide"),
                     filters = "refseq_peptide",
                     values = replace,
                     mart = Gg)

nouvelledftoEntrez <- df %>% 
  left_join(resXPtoEntrez, by = c(V1 = "refseq_peptide_predicted")) %>%
  left_join(resNPtoEntrez, by = c(V1 = "refseq_peptide")) %>%
  relocate(Orthogroup, V1, entrezgene_id.x, entrezgene_id.y)

dftoEntrez <- nouvelledftoEntrez %>% dplyr::select(c(Orthogroup, V1, entrezgene_id.x, entrezgene_id.y))
dftoEntrez[is.na(dftoEntrez)] <- ""
dftoEntrez <- dftoEntrez %>% unite(Entrez, entrezgene_id.x, entrezgene_id.y, sep = "") %>%
  filter(Entrez != "")

InterestedtoEntrez <- Interested %>% left_join(dftoEntrez)
InterestedtoEntrez$Entrez

inputtoEntrez <- factor(as.integer(dftoEntrez$Entrez %in% InterestedtoEntrez$Entrez))
names(inputtoEntrez) <- dftoEntrez$Entrez

head(inputtoEntrez)
table(inputtoEntrez)

topgo <- new("topGOdata",
             ontology = "BP",
             allGenes = inputtoEntrez,
             geneSel = function(x)(x == 1),
             nodeSize = 5,
             annot = annFUN.org,
             mapping = "org.Gg.eg.db",
             ID = "ENTREZ")

##Heavy Watal code
library(org.Hs.eg.db)
entrez_ids = mappedkeys(org.Hs.egGO)
scores = runif(length(entrez_ids), 0, 1)
named_scores = setNames(scores, entrez_ids)

class(named_scores)
named_scores
##OG内の代表遺伝子抽出
