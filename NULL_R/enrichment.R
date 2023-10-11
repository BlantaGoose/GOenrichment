##topGOでenrichment analysisしたいよね
library(biomaRt)
library(tidyverse)
library(topGO)
library(org.Gg.eg.db)

##全種と興味ある種とでOGs-GOtermを作成する。
##0.それぞれのOGsにどんなgeneが含まれているかのデータを用意。
orthogroups <- read_tsv("input/crowsniwatori/Orthogroups.tsv")

##加工によってすべてのOGsを取り出す。
##fams.csvファイルにGUIですることは以下の通り。
##newickが含まれている一行目を削除。
##各行の種名をOrthogroups.tsvと共通に。
moto <- read_csv("../cafe/summary/0720summary/0720summary_fams.csv", col_names = FALSE) %>%
  separate(X1, into = c("Node","X1"), sep = ":")

##後々biomaRtでリファレンスにする種の行を取り出す。
outgroupmoto <- moto %>% filter(Node %in% c("F.alb<10>",
                                            "P.maj<14>", "T.gut<12>",
                                            "G.gal<0>"))
##NAの列がうざったいので、確認してselectで消去しよう
outgroupmoto1 <- outgroupmoto %>% dplyr::select(Node:X17) %>%
  t()
outgroupmoto2 <- outgroupmoto1 %>% str_remove("\\[.+\\]")
outgroupmoto2[2] <- outgroupmoto2[2] %>% str_remove("\\\t")

##一つ以上のノードで急進化したGFをまとめたdata.frameを作ろう。
Allmoto <- moto %>% filter(Node == "Overall rapid ") %>% dplyr::select(!1)
Allmoto[1] <- Allmoto[1] %>% str_remove("\\\t")
All1 <- Allmoto %>% t() %>% as.data.frame()
All2 <- All1 %>% left_join(orthogroups, by = c(V1 =  "Orthogroup"))
All3 <- All2 %>% pivot_longer(2:ncol(All2), names_to = "species", values_to = "refseq") %>%
  na.omit()
##All2では、全急進化GFに、各生物種でどんなgeneが入っているのか紐付けした。


##refseqのrefseq_peptide_predictedのお尻(XPaaaaaa.1→XPaaaaaa)を削る
##まず、三列目を分解して分ける
All3_2 <- All3 %>% 
  separate(refseq, 
           into = c("1","2","3","4","5","6","7","8","9","10","11", "12", "13",
                    "14","15","16","17","18","19","20","21","22","23","24","25",
                    "26","27","28","29","30","31","32","33","34","35","36","37",
                    "38","39","40","41","42","43","44","45", "46","47","48","49",
                    "50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68"),#"69","70"), 
           sep = ",")
##すべての行でNAになるのが何列目からか、上を走らせて確認してから、NAの列を削除するのが良さそう

##すべての列について、お尻削らせる
All3_3 <- All3_2 %>% 
  mutate(across(.cols = !c(V1, species),
                .fns = str_replace, 
                pattern = "\\.\\d",
                replacement = ""))
All4 <- All3_3 %>% 
  pivot_longer(3:ncol(All3_3), names_to = "NUMBER", values_to = "refseq") %>%
  na.omit() %>% 
  dplyr::select(!NUMBER)
##mutate(across)は、指定した列に対して一括で操作を行う。strをするときは、上のようにバラさないといけないっぽい
##https://community.rstudio.com/t/across-and-stringr-together-with-mutate/93439


##ref$speciesを三種それぞれに分解して、三種それぞれのrefseq情報を作成。
Falb <- All4[All4$species == "Ficedula_albicollis",]
Pmaj <- All4[All4$species == "Parus_major",]
Tgut <- All4[All4$species == "Taeniopygia_guttata",]
Ggal <- All4[All4$species == "Gallus_gallus",]
##急進化が確認されたOG全体を背景とするために、外群3種のEntrezgene_id紐付け情報を作成する

##2.biomaRtで参照するデータの準備
listMarts() ##どんなmartがあるのか調べる。
db <- useMart("ENSEMBL_MART_ENSEMBL") ##使いたいbiomartを指定
listDatasets(db) ##どんなデータがストックされているか調べる
fa <- useDataset("falbicollis_gene_ensembl", mart = db) ##martはデータベースの名前
pm <- useDataset("pmajor_gene_ensembl", mart = db)
tg <- useDataset("tguttata_gene_ensembl", mart = db)
gg <- useDataset("ggallus_gene_ensembl", mart = db)

##変換元の属性名を検索
filters <- listFilters(gg) ##引数はfaなど
filters[grep("refseq", filters[,1]),] ##input(すなわち、REFSEQ)はNCBIから撮ってきたXP_xxxxxx(タンパク質のAA配列)
##変換先は？遺伝子名をEnsembl outputに直したいのでensembl_peptide_id。元のrefseq_peptide(_predicted)も入れておこう。
attributes <- listAttributes(gg)
attributes[grep("", attributes[,1]),] ##ensembl,go_id,seqで検索

##3.それじゃあ変換しよう
##念の為三種のアノテーション情報を用いる。entrezgene_idよりもgo_idが欲しかったら、attributesを変更してね
resFa <- getBM(attributes =c("refseq_peptide_predicted", "go_id"),  #attributesに入れる値は、listAttributes()で調べたやつ
               filters = "refseq_peptide_predicted", #filters に入れる値は、listFilters()で探したやつ
               values = Falb$refseq, #valuesには、変換元のデータが入っている列を指定。
               mart = fa)

resPm <- getBM(attributes =c("refseq_peptide_predicted", "go_id"),  #attributesに入れる値は、listAttributes()で調べたやつ
               filters = "refseq_peptide_predicted", #filters に入れる値は、listFilters()で探したやつ
               values = Pmaj$refseq, #valuesには、変換元のデータが入っている列を指定。
               mart = pm)

resTg <- getBM(attributes =c("refseq_peptide_predicted", "go_id"),  #attributesに入れる値は、listAttributes()で調べたやつ
               filters = "refseq_peptide_predicted", #filters に入れる値は、listFilters()で探したやつ
               values = Tgut$refseq, #valuesには、変換元のデータが入っている列を指定。
               mart = tg)

resGg <- getBM(attributes =c("refseq_peptide_predicted", "go_id"),  #attributesに入れる値は、listAttributes()で調べたやつ
               filters = "refseq_peptide_predicted", #filters に入れる値は、listFilters()で探したやつ
               values = Ggal$refseq, #valuesには、変換元のデータが入っている列を指定。
               mart = gg)


##これで、遺伝子にEntrezgene_idを紐付けできた。
res1 <- bind_rows(resFa, resPm)
res2 <- bind_rows(res1, resTg)
##res <- bind_rows(res2, resGg)
##resはアノテーションできた外群全部のせ。重複行を削除すれば、どんなGOが出せたかわかる。
uniqueGO <- res2 %>% distinct(go_id) %>%
  filter(str_detect(go_id, "GO"))

uniqueGO %>% write.table("biomaRt/crowsniwatori/niwatoriGO.txt",
                           quote = F, col.names = F, row.names = F)

#class(All4$refseq)
background <- full_join(All4, res2, by = c(refseq = "refseq_peptide_predicted")) %>%
  na.omit() %>%
  dplyr::select("V1", "go_id")
##backgroundはOrthoFinderで得られた全てのOrthogroupについて、外群三種のgo_idが紐付けされたやつ

Kyoumi1 <- read_tsv("../cafe/summary/0720summary/CommonGF.txt", col_names = "Orthogroup") %>%
  separate(Orthogroup, into=c("number", "Orthogroup")) %>% 
  dplyr::select(Orthogroup) %>%
  filter(Orthogroup != "object")

Kyoumi2 <- Kyoumi1 %>% left_join(background, by = c(Orthogroup = "V1"))

Kyoumi3 <- background %>% left_join(Kyoumi2, by = "go_id") %>%
  mutate(bool = is.na(Orthogroup)) %>%
  filter(bool == FALSE) %>%
  dplyr::select(!c(bool, V1)) %>%
  relocate(Orthogroup) %>%
  distinct(go_id, .keep_all = TRUE)

Kyoumi3 %>% write.csv("biomaRt/crowsniwatori/niwatoriGoterm.csv",
                      row.names = F)

##Kyoumi3は興味のあるクレード、生物種で共通(CommonGF.txtの場合)または特異的(SpecificGF.txtの場合)

##NULL次の目標は、background-Kyoumi2に1、Kyoumi2に0をふって、topGOdataを作成すること
#nouvelle <- background %>% left_join(Kyoumi2, by = "entrezgene_id") #%>%
#  mutate(bool = is.na(Orthogroup)) %>%
#  mutate(detection = if_else(bool, 1, 0))

#nouvelle1 <- nouvelle %>% dplyr::select(entrezgene_id) %>%
#  t() %>% as.character()
#nouvelle2 <- nouvelle %>% dplyr::select(detection) %>%
#  t() %>% as.numeric()
#Genes <- setNames(nouvelle2, nouvelle1)

topgo <- new("topGOdata",
             ontology = "BP",           # オントロジー
             allGenes = Genes,            # FDR 値
             geneSel = function(x) {x < 1},
             annot = annFUN.org,        # アノテーション
             mapping = "org.Gg.eg.db",  # アノテーションデータ
             ID = "entrez")

##otamesi
#x <- c(1,2,3,4,5)
#y <- c(101819515, 100222674, 101817172, 100217483, 100220410)
#xy <- setNames(x, y)

##おそらく、mappingに入れる生物種とallGenesに入っているEntrezgene_idを使用した生物種が一致している必要がある。
##次は、Entrezgene_idを"org.Gg.eg.db"にするために、Gallus_gallusで置き換える→だめでした
##Heavy Watalさんのブログのやつではなぜ成功するのだ廊下。