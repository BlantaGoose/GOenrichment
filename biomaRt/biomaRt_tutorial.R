#biomaRtから利用可能な情報を調べて取り出そう
BiocManager::install("biomaRt")
library(biomaRt)

listMarts() #一列目のbiomartの部分を指定する。
db <- useMart("ensembl") #useMart()で使いたいプラットフォームの選択
listDatasets(db) #martで利用可能なデータセットを列挙。
hg <- useDataset("hsapiens_gene_ensembl", mart = db) #ヒトデータを抽出

#検索式
#まず、FiltersとAttributesにどんな属性を指定できるか調べる
listFilters(hg)
listAttributes(hg)

#Case1. 実際にいくつかのEnsembl IDを与え、GO_idを検索
ensid <- c("ENSG00000000003", "ENSG00000000457", "ENSG00000000938",
           "ENSG00000006327", "ENSG00000011405", "ENSG00000001497")
ensid
class(ensid)
db <- useMart("ensembl")
hg <- useDataset("hsapiens_gene_ensembl", mart = db)
res <- getBM(attributes = c("ensembl_gene_id", "go_id"), #attributes()でoutputで欲しい属性を指定
             filters = "ensembl_gene_id", #検索時に条件を指定。valueで指定したキーワードの属性
             values = ensid, #検索キーワード
             mart = hg) #useMartで作成したオブジェクト
head(res)
tail(res)

#Otameshi カラスクレードで実際に確認できた遺伝子をいくつか、ensidみたいな形にしてアノテーションできるかな
gg <- c("XP_420402")
gg
class(gg)
listMarts()
db <- useMart("ENSEMBL_MART_ENSEMBL")
listDatasets(db)
hg <- useDataset("ggallus_gene_ensembl", mart = db)

filters <- listFilters(hg)
filters[grep("seq", filters[,1]),]
re <- getBM(attributes = c("ensembl_gene_id"),#"refseq_peptide_predicted", "entrezgene_id", 
                           #"ensembl_gene_id", "go_id"),
            filters = "refseq_peptide_predicted",
            values = gg,
            mart = hg)
re

#Case2. Ensembl ID→Entrez IDを検索
res <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "go_id"),
             filters = "ensembl_gene_id", 
             values = ensid,
             mart = hg)
head(res)
tail(res)