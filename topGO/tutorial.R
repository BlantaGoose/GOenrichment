install.packages("BiocManager")
BiocManager::install(c("Biostrings", "GenomicRanges", "rtracklayer"))
BiocManager::install("topGO")

#To view documentation for the version of this package
browseVignettes("topGO")

#1. Data preparation。R objectの中にgene-to-GO annotationなど色々集めたりする
#CAFEデータを準備し、topGOdataの作成
tg_data = new("topGOdata",
              ontology = "", #3 term input
              allGenes = , #named vector
              geneSelectionFun = function(x) {x < 0.01},
              nodeSize = 10, #filthering with minimum size of gene
              annotationFun = annFUN.gene2GO,
              ID = "ensembl")