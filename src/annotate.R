# script to add some annotations to results files.

library(AnnotationDbi)
library(AnnotationHub)

hub <- AnnotationHub()
query(hub, "Glycine max")
gm <- hub[['AH66196']]

columns(gm)


entrez <- read_tsv('~/depot/projects/Kovinich/Kovinich_2019_01/Data/entrez.tsv')
entrez <- entrez %>% dplyr::select(GeneID=gene_stable_id, Entrez=xref)
entrez$Entrez <- as.character(entrez$Entrez)


res <- read_tsv('~/depot/projects/Kovinich/Kovinich_2019_01/deseq/H2O_All_results.txt')

res <- inner_join(res, entrez)


my_concat <- function(x){paste(x, sep="|", collapse="|")}


# Get gene names
res$Symbol <- mapIds(gm, keys=res$Entrez, column='SYMBOL', keytype='ENTREZID', multiVals=my_concat)
res$GeneName <- mapIds(gm, keys=res$Entrez, column='GENENAME', keytype='ENTREZID', multiVals=my_concat)
res$Ontology <- mapIds(gm, keys=res$Entrez, column='ONTOLOGY', keytype='ENTREZID', multiVals=my_concat)
res$OntologyAll <- mapIds(gm, keys=res$Entrez, column='ONTOLOGYALL', keytype='ENTREZID', multiVals=my_concat)
res$GO <- mapIds(gm, keys=res$Entrez, column='GO', keytype='ENTREZID', multiVals=my_concat)
res$GOAll <- mapIds(gm, keys=res$Entrez, column='GOALL', keytype='ENTREZID', multiVals=my_concat)



write_tsv(res, '~/depot/projects/Kovinich/Kovinich_2019_01/deseq/H2O_All_results_annot.txt')
