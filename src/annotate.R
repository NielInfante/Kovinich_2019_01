# script to add some annotations to results files.

library(tidyverse)
library(AnnotationDbi)
library(AnnotationHub)

hub <- AnnotationHub()
query(hub, "Glycine max")
gm <- hub[['AH66196']]

#columns(gm)


entrez <- read_tsv('~/depot/projects/Kovinich/Kovinich_2019_01/Data/entrez.tsv')
entrez <- entrez %>% dplyr::select(GeneID=gene_stable_id, Entrez=xref)
entrez$Entrez <- as.character(entrez$Entrez)

my_concat <- function(x){paste(x, sep="|", collapse="|")}


doTheAnnotation('~/depot/projects/Kovinich/Kovinich_2019_01/deseq/All_results.txt')


doTheAnnotation('~/depot/projects/Kovinich/Kovinich_2019_01/deseq/H2O_MYB29A2_vs_pGWB2_results.txt') 
doTheAnnotation('~/depot/projects/Kovinich/Kovinich_2019_01/deseq/WGE_MYB29A2_vs_pGWB2_results.txt')
doTheAnnotation('~/depot/projects/Kovinich/Kovinich_2019_01/deseq/H2O_All_results.txt')  
doTheAnnotation('~/depot/projects/Kovinich/Kovinich_2019_01/deseq/H2O_NAC42_vs_pGWB2_results.txt')    
doTheAnnotation('~/depot/projects/Kovinich/Kovinich_2019_01/deseq/WGE_NAC42_vs_pGWB2_results.txt')

# Overlap Files
doTheAnnotation('~/depot/projects/Kovinich/Kovinich_2019_01/overlap/H2O_WGE_MYB29A2.txt')
doTheAnnotation('~/depot/projects/Kovinich/Kovinich_2019_01/overlap/H2O_WGE_NAC42.txt')
doTheAnnotation('~/depot/projects/Kovinich/Kovinich_2019_01/overlap/NAC42_MYB29A2_H2O.txt')
doTheAnnotation('~/depot/projects/Kovinich/Kovinich_2019_01/overlap/NAC42_MYB29A2_WGE.txt')
doTheAnnotation('~/depot/projects/Kovinich/Kovinich_2019_01/overlap/wge_MYB29A2_old.txt')
doTheAnnotation('~/depot/projects/Kovinich/Kovinich_2019_01/overlap/wge_old.txt')									

doTheAnnotation('~/depot/projects/Kovinich/Kovinich_2019_01/overlap/old_Up_wge_MYB29A2_not.txt')									
doTheAnnotation('~/depot/projects/Kovinich/Kovinich_2019_01/overlap/old_Up_wge_MyvpG_not.txt')									
doTheAnnotation('~/depot/projects/Kovinich/Kovinich_2019_01/overlap/wge_MYB29A2_Up_old_not.txt')									
doTheAnnotation('~/depot/projects/Kovinich/Kovinich_2019_01/overlap/wge_MyvpG_Up_old_not.txt')									



# Define the Function
doTheAnnotation <- function(fileName){

	res <- read_tsv(fileName)

	res <- inner_join(res, entrez)

	# Get gene names
	res$Symbol <- mapIds(gm, keys=res$Entrez, column='SYMBOL', keytype='ENTREZID', multiVals=my_concat)
	res$GeneName <- mapIds(gm, keys=res$Entrez, column='GENENAME', keytype='ENTREZID', multiVals=my_concat)
	res$Ontology <- mapIds(gm, keys=res$Entrez, column='ONTOLOGY', keytype='ENTREZID', multiVals=my_concat)
	res$OntologyAll <- mapIds(gm, keys=res$Entrez, column='ONTOLOGYALL', keytype='ENTREZID', multiVals=my_concat)
	res$GO <- mapIds(gm, keys=res$Entrez, column='GO', keytype='ENTREZID', multiVals=my_concat)
	res$GOAll <- mapIds(gm, keys=res$Entrez, column='GOALL', keytype='ENTREZID', multiVals=my_concat)

	
	write_tsv(res, paste0(sub('\\.txt$', '', fileName), '_annot.txt'))
}
