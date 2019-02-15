library(tidyverse)

# Script to find overlapping gene sets


# Want up in NAC42  vs pGWB2 in H2O
#  Overlap with
# Up in MYB29A2 vs pGWB2 in H2O


nac <- read_tsv('~/depot/projects/Kovinich/Kovinich_2019_01/deseq/H2O_NAC42_vs_pGWB2_results.txt')
myb <- read_tsv('~/depot/projects/Kovinich/Kovinich_2019_01/deseq/H2O_MYB29A2_vs_pGWB2_results.txt')


# nac - Positive fold change means up in NAC42
# myb - Positive fold change means up in MYB29A2


nac <- nac %>% select(Gene, NAC42_FC=log2FoldChange, NAC42_p=padj)
myb <- myb %>% select(Gene, MYB29A2_FC=log2FoldChange, MYB29A2_p=padj)

overlap <- inner_join(nac, myb)

overlap_out <- overlap %>% filter(NAC42_FC > 1 &
															NAC42_p < 0.05 &
															MYB29A2_FC > 1 &
															MYB29A2_p < 0.05)

write_tsv(overlap_out, '~/depot/projects/Kovinich/Kovinich_2019_01/overlap/NAC42_MYB29A2_H2O.txt')



# Want up in NAC42  vs pGWB2 in WGE
#  Overlap with
# Up in MYB29A2 vs pGWB2 in WGE


nac <- read_tsv('~/depot/projects/Kovinich/Kovinich_2019_01/deseq/WGE_NAC42_vs_pGWB2_results.txt')
myb <- read_tsv('~/depot/projects/Kovinich/Kovinich_2019_01/deseq/WGE_MYB29A2_vs_pGWB2_results.txt')

# nac - Positive fold change means up in NAC42
# myb - Positive fold change means up in MYB29A2

nac <- nac %>% select(Gene, NAC42_FC=log2FoldChange, NAC42_p=padj)
myb <- myb %>% select(Gene, MYB29A2_FC=log2FoldChange, MYB29A2_p=padj)

overlap <- inner_join(nac, myb)

overlap_out <- overlap %>% filter(NAC42_FC > 1 &
																		NAC42_p < 0.05 &
																		MYB29A2_FC > 1 &
																		MYB29A2_p < 0.05)

write_tsv(overlap_out, '~/depot/projects/Kovinich/Kovinich_2019_01/overlap/NAC42_MYB29A2_WGE.txt')




# Want up in NAC42  vs pGWB2 in H2O
#  Overlap with
# Want up in NAC42  vs pGWB2 in WGE


h2o <- read_tsv('~/depot/projects/Kovinich/Kovinich_2019_01/deseq/H2O_NAC42_vs_pGWB2_results.txt')
wge <- read_tsv('~/depot/projects/Kovinich/Kovinich_2019_01/deseq/WGE_NAC42_vs_pGWB2_results.txt')

# h20 - Positive fold change means up in NAC42
# wge - Positive fold change means up in NAC42

h2o <- h2o %>% select(Gene, H2O_FC=log2FoldChange, H2O_p=padj)
wge <- wge %>% select(Gene, WGE_FC=log2FoldChange, WGE_p=padj)

overlap <- inner_join(h2o, wge)

overlap_out <- overlap %>% filter(H2O_FC > 1 &
																		H2O_p < 0.05 &
																		WGE_FC > 1 &
																		WGE_p < 0.05)

write_tsv(overlap_out, '~/depot/projects/Kovinich/Kovinich_2019_01/overlap/H2O_WGE_NAC42.txt')





# Want up in NAC42  vs pGWB2 in H2O
#  Overlap with
# Want up in NAC42  vs pGWB2 in WGE


h2o <- read_tsv('~/depot/projects/Kovinich/Kovinich_2019_01/deseq/H2O_NAC42_vs_pGWB2_results.txt')
wge <- read_tsv('~/depot/projects/Kovinich/Kovinich_2019_01/deseq/WGE_NAC42_vs_pGWB2_results.txt')

# h20 - Positive fold change means up in NAC42
# wge - Positive fold change means up in NAC42

h2o <- h2o %>% select(Gene, H2O_FC=log2FoldChange, H2O_p=padj)
wge <- wge %>% select(Gene, WGE_FC=log2FoldChange, WGE_p=padj)

overlap <- inner_join(h2o, wge)

overlap_out <- overlap %>% filter(H2O_FC > 1 &
																		H2O_p < 0.05 &
																		WGE_FC > 1 &
																		WGE_p < 0.05)

write_tsv(overlap_out, '~/depot/projects/Kovinich/Kovinich_2019_01/overlap/H2O_WGE_NAC42.txt')



# Want up in MYB29A2  vs pGWB2 in H2O
#  Overlap with
# Want up in MYB29A2  vs pGWB2 in WGE


h2o <- read_tsv('~/depot/projects/Kovinich/Kovinich_2019_01/deseq/H2O_MYB29A2_vs_pGWB2_results.txt')
wge <- read_tsv('~/depot/projects/Kovinich/Kovinich_2019_01/deseq/WGE_MYB29A2_vs_pGWB2_results.txt')

# h20 - Positive fold change means up in MYB29A2
# wge - Positive fold change means up in MYB29A2

h2o <- h2o %>% select(Gene, H2O_FC=log2FoldChange, H2O_p=padj)
wge <- wge %>% select(Gene, WGE_FC=log2FoldChange, WGE_p=padj)

overlap <- inner_join(h2o, wge)

overlap_out <- overlap %>% filter(H2O_FC > 1 &
																		H2O_p < 0.05 &
																		WGE_FC > 1 &
																		WGE_p < 0.05)

write_tsv(overlap_out, '~/depot/projects/Kovinich/Kovinich_2019_01/overlap/H2O_WGE_MYB29A2.txt')



# Overlap and non-overlap of gene up in old WGE 
# and WGE NAC42 vs pGWB2

wge <- read_tsv('~/depot/projects/Kovinich/Kovinich_2019_01/deseq/WGE_NAC42_vs_pGWB2_results.txt')
old <- read_tsv('~/depot/projects/Kovinich/Kovinich_2019_01/previousWork/withID_WGEvsH2O.txt')

wge <- wge %>% select(Gene, WGE_FC=log2FoldChange, WGE_p=padj)
old <- old %>% select(Gene=Glyma2.0, OLD_FC=log2FoldChange, OLD_p=padj)

overlap <- inner_join(old, wge)

overlap_out <- overlap %>% filter(OLD_FC > 0 &
																		OLD_p < 0.05 &
																		WGE_FC > 0 &
																		WGE_p < 0.05)

write_tsv(overlap_out, '~/depot/projects/Kovinich/Kovinich_2019_01/overlap/')


ol <- ol %>% mutate(Up= ifelse()    )

# Venn style

ol <- overlap %>% filter(OLD_p < 0.05 | WGE_p < 0.05)

wge_up <- ol %>% filter(WGE_FC > 0) %>% select(Gene)
old_up <- ol %>% filter(OLD_FC > 0) %>% select(Gene)


library(VennDiagram)
png('pics/hallmark_venn.png')

vd <- venn.diagram(list("WGE Up"=wge_up$Gene, "Old Up"=old_up$Gene), 
									 fill=3:4, alpha=0.4, filename=NULL, cex=1.3, cat.cex=1.4,
									 main="")
grid.newpage()  ;  grid.draw(vd)

dev.off()








