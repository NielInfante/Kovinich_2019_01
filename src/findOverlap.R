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


