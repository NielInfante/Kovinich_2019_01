library(tidyverse)

d <- read.table('~/depot/projects/Kovinich/Kovinich_2019_01/previousWork/WGEvsH2O.txt',
								header=T, sep="\t", stringsAsFactors = F)


d$GLYMA_1.1 <- substr(d$GLYMA_1.1, 7, 14)
head(d)

n <- read.table('~/depot/projects/Kovinich/Kovinich_2019_01/deseq/WGE_NAC42_vs_pGWB2_results.txt',
								header=T, sep="\t", stringsAsFactors = F)



head(n)

trans <- read.table('~/depot/projects/Kovinich/Kovinich_2019_01/Data/Glyma_11_to_Glyma_20_Correspondence_Full.csv',
									header=T, sep="\t", skip=2, stringsAsFactors = F)


head(trans)
#trans$Glyma.1.1 <- toupper(trans$Glyma.1.1)
#trans$Glyma.1.1 <- sub('a','a.', trans$Glyma.1.1)
#trans$Glyma.1.1 <- sub('g','G', trans$Glyma.1.1)
#trans$Glyma.1.1 <- sub('s','S', trans$Glyma.1.1)
trans$Glyma.1.1 <- substr(trans$Glyma.1.1,6,99)
trans$Glyma.1.1 <- toupper(trans$Glyma.1.1)

trans$Glyma2.0 <- toupper(trans$Glyma2.0)
trans$Glyma2.0 <- sub('[.]','_',trans$Glyma2.0)

head(trans)
substr(d[1:5,1], 1,15)


j <- inner_join(d, trans, by= c('GLYMA_1.1'='Glyma.1.1'))



head(j)

dim(j)
