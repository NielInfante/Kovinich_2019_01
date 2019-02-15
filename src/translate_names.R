library(tidyverse)

# Read Previously done DE
d <- read.table('~/depot/projects/Kovinich/Kovinich_2019_01/previousWork/WGEvsH2O.txt',
								header=T, sep="\t", stringsAsFactors = F)


# Take only a portion of the glyma 1.1 ID
d$GLYMA_1.1 <- substr(d$GLYMA_1.1, 7, 14)

# Read the translation table
trans <- read.table('~/depot/projects/Kovinich/Kovinich_2019_01/Data/Glyma_11_to_Glyma_20_Correspondence_Full.csv',
									header=T, sep="\t", skip=2, stringsAsFactors = F)


#head(trans)

# get IDs into usable form
trans$Glyma.1.1 <- substr(trans$Glyma.1.1,6,99)
trans$Glyma.1.1 <- toupper(trans$Glyma.1.1)

trans$Glyma2.0 <- toupper(trans$Glyma2.0)
trans$Glyma2.0 <- sub('[.]','_',trans$Glyma2.0)


# Join the IDs with the data table
j <- inner_join(d, trans, by= c('GLYMA_1.1'='Glyma.1.1'))


write_tsv(j,'~/depot/projects/Kovinich/Kovinich_2019_01/previousWork/withID_WGEvsH2O.txt')




