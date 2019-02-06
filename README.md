# RNAseq
RNA seq of Soya, looking at NAC42-1 and MYB29A2 overexpressions.

From request:
can you please determine what genes are sig down and what is sig up in:
10) NAC42-pGWB2-H2Os / pGWB2-H2Os
11) NAC42-pGWB2-WGEs / pGWB2-WGEs
12) MYB29A2-pGWB2-H2Os / pGWB2-H2Os
13) MYB29A2-pGWB2-WGEs / pGWB2-WGEs
Also,
15) Overlap of upreg genes in 11 and 13
16) Overlap of upreg genes in 10 and 11
17) Overlap of upreg genes in 12 and 13
18) Nonoverlap and overlap of genes upreg in WGE hairy roots (old data) and 11
19) Nonoverlap and overlap of genes upreg in WGE hairy roots (old data) and 13




## Steps:

###  :white_check_mark: Quality Control
In the reads folder, create a fastqc folder. Then run:

```
for f in *.gz; do echo $f; fastq -o fastqc -t 40 $f; done
```

This makes QC reports for each file. To get a single report, cd to the fastqc folder and use:
```
multiqc -o multiqc .
```

Quality here is great, no need to do anything with it.


###  :white_check_mark: Quantification

Download the appropriate transcripts and build a salmon index

```
# Ensembl plant release 42
curl ftp://ftp.ensemblgenomes.org/pub/release-42/plants/fasta/glycine_max/cdna/Glycine_max.Glycine_max_v2.1.cdna.all.fa.gz -o GM_2.1_cdna.fa.gz

curl ftp://ftp.ensemblgenomes.org/pub/release-42/plants/fasta/glycine_max/ncrna/Glycine_max.Glycine_max_v2.1.ncrna.fa.gz -o GM_2.1_nc.fa.gz

gunzip *.gz
cat GM_2.1_cdna.fa GM_2.1_nc.fa > GM_2.1.fa

# activate salmon
conda activate salmon

# Build index
salmon index -t GM_2.1.fa -i GM_2.1 --type quasi -k 31

# Do the quantification and pretty up the names
for f in *R1*; do echo $f; out=${f#Nik-}; out=${out%_S*}; out=${out//-/_}; echo $out; salmon quant -i ../Data/GM_2.1 --libType A --gcBias -p 30 --numBootstraps 50 -o ../salmon/$out -1 $f -2 ${f/R1/R2}; done

conda deactivate
```




###  :white_check_mark: Differential expression

Set up the stats, and use DESeq2 to find differential expression.

##### H2O_NAC42_vs_pGWB2

```
PCA_Group <- 'Group'
design =~ Group
contrast <- c('Group','NAC42','pGWB2')
outPrefix <- 'H2O_NAC42_vs_pGWB2'

meta <- metadata %>% dplyr::filter(Group %in% c('NAC42', 'pGWB2'), Treatment=='H2O')
```


##### WGE_NAC42_vs_pGWB2
```
PCA_Group <- 'Group'
design =~ Group
contrast <- c('Group','NAC42','pGWB2')
outPrefix <- 'WGE_NAC42_vs_pGWB2'

meta <- metadata %>% dplyr::filter(Group %in% c('NAC42', 'pGWB2'), Treatment=='WGE')
```

##### H2O_MYB29A2_vs_pGWB2
```
PCA_Group <- 'Group'
design =~ Group
contrast <- c('Group','MYB29A2','pGWB2')
outPrefix <- 'WGE_MYB29A2_vs_pGWB2'

meta <- metadata %>% dplyr::filter(Group %in% c('MYB29A2', 'pGWB2'), Treatment=='H2O')
```

##### WGE_MYB29A2_vs_pGWB2
```
PCA_Group <- 'Group'
design =~ Group
contrast <- c('Group','MYB29A2','pGWB2')
outPrefix <- 'WGE_MYB29A2_vs_pGWB2'

meta <- metadata %>% dplyr::filter(Group %in% c('MYB29A2', 'pGWB2'), Treatment=='WGE')
```

###  :white_check_mark: Emperor
Make an interactive PCA to explore data.

### :large_orange_diamond:  <!--:white_check_mark:-->  Overlaps of Up regulated genes

Use findOverlaps.R.
Restrict to genes that are at least 2 fold upregulated, with an adjusted p-value of less than 0.05.



### :x: <!--:large_orange_diamond: :white_check_mark:--> GO Analysis

### :x: <!--:large_orange_diamond: :white_check_mark:-->  Pathway Analysis
