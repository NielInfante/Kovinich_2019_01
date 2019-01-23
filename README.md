# RNAseq
RNA seq of Soya, looking at NAC42-1 and MYB29A2 overexpressions.


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



In the reads folder, use salmon to quantify each read


###  :large_orange_diamond: <!--:white_check_mark:--> Differential expression

Use the R script doDESeq.R to find differentially expressed genes, and produce some figures

### :x: <!--:large_orange_diamond: :white_check_mark:--> Emperor
Make an interactive PCA to explore data.

### :x: <!--:large_orange_diamond: :white_check_mark:--> GO Analysis

### :x: <!--:large_orange_diamond: :white_check_mark:-->  Pathway Analysis
