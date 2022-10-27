## Downloads geneome sequence
For P. japonicus
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/017/312/705/GCF_017312705.1_Mj_TUMSAT_v1.0/GCF_017312705.1_Mj_TUMSAT_v1.0_genomic.fna.gz
For WSSV
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/848/085/GCF_000848085.2_ViralProj14616/GCF_000848085.2_ViralProj14616_genomic.fna.gz
combined two genome fasta files into one file.

## Creation of STAR index
```
STAR \
--runThreadN 40 \
--runMode genomeGenerate \
--genomeDir Dir_to_genome \
--genomeFastaFiles file_of_genome \
--sjdbGTFfeatureExon exon \
--sjdbGTFfile PjWSSV.gtf
```

## Mapping against genome sequence and count UMIs/Genes by STARsolo
```
STAR \
--runThreadN X \
--runMode alignReads \
--outSAMtype BAM SortedByCoordinate \
--sysShell /bin/bash \
--genomeDir Dir_to_index \
--readFilesIn Read2.fastq.gz Read1.fastq.gz \
--readFilesCommand 'gzip -c -d' \
--soloCBwhitelist None \
--soloType CB_UMI_Simple \
--soloCBmatchWLtype 1MM_multi \
--soloCBstart  1 \
--soloCBlen 12 \
--soloUMIstart 13 \
--soloUMIlen 8 \
--soloBarcodeReadLength 0 \
--sjdbGTFfile Dir_to_gtf \
--soloCellFilter CellRanger2.2 1500 0.99 10 \
--outFilterMultimapNmax 1 \
--outSAMattributes NH HI AS nM CB UB CR CY UR UY \
--outFileNamePrefix Dir_to_out \
--outReadsUnmapped Fastx \
--quantMode GeneCounts \
--bamRemoveDuplicatesType UniqueIdentical \
--outTmpDir Dir_to_tmp \
--soloFeatures Gene
```
