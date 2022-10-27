## Creation of Sequence Dictionary

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
