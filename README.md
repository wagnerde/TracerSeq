TracerSeq
=========

This pipeline facilitates processing of scRNA-seq TracerSeq FASTQ files to produce a TracerClones x Cells counts table.  

These functions require that input FASTQ files have already been demultiplexed by sample (e.g. steps 1-3 of the inDrops.py pipeline). Inputs are expected to match the format of inDrops.py "sorted FASTQ" files.  

This code was written in Matlab (2017a) and requires the 'Statistics and Machine Learning Toolbox'.

TracerSeq was developed and used to analyze zebrafish embryonic development in:  
**Single-cell mapping of gene expression landscapes and lineage in the zebrafish embryo.**  Wagner DE, Weinreb C, Collins ZM, Briggs JA, Megason SG, Klein AM. Science 26 Apr 2018. [doi:10.1126/science.aar4362](http://science.sciencemag.org/content/early/2018/04/25/science.aar4362)

## Usage ##

1. First run the parse_TracerFastQ.m function to process raw FASTQ files. This function will filter abundant inDrops cell barcodes and UMIs, perform UMI error-validation, and then determine a consensus sequence for each unique TracerSeq barcode detected in each cell. Consensus sequences for each barcode are then writted to a csv file each tagged with its associated inDrops cell barcode.
2. Next run the parse_TracerClones.m function to perform TracerSeq barcode correction, assign barcodes to clones, and save a TracerSeq Clones x Cells counts table.  In this final table, each TracerSeq mRNA is a row with the following columns:
	column1: inDrops cell barcode 
	column2: TracerSeq clone # assignment
	column3: TracerSeq barcode sequence (corrected)
	column4: UMI counts

### Inputs ###

Input sequencing data must be sample demultiplexed such that each FASTQ file corresponds to a single inDrops library.  FASTQ formats should match that of the inDrops.py filtered/sorted FASTQ output (see below).  The second line contains the biological read, and the headeris formatted as follows: 
@:inDropsCellBarcodePart1-inDropsCellBarcodePart2:UMI:etc

Example:
'''
@CCAGACAG-CACAAGGC:CCAAAT:NS500422_445_H77LWBGX2_1_11101_12233_1072
GGCATNGATGAGCTCTACAAATAAGAGCCGAAAATAAAGATATAATCATACGTATCCGGAA
+
AAAAA#EEEEEEEEEEEEEEEEEEEEE6EEAEEEE/EAEAEAEEEAEEEEEEAEEAEEAEE
'''

### Settings ###

### Parse TracerSeq FASTQ Files ###











