TracerSeq
=========

This pipeline facilitates processing of scRNA-seq TracerSeq FASTQ files to produce a TracerClones x Cells counts table.  

These functions require that input FASTQ files have already been demultiplexed by sample (e.g. steps 1-3 of the inDrops.py pipeline). Inputs are expected to match the format of inDrops.py "sorted FASTQ" files.  

This code was written in Matlab (2017a) and requires the 'Statistics and Machine Learning Toolbox'.

TracerSeq was developed and used to analyze zebrafish embryonic development in:  
**Single-cell mapping of gene expression landscapes and lineage in the zebrafish embryo.**  Wagner DE, Weinreb C, Collins ZM, Briggs JA, Megason SG, Klein AM. Science 26 Apr 2018. [doi:10.1126/science.aar4362](http://science.sciencemag.org/content/early/2018/04/25/science.aar4362)

### Inputs ###

Input sequencing data must be sample demultiplexed such that each FASTQ file corresponds to a single inDrops library. FASTQ formats should match that of the inDrops.py filtered/sorted FASTQ output (see below).  The second line contains the biological read, and the header is formatted as follows: ```@:inDropsCellBarcodePart1-inDropsCellBarcodePart2:UMI:AdditionalInfo```

Example:
```
@CCAGACAG-CACAAGGC:CCAAAT:NS500422_445_H77LWBGX2_1_11101_12233_1072
GGCATNGATGAGCTCTACAAATAAGAGCCGAAAATAAAGATATAATCATACGTATCCGGAA
+
AAAAA#EEEEEEEEEEEEEEEEEEEEE6EEAEEEE/EAEAEAEEEAEEEEEEAEEAEEAEE
```

### Usage ###

1. First run the parse_TracerFastQ.m function to process raw FASTQ files. This function will filter abundant inDrops cell barcodes and UMIs, perform UMI error-validation, and then determine a consensus sequence for each unique TracerSeq barcode detected in each cell. Consensus sequences for each barcode are then writted to a csv file each tagged with its associated inDrops cell barcode.
2. Next run the parse_TracerClones.m function to perform TracerSeq barcode correction, assign barcodes to clones, and save a TracerSeq Clones x Cells counts table.  In this final table, each TracerSeq mRNA is a row with the following columns:
	column1: inDrops cell barcode 
	column2: TracerSeq clone # assignment
	column3: TracerSeq barcode sequence (corrected)
	column4: UMI counts

The two main functions have both required and optional inputs. Each function generates outputs text files and plots.

parse_TracerFastQ.m   

```
 REQUIRED INPUTS:
 filename      Full path to an inDrops.py sorted FASTQ file
 libname       String identifier for this library

 OPTIONAL INPUT NAME/VALUE PAIRS:
 'thresh_cell'
           Minimum number of reads required to keep a cell barcode

 'thresh_UMI'
           Minimum number of reads required within a cell to keep
           a particular UMI

 OUPUTS:
 CellBC_WeightedHist.png
           Weighted histogram of reads per cell barcode and cutoff; 
           Use to inspect 'thresh_cell'.

 UMI_WeightedHist.png
           Weighted histogram of reads per UMI and cutoff; 
           Use to inspect 'thresh_UMI'.
 
 tracerSeqs.csv
           Table of identified barcodes, each TracerSeq mRNA is a row
           column1: inDrops cell barcode 
           column2: TracerSeq barcode sequence
```
parse_TracerClones.m   

```
 REQUIRED INPUTS:
 library_set      Cell array of libnames to be processed; each libname
                  corresponds to a 'libname'_tracerSeqs.csv file
 set_name         String identifier for this set of libraries 

 OPTIONAL INPUT NAME/VALUE PAIRS:
 'edit_dist_method'
           Method for measuring sequence distance between barcodes; argument 
           is passed to strdist. Acceptable inputs: 'edit' or 'levenshtein'.
           Edit: substitutions=2,insertions/deletions=1
           Levenshtein: substitutions=1,insertions/deletions=1
           (default='levenshtein')

 'edit_dist_thresh'
           Sequence distance threshold required for calling an an edit 
           relationship between two barcodes.
           (default=2)
 
 'edit_abund_thresh'
           Counts ratio threshold required for calling an an edit 
           relationship between two barcodes.
           (default=2)
 
 'plot_distances'
           Boolean flag for whether or not to plot sequence distances 
           between Tracer barcodes. 
           (default=true)

 'plot_distances_matrix_size'
           If plotting sequence distances, this parameter sets the number of
           sequences to use for the distance matrix.
           (default=1000)

 OUTPUTS:
 All files are written to the current working directory
 CellTracerCounts.csv 
           TracerSeq x Cells counts table, each TracerSeq mRNA is a row
           column1: inDrops cell barcode 
           column2: TracerSeq clone # assignment
           column3: TracerSeq barcode sequence (corrected)
           column4: UMI counts
 
 TracerNetworkGraph.png
           Plot of the directed network graph used to collapse barcode edits.
           Each node of the graph is one of the *original* TracerSeq barcodes.
           Edges connect barcodes that are within edit_dist_thresh and whose
           abundances differ by at least factor a factor of edit_abund_thresh.
 
 BarcodeEdits.csv
           A table summarizing the barcode edit search process. See documentation 
           of the 'get_barcode_edits' function for further details.

 BarcodeCorrectionTable.csv
           The lookup table used to update original TracerSeq barcodes to their 
           corrected sequences.
 
 EditDistHeatmap.png
           A clustered heatmap of a pairwise edit distances between TracerSeq 
           barcodes.
 
 EditDistHist.png
           A histogram of pairwise edit distances between TracerSeq barcodes.

```


### Running via command line ###











