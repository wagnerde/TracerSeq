TracerSeq
=========

This pipeline facilitates processing of TracerSeq FASTQ files to produce a Tracer Barcodes x inDrops Cell Barcodes counts table.  Cell barcodes can then be used to cross reference TracerSeq data to accompanying transcriptome data for the same cells. As an example, we demonstrate how TracerSeq can be combined with transcriptome-derived cluster annotations to reveal state-lineage relationships.

This code was written in Matlab (2017a) and requires the 'Statistics and Machine Learning Toolbox'.  The [strdist](https://www.mathworks.com/matlabcentral/fileexchange/17585-calculation-of-distance-between-strings?focused=5094987&tab=function) function (Eduard Polityko) is used to calculate 'edit' distances between nucleotide barcode sequences. 

TracerSeq was developed and used to analyze zebrafish embryonic development in:  
**Single-cell mapping of gene expression landscapes and lineage in the zebrafish embryo.**  Wagner DE, Weinreb C, Collins ZM, Briggs JA, Megason SG, Klein AM. Science 26 Apr 2018. [doi:10.1126/science.aar4362](http://science.sciencemag.org/content/early/2018/04/25/science.aar4362)


## Preprocessing ##

### Sequencing File Formats ###
Input sequencing data must first be sample-demultiplexed such that each FASTQ file corresponds to a single inDrops library (e.g. steps 1-3 of the [inDrops.py](https://github.com/indrops/indrops) pipeline). FASTQ formats should match that of the inDrops.py filtered/sorted FASTQ output (see below). The second line contains the biological read, and the header is formatted as follows:    ```@:inDropsCellBarcodePart1-inDropsCellBarcodePart2:UMI:AdditionalInfo```

Example:
```
@CCAGACAG-CACAAGGC:CCAAAT:NS500422_445_H77LWBGX2_1_11101_12233_1072
GGCATNGATGAGCTCTACAAATAAGAGCCGAAAATAAAGATATAATCATACGTATCCGGAA
+
AAAAA#EEEEEEEEEEEEEEEEEEEEE6EEAEEEE/EAEAEAEEEAEEEEEEAEEAEEAEE
```

### Parsing FASTQ Files ###
1. First run **parse_TracerFastQ.m** to process raw FASTQ files. This function will filter abundant inDrops cell barcodes and UMIs, perform UMI error-validation, and then determine a consensus sequence for each unique TracerSeq barcode detected in each cell. Consensus sequences for each barcode are then writted to a csv file, each tagged with its associated inDrops cell barcode.
2. Next run **parse_TracerClones.m** to perform TracerSeq barcode correction, assign barcodes to clones, and save a TracerSeq Barcodes x Cell Barcodes counts table.  In this final table, each original TracerSeq mRNA is a row with the following fields:   
  column1: inDrops cell barcode   
  column2: TracerSeq clone # assignment   
  column3: TracerSeq barcode sequence (error-corrected)   
  column4: UMI counts   
   
The two functions have both required and optional inputs. Each function writes output text files and diagnostic plots to the working directory.

**parse_TracerFastQ.m**

```
 Usage: parse_TracerFastQ(filename, libname, varargin)

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
**parse_TracerClones.m**

```
 Usage: parse_TracerClones(library_set, set_name, varargin)

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

### Preprocess via command line ###

**parse_TracerFastQ.m** (run once per FASTQ file)
```
matlab -nodesktop -nodisplay -r "parse_TracerFastQ('/full_path_to.fastq','library_name','thresh_cell',min_reads_per_cell,'thresh_UMI',min_reads_per_UMI)"
```

**parse_TracerClones.m** (run once per sample, if necessary merge multiple libraries)
```
matlab -nodesktop -nodisplay -r "parse_TracerClones("{'libname_1' 'libname_2' 'libname_3' ... }",'set_name')"
```


## Calculate State-Lineage Couplings ##

We provide an example analysis pipeline: 'script_runTRACERSEQ_Wagner2018.m'. Execute it by typing the following into the Matlab command line:     
  ```
  run('script_runTRACERSEQ_Wagner2018.m')
  ```

This pipeline contains initial steps for downloading pre-processed
single-cell transcriptome annotations and CellTracerCounts.csv files from NCBI-GEO. CellTracerCounts files can also be generated locally from FASTQ files, as described above.  These data are then loaded into a 'DataSet' structure array with multiple fields in which each record corresponds to a biological sample (e.g. 'TracerSeq embryo 1'). The pipeline next merges Tracer counts with transcriptome annotations for individual inDrops cell barcodes by calling 'merge_tracer_data'. Finally, state-lineage couplings are calculated and plotted as a heatmap by calling 'get_tracer_couplings'.

Parameter settings for 'get_tracer_couplings' are specified using optional name/value pairs. Default behavior implements settings used in [Wagner et. al. 2018](http://science.sciencemag.org/content/early/2018/04/25/science.aar4362).

```
 'thresh_UMI'
       Minimum number of UMI counts required to call an individual cell as
       'positive' for a given Tracer clone
       (default=1)
 
 'thresh_min_cells_per_clone
       Discard any Tracer clones with fewer than this number of cells 
       (default=2)

 'thresh_max_cells_per_clone'
       Discard any Tracer clones with greater than this number of cells 
       (default=Inf)

 'thresh_min_cells_per_hit'
       Minimum number of positive cells required for a given TracerClone
       to be considered as 'hitting' a given state
       (default=1)

 'nRandTrials' 
       Number of data permutations used to generate mean and st. dev
       expectations
       (default=100)

 'heatmap_distance_metric'
       Distance metric argument passed to 'linkage' for generating final
       clustered heatmap
       (default='correlation')
 
 'heatmap_linkage_method'
       Method passed to 'linkage' for generating final clustered heatmap
       (default='average')

```
## Run TRACERSEQ ##  

