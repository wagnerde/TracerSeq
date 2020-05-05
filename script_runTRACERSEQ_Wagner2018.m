%% SCRIPT: Calculate TracerSeq State-Lineage Couplings
%
% This script contains all the necessary steps for downloading
% pre-processed counts matrices, and analyzing and plotting state-lineage
% couplings from the Wagner 2018 TracerSeq datasets.

%% Set up paths
% Create download directory
if ~exist('Downloads/', 'dir'); mkdir 'Downloads'; end

% Add paths to functions
addpath('private')
addpath('othercolor')

%% Download Wagner 2018 TracerSeq data from NCBI-GEO
% Download transcriptome counts files
websave('Downloads/TracerSeq1.csv.gz','https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3067nnn/GSM3067196/suppl/GSM3067196_TracerSeq1.csv.gz');
websave('Downloads/TracerSeq2.csv.gz','https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3067nnn/GSM3067197/suppl/GSM3067197_TracerSeq2.csv.gz');
websave('Downloads/TracerSeq3.csv.gz','https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3067nnn/GSM3067198/suppl/GSM3067198_TracerSeq3.csv.gz');
websave('Downloads/TracerSeq4.csv.gz','https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3067nnn/GSM3067199/suppl/GSM3067199_TracerSeq4.csv.gz');
websave('Downloads/TracerSeq5.csv.gz','https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3067nnn/GSM3067200/suppl/GSM3067200_TracerSeq5.csv.gz');

% Download clusterID files
websave('Downloads/TracerSeq1_clustID.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3067nnn/GSM3067196/suppl/GSM3067196_TracerSeq1_clustID.txt.gz');
websave('Downloads/TracerSeq2_clustID.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3067nnn/GSM3067197/suppl/GSM3067197_TracerSeq2_clustID.txt.gz');
websave('Downloads/TracerSeq3_clustID.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3067nnn/GSM3067198/suppl/GSM3067198_TracerSeq3_clustID.txt.gz');
websave('Downloads/TracerSeq4_clustID.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3067nnn/GSM3067199/suppl/GSM3067199_TracerSeq4_clustID.txt.gz');
websave('Downloads/TracerSeq5_clustID.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3067nnn/GSM3067200/suppl/GSM3067200_TracerSeq5_clustID.txt.gz');

% Download TracerSeq counts files
websave('Downloads/TracerSeq1_CellTracerCounts.csv.gz','https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3067nnn/GSM3067196/suppl/GSM3067196_TracerSeq1_CellTracerCounts.csv.gz');
websave('Downloads/TracerSeq2_CellTracerCounts.csv.gz','https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3067nnn/GSM3067197/suppl/GSM3067197_TracerSeq2_CellTracerCounts.csv.gz');
websave('Downloads/TracerSeq3_CellTracerCounts.csv.gz','https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3067nnn/GSM3067198/suppl/GSM3067198_TracerSeq3_CellTracerCounts.csv.gz');
websave('Downloads/TracerSeq4_CellTracerCounts.csv.gz','https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3067nnn/GSM3067199/suppl/GSM3067199_TracerSeq4_CellTracerCounts.csv.gz');
websave('Downloads/TracerSeq5_CellTracerCounts.csv.gz','https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3067nnn/GSM3067200/suppl/GSM3067200_TracerSeq5_CellTracerCounts.csv.gz');

% Download transcriptome cluster ID/name key
websave('Downloads/ClusterNames.csv.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE112nnn/GSE112294/suppl/GSE112294_ClusterNames.csv.gz');

% Unzip files
%gunzip('Downloads/*.gz');
system('gunzip -f Downloads/*.gz');

%% Load data
% Create the 'DataSet' structure.  Each record in this data structure
% corresponds to a single biological replicate (an individual embryo) for
% which both transcriptome and TracerSeq-targeted inDrops libraries have 
% been sequenced. We will import cluster IDs from the transcriptome
% dataset, and a counts matrix from the TracerSeq dataset. inDrops cell
% barcodes will also be imported for each.

% specify DataSet names
DataSet(1).name = 'TracerSeq1';
DataSet(2).name = 'TracerSeq2';
DataSet(3).name = 'TracerSeq3';
DataSet(4).name = 'TracerSeq4';
DataSet(5).name = 'TracerSeq5';

% read tracer and txome data from files
for j = 1:length(DataSet)
    % import custom TracerSeq file 'CellTracerCounts.csv'
    [DataSet(j).imported.tracer_cell_barcodes, DataSet(j).imported.tracer_counts_matrix] = load_CellTracerCounts(['Downloads/' DataSet(j).name '_CellTracerCounts.csv']);
    % cell barcodes from the transcriptome dataset can be imported as column 
    % headers from the genes x cells counts matrix
    DataSet(j).imported.txome_cell_barcodes = import_csv_header(['Downloads/' DataSet(j).name '.csv']);
    % import clusterIDs for the transcriptome dataset - clusterIDs 
    % correspond to columns from the above genes x cells counts matrix
    DataSet(j).imported.txome_clustID = importdata(['Downloads/' DataSet(j).name '_clustID.txt']);
end

%% Merge tracer and txome data based on shared inDrops cell barcodes
% Use inDrops cell barcodes to index individual cells captured in both the
% TracerSeq and transcriptome inDrops libraries
DataSet = merge_tracer_data(DataSet);

%% Plot lineage couplings between cell states
tmp.thresh_UMI = 1; % minimum number of Tracer UMIs per cell
tmp.min_cells_per_clone = 5;  % lower size limit for clones
tmp.max_cells_per_clone = Inf; % no upper size limit for clones
tmp.min_cells_per_hit = 2; % require at least 2 cells to couple to a state
tmp.nRandTrials = 1000; 

get_tracer_couplings(DataSet, readtable('Downloads/ClusterNames.csv'), 'thresh_UMI', tmp.thresh_UMI, 'thresh_min_cells_per_clone', tmp.min_cells_per_clone, 'thresh_max_cells_per_clone', tmp.max_cells_per_clone, 'thresh_min_cells_per_hit', tmp.min_cells_per_hit, 'nRandTrials', tmp.nRandTrials);

