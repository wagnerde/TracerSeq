function parse_TracerClones(library_set, set_name, varargin)
% Usage: parse_TracerClones(library_set, set_name, varargin)
%
% This function performs TracerSeq barcode correction, assigns TracerSeq
% mRNAs to clones, and saves a TracerSeq x Cells counts table. Run this 
% function after first parsing the raw FASTQ files into Tracer consensus 
% sequences (i.e. 'libname'_tracerSeqs.csv). Multiple consensus files can 
% be merged.
% 
% Matlab must be run in the same working directory as the input csv files
% The Matlab PATH must include the following functions:
%   get_barcode_edits.m
%   strdist.m
%   plot_edit_dist_heatmap.m
%   subtitle.m
%
%
% REQUIRED INPUTS:
% library_set      Cell array of libnames to be processed; each libname
%                  corresponds to a 'libname'_tracerSeqs.csv file
% set_name         String identifier for this set of libraries 
%
% OPTIONAL INPUT NAME/VALUE PAIRS:
% 'edit_dist_method'
%           Method for measuring sequence distance between barcodes; argument 
%           is passed to strdist. Acceptable inputs: 'edit' or 'levenshtein'.
%           Edit: substitutions=2,insertions/deletions=1
%           Levenshtein: substitutions=1,insertions/deletions=1
%           (default='levenshtein')
%
% 'edit_dist_thresh'
%           Sequence distance threshold required for calling an an edit 
%           relationship between two barcodes.
%           (default=2)
% 
% 'edit_abund_thresh'
%           Counts ratio threshold required for calling an an edit 
%           relationship between two barcodes.
%           (default=2)
% 
% 'plot_distances'
%           Boolean flag for whether or not to plot sequence distances 
%           between Tracer barcodes. 
%           (default=true)
%
% 'plot_distances_matrix_size'
%           If plotting sequence distances, this parameter sets the number of
%           sequences to use for the distance matrix.
%           (default=1000)
%
% OUTPUTS:
% All files are written to the current working directory
% CellTracerCounts.csv 
%           TracerSeq x Cells counts table, each TracerSeq mRNA is a row
%           column1: inDrops cell barcode 
%           column2: TracerSeq clone # assignment
%           column3: TracerSeq barcode sequence (corrected)
%           column4: UMI counts
% 
% TracerNetworkGraph.png
%           Plot of the directed network graph used to collapse barcode edits.
%           Each node of the graph is one of the *original* TracerSeq barcodes.
%           Edges connect barcodes that are within edit_dist_thresh and whose
%           abundances differ by at least factor a factor of edit_abund_thresh.
% 
% BarcodeEdits.csv
%           A table summarizing the barcode edit search process. See documentation 
%           of the 'get_barcode_edits' function for further details.
%
% BarcodeCorrectionTable.csv
%           The lookup table used to update original TracerSeq barcodes to their 
%           corrected sequences.
% 
% EditDistHeatmap.png
%           A clustered heatmap of a pairwise edit distances between TracerSeq 
%           barcodes.
% 
% EditDistHist.png
%           A histogram of pairwise edit distances between TracerSeq barcodes.
%

%% PARAMETER SETTINGS
% Set defaults
def.edit_dist_thresh = 2;
def.edit_dist_method = 'levenshtein';
def.edit_abund_thresh = 2;
def.plot_distances = true;
def.plot_distances_matrix_size = 1000;
% Create parser object
parserObj = inputParser;
parserObj.FunctionName = 'parse_TracerClones';
parserObj.StructExpand = false; 
parserObj.addOptional('edit_dist_thresh',def.edit_dist_thresh);
parserObj.addOptional('edit_dist_method',def.edit_dist_method);
parserObj.addOptional('edit_abund_thresh',def.edit_abund_thresh);
parserObj.addOptional('plot_distances',def.plot_distances);
parserObj.addOptional('plot_distances_matrix_size',def.plot_distances_matrix_size);
% Parse input options
parserObj.parse(varargin{:});
settings = parserObj.Results;

%% Set up paths
if ~exist('plots', 'dir'); mkdir 'plots'; end

%% Display pipeline targets:
disp(['Identifying TracerSeq clones for ' set_name])
disp([set_name ' libraries:']); disp(library_set);

%% Load data:
original_tracerSeqs = {};
original_cellBcodes = {};
for j = 1:length(library_set)
    fileID = fopen([library_set{j} '_tracerSeqs.csv']);
    next_file = textscan(fileID,'%s %s','Delimiter',',');
    original_tracerSeqs = [original_tracerSeqs; next_file{2}];
    % append library name as a prefix to the cell barcode
    library_prefices_next = repmat({[library_set{j} '-']},[size(next_file{1},1),1]);
    cellBcodes_next = cellstr(cell2mat([library_prefices_next, next_file{1}]));
    original_cellBcodes = [original_cellBcodes; cellBcodes_next];
    fclose(fileID);    
end
disp('done loading data')

%% Eliminate barcodes less than 18bp long
short = (cellfun(@length, original_tracerSeqs)<18);
original_tracerSeqs(short)=[]; 
original_cellBcodes(short)=[];

%% Find barcodes with exact matches
[unique_tracerSeqs, ~, unique_tracerSeqs_ic] = unique(original_tracerSeqs);
unique_tracerSeq_ind = 1:length(unique_tracerSeqs);
unique_tracerSeq_counts = histc(unique_tracerSeqs_ic,unique_tracerSeq_ind);
nUniqueTracerSeqs = length(unique_tracerSeqs);

%% Find Tracer barcodes that are likely to be an 'edited' version of 
% another more abundant barcode sequence
[~, BarcodeEditTable] = get_barcode_edits(unique_tracerSeqs, unique_tracerSeq_counts, settings.edit_abund_thresh, settings.edit_dist_thresh, settings.edit_dist_method);
writetable(BarcodeEditTable, [set_name '_BarcodeEdits.csv'],'WriteVariableNames',1);
disp('calculated edit distances')

%% Use a network graph to group together 'neighborhoods' of edited barcodes
% a single original barcode could undergo multiple rounds of editing
% each node is a Tracer barcode
% directed edges connect parent barcodes to edited versions
G = digraph(false(nUniqueTracerSeqs,nUniqueTracerSeqs));
S_edit = BarcodeEditTable.original_ind(~BarcodeEditTable.valid);
T_edit = BarcodeEditTable.parent_original_ind(~BarcodeEditTable.valid);
for j = 1:length(S_edit)
    if findedge(G, S_edit(j), T_edit(j)) == 0
        G = addedge(G, S_edit(j), T_edit(j));
    end
end

%% Assemble a barcode error 'correction' table for updating sequences;
% barcodes belonging to the same graph connected component are clones;
% within each clone, the most abundant barcode is the 'corrected' 
% nucleotide sequence for all members of the clone
CorrectionTable = table(unique_tracerSeq_ind', unique_tracerSeq_ind', unique_tracerSeqs, unique_tracerSeq_counts, cell(nUniqueTracerSeqs,1), 'VariableNames', {'OriginalID', 'CloneID', 'TracerSeq_orig', 'nCountsAllCells', 'TracerSeq_cor'});
CorrectionTable.CloneID = conncomp(G, 'Type', 'Weak')';
nClones = length(unique(CorrectionTable.CloneID));
for j = 1:nClones
    rows_for_this_clone = find(CorrectionTable.CloneID==j);
    [~, max_count_ind] = max(CorrectionTable.nCountsAllCells(rows_for_this_clone));
    CorrectionTable.TracerSeq_cor(rows_for_this_clone) = CorrectionTable.TracerSeq_orig(rows_for_this_clone(max_count_ind));
end
writetable(CorrectionTable, [set_name '_BarcodeCorrectionTable.csv'],'WriteVariableNames',1);

% Use the above table to perform barcode correction on all original TracerSeq mRNAs
% Final UMI counts are the number of times each *corrected* clone-cell pair was observed
for j = 1:length(original_tracerSeqs)
    loc_in_correction_table = strcmp(original_tracerSeqs(j), CorrectionTable.TracerSeq_orig);
    original_CloneID(j,1) = CorrectionTable.CloneID(loc_in_correction_table);    
    original_TracerSeq_cor(j,1) = CorrectionTable.TracerSeq_cor(loc_in_correction_table);    
end
[FinalTable_CloneCellPairs, ~, FinalTable_CloneCellPairs_ic] = unique(table(original_cellBcodes, original_CloneID, original_TracerSeq_cor));
FinalTable_CloneCellPairs_ind = 1:height(FinalTable_CloneCellPairs);
CloneCellPair_UMIcounts = histc(FinalTable_CloneCellPairs_ic,FinalTable_CloneCellPairs_ind);
FinalTable_CloneCellPairs = [FinalTable_CloneCellPairs table(CloneCellPair_UMIcounts)];
FinalTable_CloneCellPairs.Properties.VariableNames = {'Cell_Barcode', 'TracerSeq_CloneID', 'TracerSeq_CorrectedBarcode', 'UMI_Counts'};
writetable(FinalTable_CloneCellPairs, [set_name '_CellTracerCounts.csv'],'WriteVariableNames',0);
disp('done parsing clones')

%% Plot network graph
figure('Position', [100 600 600 600]);
ax1 = gca; ax1.Position = [0 0 1 1];
f1 = gcf; f1.Color = [0.95 0.95 0.95 0.95];
P = plot(G, 'Layout', 'force','Iterations',100, 'LineWidth', 2);
P.EdgeColor = 'black';
axis off
% node size is proportional to barcode abundance
P.MarkerSize = (log10(BarcodeEditTable.original_counts)+0.1)*2;
% print stats to plot
nEdges_total = height(G.Edges);
subtitle({['nTracers: ' num2str(nUniqueTracerSeqs)],...
          ['nClones: ' num2str(nClones)],...
          ['nEdges: ' num2str(nEdges_total)]},...
           'TopLeft',[0.02 0.01]);
saveas(gcf, ['plots/' set_name '_TracerNetworkGraph.png']) 
  
%% Plot histograms: # clones per cell 
[~, ~, each_cell_ic] = unique(original_cellBcodes);
each_cell = unique(each_cell_ic);
for j = 1:length(each_cell)
    original_rows = (each_cell(j) == each_cell_ic);
    clones_in_this_cell{j,1} = unique(original_CloneID(original_rows));
    nClones_in_this_cell(j) = length(clones_in_this_cell{j,1});
end

figure;
[counts, centers] = hist(nClones_in_this_cell,max(nClones_in_this_cell));
bar(centers,log10(counts))
set(gca,'fontsize',12)
xlabel('# Clones per Cell')
ylabel('# Cells (log10)')
title([set_name ': Clones per Cell']);
st = subtitle({['Mean clones per cell = ', num2str(mean(nClones_in_this_cell))],...
               ['Median clones per cell = ', num2str(median(nClones_in_this_cell))],...
               ['Max clones per cell = ' num2str(max(nClones_in_this_cell))]},'TopRight',[0.02 0.01]);
saveas(gcf, ['plots/' set_name '_ClonesPerCellLog.png']) 

figure;
bar(centers, counts)
set(gca,'fontsize',12)
xlabel('# Clones per Cell')
ylabel('# Cells')
title([set_name ': Clones per Cell']);
st = subtitle({['Mean clones per cell = ', num2str(mean(nClones_in_this_cell))],...
               ['Median clones per cell = ', num2str(median(nClones_in_this_cell))],...
               ['Max clones per cell = ' num2str(max(nClones_in_this_cell))]},'TopRight',[0.02 0.01]);
saveas(gcf, ['plots/' set_name '_ClonesPerCell.png']) 


%% Plot histograms: # cells per clone 
[~, ~, each_clone_ic] = unique(original_CloneID);
each_clone = unique(each_clone_ic);
for j = 1:length(each_clone)
    original_rows = (each_clone(j) == each_clone_ic);
    cells_in_this_clone{j,1} = unique(original_cellBcodes(original_rows));
    nCells_in_this_clone(j) = length(cells_in_this_clone{j,1});
end

figure
[counts, centers] = hist(nCells_in_this_clone,max(nCells_in_this_clone));
bar(centers, log10(counts))
set(gca,'fontsize',12)
xlabel('# Cells per Clone')
ylabel('# Clones (log10)')
title([set_name ': Cells per Clone']);
st = subtitle({['Mean cells per clone = ', num2str(mean(nCells_in_this_clone))],...
               ['Median cells per clone = ', num2str(median(nCells_in_this_clone))],...
               ['Max cells per clone = ' num2str(max(nCells_in_this_clone))]},'TopRight',[0.02 0.01]);
saveas(gcf, ['plots/' set_name '_CellsPerCloneLog.png']) 

figure
bar(centers, counts)
set(gca,'fontsize',12)
xlabel('# Cells per Clone')
ylabel('# Clones')
title([set_name ': Cells per Clone']);
st = subtitle({['Mean cells per clone = ', num2str(mean(nCells_in_this_clone))],...
               ['Median cells per clone = ', num2str(median(nCells_in_this_clone))],...
               ['Max cells per clone = ' num2str(max(nCells_in_this_clone))]},'TopRight',[0.02 0.01]);
saveas(gcf, ['plots/' set_name '_CellsPerClone.png']) 

%% Plot heatmap of Tracer-Tracer edit distances
% this is computationally expensive, but makes for a nice visual
% construct distance matrix for a random subset of the original Tracer 
% barcodes- illustrates how well separated the original Tracers are in 
% sequence space, even prior to any barcode correction.
% then repeat for the corrected barcodes
if settings.plot_distances
    
    % get a random subset of the original Tracer sequences
    rand_ind = randperm(length(original_tracerSeqs), settings.plot_distances_matrix_size);
    original_tracerSeqs_subset = original_tracerSeqs(rand_ind);
    
    % plot heatmap for the original sequences
    seqdist_matrix_original = plot_edit_dist_heatmap(original_tracerSeqs_subset, [set_name ': Raw TracerSeq Edit Distances']);
    saveas(gcf, ['plots/' set_name '_EditDistHeatmap.png']) 
    
    % plot heatmap for the corrected sequences
    for j = 1:settings.plot_distances_matrix_size
        correction_ind = find(strcmp(original_tracerSeqs_subset{j},CorrectionTable.TracerSeq_orig));
        corrected_tracerSeqs_subset{j} = CorrectionTable.TracerSeq_cor{correction_ind};
    end
    seqdist_matrix_corrected = plot_edit_dist_heatmap(corrected_tracerSeqs_subset, [set_name ': Corrected TracerSeq Edit Distances']);
    saveas(gcf, ['plots/' set_name '_EditDistHeatmapCorrected.png']) 
    
end

%% Plot histogram of Tracer-Tracer edit distances
% histogram of pairwise distances from the matrices calculated above
if settings.plot_distances
    
    % only look at the upper triangle
    orig_seqdist_matrix_upper_only = seqdist_matrix_original(~tril(ones(size(seqdist_matrix_original))));
    corrected_seqdist_matrix_upper_only = seqdist_matrix_corrected(~tril(ones(size(seqdist_matrix_corrected))));

    % plot for original sequences
    figure;
    [counts, centers] = hist(orig_seqdist_matrix_upper_only,0:1:max(orig_seqdist_matrix_upper_only));
    bar(centers, counts)
    set(gca,'fontsize',12)
    xlim([-1 max(orig_seqdist_matrix_upper_only(:))+1]);
    xlabel('Edit Distance')
    ylabel('# Tracer Pairs')
    title([set_name ': Edit Distance Histogram']);
    saveas(gcf, ['plots/' set_name '_EditDistHist.png']) 
    
    % plot original (log10)
    figure;
    bar(centers, log10(counts))
    set(gca,'fontsize',12)
    xlim([-1 max(orig_seqdist_matrix_upper_only(:))+1]);
    xlabel('Edit Distance')
    ylabel('# Tracer Pairs (log10)')
    title([set_name ': Edit Distance Histogram (Log10)']);
    saveas(gcf, ['plots/' set_name '_EditDistHistLog.png']) 

    % plot for the corrected sequences
    figure;
    [counts, centers] = hist(corrected_seqdist_matrix_upper_only,0:1:max(corrected_seqdist_matrix_upper_only));
    bar(centers, counts)
    set(gca,'fontsize',12)
    xlim([-1 max(corrected_seqdist_matrix_upper_only(:))+1]);
    xlabel('Edit Distance')
    ylabel('# Tracer Pairs')
    title([set_name ': Edit Distance Histogram Corrected']);
    saveas(gcf, ['plots/' set_name '_EditDistHistCorrected.png']) 
  
    % plot corrected (log10)
    figure;
    bar(centers, log10(counts))
    set(gca,'fontsize',12)
    xlim([-1 max(corrected_seqdist_matrix_upper_only(:))+1]);
    xlabel('Edit Distance')
    ylabel('# Tracer Pairs (log10)')
    title([set_name ': Edit Distance Histogram Corrected (Log10)']);
    saveas(gcf, ['plots/' set_name '_EditDistHistCorrectedLog.png'])
      
end

disp('done plotting sequence distances')
