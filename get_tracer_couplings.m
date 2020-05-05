function  Couplings = get_tracer_couplings(DataSet, ClusterNames, varargin)
%% Usage: Couplings = get_tracer_couplings(DataSet, ClusterNames, varargin)
% 
% This function
% 
% Couplings are based on a z-score vs random trials (cell state assignments 
% are shuffled across all cells from all datasets).
%
%
% REQUIRED INPUTS:
% 'DataSet'
%       A structure array with multiple fields, and one record for each 
%       biological sample. 'Merged' fields required must be first generated
%       by the 'merge_tracer_data' function.
%       DataSet.name
%               Short sample name, e.g. 'TracerSeq1'
%
%       DataSet.merged.tracer_counts_matrix
%               Matrix of TracerSeq UMI counts for cell barcodes that were 
%               also present in the transcriptome dataset.
%
%       DataSet.merged.tracer_cell_barcodes
%               List of inDrops cell barcodes, one barcode for each column
%               in the merged TracerSeq counts matrix
%
%       DataSet.merged.clustID
%               A numerical array of cell-clusterID annotations, one for 
%               each column in the merged TracerSeq counts matrix
%
% 'ClusterNames'
%       A table consisting of the following two columns, used to decode
%       clusterID assignments into cluster names:       
%       ClusterNames.ClusterID
%               A unique integer for each clusterID assignment in the
%               transcriptome data analysis
%
%       ClusterNames.ClusterName
%               Short names/descriptions corresponding to each clusterID.
% 
% Optional input name/value pairs: 
% 'thresh_UMI'
%       Minimum number of UMI counts required to call an individual cell as
%       'positive' for a given Tracer clone
%       (default=1)
% 
% 'thresh_min_cells_per_clone
%       Discard any Tracer clones with fewer than this number of cells 
%       (default=5)
%
% 'thresh_max_cells_per_clone'
%       Discard any Tracer clones with greater than this number of cells 
%       (default=Inf)
%
% 'thresh_min_cells_per_hit'
%       Minimum number of positive cells required for a given TracerClone
%       to be considered as 'hitting' a given state
%       (default=2)
%
% 'nRandTrials' 
%       Number of data permutations used to generate mean and st. dev
%       expectations
%       (default=1000)
%
% 'heatmap_distance_metric'
%       Distance metric argument passed to 'linkage' for generating final
%       clustered heatmap
%       (default='correlation')
% 
% 'heatmap_linkage_method'
%       Method passed to 'linkage' for generating final clustered heatmap
%       (default='average')
%
% OUTPUTS:
% 'Couplings'
%       A structure array with the following fields
%       Couplings.matrix
%               Values from the correlation ('c') or zscore ('z') heatmap 
%               matrices, reordered based on the plotted dendrogram
%       
%       Couplings.state_names_ord
%               Cluster names, reordered based on the plotted dendrogram
% 
%       Couplings.settings
%               A list of all parameters / values
%

%% PARAMETER SETTINGS

% set defaults
def.thresh_UMI = 1;
def.thresh_min_cells_per_clone = 5; 
def.thresh_max_cells_per_clone = Inf; 
def.thresh_min_cells_per_hit = 2; 
def.nRandTrials = 1000; 
def.heatmap_distance_metric = 'correlation';
def.heatmap_linkage_method = 'average';

% create parser object
parserObj = inputParser;
parserObj.FunctionName = 'plot_tracer_couplings';
parserObj.StructExpand = false; 
parserObj.addOptional('thresh_UMI',def.thresh_UMI);
parserObj.addOptional('thresh_min_cells_per_clone',def.thresh_min_cells_per_clone);
parserObj.addOptional('thresh_max_cells_per_clone',def.thresh_max_cells_per_clone);
parserObj.addOptional('thresh_min_cells_per_hit',def.thresh_min_cells_per_hit);
parserObj.addOptional('nRandTrials',def.nRandTrials);
parserObj.addOptional('heatmap_distance_metric',def.heatmap_distance_metric);
parserObj.addOptional('heatmap_linkage_method',def.heatmap_linkage_method);

% parse input options
parserObj.parse(varargin{:});
settings = parserObj.Results;

%% Load data

nDataSets = length(DataSet);

% build the first cell-x-clone matrix
for f = 1
    X = DataSet(f).merged.tracer_counts_matrix;
    cell_stateID = DataSet(f).merged.clustID';
end

% append additional cell-x-clone matrices as block diagonal
for f = 2:nDataSets   
    
    % build the next cell-x-clone matrix    
    X_next = DataSet(f).merged.tracer_counts_matrix;
    cell_stateID_next = DataSet(f).merged.clustID';

    % append 
    X = blkdiag(X, X_next);
    cell_stateID = [cell_stateID cell_stateID_next];   
end

%% Filter cells and clones

% convert counts matrix to boolean based on UMI threshold
X = +(X >= settings.thresh_UMI);

% remove cells that lack a clusterID state assignment
% remove clones that now have zero cells
clone_filt = sum(X,2) < settings.thresh_min_cells_per_clone | sum(X,2) > settings.thresh_max_cells_per_clone;
cell_filt = sum(X,1)==0;
cell_filt_tally = sum(cell_filt);
clone_filt_tally = sum(clone_filt);
while sum(cell_filt) > 0 || sum(clone_filt) > 0
    % perform one round of filtering
    X(clone_filt,:)=[];
    X(:,cell_filt)=[]; 
    cell_stateID(cell_filt)=[];
    % recalculate empty rows and columns 
    clone_filt = sum(X,2) < settings.thresh_min_cells_per_clone | sum(X,2) > settings.thresh_max_cells_per_clone;
    cell_filt = sum(X,1)==0;
    cell_filt_tally = cell_filt_tally + sum(cell_filt);
    clone_filt_tally = clone_filt_tally+ sum(clone_filt);
end
if cell_filt_tally > 0 || clone_filt_tally > 0
    disp(['Filtering ' num2str(cell_filt_tally) ' cells and ' num2str(clone_filt_tally) ' clones based on set thresholds'])
end

%% Filter out states with zero hits

% >1 values of the parameter: thresh_min_cells_per_hit will cause some
% states to no longer meet the 'hit' threshold. 

% get a list of all the states in the tracer dataset
unique_state_IDs = unique(cell_stateID);
nStates = length(unique_state_IDs);

% determine the number of 'hits' for each state
for j = 1:nStates
    % index the cells assigned to this particular j state
    cells_in_state_j = (unique_state_IDs(j) == cell_stateID);
    clone_hits_in_state_j(j) = sum((sum(X(:,cells_in_state_j),2) >= settings.thresh_min_cells_per_hit));
end

% omit states with zero hits from the analysis
nStates_to_remove = sum(clone_hits_in_state_j==0);
if nStates_to_remove > 0
    unique_state_IDs(clone_hits_in_state_j==0) = [];
    nStates = length(unique_state_IDs);
    disp(['Omitting ' num2str(nStates_to_remove) ' states with zero hits'])
end

%% Calculate state-state couplings 

disp(['Analyzing couplings across ' num2str(size(X, 2)) ' cells, ', num2str(size(X, 1)) ' clones, and ', num2str(nStates) ' state pairs'])

couplings = zeros(nStates,nStates);

initiate_timer = tic;
for j = 1:nStates
    
    % index the cells assigned to this particular j state
    cells_in_state_j = (unique_state_IDs(j) == cell_stateID);
    clone_hits_in_state_j = (sum(X(:,cells_in_state_j),2) >= settings.thresh_min_cells_per_hit);
    
    for k = j:nStates

        % index the cells assigned to this particular k state
        cells_in_state_k = (unique_state_IDs(k) == cell_stateID);
        clone_hits_in_state_k = sum(X(:,cells_in_state_k),2) >= settings.thresh_min_cells_per_hit;    

        % state couplings definition:
        % sum the number of times any clone hit both state j and state k
        couplings(j,k) = sum(clone_hits_in_state_j & clone_hits_in_state_k);

    end
end

% reconstruct a symmetric matrix
couplings = symmetrize_matrix(couplings);

disp('Done calculating state couplings')

%% Calculate state-state couplings on randomized data

cycle_calc_time = toc(initiate_timer);
total_time_needed = cycle_calc_time * settings.nRandTrials;
if total_time_needed < 60
    disp(['Estimating ' num2str(round(total_time_needed)) ' sec will be needed to perform ' num2str(settings.nRandTrials) ' random trials...'])
else
    disp(['Estimating ' num2str(round(total_time_needed/60, 1)) ' min will be needed to perform ' num2str(settings.nRandTrials) ' random trials...'])
end

couplings_rand = zeros(nStates,nStates,settings.nRandTrials);

for r = 1:settings.nRandTrials

    % randomize cell state assignments
    rand_ind = randperm(length(cell_stateID));
    cell_stateID_rand = cell_stateID(rand_ind); 

    for j = 1:nStates

        % index the cells assigned to this particular j state
        cells_in_state_j = (unique_state_IDs(j) == cell_stateID_rand);
        clone_hits_in_state_j = (sum(X(:,cells_in_state_j),2) >= settings.thresh_min_cells_per_hit);

        for k = j:nStates 

                % index the cells assigned to this particular k state
                cells_in_state_k = (unique_state_IDs(k) == cell_stateID_rand);
                clone_hits_in_state_k = sum(X(:,cells_in_state_k),2) >= settings.thresh_min_cells_per_hit;    

                % state couplings definition:
                % sum the number of times any clone hit both state j and state k
                couplings_rand(j,k,r) = sum(clone_hits_in_state_j & clone_hits_in_state_k);

        end 
    end

    % reconstruct a symmetric matrix
    couplings_rand(:,:,r) = symmetrize_matrix(couplings_rand(:,:,r));

end

disp('Done calculating randomized couplings')

%% Calculate z-scores and correlations

% calculate coupling z-scores
couplings_rand_mean = mean(couplings_rand,3);
couplings_rand_stdev = std(couplings_rand,0,3);
couplings_zscore = (couplings - couplings_rand_mean) ./ couplings_rand_stdev;
couplings_zscore(isnan(couplings_zscore)) = 0;
couplings_zscore(isinf(couplings_zscore)) = 0;

% calculate correlations
couplings_correlation = corrcoef(couplings_zscore, 'Rows', 'pairwise');
couplings_correlation(isnan(couplings_correlation)) = 0;

%% Plot correlations heatmap

matrix_for_heatmap = couplings_correlation;

f = figure('position',[100 100 800 800]);
axis off;

% draw left dendrogram
ax2 = axes('position', [0.02, 0.05, 0.07, 0.8]);
Z = linkage(matrix_for_heatmap, settings.heatmap_linkage_method, settings.heatmap_distance_metric);
[H, ~, ~] =  dendrogram(Z,size(matrix_for_heatmap,1),'Orientation', 'left');
set(H, 'Color',[0.2 0.2 0.2]);
set(ax2,'YDir','reverse'),
axis off;

% draw top dendrogram
ax3= axes('position', [0.1 0.855 0.6 0.06]); 
[H, ~, outperm_c] =  dendrogram(Z,size(matrix_for_heatmap,1));
set(H, 'Color',[0.2 0.2 0.2]);
axis off;
title('State Coupling Correlations')

% import and reorder cluster/state names 
[~, locs] = ismember(unique_state_IDs(outperm_c),ClusterNames.ClusterID);
reordered_state_names_c = ClusterNames.ClusterName(locs);

% plot heatmap
ax1 = axes('position', [0.1 0.05 0.6 0.8], 'YAxisLocation', 'right');
nSharedClones_reorder = imresize(matrix_for_heatmap(outperm_c,outperm_c),1,'nearest');
imagesc(nSharedClones_reorder);
set(gca,'yaxislocation','right');
set(gca,'XTick',1:length(outperm_c),'XTickLabel', [])
set(gca,'YTick',1:length(outperm_c),'YTickLabel', reordered_state_names_c,'TickLength',[0 0])
box off

% adjust colorbar
cb = colorbar('eastoutside', 'position',[0.92, 0.92, 0.02, 0.06]);
cb.Label.HorizontalAlignment = 'Right';

% link axes for zooming
linkaxes([ax1, ax3],'x');
linkaxes([ax1, ax2],'y');

try
    addpath('othercolor')
    colormap((othercolor('RdYlBu10',1000)));
    f.Children(2).CLim = [-1 1];
end

%% Plot z-score heatmap

matrix_for_heatmap = couplings_zscore;

f = figure('position',[100 100 800 800]);
axis off;

% draw left dendrogram
ax2 = axes('position', [0.02, 0.05, 0.07, 0.8]);
Z = linkage(matrix_for_heatmap, settings.heatmap_linkage_method, settings.heatmap_distance_metric);
[H, ~, ~] =  dendrogram(Z,size(matrix_for_heatmap,1),'Orientation', 'left');
set(H, 'Color',[0.2 0.2 0.2]);
set(ax2,'YDir','reverse'),
axis off;

% draw top dendrogram
ax3= axes('position', [0.1 0.855 0.6 0.06]); 
[H, ~, outperm_z] =  dendrogram(Z,size(matrix_for_heatmap,1));
set(H, 'Color',[0.2 0.2 0.2]);
axis off;
title('State Coupling Z-Scores')

% import and reorder cluster/state names 
[~, locs] = ismember(unique_state_IDs(outperm_z),ClusterNames.ClusterID);
reordered_state_names_z = ClusterNames.ClusterName(locs);

% plot heatmap
ax1 = axes('position', [0.1 0.05 0.6 0.8], 'YAxisLocation', 'right');
nSharedClones_reorder = imresize(matrix_for_heatmap(outperm_z,outperm_z),1,'nearest');
imagesc(nSharedClones_reorder);
set(gca,'yaxislocation','right');
set(gca,'XTick',1:length(outperm_z),'XTickLabel', [])
set(gca,'YTick',1:length(outperm_z),'YTickLabel', reordered_state_names_z,'TickLength',[0 0])
box off

% adjust colorbar
cb = colorbar('eastoutside', 'position',[0.92, 0.92, 0.02, 0.06]);
cb.Label.HorizontalAlignment = 'Right';

% link axes for zooming
linkaxes([ax1, ax3],'x');
linkaxes([ax1, ax2],'y');

% adjust colormap
try
    addpath('othercolor')
    colormap((othercolor('RdYlBu10',1000)));
    f.Children(2).CLim = [-5 5];
end

%% Outputs:

% export couplings matrices
Couplings.matrix_c = couplings_correlation(outperm_c,outperm_c);
Couplings.matrix_z = couplings_zscore(outperm_z,outperm_z);

% etc
Couplings.state_names_ord_c = reordered_state_names_c';
Couplings.state_names_ord_z = reordered_state_names_z';
Couplings.settings = settings;

