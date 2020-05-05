function DataSet = merge_tracer_data(DataSet)
%Usage:  DataSet = merge_tracer_data(DataSet)
%
% This function analyzes inDrops cell barcode lists to identify individual 
% cells that were captured in both the Tracer and transcriptome datasets. 
% The DataSet structure is then updated with a new set of 'merged' fields 
% in which clusterIDs derived from the transcriptome analysis are indexed
% to the tracer counts matrix.
% 
% Cell barcodes present in either the transcriptome or tracer datasets (but
% not both) are discarded.
% 
% After merging, some cells will have been removed. Any clones that now
% lack cells are also discarded.
%
% REQUIRED INPUTS:
% Inputs are supplied via 'DataSet', a structure array with multiple fields, 
% and one record for each biological sample
%       DataSet.name
%               Short sample name, e.g. 'TracerSeq1'
%
%       DataSet.imported.tracer_counts_matrix
%               Matrix of TracerSeq UMI counts. Each row corresponds to a
%               unique TracerSeq clone, each column corresponds to a single
%               cell
%
%       DataSet.imported.tracer_cell_barcodes 
%               List of inDrops cell barcodes, one barcode for each column
%               in the TracerSeq counts matrix
%
%       DataSet.imported.txome_clustID
%               A numerical array of cell-clusterID annotations from the
%               transcriptome analysis
%
%       DataSet.imported.txome_cell_barcodes 
%               List of inDrops cell barcodes, one barcode for each
%               clusterID assignment in the above array
%
% OUTPUTS: 
% The following fields will be outputted for each record in DataSet
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

%% CODE:

for j = 1:length(DataSet)
    
    disp(['Merging Tracer Dataset #' num2str(j)])

    % flag cell barcodes that are present in both the txome and tracer datasets
    [~, tmp.tracer_ind, tmp.txome_ind] = intersect(DataSet(j).imported.tracer_cell_barcodes, DataSet(j).imported.txome_cell_barcodes);
    disp(['Found n=' num2str(length(tmp.tracer_ind)) ' shared cell barcodes between the Tracer and transcriptome datasets'])
    
    % merge tracer and txome datasets
    DataSet(j).merged.tracer_counts_matrix = DataSet(j).imported.tracer_counts_matrix(:,tmp.tracer_ind);
    DataSet(j).merged.tracer_cell_barcodes = DataSet(j).imported.tracer_cell_barcodes(tmp.tracer_ind);
    DataSet(j).merged.clustID = DataSet(j).imported.txome_clustID(tmp.txome_ind);
    
    % filter cells & clones
    % remove cells that lack a clusterID state assignment
    % remove clones that now have zero cells after applying the cell filter
    tmp.clone_filt = sum(DataSet(j).merged.tracer_counts_matrix,2)==0;
    tmp.cell_filt = sum(DataSet(j).merged.tracer_counts_matrix,1)==0 | isnan(DataSet(j).merged.clustID)';    
    while sum(tmp.cell_filt) > 0 || sum(tmp.clone_filt) > 0
        % perform one round of filtering
        disp(['Filtering ' num2str(sum(tmp.cell_filt)) ' cells and ' num2str(sum(tmp.clone_filt)) ' clones'])
        DataSet(j).merged.tracer_counts_matrix(tmp.clone_filt,:)=[];
        DataSet(j).merged.tracer_counts_matrix(:,tmp.cell_filt)=[];
        DataSet(j).merged.tracer_cell_barcodes(tmp.cell_filt) = [];    
        DataSet(j).merged.clustID(tmp.cell_filt) = [];
        % recalculate empty rows and columns 
        tmp.clone_filt = sum(DataSet(j).merged.tracer_counts_matrix,2)==0;
        tmp.cell_filt = sum(DataSet(j).merged.tracer_counts_matrix,1)==0;
    end
    
    disp('Done')
end