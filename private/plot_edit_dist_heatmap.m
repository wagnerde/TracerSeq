function  dist_matrix = plot_edit_dist_heatmap(sequences, plot_title)
%% Usage: dist_matrix = plot_edit_dist_heatmap(sequences, plot_title)
%
% This function generates a clustered heatmap of all pairwise 'Levenshtein' 
% distances, calculated by strdist (Eduard Polityko), between a set of 
% character strings.
%
% INPUTS:
% sequences 
% 		A cell array of strings (e.g. nucleotide sequences)
% plot_title
% 		Name for the plot
%
% OUTPUTS:
% dist_matrix
%		A matrix of all pairwise levenshtein distances
% 

%% CODE:

nSeq = length(sequences);
seqdist_matrix_upper = zeros(nSeq,nSeq);

% get edit distances just for the upper triangle of the matrix
for j = 1:nSeq
    for k = 1:nSeq
        if k >= j % upper triangle 
            seqdist_matrix_upper(j,k) = strdist(sequences{j}, sequences{k});
        end
    end
end

% reconstruct complete matrix
dist_matrix = seqdist_matrix_upper + seqdist_matrix_upper.';

% generate figure
fg = figure('position',[100 100 600 600]);
axis off;

% add plot title
title(plot_title ,'Interpreter', 'none');

% draw left dendrogram
ax2 = axes('position', [0.02, 0.05, 0.07, 0.8]);
Z = linkage(dist_matrix,'average','correlation');
[H, ~, outperm] =  dendrogram(Z,size(dist_matrix,1),'Orientation', 'left');
set(H, 'Color',[0.2 0.2 0.2]);
set(ax2,'YDir','reverse'),
axis off;

% draw top dendrogram
ax3 = axes('position', [0.1 0.855 0.8 0.06]); 
[H, ~, outperm] =  dendrogram(Z,size(dist_matrix,1));
set(H, 'Color',[0.2 0.2 0.2]);
axis off;

% draw heatmap
ax1 = axes('position', [0.1 0.05 0.8 0.8]);
seqdist_matrix_reordered = imresize(dist_matrix(outperm,outperm),1,'nearest');
h = imagesc(seqdist_matrix_reordered);
axis square
axis off
colormap(flipud(parula))
%colormap(othercolor('GrMg_16'))
%set(gca,'clim', [0 10])

% add color bar
cb = colorbar('eastoutside', 'position',[0.92, 0.92, 0.02, 0.06]);
cb.Label.HorizontalAlignment = 'Right'; 
