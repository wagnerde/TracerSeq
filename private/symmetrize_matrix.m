function  B = symmetrize_matrix(A)
%% Usage: B = symmetrize_matrix(A)
%
% Symmetrizes a trianglular matrix 
%
%% CODE:

[n,~] = size(A);
B = A' + A;
B(1:n+1:end) = diag(A);

