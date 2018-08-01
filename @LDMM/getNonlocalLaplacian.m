function L = getNonlocalLaplacian(data, graphLaplaicanOptions)
% getNonlocalLaplacian: construct nonlocal Laplacian from a data matrix
%     data:    the matrix of dimension N-by-D, where N is the number of
%              data points and D is the dimension of the data
%     graphLaplacianOptions: 
%          --- symLaplacian  { [false] | true }
%          --- NN            { [50] | (positive integer < N) }
%          --- SelfTuningNN  { [20] | (positive integer < NN )}
%
%
%   last modified by Tingran Gao (trgao10@math.duke.edu)
%   last modified on June 24, 2016

symLaplacian = LDMM.getoptions(graphLaplaicanOptions, 'symLaplacian', false);
NN = LDMM.getoptions(graphLaplaicanOptions, 'NN', 50);
SelfTuningNN = LDMM.getoptions(graphLaplaicanOptions, 'SelfTuningNN', 20);
MaxNumCompUpperBound = LDMM.getoptions(graphLaplaicanOptions, 'MaxNumCompUpperBound', 2^10);

[N, ~] = size(data);

kdtree = vl_kdtreebuild(data');
[idx, dist] = vl_kdtreequery(kdtree, data', data', 'NumNeighbors', NN, 'MaxComparisons', min(N,MaxNumCompUpperBound));

sigma = sparse(1:N, 1:N, 1./max(dist(SelfTuningNN+1,:), 1e-2), N, N);
% sigma = sparse(1:N, 1:N, 1./max(median(dist), 1e-2), N, N);

if ~symLaplacian
    w = exp(-(dist*sigma).^2);
else
%     sigma_row = 1./max(median(dist), 1e-2);
    sigma_row = 1./max(dist(SelfTuningNN+1,:), 1e-2);
    sigma_symmetrize = sigma_row(idx);
    w = exp(-(sqrt(sigma_symmetrize).*dist*sqrt(sigma)).^2);
end

L = sparse(repmat(1:N, NN, 1), double(idx), w, N, N);

end
