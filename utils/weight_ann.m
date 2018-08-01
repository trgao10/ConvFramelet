function y=weight_ann(data)
% the function to compute the weight matrix using kdtree

[m,n]=size(data);

num_s=50; % number of neighbors

kdtree = vl_kdtreebuild(data);
[idx, dist] = vl_kdtreequery(kdtree, data, data, 'NumNeighbors', num_s, 'MaxComparisons', min(n,2^10));

id_row=repmat([1:n],num_s,1);
id_col=double(idx);

sigma=sparse([1:n],[1:n],1./max(dist(21,:),1e-2),n,n);
w=exp(-(dist*sigma).^2); %% old

%%% START --- try new and symmetric weights
% sigma_row = 1./max(dist(21,:),1e-2);
% sigma_symmetrize = sigma_row(idx);
% w=exp(-(sqrt(sigma_symmetrize).*dist*sqrt(sigma)).^2);
%%% END --- if not working, get back to the original weights

y=sparse(id_row,id_col,w,n,n);

%%% START --- try extra symmetrization
% y = max(y, y');
% y = min(y, y');
% y = (y+y')/2;
%%% END --- if not working revert by commenting the line above

end    