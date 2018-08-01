function w = gen_weight(F, V)

sym = 0;

% generate weight for fixed dictionaries
w = sqrt(sum(abs(F*V).^2,1));
w = w/max(w);
% (optional) symmetrization, assume V is ordered as vectorized 2d matrix
if sym
    p = uint(sqrt(size(V,1)));
    w = reshape(w, [p, p]);
    w = (w+w')/2;
    w = w(:);
end