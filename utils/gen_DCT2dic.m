function Dct = gen_DCT2dic(px,py)
% generate dct of size px-by-py

Dct = zeros(px*py, px*py);

c = zeros(px, py);
for i = 1:px*py
    c(i) = 1;
    Dct(:,i) = reshape(idct2(c),[],1);
    c(i) = 0;
end