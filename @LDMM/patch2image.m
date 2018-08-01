function y = patch2image(f,x1,x2,px,py,mi,ni)

fe=zeros(px+mi,py+ni);
w=zeros(px+mi,py+ni);


me=px+mi;
ne=py+ni;

for j1=1:px
    for j2=1:py
        f_tmp=zeros(me,ne);
        f_tmp((x2+j2-2)*me+x1+j1-1)=f(:,(j1-1)*py+j2);
        fe=fe+f_tmp;
        f_tmp=zeros(me,ne);
        f_tmp((x2+j2-2)*me+x1+j1-1)=1;
        w=w+f_tmp;
    end
end

w=w(1:mi,1:ni);

if(sum(sum(w<1/2))>0)
    fprintf('patches do not cover the image.\n')
    y=0;
    return;
end

y=fe(1:mi,1:ni)./w;

end