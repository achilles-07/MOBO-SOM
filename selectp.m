function p=selectp(rep,p1)
tempp=zeros(p1,1);
for i=1:p1
tempp(i)=rep(i).GridIndex;
end
[~,b]=sort(tempp);
if(rand<=0.5)
p=b(1);
else
p=b(end);
end
end

