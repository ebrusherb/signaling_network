function L = dividecontours(C)
L=cell(1);
q=size(C,2);
location=1;
i=1;
while location<=q
    go=C(2,location);
    L{i}=C(:,location+(1:go));
    location=location+1+go;
    i=i+1;
end
    