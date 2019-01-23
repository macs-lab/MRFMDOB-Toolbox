function Mr = Mr_prd(Apara,n,L,r)
    Mr = zeros(2*n*L);
    for i=1:2*n*(L-1)
        for j=1:length(Apara)
            Mr(j+i-1,i) = Apara(j);
        end
    end
    for i=2*n*(L-1)+1:2*n*L
        Mr(r+L*(i-(2*n*(L-1)+1)),i)=1;
    end
end