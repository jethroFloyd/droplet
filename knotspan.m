function[knots]=knotspan(M,N,degree)
for j=1:degree+1
    knotsu(j)=0;
end

for j=degree+2:(N)
    knotsu(j)=knotsu(j-1)+1/(N-degree);
end
for j=(N+1):(N+degree+1)
    knotsu(j)=1;
end

for j=1:degree+1
    knotsv(j)=0;
end

for j=degree+2:(M)
    knotsv(j)=knotsv(j-1)+1/(M-degree);
end
for j=(M+1):(M+degree+1)
    knotsv(j)=1;
end
knots={knotsu, knotsv};
end