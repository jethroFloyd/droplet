function[A]=basis_matrix(M,N,uu,degree,knots)
for i=2:(M-1) 
span = findspan(M-1,degree,uu(i),knots); 
NN = basisfun(span,uu(i),degree,knots); 
l = span-degree; 
for k=1:degree 
A(i,l+k) = NN(k); 
end 
end 
% Set first and last Matrix Value 
A(1,1) = 1.0; 
A(M,N) = 1.0; 
end