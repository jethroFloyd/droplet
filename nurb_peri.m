function[surface_area]=nurb_peri(p,M,N)
for i=1:M
    peri=0;
for j=1:N-1
    peri=peri+sqrt((p(1,j+1,i)-p(1,j,i))^2+(p(3,j+1,i)-p(3,j,i))^2);
end
slice(i)=peri;
end
surface_area=0;
for i=1:M-1
   surface_area=surface_area+abs((p(2,1,i+1)- p(2,1,i)))*slice(i);
end
end