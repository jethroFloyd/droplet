function[vol]=nurb_volume(p,M,N)
 for j=1:N
     area_arr(j)=trapz(p(1,:,j),p(3,:,j));
    end   
    
vol=0;
for i=1:M-1
   vol=vol+abs((p(2,1,i+1)- p(2,1,i)))*area_arr(i);
end
end


