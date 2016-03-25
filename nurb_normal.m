function[normal,norm_dir]=nurb_normal(jac,M,N)
for j=1:N
    for i=1:M
        normal(:,i,j)=cross(jac{1}(:,i,j),jac{2}(:,i,j));
    end
end
%     quiver3(p(1,1:99,1:99),p(2,1:99,1:99),p(3,1:99,1:99),norm(1,1:99,1:99),norm(2,1:99,1:99),norm(3,1:99,1:99));        
%% direction of normal   
  for i=1:N
      for j=1:M
          for k=1:3
         norm_dir(k,j,i)=normal(k,j,i)/sqrt(normal(1,j,i)^2+normal(2,j,i)^2+normal(3,j,i)^2);
          end
      end
  end
end