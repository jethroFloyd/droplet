function[R]=nurb_curvature(jac,hess,norm,M,N)
for j=1:M
    for i=1:N
        A(i,j)=dot(cross(jac{1}(:,i,j),jac{2}(:,i,j)),hess{1,1}(:,i,j));
        B(i,j)=dot(cross(jac{1}(:,i,j),jac{2}(:,i,j)),hess{1,2}(:,i,j));
        C(i,j)=dot(cross(jac{1}(:,i,j),jac{2}(:,i,j)),hess{2,2}(:,i,j));
    end
end
 for j=1:N
      for i=1:M
          
             D(i,j)= sqrt(norm(1,i,j)^2+norm(2,i,j)^2+norm(3,i,j)^2);
      end
 end
for i=1:M
    for j=1:N
        R(i,j)=(A(i,j)*C(i,j)-B(i,j)^2)/D(i,j)^4;
    end
end
end