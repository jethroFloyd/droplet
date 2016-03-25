function[RR]=gaussian_nurb_curvature(jac,hess,normal,M,N)
for j=1:M
    for i=1:N
        A(i,j)=dot(cross(jac{1}(:,i,j),jac{2}(:,i,j)),hess{1,1}(:,i,j));
        B(i,j)=dot(cross(jac{1}(:,i,j),jac{2}(:,i,j)),hess{1,2}(:,i,j));
        C(i,j)=dot(cross(jac{1}(:,i,j),jac{2}(:,i,j)),hess{2,2}(:,i,j));
    end
end
 for j=1:N
      for i=1:M
          
             D(i,j)= sqrt(normal(1,i,j)^2+normal(2,i,j)^2+normal(3,i,j)^2);
      end
 end
% for i=1:M
%     for j=1:N
%         R(i,j)=(A(i,j)*C(i,j)-B(i,j)^2)/D(i,j)^4;
%     end
% end
R=-(A.*C-B.^2)./(D.^4);



%% mean curvature %%
for j=1:N
    for i=1:M
        magQ_v(i,j)=norm(jac{1,1}(:,i,j));
        magQ_u(i,j)=norm(hess{1,2}(:,i,j));
    end
end
for j=1:N
    for i=1:M
        mid(i,j)=2*B(i,j)*dot(jac{1,1}(:,i,j),jac{1,2}(:,i,j));
    end
end

for j=1:N
    for i=1:M
        RR(i,j)=(A(i,j)*magQ_v(i,j)^2-mid(i,j)+C(i,j)*magQ_u(i,j))/2*D(i,j)^3;
    end
end
%      