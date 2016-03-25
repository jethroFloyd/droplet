clear all;
close all;
clc
addpath('C:\Program Files\MATLAB\R2011a\toolbox\nurbs_toolbox');
radius_major=5;
 
fullangle_theta1=0;
fullangle_theta2=(pi/3);
fullangle_phi=2*pi;
M=100;
N=100;
degree=4;
theta=linspace(0,(pi/3),N);
phi=linspace(0,2*pi,M);
%% defining the weights of the nurb %%
for i=1:N
    for j=1:M
    w(j,i)=1;
    end
end
%w(:,N)=1.5;
for i=1:N
    
        for k=1:M
            pp(1,k,i)=radius_major*sin(theta(i))*cos(phi(k))*w(k,i);
        end
        for k=1:N
            pp(2,k,i)=radius_major*sin(theta(i))*sin(phi(k))*w(k,i);
        end
        for k=1:N
            pp(3,k,i)=radius_major*cos(theta(i))*w(k,i);
        end
        for k=1:N
            pp(4,k,i)=w(k,i);
        end
end
degree=4;
knots=knotspan(M,N,degree);
%theta_norm=theta/max(theta);
%phi_norm=phi/max(phi);
theta_norm=linspace(0,1,N);
phi_norm=linspace(0,1,M);
srf=nrbmak(pp,knots);
[p,basis]=nrbeval(srf,{theta_norm,phi_norm});
plot3(p(1,:),p(2,:),p(3,:))

[dsurf,d2surf]=nrbderiv_modified(srf);
[pnt, jac, hess] = nrbdeval_modified(srf, dsurf, d2surf, {theta_norm, phi_norm});
[normal,norm_dir]=nurb_normal(jac,M,N);

[radius]=gaussian_nurb_curvature(jac,hess,normal,M,N);


 for k=1:21
      for j=1:21
          for i=1:3
          
              pp(i,j,k)=(pp(i,j,k)-norm_dir(i,j,k)*.2);
              
          end
          
      end
 end
 srf=nrbmak(pp,knots);
 
%  w(:)=0.00000000001;
%  %for k=1:21
%      %for j=1:21
%          w(21,21)=10;
%      %end
%  %end
%  for i=1:N
%     
%         for k=1:M
%             pp(1,k,i)=p(1,k,i)*w(k,i);
%         end
%         for k=1:M
%             pp(2,k,i)=p(2,k,i)*w(k,i);
%         end
%         for k=1:M
%             pp(3,k,i)=p(3,k,i)*w(k,i);
%         end
%         for k=1:M
%             pp(4,k,i)=w(k,i);
%         end
% end
%  
%  srf=nrbmak(pp,knots)
%  p=nrbeval(srf,{theta_norm,phi_norm});
%  
% 
%  