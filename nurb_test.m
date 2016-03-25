clear all;
close all;
addpath('C:\Program Files\MATLAB\R2011a\toolbox\nurbs_toolbox')

R=5;
contact_angle=120;
base_radius=sind(contact_angle)*R;
fullangle=360;
N=100;
M=100;
r(1)=0;
for i=2:N
    r(i)=r(i-1)+base_radius/(N-1);
end
theta(1)=0;
for i=2:M
    theta(i)=theta(i-1)+fullangle/(M-1);
end
for i=1:M
    
        for k=1:N
            pp(1,k,i)=r(k)*cosd(theta(i));
        end
        for k=1:N
            pp(2,k,i)=r(k)*sind(theta(i));
        end
        for k=1:N
            pp(3,k,i)=sqrt(R^2-r(k)^2)-R*cos(contact_angle);
        end
end
degree=4;
for j=1:degree
    knotsu(j)=0;
end
for j=degree+1:(N+degree-5)
    knotsu(j)=knotsu(j-1)+1/N;
end
for j=(N+degree-4):(degree+N-1)
    knotsu(j)=1;
end

for j=1:degree
    knotsv(j)=0;
end
for j=degree+1:(M+degree-5)
    knotsv(j)=knotsv(j-1)+1/M;
end
for j=(M+degree-4):(degree+M-1)
    knotsv(j)=1;
end


knots={knotsu, knotsv};
srf=nrbmak(pp,knots);
nrbplot(srf,[100 100])
hold off


dsurf=nrbderiv(srf);
d2surf_u=nrbderiv(dsurf{1});
d2surf_v=nrbderiv(dsurf{2});
r_norm=r./max(r);
theta_norm=theta./max(theta);

% for i=1:N+1
%     for j=1:N+1
%     pnt_x(i,j)=r(j)*cosd(theta(i));
%     end
%     
% end
% for i=1:N+1
%     for j=1:N+1
%     pnt_y(i,j)=r(j)*sind(theta(i));
%     end
%     
% end


%  p = nrbeval(srf,{linspace(-5.0,5.0,20) linspace(-5.0,5.0,20)});
%  h = surf(squeeze(p(1,:,:)),squeeze(p(2,:,:)),squeeze(p(3,:,:)));
%  set(h,'FaceColor','blue','EdgeColor','blue');
%  title('First derivatives over a test surface.');
% 
%  npts = 5;
%  tt = linspace(0.0,1.0,npts);
%  dsrf = nrbderiv(srf);
% tt={pnt_x,pnt_y};
 [p1, dp] = nrbdeval(srf, dsurf, {r_norm theta_norm});
 [p2,dp_u]=nrbdeval(dsurf{1},d2surf_u,{r_norm,theta_norm});
 [p3,dp_v]=nrbdeval(dsurf{2},d2surf_v,{r_norm,theta_norm});

 for i=1:N-1
    
         for k=1:M-1
             zeta1=dp{1}(:,1:99,i);
             zeta2=dp{2}(:,1:99,i);
             zeta_norm(:,k)=cross(zeta1(:,k),zeta2(:,k));
         end
          norm(:,:,i)=zeta_norm;
 end
 
  for i=1:N-1
      for j=1:M-1
          for k=1:3
         norm_dir(k,j,i)=norm(k,j,i)/sqrt(norm(1,j,i)^2+norm(2,j,i)^2+norm(3,j,i)^2);
          end
      end
  end
  
  %%  after perturbation %%
  for i=1:2:21
      for j=1:2:21
          for k=1:3
          pp(k,j,i)=pp(k,j,i)-norm_dir(k,j,i)*1.5;
          end
      end
  end
  srf=nrbmak(pp,knots);
  nrbplot(srf,[100 100])
  hold on
  dsurf=nrbderiv(srf);
d2surf_u=nrbderiv(dsurf{1});
d2surf_v=nrbderiv(dsurf{2});
[p1, dp] = nrbdeval(srf, dsurf, {r_norm theta_norm});
 [p2,dp_u]=nrbdeval(dsurf{1},d2surf_u,{r_norm,theta_norm});
 [p3,dp_v]=nrbdeval(dsurf{2},d2surf_v,{r_norm,theta_norm});

 for i=1:N-1
    
         for k=1:M-1
             zeta1=dp{1}(:,1:99,i);
             zeta2=dp{2}(:,1:99,i);
             zeta_norm(:,k)=cross(zeta1(:,k),zeta2(:,k));
         end
          norm(:,:,i)=zeta_norm;
 end
  
   quiver3(p1(1,1:99,1:99),p1(2,1:99,1:99),p1(3,1:99,1:99),norm(1,1:99,1:99),norm(2,1:99,1:99),norm(3,1:99,1:99));
  %quiver3(p1(1,:,:,1:99),p1(2,:,:,1:99),p1(3,:,:,1:99),norm(1,:,:,:),norm(2,:,:,:),norm(3,:,:,:))
%  up2 = vecnorm(dp{1});
%  vp2 = vecnorm(dp{2});
%  
%  hold on;
%  plot3(p1(1,:),p1(2,:),p1(3,:),'ro');
%  h1 = quiver3(p1(1,:),p1(2,:),p1(3,:),up2(1,:),up2(2,:),up2(3,:));
%  h2 = quiver3(p1(1,:),p1(2,:),p1(3,:),vp2(1,:),vp2(2,:),vp2(3,:));
%  set(h1,'Color','black');
%  set(h2,'Color','black');
% 
%
% nrbplot(dsurf,[100 100])