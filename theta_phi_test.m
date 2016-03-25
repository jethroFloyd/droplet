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
theta(1)=0;
degree=4;

% for i=2:N
%     theta(i)=theta(i-1)+(fullangle_theta2-fullangle_theta1)/(N);
% end
% phi(1)=0;
% for i=2:M
%     phi(i)=phi(i-1)+fullangle_phi/(M);
% end
theta=linspace(0,(pi/3),N);
phi=linspace(0,2*pi,M);
% theta=0:(1/N-degree-2):pi/2;
% phi=0:(1/M-degree-2):2*pi;
% N=length(theta);
% M=length(phi);
%% building nurb data structure %%
for i=1:M
    
        for k=1:N
            pp(1,k,i)=radius_major*sin(theta(i))*cos(phi(k));
        end
        for k=1:N
            pp(2,k,i)=radius_major*sin(theta(i))*sin(phi(k));
        end
        for k=1:N
            pp(3,k,i)=radius_major*cos(theta(i));
        end
end
degree=4;
for j=1:degree+1
    knotsu(j)=0;
end

for j=degree+2:(N+degree-6)
    knotsu(j)=knotsu(j-1)+1/(N-degree-2);
end
for j=(N+degree-5):(N+degree-1)
    knotsu(j)=1;
end

for j=1:degree+1
    knotsv(j)=0;
end

for j=degree+2:(M+degree-6)
    knotsv(j)=knotsv(j-1)+1/(M-degree-2);
end
for j=(M+degree-5):(M+degree-1)
    knotsv(j)=1;
end
knots={knotsu, knotsv};
theta_norm=theta/max(theta);
phi_norm=phi/max(phi);




srf=nrbmak(pp,knots);
 nrbplot(srf,{theta_norm,phi_norm});
% hold on
%% evaluating nurb data structure 
[p,basis]=nrbeval(srf,{theta_norm,phi_norm});

for i=1:M
    Z(i)=p(3,1,i);
end
minimum_z=find(Z==min(Z));
vol=nurb_volume(p,M,N);

%% evaluating derivatives 
[dsurf,d2surf]=nrbderiv_modified(srf);
[pnt, jac, hess] = nrbdeval_modified(srf, dsurf, d2surf, {theta_norm, phi_norm});
[normal,norm_dir]=nurb_normal(jac,M,N);
[radius]=gaussian_nurb_curvature(jac,hess,normal,M,N);
%  plot3(p(1,:),p(2,:),p(3,:))
% 
% hold on
%   quiver3(p(1,:,:),p(2,:,:),p(3,:,:),normal(1,:,:),normal(2,:,:),normal(3,:,:))
%% evaluating normals



% for j=1:N
%     for i=1:N
%         norm(:,i,j)=cross(jac{1}(:,i,j),jac{2}(:,i,j));
%     end
% end
% %     quiver3(p(1,1:99,1:99),p(2,1:99,1:99),p(3,1:99,1:99),norm(1,1:99,1:99),norm(2,1:99,1:99),norm(3,1:99,1:99));        
% %% direction of normal   
%   for i=1:N
%       for j=1:M
%           for k=1:3
%          norm_dir(k,j,i)=norm(k,j,i)/sqrt(norm(1,j,i)^2+norm(2,j,i)^2+norm(3,j,i)^2);
%           end
%       end
%   end
  %% finding the ground value of Z %%
% for i=1:M
%     Z(i)=p(3,1,i);
% end
% min=find(Z==min(Z));
  
  %% test perturbation %%
  %for k=1:21
  for k=45:48
      for j=1:51
          for i=1:3
              
          if k~=minimum_z
%           for j=1:M
%               for k=1:3
              p(i,j,k)=p(i,j,k)+norm_dir(i,j,k)*.2;
          end
          end
      end
  end
  %%
srf=nrbmak(p,knots);
p=nrbeval(srf,{theta_norm,phi_norm});

[dsurf,d2surf]=nrbderiv_modified(srf);
[pnt, jac, hess] = nrbdeval_modified(srf, dsurf, d2surf, {theta_norm, phi_norm});
[normal,norm_dir]=nurb_normal(jac,M,N);

[radius]=gaussian_nurb_curvature(jac,hess,normal,M,N);
vol=nurb_volume(p,M,N);
 plot3(p(1,:),p(2,:),p(3,:))

hold on
  quiver3(p(1,:,:),p(2,:,:),p(3,:,:),normal(1,:,:),normal(2,:,:),normal(3,:,:))
  daspect([1 1 1])
 
 
%for loops=1:20
    %% calculating minimum radius value
    [m_max,n_max,m_min,n_min]=opt_radius(radius,M);
    %% push-pop algo 
    for loop=1:20
        for k=1:3
        p(k,m_min,n_min)=p(k,m_min,n_min)+norm_dir(k,m_min,n_min)*.1;
        end
   %      vol2=nurb_volume(p,M,N);
        ashesh=1;
        kk=.01
        count=1;
        while(ashesh==1||(vol_new-vol)>.1)
            ashesh=2;
        for k=1:3
        p(k,m_max,n_max)=p(k,m_max,n_max)+norm_dir(k,m_max,n_max)*kk;
        end
        vol_new=nurb_volume(p,M,N)
        if((vol_new-vol)>.1)
           for k=1:3 
              p(k,m_max,n_max)=p(k,m_max,n_max)-norm_dir(k,m_max,n_max)*kk;
           end
        end
        kk=kk+(vol_new-vol)/10;
         figure(1)
         plot3(p(1,:),p(2,:),p(3,:))
view(120,30)
         daspect([1 1 1])
        hold off
        count=count+1;
        end
        srf=nrbmak(p,knots);
        p=nrbeval(srf,{theta_norm,phi_norm});
       [dsurf,d2surf]=nrbderiv_modified(srf);
       [pnt, jac, hess] = nrbdeval_modified(srf, dsurf, d2surf, {theta_norm, phi_norm});
       [normal,norm_dir]=nurb_normal(jac,M,N);
       radius=nurb_curvature(jac,hess,normal,M,N);
        [m_max,n_max,m_min,n_min]=opt_radius(radius,M);
%         figure(1)
%         
%         daspect([1 1 1])
%         hold off
        
    end
%     
% plot3(p(1,:),p(2,:),p(3,:))
        
    
      

%plot3(p(1,:),p(2,:),p(3,:))
%hold on
%quiver3(p(1,1:99,1:99),p(2,1:99,1:99),p(3,1:99,1:99),norm(1,1:99,1:99),norm(2,1:99,1:99),norm(3,1:99,1:99));        
% %% direction of normal   
  %%
  
%% radius of curvature %%
% for j=1:M
%     for i=1:N
%         A(i,j)=dot(cross(jac{1}(:,i,j),jac{2}(:,i,j)),hess{1,1}(:,i,j));
%         B(i,j)=dot(cross(jac{1}(:,i,j),jac{2}(:,i,j)),hess{1,2}(:,i,j));
%         C(i,j)=dot(cross(jac{1}(:,i,j),jac{2}(:,i,j)),hess{2,2}(:,i,j));
%     end
% end
%  for j=1:N
%       for i=1:M
%           
%              D(i,j)= sqrt(norm(1,i,j)^2+norm(2,i,j)^2+norm(3,i,j)^2);
%       end
%  end
% for i=1:M
%     for j=1:N
%         R(i,j)=(A(i,j)*C(i,j)-B(i,j)^2)/D(i,j)^4;
%     end
% end
%% finding the minimum of R %%
% for i=1:M
%     Z(i)=p(3,1,i);
% end
% min=find(Z==min(Z));


% %tri=delaunay(xdata',ydata',zdata');
% %  xsort=sort(xdata);
% %  for i=1:length(xsort)
% %      index=find(xdata==xsort(i));
% %      ysort(i)=ydata(index);
% %      zsort(i)=zdata(index);
% % end
%  % V = dblquad(@(xsort,ysort) zsort,min(xsort),max(xsort),min(ysort),max(ysort))
% %  col(:,1)=xsort;
% %  col(:,2)=ysort;
% %  col(:,3)=zsort;
% 
% % F2 = TriScatteredInterp(xdata',ydata',zdata');
% % %q1 = quad2d(@(xdata,ydata) abs(F2(xdata,ydata)),min(xdata),max(xdata),min(ydata),max(ydata),'AbsTol',0.01)
% % Z= zdata- min(zdata);
% % X=xdata;
% % Y=ydata;
% % TRI= DelaunayTri(X',Y',Z');
% % vol = 0;
% % area = 0;
% % nTri = size(TRI,1); % Number of triangles
% % if (nTri > 3) % Need at least 4 triangles to form a closed surface
% % for i=1:nTri
% % U = [X(TRI(i,1)) Y(TRI(i,1)) Z(TRI(i,1))]; % First point of triangle
% % V = [X(TRI(i,2)) Y(TRI(i,2)) Z(TRI(i,2))]; % Second point of triangle
% % W = [X(TRI(i,3)) Y(TRI(i,3)) Z(TRI(i,3))]; % Third point of triangle
% % A = V - U; % One side of triangle (from U to V)
% % B = W - U; % Another side of triangle (from U to W)
% % C = cross(A, B); % Length of C equals to the area of the parallelogram [A,B]
% % normC = norm(C);
% % a = 0.5 * normC; % Area of triangle [U,V,W]
% % P = (Z(TRI(i,1)) + Z(TRI(i,2)) + Z(TRI(i,3)) ) / 3; % Middle of triangle
% % vol = vol + P *a;
% % area = area + a;
% % end
% % end
% vol = abs(vol);



%     for j=1:N
%      area_arr(j)=polyarea(p(1,:,j),p(2,:,j));
%     end   
%     
% vol=0;
% for i=1:M-1
%    vol=vol+abs((p(3,1,i+1)- p(3,1,i)))*area_arr(i);
% end
