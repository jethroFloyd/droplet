clear all;
close all;
clc
addpath('C:\Program Files\MATLAB\R2011a\toolbox\nurbs_toolbox');
radius_major=5;
 
fullangle_theta1=0;
fullangle_theta2=(pi/2);
fullangle_phi=2*pi;
M=100;
N=100;
theta(1)=0;
degree=99;

% for i=2:N
%     theta(i)=theta(i-1)+(fullangle_theta2-fullangle_theta1)/(N);
% end
% phi(1)=0;
% for i=2:M
%     phi(i)=phi(i-1)+fullangle_phi/(M);
% end
theta=linspace(0,(pi/3),N);
phi=linspace(0,2*pi,M);

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

% for j=1:degree+1
%     knotsu(j)=0;
% end
% 
% for j=degree+2:(N-2)
%     knotsu(j)=knotsu(j-1)+1/(N-degree-2);
% end
% for j=(N-1):(N+degree-1)
%     knotsu(j)=1;
% end
% 
% for j=1:degree+1
%     knotsv(j)=0;
% end
% 
% for j=degree+2:(M-2)
%     knotsv(j)=knotsv(j-1)+1/(M-degree-2);
% end
% for j=(M-1):(M+degree-1)
%     knotsv(j)=1;
% end
% knots={knotsu, knotsv};
knots=knotspan(M,N,degree);
theta_norm=theta/max(theta);
phi_norm=phi/max(phi);




srf=nrbmak(pp,knots);

 nrbplot(srf,[M,N]);
 [p,basis]=nrbeval_octave(srf,{theta_norm,phi_norm});
 
% hold on
%% evaluating nurb data structure 
%[p,basis]=nrbeval(srf,{theta_norm,phi_norm});
%min_p=pp(:,:,N);
%p(:,:,N)=min_p;

% for i=1:M
%     Z(i)=p(3,1,i);
% end
% minimum_z=find(Z==min(Z));
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

  %% test perturbation %%
  %for k=1:21
  for k=45:55
      for j=30:35
          for i=1:3
              
          %if k~=minimum_z
%           for j=1:M
%               for k=1:3
              p(i,j,k)=p(i,j,k)+norm_dir(i,j,k)*.1;
          %end
          end
      end
  end
%%
srf=nrbmak(p,knots);
nrbplot(srf,[M,N])
[p,basis]=nrbeval_octave(srf,{theta_norm,phi_norm});


[dsurf,d2surf]=nrbderiv_modified(srf);
[pnt, jac, hess] = nrbdeval_modified(srf, dsurf, d2surf, {theta_norm, phi_norm});
[normal,norm_dir]=nurb_normal(jac,M,N);

[radius]=gaussian_nurb_curvature(jac,hess,normal,M,N);


vol=nurb_volume(p,M,N);
 


 
 
%for loops=1:20
    %% calculating minimum radius value
    [m_max,n_max,m_min,n_min]=opt_radius(radius,M);
    %% push-pop algo 
    for loop=1:100
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
        p(k,m_max,n_max)=p(k,m_max,n_max)-norm_dir(k,m_max,n_max)*kk;
        end
        srf=nrbmak(p,knots);
        [p1,base]=nrbeval_octave(srf,{theta_norm,phi_norm})
        vol_new=nurb_volume(p1,M,N);
        if((vol_new-vol)>.1)
           for k=1:3 
              p(k,m_max,n_max)=p(k,m_max,n_max)+norm_dir(k,m_max,n_max)*kk;
           end
        end
        kk=kk+(vol_new-vol)/100;
        srf=nrbmak(p,knots);
       h=figure(1)
        nrbplot(srf,[M,N])
        hold off
        view(180,30)
        %axis([0 1 0 1 0 1])
        hold off
        count=count+1;
        end
        area_array(loop)=nurb_peri(p,M,N);
        vol_arr(loop)=vol_new;
%         figure(1)
%         plot(area_array);
%         hold on
       srf=nrbmak(p,knots);
       [p,basis]=nrbeval_octave(srf,{theta_norm,phi_norm})
       [dsurf,d2surf]=nrbderiv_modified(srf);
       [pnt, jac, hess] = nrbdeval_modified(srf, dsurf, d2surf, {theta_norm,phi_norm});
       [normal,norm_dir]=nurb_normal(jac,M,N);
       radius=nurb_curvature(jac,hess,normal,M,N);
       [m_max,n_max,m_min,n_min]=opt_radius(radius,M);
%           F = getframe(h);

 end