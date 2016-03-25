clear all
close all
clc
addpath('C:\Program Files\MATLAB\R2011a\toolbox\nurbs_toolbox');
r=1;
M=100;
N=100;
u=linspace(-1,1,M);
v=linspace(-1,1,N);
 
degree=3;
for i=1:N
    for j=1:M
        
            pp(3,j,i)=(1-(u(i))^2)*(1-(v(j))^2)/((1+(u(i))^2)*(1+(v(j))^2));
            pp(1,j,i)=2*v(j)*(1-(u(i))^2)/((1+(u(i))^2)*(1+(v(j))^2));
            pp(2,j,i)=2*u(i)*(1+(v(j))^2)/((1+(u(i))^2)*(1+(v(j))^2));
    end
end

%% dent %%

  for k=1:N
      
      for j=1:M
       r=sqrt(pp(1,j,k)^2+pp(2,j,k)^2);   

          pp(3,j,k)=pp(3,j,k)*(1/(1+exp((-30*r))));
          
      end
  end
%   for k=35:45
%       for j=20:35
%           for i=1:3
%           p(i,j,k)=p(i,j,k)-norm_dir(i,j,k)*.1;
%           
%           end
%       end
%   end


% for i=1:3
%     p(i,31,31)= p(i,31,31)+norm_dir(i,31,31)*.1;
% end
% for i=1:3
%     p(i,41,41)= p(i,41,41)-norm_dir(i,41,41)*.1;
% end
vol_a=nurb_volume(pp,M,N);
[uu,vv]=centripetal_param(pp,M,N);
knotsu=KVcent(M,degree,uu);
knotsv=KVcent(N,degree,vv);
knots={knotsu,knotsv};
srf=nrbmak(pp,knots);
nrbplot(srf,[M,N])
p=nrbeval(srf,{uu,vv});
[dsurf,d2surf]=nrbderiv_modified(srf);
[pnt, jac, hess] = nrbdeval_modified(srf, dsurf, d2surf, {uu,vv});
[normal,norm_dir]=nurb_normal(jac,M,N);
[radius]=gaussian_nurb_curvature(jac,hess,normal,M,N);
%vol=nurb_volume(p,M,N);
[m_max,n_max,m_min,n_min]=opt_radius(radius,M);
[p,vol]=correction(p,radius,norm_dir,M,N,degree,vol_a);
counter=1;
for loop=1:90
        for k=1:3
        p(k,m_max,n_max)=p(k,m_max,n_max)+norm_dir(k,m_max,n_max)*.01;
        end
        %% Kp parameter definition %%
%         [uu_temp,vv_temp]=centripetal_param(p,M,N);
%         knotsu_temp=KVcent(M,degree,uu);
%         knotsv_temp=KVcent(N,degree,vv);
%         knots_temp={knotsu,knotsv};
%         srf_temp=nrbmak(p,knots);
%         p1=nrbeval(srf_temp,{uu_temp,vv_temp});
%         vol_temp=nurb_volume(p1,M,N);
%         del_vol=abs(vol-vol_temp);
%         P = 100*sqrt(abs(del_vol));
       %%
        ashesh=1;
        kk=.001;
        count=1;
        e=[];
        while(ashesh==1||abs(vol-vol_new)>vol/100)
            ashesh=2;
        for k=1:3
        p(k,m_min,n_min)=p(k,m_min,n_min)-norm_dir(k,m_min,n_min)*kk;
        end
        [uu_temp,vv_temp]=centripetal_param(p,M,N);
        knotsu_temp=KVcent(M,degree,uu_temp);       % i added temp%
        knotsv_temp=KVcent(N,degree,vv_temp);
        knots_temp={knotsu,knotsv};
        srf_temp=nrbmak(p,knots);
        p1=nrbeval(srf_temp,{uu_temp,vv_temp});
       
%         [dsurf,d2surf]=nrbderiv_modified(srf_temp);
%        [pnt, jac, hess] = nrbdeval_modified(srf_temp, dsurf, d2surf, {uu,vv});
%        [normal,norm_dir]=nurb_normal(jac,M,N);
%        [radius]=gaussian_nurb_curvature(jac,hess,normal,M,N);
%        [m_max,n_max,m_min,n_min]=opt_radius(radius,M);
       % vol_new=nurb_volume(p1,M,N)+.0022*count;
       [p,vol_new]=correction(p,radius,norm_dir,M,N,degree,vol);
       vol_new
        kk=kk+(vol-vol_new)*.01
        count=count+1;
        counter=counter+1;
        e(count)=vol-vol_new;
        [dsurf,d2surf]=nrbderiv_modified(srf_temp);
[pnt, jac, hess] = nrbdeval_modified(srf_temp, dsurf, d2surf, {uu,vv});
[normal,norm_dir]=nurb_normal(jac,M,N);
[radius]=gaussian_nurb_curvature(jac,hess,normal,M,N);
 [m_max,n_max,m_min,n_min]=opt_radius(radius,M);
        end
%          if((vol-vol_new)>.1)
%            for k=1:3 
%               p(k,m_min,n_min)=p(k,m_min,n_min)-norm_dir(k,m_min,n_min)*kk;
%            end
%         end
%         kk=kk+(vol-vol_new)/100;
%         count=count+1;
%         end
        [uu,vv]=centripetal_param(p,M,N);
        knotsu=KVcent(M,degree,uu);
        knotsv=KVcent(N,degree,vv);
        knots={knotsu,knotsv};
        srf=nrbmak(p,knots);
        p=nrbeval(srf,{uu,vv});
        figure(1)
      
        nrbplot(srf,[M,N]);
        daspect([1 1 1])
        view(-60,30)
        
%         plot(radius(50,2:99))
%         axis([0 100 -2 2])
% plot3(p(1,:),p(2,:),p(3,:));
       % hold on
%         plot3(p(1,m_max,n_max),p(2,m_max,n_max),p(3,m_max,n_max),'r*');
%         plot3(p(1,m_min,n_min),p(2,m_min,n_min),p(3,m_min,n_min),'r*');
%          quiver3(pnt(1,m_min,n_min),pnt(2,m_min,n_min),pnt(3,m_min,n_min),norm_dir(1,m_min,n_min),norm_dir(2,m_min,n_min),norm_dir(3,m_min,n_min));
        
       
       hold off
       
       vol_arr(loop)=vol_new;
       area_arr(loop)=nurb_peri(p,M,N);
%        CA=contact_angle(p,M,N);
%        contact_angle(loop)=CA;
       %contact_angle_circle(loop)=CA_circle;
[dsurf,d2surf]=nrbderiv_modified(srf);
[pnt, jac, hess] = nrbdeval_modified(srf, dsurf, d2surf, {uu,vv});
[normal,norm_dir]=nurb_normal(jac,M,N);
[radius]=gaussian_nurb_curvature(jac,hess,normal,M,N);
 [m_max,n_max,m_min,n_min]=opt_radius(radius,M);
end
