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
[uu,vv]=centripetal_param(pp,M,N);
knotsu=KVcent(M,degree,uu);
knotsv=KVcent(N,degree,vv);
knots={knotsu,knotsv};
srf=nrbmak(pp,knots);
p=nrbeval(srf,{uu,vv});
p=nrbeval(srf,{uu,vv});
[dsurf,d2surf]=nrbderiv_modified(srf);
[pnt, jac, hess] = nrbdeval_modified(srf, dsurf, d2surf, {uu,vv});
[normal,norm_dir]=nurb_normal(jac,M,N);
[radius]=gaussian_nurb_curvature(jac,hess,normal,M,N);
[radius]=gaussian_nurb_curvature(jac,hess,normal,M,N);
%% dent %%
mid=50;
count=0;
  for k=45:50
      x=(mid-count)/10
      for j=mid-count:mid+count
          
%           for i=1:3
%           p(i,j,k)=p(i,j,k)+norm_dir(i,j,k)*.1;
          
%           end
          p(3,j,k)=p(3,j,k)*(1/(1+exp((-x))));
          if j>=50
              x=x+1;
          else
              x=x-1;
      end
      
      end
      count=count+1;
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
[uu,vv]=centripetal_param(p,M,N);
knotsu=KVcent(M,degree,uu);
knotsv=KVcent(N,degree,vv);
knots={knotsu,knotsv};
srf=nrbmak(p,knots);
nrbplot(srf,[M,N])
p=nrbeval(srf,{uu,vv});
[dsurf,d2surf]=nrbderiv_modified(srf);
[pnt, jac, hess] = nrbdeval_modified(srf, dsurf, d2surf, {uu,vv});
[normal,norm_dir]=nurb_normal(jac,M,N);
[radius]=gaussian_nurb_curvature(jac,hess,normal,M,N);
vol=nurb_volume(p,M,N);
[m_max,n_max,m_min,n_min]=opt_radius(radius,M);
for loop=1:90
        for k=1:3
        p(k,m_max,n_max)=p(k,m_max,n_max)-norm_dir(k,m_max,n_max)*.01;
        end
         ashesh=1;
        kk=.001
        count=1;
        while(ashesh==1||(vol-vol_new)>.1)
            ashesh=2;
        for k=1:3
        p(k,m_min,n_min)=p(k,m_min,n_min)+norm_dir(k,m_min,n_min)*kk;
        end
        [uu_temp,vv_temp]=centripetal_param(p,M,N);
        knotsu_temp=KVcent(M,degree,uu);
        knotsv_temp=KVcent(N,degree,vv);
        knots_temp={knotsu,knotsv};
        srf_temp=nrbmak(p,knots);
        p1=nrbeval(srf_temp,{uu_temp,vv_temp});
        vol_new=nurb_volume(p1,M,N);
         if((vol-vol_new)>.1)
           for k=1:3 
              p(k,m_min,n_min)=p(k,m_min,n_min)-norm_dir(k,m_min,n_min)*kk;
           end
        end
        kk=kk+(vol-vol_new)/100;
        count=count+1;
        end
        [uu,vv]=centripetal_param(p,M,N);
        knotsu=KVcent(M,degree,uu);
        knotsv=KVcent(N,degree,vv);
        knots={knotsu,knotsv};
        srf=nrbmak(p,knots);
        p=nrbeval(srf,{uu,vv});
        figure(1)
      
        nrbplot(srf,[M,N])
        daspect([1 1 1])
        view(-60,-60)
       
       hold off
       vol_arr(loop)=vol_new;
       area_arr(loop)=nurb_peri(p,M,N);
       [CA,CA_circle]=contact_angle(p,M,N);
       contact_angle(loop)=CA;
       contact_angle_circle(loop)=CA_circle;
[dsurf,d2surf]=nrbderiv_modified(srf);
[pnt, jac, hess] = nrbdeval_modified(srf, dsurf, d2surf, {uu,vv});
[normal,norm_dir]=nurb_normal(jac,M,N);
[radius]=gaussian_nurb_curvature(jac,hess,normal,M,N);
 [m_max,n_max,m_min,n_min]=opt_radius(radius,M);
end
