clear all;
close all;
clc
addpath('C:\Program Files\MATLAB\R2011a\toolbox\nurbs_toolbox');
r=1;
M=100;
N=100;
u=linspace(-1,1,M);
uu=linspace(0,1,M);
vv=linspace(0,1,N);
v=linspace(-1,1,N);

% for i=1:N
%     for j=1:M
%         
%             pp(1,j,i)=(1-(u(i))^2)*(1-(v(j))^2)/((1+(u(i))^2)*(1+(v(j))^2));
%             pp(2,j,i)=2*v(j)*(1-(u(i))^2)/((1+(u(i))^2)*(1+(v(j))^2));
%             pp(3,j,i)=2*u(i)*(1+(v(j))^2)/((1+(u(i))^2)*(1+(v(j))^2));
%     end
% end
for i=1:N
    for j=1:M
        
            pp(3,j,i)=(1-(v(j))^2)*(1-(u(i))^2)/((1+(v(j))^2)*(1+(u(i))^2));
            pp(1,j,i)=2*u(i)*(1-(v(j))^2)/((1+(v(j))^2)*(1+(u(i))^2));
            pp(2,j,i)=2*v(j)*(1+(u(i))^2)/((1+(v(j))^2)*(1+(u(i))^2));
    end
end
degree=5;
knots=knotspan(M,N,degree);
srf=nrbmak(pp,knots);
%p=nrbeval(srf,{uu,vv});
nrbplot(srf,[M,N])
p=nrbeval(srf,{linspace(0,1,M),linspace(0,1,N)});
[dsurf,d2surf]=nrbderiv_modified(srf);
[pnt, jac, hess] = nrbdeval_modified(srf, dsurf, d2surf, {u,v});
[normal,norm_dir]=nurb_normal(jac,M,N);
[radius]=gaussian_nurb_curvature(jac,hess,normal,M,N);
%% dent %%
  for k=45:55
      for j=30:35
          for i=1:3
          p(i,j,k)=p(i,j,k)+norm_dir(i,j,k)*.2;
          
          end
      end
  end
  for k=35:40
      for j=30:35
          for i=1:3
          p(i,j,k)=p(i,j,k)-norm_dir(i,j,k)*.2;
          
          end
      end
  end
  
  srf=nrbmak(p,knots);
  nrbplot(srf,[M,N])
  [dsurf,d2surf]=nrbderiv_modified(srf);
  [pnt, jac, hess] = nrbdeval_modified(srf, dsurf, d2surf, {u,v});
  [normal,norm_dir]=nurb_normal(jac,M,N);
  [radius]=gaussian_nurb_curvature(jac,hess,normal,M,N);
  vol=nurb_volume(p,M,N);
 [m_max,n_max,m_min,n_min]=opt_radius(radius,M);
  filename= [ 'optimization','.avi'];
  vidObj = VideoWriter(filename);
 open(vidObj);
 for loop=1:41
        for k=1:3
        p(k,m_min,n_min)=p(k,m_min,n_min)-norm_dir(k,m_min,n_min)*.1;
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
        srf=nrbmak(pp,knots);
        p1=nrbeval(srf,{u,v})
        vol_new=nurb_volume(p1,M,N);
        if((vol_new-vol)>.1)
           for k=1:3 
              p(k,m_max,n_max)=p(k,m_max,n_max)-norm_dir(k,m_max,n_max)*kk;
           end
        end
        kk=kk+(vol_new-vol)/100;
        srf=nrbmak(p,knots);
       h=figure(1)
        nrbplot(srf,[2*M,2*N])
        hold off
        view(60,30)
        %axis([0 1 0 1 0 1])
        hold off
        count=count+1;
        end
       srf=nrbmak(p,knots);
       p=nrbeval(srf,{u,v});
       [dsurf,d2surf]=nrbderiv_modified(srf);
       [pnt, jac, hess] = nrbdeval_modified(srf, dsurf, d2surf, {u,v});
       [normal,norm_dir]=nurb_normal(jac,M,N);
       radius=nurb_curvature(jac,hess,normal,M,N);
       [m_max,n_max,m_min,n_min]=opt_radius(radius,M);
          F = getframe(h);

 end
     writeVideo(vidObj,F);