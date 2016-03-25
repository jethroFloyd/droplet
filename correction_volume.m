clc
clear;
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

%% introduce bump %%
 for k=1:N
      
      for j=1:M
       r=sqrt(pp(1,j,k)^2+pp(2,j,k)^2);   

          pp(3,j,k)=pp(3,j,k)*(1/(1+exp((-30*r))));
          
      end
  end

vol_a=nurb_volume(pp,M,N);

[uu,vv]=centripetal_param(pp,M,N);
knotsu=KVcent(M,degree,uu);
knotsv=KVcent(N,degree,vv);
knots={knotsu,knotsv};
srf=nrbmak(pp,knots);
p=nrbeval(srf,{uu,vv});
vol_new=nurb_volume(p,M,N);

%p=nrbeval(srf,{uu,vv});
[dsurf,d2surf]=nrbderiv_modified(srf);
[pnt, jac, hess] = nrbdeval_modified(srf, dsurf, d2surf, {uu,vv});
[normal,norm_dir]=nurb_normal(jac,M,N);
[radius]=gaussian_nurb_curvature(jac,hess,normal,M,N);
ashesh=1;
count=1;
while(ashesh==1||(vol_a-vol_new)>vol_a/10000)
    ashesh=2;
    P=.001;
    for i=2:N-1
        for j=2:M-1
            for k=1:3
                p(k,j,i)=p(k,j,i)+norm_dir(k,j,i)*(-sign(radius(k,i)))*P;
            end
        end
    end
    
[uu,vv]=centripetal_param(p,M,N);
knotsu=KVcent(M,degree,uu);
knotsv=KVcent(N,degree,vv);
knots={knotsu,knotsv};
srf=nrbmak(p,knots);
p=nrbeval(srf,{uu,vv});
vol_new=nurb_volume(p,M,N);
P=(vol_a-vol_new)/1000;
e(count)=vol_a-vol_new;

% [dsurf,d2surf]=nrbderiv_modified(srf);
% [pnt, jac, hess] = nrbdeval_modified(srf, dsurf, d2surf, {uu,vv});
% [normal,norm_dir]=nurb_normal(jac,M,N);
% [radius]=gaussian_nurb_curvature(jac,hess,normal,M,N);
count=count+1;
figure(1)
plot(p(1,:,50),p(3,:,50));


hold off
end
%[radius]=gaussian_nurb_curvature(jac,hess,normal,M,N);