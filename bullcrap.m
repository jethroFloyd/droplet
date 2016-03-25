clear all
close all
clc
addpath('C:\Program Files\MATLAB\R2011a\toolbox\nurbs_toolbox');
r=1;
M=200;
N=200;
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

% for i=1:M
% crv{i}=nrbmak(pp(:,:,i),knots{1});
% end
srf=nrbmak(pp,knots);
 p=nrbeval(srf,{uu,vv});
for j=1:100
    [uu,vv]=centripetal_param(p,M,N);
   knotsu=KVcent(M,degree,uu);
   knotsv=KVcent(N,degree,vv);
    knots={knotsu,knotsv};
    
%     for i=1:M
%         crv{i}=nrbmak(p(:,:,i),knots{1});
%     end
    srf=nrbmak(p,knots);
    p=nrbeval(srf,{uu,vv});
    
    vol(j)=nurb_volume(p,M,N);
    figure(1)
    nrbplot(srf,[M,N]);
    hold off
    
%     pause(.5)
end