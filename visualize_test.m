clear all;
close all;
addpath('C:\Program Files\MATLAB\R2011a\toolbox\nurbs_toolbox')

base_radius=5;


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
r_norm=r./max(r);
theta_norm=theta./max(theta);
%% algorithm %%
for i=1:M
    
        for k=1:N
            pp(1,k,i)=r(k);
        end
        for k=1:N
            pp(2,k,i)=sind(theta(i));
        end
        for k=1:N
            pp(3,k,i)=r(k)*sind(theta(i));
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
p=nrbeval(srf,{r_norm,theta_norm});

nrbplot(srf,[100 100])