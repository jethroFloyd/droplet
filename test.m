close all;
clear all;
% 
% theta=360;
% N=30;
% R=1;
% 
% alpha=90-theta:2*theta/N:90+theta;
% 
% pp = [R*cosd(alpha') R*sind(alpha')-R*cosd(theta)];
% x=pp(:,1);
% y=pp(:,2);
% count=1;
% for i=1:length(pp(:,1))
%     for j=1:length(pp(:,2))
% z(i,j)=sqrt(x(i)^2+y(j)^2);
% end
% end
% figure(1)
% surf(x,y,z)
% % for i=5:11
% %     for j=5:11
% %         z(i,j)=z(i,j)+2;
% %     end
% % end
% % figure(2)
% % surf(x,y,z)
% 
% R = 7;
% [X,Y] = meshgrid(-10:.1:10);
% Z = sqrt(R.^2 - X.^2 - Y.^2);
% Z(imag(Z) ~= 0) = 0;
% 
% surfnorm(Z);
R=5;
x=-R:.1:R;
y=-R:.1:R;
pp=[];

count=1;

for i=1:length(y)
    for j=1:length(x)
        if((R^2-x(j)^2-y(i)^2)>=0)
            pp(count,1)=x(j);
            pp(count,2)=y(i);
            count=count+1;
        end
    end
end
for i=1:length(pp(:,1))
    z(i)=sqrt(R^2-pp(i,1)^2-pp(i,2)^2);
end
pp(:,3)=z(1,:);
%% using de-triang method %
DT =  delaunayTri(pp);
tri=delaunay(pp(:,1),pp(:,2),pp(:,3));
%% using meshgrids %%
[xx,yy]=meshgrid(min(pp(:,1)):.1:max(pp(:,1)),min(pp(:,2)):.1:max(pp(:,2)));
zz = griddata(pp(:,1),pp(:,2),pp(:,3),xx,yy);
%% using triplets %%
% for i=1:length(pp(:,1))
%       Z(i,i)=sqrt(R^2-pp(i,1)^2-pp(i,2)^2);
% end
% surf(pp(:,1).pp(:,2),Z);

% trisurf(tri,pp(:,1),pp(:,2),pp(:,3))
surf(xx,yy,zz);
%% tri-scatter %%
% F = TriScatteredInterp(pp(:,1), pp(:,2), pp(:,3));
% zz=F(xx,yy);


%        surf(xx,yy,zz)  
% out= tricurv_v01(tri,pp);

% [nx ny nz]=surfnorm(xx,yy,zz);
 [Cmean,Cgaussian,Dir1,Dir2,Lambda1,Lambda2]=curv_test(DT,3)
