clear all;
close all;
N=20;
% R=3;
% theta=0:.1:360;
% x=R*cosd(theta);
% y=R*sind(theta);
% [X,Y]=meshgrid(min(x):.1:max(x),min(y):.1:max(y));
% Z=sqrt(R^2-X.^2-Y.^2);
% surf(X,Y,Z)
[x y z]=sphere(N)
xx=x((N/2+1):end,:);
yy=y((N/2+1):end,:);
zz=z((N/2+1):end,:);


[nx,ny,nz]=surfnorm(xx,yy,zz);
[ii,jj]=find(zz==0);
for i=1:length(ii)
    x_pinned(i)=xx(ii(i),jj(i));
    y_pinned(i)=yy(ii(i),jj(i));
end
% for i=1:length(ii)
%     xx(ii(i),jj(i))=xx(ii(i),jj(i))+.01*nx(ii(i),jj(i))/sqrt(nx(ii(i),jj(i)).^2+nx(ii(i),jj(i)).^2+nx(ii(i),jj(i)).^2);
%     yy(ii(i),jj(i))=yy(ii(i),jj(i))+.01*ny(ii(i),jj(i))/sqrt(nx(ii(i),jj(i)).^2+nx(ii(i),jj(i)).^2+nx(ii(i),jj(i)).^2);
%     zz(ii(i),jj(i))=zz(ii(i),jj(i))+.01*nz(ii(i),jj(i))/sqrt(nx(ii(i),jj(i)).^2+nx(ii(i),jj(i)).^2+nx(ii(i),jj(i)).^2);
% end

% for i=1:3:N/2
%     for j=1:3:N/2
%        if zz(i,j)==0
%            zz(i,j)=zz(i,j)+0;
%        else
%         zz(i,j)=zz(i,j)+.3;
%        end
%     end
% end
[K,H,p1,p2]=surfature_new(xx,yy,zz);
surf(xx,yy,zz)
hold on
quiver3(xx,yy,zz,nx,ny,nz)
daspect([1 1 1])