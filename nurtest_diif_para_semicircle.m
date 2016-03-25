clear all;
close all;
clc
addpath('C:\Program Files\MATLAB\R2011a\toolbox\nurbs_toolbox');
r=1;
M=100;
N=100;
u=linspace(-1,1,M);
v=linspace(-1,1,N);
uu=linspace(0,1,M);
vv=linspace(0,1,N);
% xx=linspace(0,1,2*M);
% yy=linspace(0,1,2*N);


for i=1:N
    for j=1:M
        
            pp(3,j,i)=(1-(u(i))^2)*(1-(v(j))^2)/((1+(u(i))^2)*(1+(v(j))^2));
            pp(1,j,i)=2*v(j)*(1-(u(i))^2)/((1+(u(i))^2)*(1+(v(j))^2));
            pp(2,j,i)=2*u(i)*(1+(v(j))^2)/((1+(u(i))^2)*(1+(v(j))^2));
    end
end
% for i=1:N
%     for j=1:M
%         
%             pp(3,j,i)=(1-(v(j))^2)*(1-(u(i))^2)/((1+(v(j))^2)*(1+(u(i))^2));
%             pp(1,j,i)=2*u(i)*(1-(v(j))^2)/((1+(v(j))^2)*(1+(u(i))^2));
%             pp(2,j,i)=2*v(j)*(1+(u(i))^2)/((1+(v(j))^2)*(1+(u(i))^2));
%     end
% end
degree=3;
knots=knotspan(M,N,degree);

% basis=basis_matrix(M,N,uu,degree,knots{1});
% Q=(reshape(pp,[3,M*N]))';
% basis2=kron(basis,basis);
% %basis2(1,1)=10^5;
% %basis2(length(basis2),length(basis2))=10^5;
% [L,U]=lu(basis2)
% % P=inv(basis2)*Q
% P=U\(L\Q);
%p_prime=reshape(P',[3,M,N]);

% uu(1,:)=0;
% uu(N,:)=0;
% 
% for i=1:N
%     x(1,:,i)=pp(1,:,i);
%     x(2,:,i)=pp(3,:,i);
%     uu(i,:)=centripetal_param(x);
% end

% for i=1:M
%  
%     P=(basis)\x(:,:,i);
%     control(1,:,i)=P(:,1);
%     control(3,:,i)=P(:,2);
%     control(2,:,i)=pp(2,:,i);
% end

%control=interpolation(pp,basis);
for i=1:M
crv{i}=nrbmak(pp(:,:,i),knots{1});
end
srf=nrbloft(crv);
%[uu]=centripetal_param(pp,M,N);
p=nrbeval(srf,{uu,vv});

surfl(squeeze(p(1,:,:)),squeeze(p(2,:,:)),squeeze(p(3,:,:)))
% for i=1:N
%     p(:,:,i)=bspeval(degree,crv{i}.coefs,knots{1},uu(i,:))
% end
% for i=1:10
%     p=nrbeval(srf,{uu,vv});
%     for z=1:M
%     crv{z}=nrbmak(p(:,:,z),knots{1});
%     end 
% srf=nrbloft(crv);
% nrbplot(srf,[2*M,2*N]);
% p1=nrbeval(srf,{linspace(0,1,2*M),linspace(0,1,2*N)})
% vol(i)=nurb_volume(p1,2*M,2*N)
% pause(.1)
% end

% nrbplot(srf,[M,N]);



%[p,N]=nrbeval(srf,{uu,vv});
% nrbplot(srf,[M,N])
% daspect([1 1 1])
[p,Basis]=nrbeval_octave(srf,{uu,vv});
[dsurf,d2surf]=nrbderiv_modified(srf);
[pnt, jac, hess] = nrbdeval_modified(srf, dsurf, d2surf, {uu,vv});
[normal,norm_dir]=nurb_normal(jac,M,N);
[radius]=gaussian_nurb_curvature(jac,hess,normal,M,N);
%% dent %%
%   for k=35:45
%       for j=20:45
%           for i=1:3
%           p(i,j,k)=p(i,j,k)+norm_dir(i,j,k)*.1;
%           
%           end
%       end
%   end
%   for k=35:45
%       for j=20:35
%           for i=1:3
%           p(i,j,k)=p(i,j,k)-norm_dir(i,j,k)*.1;
%           
%           end
%       end
%   end
for i=1:3
    p(i,31,31)= p(i,31,31)+norm_dir(i,31,31)*.1;
end
for i=1:3
    p(i,41,41)= p(i,41,41)-norm_dir(i,41,41)*.1;
end
%%
  for i=1:M
crv{i}=nrbmak(p(:,:,i),knots{1});
end
srf=nrbloft(crv);
%[uu]=centripetal_param(p,M,N)
p=nrbeval(srf,{uu,vv});


% nrbplot(srf,[M,N]);
%   Q=(reshape(p,[3,M*N]))';
% basis2=kron(basis,basis);
% [L,U]=lu(basis2);
% %basis2(1,1)=10^5;
% %basis2(length(basis2),length(basis2))=10^5;
% P=U\(L\Q)
% p_prime=reshape(P',[3,M,N]);
%   %control=interpolation(p,basis)
%   srf=nrbmak(p_prime,knots);
 h= figure(1)
%   nrbplot(srf,[M,N]);
surfl(squeeze(p(1,:,:)),squeeze(p(2,:,:)),squeeze(p(3,:,:)))
  %view(-150,30)
 
  [dsurf,d2surf]=nrbderiv_modified(srf);
  [pnt, jac, hess] = nrbdeval_modified(srf, dsurf, d2surf, {uu,vv});
  [normal,norm_dir]=nurb_normal(jac,M,N);
  [radius]=gaussian_nurb_curvature(jac,hess,normal,M,N);
  p1=nrbeval(srf,{uu,vv});
  vol=nurb_volume(p1,M,N);
 [m_max,n_max,m_min,n_min]=opt_radius(radius,M);
%   filename= [ 'optimization','.avi'];
%   vidObj = VideoWriter(filename);
%  open(vidObj);
 for loop=1:3
        for k=1:3
        p(k,m_max,n_max)=p(k,m_max,n_max)+norm_dir(k,m_max,n_max)*.01;
        end
   %      vol2=nurb_volume(p,M,N);
        ashesh=1;
        kk=.01
        count=1;
        while(ashesh==1||(vol_new-vol)>.1)
            ashesh=2;
        for k=1:3
        p(k,m_min,n_min)=p(k,m_min,n_min)+norm_dir(k,m_min,n_min)*kk;
        end
        %control=interpolation(p,basis)
      for i=1:M
      crv{i}=nrbmak(p(:,:,i),knots{1});
      end
         srf=nrbloft(crv);
         %[uu]=centripetal_param(p,M,N);
        [p1,base]=nrbeval_octave(srf,{uu,vv})
        vol_new=nurb_volume(p1,M,N);
        if((vol_new-vol)>.1)
           for k=1:3 
              p(k,m_min,n_min)=p(k,m_min,n_min)-norm_dir(k,m_min,n_min)*kk;
           end
        end
        kk=kk+(vol_new-vol)/100;
        
        %control=interpolation(p,basis)
        for i=1:M
        crv{i}=nrbmak(p(:,:,i),knots{1});
        end
        srf=nrbloft(crv);
        %[uu]=centripetal_param(p,M,N);
        p=nrbeval(srf,{uu,vv});
       h=figure(1)
        nrbplot(srf,[M,N])
        hold on
        plot3(p(1,:,31),p(2,:,31),p(3,:,31),'r*')
        hold on
        daspect([1 1 1])
        %surfl(squeeze(p(1,:,:)),squeeze(p(2,:,:)),squeeze(p(3,:,:)))
        %plot3(p(1,:),p(2,:),p(3,:))
%         hold on
        plot3(p(1,m_max,n_max),p(2,m_max,n_max),p(3,m_max,n_max),'b*')
        plot3(p(1,m_min,n_min),p(2,m_min,n_min),p(3,m_min,n_min),'b*')
        
        
 view(-90,30)
        %axis([0 1 0 1 0 1])
        hold off
        count=count+1;
        end
        area_array(loop)=nurb_peri(p,M,N);
        vol_arr(loop)=vol_new;
%         figure(1)
%         plot(area_array);
%         hold on
       %control=interpolation(p,basis)
       for i=1:M
       crv{i}=nrbmak(p(:,:,i),knots{1});
      end
srf=nrbloft(crv);
%[uu]=centripetal_param(p,M,N);
       [p,basis]=nrbeval_octave(srf,{uu,vv})
       [dsurf,d2surf]=nrbderiv_modified(srf);
       [pnt, jac, hess] = nrbdeval_modified(srf, dsurf, d2surf, {uu,vv});
       [normal,norm_dir]=nurb_normal(jac,M,N);
       radius=nurb_curvature(jac,hess,normal,M,N);
       [m_max,n_max,m_min,n_min]=opt_radius(radius,M);
%           F = getframe(h);

 end
%      writeVideo(vidObj,F);