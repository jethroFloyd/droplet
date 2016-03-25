function[p,vol_new]=correction(p,radius,norm_dir,M,N,degree,vol_a)
for loop=1:200
   
    P=.0001;
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
%e(count)=vol_a-vol_new;

% [dsurf,d2surf]=nrbderiv_modified(srf);
% [pnt, jac, hess] = nrbdeval_modified(srf, dsurf, d2surf, {uu,vv});
% [normal,norm_dir]=nurb_normal(jac,M,N);
% [radius]=gaussian_nurb_curvature(jac,hess,normal,M,N);
%count=count+1;
% figure(1)
% plot(p(1,:,50),p(3,:,50));
% 
% 
% hold off
end
%[radius]=gaussian_nurb_curvature(jac,hess,normal,M,N);