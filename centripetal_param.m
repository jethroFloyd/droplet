function[u_param,v_param]=centripetal_param(pp,M,N)

uu(1,:)=0;
uu(N,:)=0;
uu(:,1)=0;
uu((2:N-1),N)=1;
d=[];
for j=2:N-1
    d(j)=0;
    for i=2:M-1
        d(j)=d(j)+sqrt(norm(pp(:,i,j)-pp(:,i-1,j)));
         %d(j)=d(j)+(norm(pp(:,i,j)-pp(:,i-1,j)));
    end
end
        
    
for j=2:N-1
    for i=2:M-1
        uu(j,i)=uu(j,i-1)+ (sqrt(norm(pp(:,i,j)-pp(:,i-1,j))))/d(j);
        %uu(j,i)=uu(j,i-1)+ ((norm(pp(:,i,j)-pp(:,i-1,j))))/d(j);
    end
end
%u_param=uu;
%% avg %%
u_param(1)=0;
u_param(M)=1;
for i=2:M-1
    u_param(i)=sum(uu(:,i))/M;
end
vv(1,:)=0;
vv(N,:)=0;
vv(:,1)=0;
vv((2:N-1),N)=1;
c=[];
for i=2:N-1
    c(i)=0;
    for j=2:M-1
        c(i)=c(i)+sqrt(norm(pp(:,i,j)-pp(:,i,j-1)));
         %c(i)=c(i)+(norm(pp(:,i,j)-pp(:,i,j-1)));
    end
end
for j=2:N-1
    for i=2:M-1
        vv(j,i)=vv(j,i-1)+ (sqrt(norm(pp(:,j,i)-pp(:,j,i-1))))/c(j);
        %vv(j,i)=vv(j,i-1)+ ((norm(pp(:,j,i)-pp(:,j,i-1))))/c(j);
    end
end
%v_param=vv;
v_param(1)=0;
v_param(M)=1;
for i=2:M-1
    v_param(i)=sum(vv(:,i))/N;
end
end