function[knotsu]=KVcent(n,degree,u)
for j=1:degree+1
    knotsu(j)=0;
end
% KV = zeros(n+degree+1,1); 

% for i=size(KV,1)-p:size(KV,1) 
% KV(i) = 1.0; 
% end 
% 
% for j=1:size(u,1)-p 
% tmp = 0.0; 
% for i=j:j+p-1 
% tmp = tmp+u(i); 
% end 
% KV(j+p) = tmp/p; 
% end 
% end
count=1;
for j=degree+2:(n)
    knotsu(j)=(1/(degree))*sum(u(count+1:count+degree));
    count=count+1;
end
for j=(n+1):(n+degree+1)
    knotsu(j)=1;
end