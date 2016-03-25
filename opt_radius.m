function[m_max,n_max,m_min,n_min]=opt_radius(radius,M)
count=1;
for i=2:M-1
temp(:,count)=radius(2:M-1,i);
count=count+1;
end
%temp(:,minimum_z)=NaN;

[max_r,ind]=max(temp(:));
[m_max,n_max]=ind2sub(size(temp),ind);
[min_r,ind]=min(temp(:));
[m_min,n_min]=ind2sub(size(temp),ind);
m_max=m_max+1;
n_max=n_max+1;
m_min=m_min+1;
n_min=n_min+1;
end