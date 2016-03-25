function[v_new]=search(v)
v_temp=v;
v_new=[];
for i=1:(length(v))
    for j=1:(length(v))
        if j~=i
    if v_temp(j,:)==v_temp(i,:);
        v_temp(j,:)=0;
    end
        end
    end
end
count=1;
for i=1:length(v_temp)
    if v_temp(i,:)~=0
        v_new(count,:)=v_temp(i,:);
        count=count+1;
    end
end
end
