function[control]=interpolation(pp,basis)

for i=1:length(pp)
    x(:,1,i)=pp(1,:,i);
    x(:,2,i)=pp(3,:,i);
end
for i=1:length(pp)
 
    P=(basis)\x(:,:,i);
    control(1,:,i)=P(:,1);
    control(3,:,i)=P(:,2);
    control(2,:,i)=pp(2,:,i);
end
end