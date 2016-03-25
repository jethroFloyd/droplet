clear all;
close all;
[v, f, n, c, stltitle]=stlread('C:\Users\BabaChattopadhyay\Desktop\3-d droplet\test_nurb_solid.STL',false);
%%[v1, f1, n1, c1, stltitle1]=stlread('C:\Users\BabaChattopadhyay\Desktop\3-d droplet\droplet.STL',false);
v_new=search(v);
% v1_new=search(v1);
vv(:,1)=v_new(:,1);
vv(:,2)=v_new(:,3);
vv(:,3)=v_new(:,2);
% alphavol(vv,100,1)
% 
% alphavol(vv,100,1);
% degree=4;
% knots=[];
% N=length(vv);
% knot_len=degree+N-1;
% for i=1:3
%     knots(i)=0;
% end
% count=1;
% %     for i=5:(knot_len-degree)
% %         knots(i)=count;
% %         count=count+1;
% %     end
% %     for i=(knot_len-degree)+1:knot_len
% %         knots(i)=count;
% %     end
% %         
% 
% %% nurbgrid_x %%
% v_temp=vv;
% v_new_nurb=[];
% for i=1:(length(vv))
%     for j=1:(length(vv))
%         if j~=i
%     if v_temp(j,1)==v_temp(i,1);
%         v_temp(j,1)=0;
%     end
%         end
%     end
% end
% count=1;
% for i=1:length(v_temp)
%     if v_temp(i,1)~=0
%         v_new_nurb(count)=v_temp(i,1);
%         count=count+1;
%     end
% end
% for i=1:length(vv)
%     if vv(j,2)<1/10000000000
%         vv(j,2)=0;
%     end
% end
% 
% %% nurbgrid_y %%
% count=1;
% 
% for i=1:length(v_new_nurb)
%     for j=1:length(vv)
%         if v_new_nurb(i)==vv(j,1);
%             xy_new_nurb(count,1)=v_new_nurb(i);
%             xy_new_nurb(count,2)=vv(j,2);
%             xy_new_nurb(count,3)=vv(j,3);
%             count=count+1;
%         end
%     end
%     
% end
% count=1;
% for i=1:length(xy_new_nurb(:,1))-1
%     if xy_new_nurb(i,1)~=xy_new_nurb(i+1,1)
%         count=count+1;
%     end
% end
% countt=1;
% counter=1;
% for i=1:length(v_new_nurb)
%     for j=1:length(xy_new_nurb(:,1))
%         if v_new_nurb(i)==xy_new_nurb(j,1);
%         final_matrix(2,:,i)=xy_new_nurb(j,2);
%         end
%     end
% end
        

% for i=1:length(v_new_nurb
% final_matrix(1,:)=v_new_nurb(:,1);
% for i=1:length(xy_nurb_new)
%     if xy_nurb_new(i,1)=v_new_nurb(i)
%         final_matrix(2,:)=

% tri=delaunay(v_new(:,1),v_new(:,2),v_new(:,3));

%  trisurf(tri,v_new(:,1),v_new(:,2),v_new(:,3));
%  hold on

 
% v_new(10,3)=v_new(10,3)*10;

 
 DT=delaunayTri(v)
tri=delaunay(v_new(:,1),v_new(:,3),v_new(:,2));
tr1=delaunay(v1_new(:,3),v1_new(:,2),v1_new(:,1));
figure(1)
subplot(2,1,1)
trisurf(tri,v_new(:,1),v_new(:,3),v_new(:,2))
subplot(2,1,2)
trisurf(tri,v1_new(:,1),v1_new(:,2),v1_new(:,3))
% % hh=alphaShape(v_new(:,1),v_new(:,3),v_new(2))
% [zgrid,xgrid,ygrid]=gridfit(v_new(:,1),v_new(:,3),v_new(:,2),f(:,1),f(:,3))
% tr=TriRep(tri,v_new(:,1),v_new(:,3),v_new(:,2));
% [ttt, Xb] = freeBoundary(tr)