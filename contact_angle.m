function[CA]=contact_angle(p,M,N)
x=p(1,:,M/2);
z=p(1,:,N/2);
CA=atand((z(2)-z(1))/(x(2)-x(1)));
%[xc,yc,Re]=circfit(x,z);
 %xfit = Re*cosd(0:.1:180)+xc; yfit = Re*sind(0:.1:180)+yc;
%         ydash = sqrt(Re^2-yc^2)/yc;
%         CA_circle=atand(ydash);
end
