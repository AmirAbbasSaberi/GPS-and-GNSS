function [I1,E,A] =  Klobuchar(sat_gps1,pos,filename,c)
wgs84 = wgs84Ellipsoid('kilometer');
x = sat_gps1(:,21);
y = sat_gps1(:,22);
z = sat_gps1(:,23);


[latp,lonp,~] = ecef2geodetic(wgs84,pos(1),pos(2),pos(3));

R = [-sind(lonp) cosd(lonp) 0
    -cosd(lonp)*sind(latp) -sind(lonp)*sind(latp) cosd(latp)
    cosd(lonp)*cosd(latp) sind(lonp)*cosd(latp) sind(latp)];
r1 = [x y z];
r1 = r1';
r0 = repmat(pos,1,size(r1,2));
dr = r1-r0;
dr_prime = R*dr;
dr_prime = dr_prime';
e_dr = dr_prime./repmat(sqrt(sum((dr_prime.^2),2)),1,3);
az = atan2d(dr_prime(:,1),dr_prime(:,2));
elv = 90-acosd(dr_prime(:,3)./sqrt(sum((dr_prime.^2),2)));
for j = 1:32
    i = find(sat_gps1(:,1) == j);
    dr1 = dr';
    dr1 = sqrt(sum((dr1.^2),2));
    dr1 = dr1(i).*cos(deg2rad(elv(i)));
    i_az = deg2rad(az(i));
    i_elv = deg2rad(elv(i));
%     polarplot(i_az,dr1,'.')
%     hold on
end
% hold off


E = deg2rad(elv);
A = deg2rad(az);
RE = 6378;
h = 350;
si = pi/2-E-asin(RE/(RE+h)*cos(E));
phi_i = asin(sind(latp).*cos(si)+cosd(latp).*sin(si).*cos(A));
landa_i = deg2rad(lonp)+(si.*sin(A))./(cos(phi_i));
phi_m = rad2deg(asin(sin(phi_i).*sin(deg2rad(78.3))+cos(phi_i).*cos(deg2rad(78.3)).*cos(landa_i-deg2rad(291.0))));

JD = juliand(sat_gps1(:,2),sat_gps1(:,3),sat_gps1(:,4),sat_gps1(:,5),sat_gps1(:,6),sat_gps1(:,7)) ;
JD0 = juliand(1980,1,6,0,0,0) ;
gps_week = floor( (JD - JD0) / 7 );
t_gps = round(((((JD-2444244.5)/7)-gps_week)*24*60*60*7)/0.5)*0.5; %% GPS secend

t = 43200*landa_i/pi+t_gps;
i = find(t > 86400);
if isempty(i) == 0
    t(i) = t(i) -86400;
end
i = find(t < 86400);
if isempty(i) == 0
    t(i) = t(i) +86400;
end

[a,b]=ion(filename);
AI = sum(a(1)+a(2)*phi_m/pi+a(3)*(phi_m/pi).^2+a(4)*(phi_m/pi).^3,2);
i = AI<0;
AI(i,1) = 0;
PI = sum(b(1)+b(2)*phi_m/pi+b(3)*(phi_m./pi).^2+b(4)*(phi_m/pi).^3,2);
i = PI<72000;
PI(i,1) = 0;
XI = (2*pi*(t-50400)./PI);
F = (1-(RE/(RE+h)*cos(E)).^2).^(-0.5);
i = find(abs(XI)<pi/2);
I1(i,1) = (5*10^-9+AI(i).*cos(XI(i))).*F(i);
i = find(abs(XI)>=pi/2)*c;
I1(i,1) = (5*10^-9)*F(i)*c;



