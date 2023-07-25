function Tr = Collins(pos,E)
wgs84 = wgs84Ellipsoid('kilometer');
[latp,~,hp] = ecef2geodetic(wgs84,pos(1),pos(2),pos(3));
k1 = 77.604;
k2 = 382000;
Rd = 287.054;
gm = 9.784;
g = 9.80665;
D = 2;

[T,P,e,B,l] = linear_interpolation(deg2rad(latp),D);
Tr_z0_d = (10^-6)*k1*Rd*P/gm;
Tr_z0_w = (10^-6)*k2*Rd*e/(((l+1)*gm-B*Rd)*T);

Tr_z_d = ((1-(B*h/T))^(g/(Rd*B)))*Tr_z0_d;
Tr_z_w = ((1-(B*h/T))^(((l+1)*g/(Rd*B))-1))*Tr_z0_w;

M = 1.001./(sqrt(0.002001+sin(E).^2));

Tr = (Tr_z_d+Tr_z_w).*M;

