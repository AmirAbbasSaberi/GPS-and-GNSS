function Tr_s = Sustamiinen(pos,E)
wgs84 = wgs84Ellipsoid('kilometer');
[latp,~,hp] = ecef2geodetic(wgs84,pos(1),pos(2),pos(3));
hp = hp/1000;
g = 1 - 0.0026*cosd(2*latp)-0.00000028*hp;
D = 2;
% pr = 1013.25;
% Tr = 18+273;
% hr = 0;
% Hr = 50/100;
% p = pr*(1-0.0000226*(hp-hr)).^5.225;
% T = Tr - 0.0065*(hp-hr);
% H = Hr*e.^(-0.0006396*(hp-hr));

M = 1.001./(sqrt(0.002001+sin(E).^2));

[T,P,e,~,~] = linear_interpolation(deg2rad(latp),D);
Tr_dw = 0.002277./(g).*(P+(1255./T+0.05).*e);
Tr_s = (Tr_dw).*M;