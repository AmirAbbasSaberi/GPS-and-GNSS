clc
clear
close all
format long g
% git remote add origin git@github.com:Arsenic2000/GPS_matlab_master.git
% git branch -M main
% git push -u origin main
%%
tic
path = 'mizu0770.15o';
[pos,sat] = mizuRinexExtraction(path);
toc
clear P1 P2 L1 L2 L5 n p  k i hour minute sec month line s S1 S2 S5 year t C1 C2 C5 number number1 number2 sat_num day end_of_header
close all
c= 299792458;
we = 7.29211514671*10^-5;
iP1 = find(strcmpi(sat(1,:) , 'P1'));
gnss_name = char(sat(:,1));
index = find(strcmpi(gnss_name(:,1),"G") );
sat_gps = sat(index,:);
gps_name = char(sat_gps(:,1));
sat_gps = [str2double(string(gps_name(:,2:3))) cell2mat(sat_gps(:,2:end))];
t_rcv = sat_gps(:,end);
P = sat_gps(:,iP1);
t_emission = t_rcv - P/c;
filename = 'brdc0770.15n';
sat_gps(:,end+1) = t_emission;
[data] = Extract_Navigation_data (filename);

sat_gps = gps_time_delay_coordiante(data,sat_gps,pos,c,we);


sat_gps (:,30) = sat_gps (:,24)+sat_gps (:,25);

i = find(sat_gps(:,1) == 5);
deltaT_sat = sat_gps(i,end-2);
time = sat_gps(i,5)*3600 + sat_gps(i,6)*60 + sat_gps(i,7);


%% Bancroft

clear aj B
n = 6;
B = [sat_gps(1:n,21:23),sat_gps(1:n,14)+sat_gps(1:n,26)+sat_gps(1:n,end-1)*c];



M = eye(4);
M(4,4) = -1;
for i = 1:n
    b = B(i,:);
    aj(i,1) = 0.5*b*M*b';
    
end


syms x

eq1 = Lorentz_dot((B'*B)^-1*B'*ones(n,1),(B'*B)^-1*B'*ones(n,1),M)*x^2 ...
    +2*(Lorentz_dot((B'*B)^-1*B'*ones(n,1),(B'*B)^-1*B'*aj,M)-1)*x ...
    +Lorentz_dot((B'*B)^-1*B'*aj,(B'*B)^-1*B'*aj,M) == 0;

lambda = solve(eq1);
lambda = max(double(lambda));


r_cdt = M*(B'*B)^-1*B'*(lambda*ones(n,1)+aj);
ground_pos0 = r_cdt(1:3);

diff_pos = ground_pos0 - pos;


sat_gps(:,end+1) = sqrt(sum(sat_gps(:,21:23).^2,2));

k = sat_gps(:,14) - sat_gps(:,28);
s = [k sat_gps(:,26)];
p = k+sat_gps(:,26);
i = find(p<0);
sat_gps(i,:) = [];
clear aj B
n = numel(sat_gps(:,1));
B = [sat_gps(1:n,21:23),sat_gps(1:n,14)+sat_gps(1:n,26)];



M = eye(4);
M(4,4) = -1;
for i = 1:n
    b = B(i,:);
    aj(i,1) = 0.5*b*M*b';
    
end


syms x

eq1 = Lorentz_dot((B'*B)^-1*B'*ones(n,1),(B'*B)^-1*B'*ones(n,1),M)*x^2 ...
    +2*(Lorentz_dot((B'*B)^-1*B'*ones(n,1),(B'*B)^-1*B'*aj,M)-1)*x ...
    +Lorentz_dot((B'*B)^-1*B'*aj,(B'*B)^-1*B'*aj,M) == 0;

lambda = solve(eq1);
lambda = max(double(lambda));


r_cdt = M*(B'*B)^-1*B'*(lambda*ones(n,1)+aj);
ground_pos0 = r_cdt(1:3);

diff_pos = ground_pos0 - pos;



%% Tr Ion

%-------------- Klobuchar Algorithm ---------------
wgs84 = wgs84Ellipsoid('kilometer');
x = sat_gps(:,21);
y = sat_gps(:,22);
z = sat_gps(:,23);
[lat,lon,h] = ecef2geodetic(wgs84,x,y,z);

[latp,lonp,hp] = ecef2geodetic(wgs84,pos(1),pos(2),pos(3));

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
    i = find(sat_gps(:,1) == j);
    dr1 = dr';
    dr1 = sqrt(sum((dr1.^2),2));
    dr1 = dr1(i).*cos(deg2rad(elv(i)));
    i_az = deg2rad(az(i));
    i_elv = deg2rad(elv(i));
    polarplot(i_az,dr1,'.')
    hold on
end
hold off


E = deg2rad(elv);
A = deg2rad(az);
RE = 6378;
h = 350;
si = pi/2-E-asin(RE/(RE+h)*cos(E));
phi_i = asin(sind(latp).*cos(si)+cosd(latp).*sin(si).*cos(A));
landa_i = deg2rad(lonp)+(si.*sin(A))./(cos(phi_i));
phi_m = rad2deg(asin(sin(phi_i).*sin(deg2rad(78.3))+cos(phi_i).*cos(deg2rad(78.3)).*cos(landa_i-deg2rad(291.0))));

JD = juliand(sat_gps(:,2),sat_gps(:,3),sat_gps(:,4),sat_gps(:,5),sat_gps(:,6),sat_gps(:,7)) ;
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
i = find(AI<0);
AI(i,1) = 0;
PI = sum(b(1)+b(2)*phi_m/pi+b(3)*(phi_m./pi).^2+b(4)*(phi_m/pi).^3,2);
i = find(PI<72000);
PI(i,1) = 0;
XI = (2*pi*(t-50400)./PI);
F = (1-(RE/(RE+h)*cos(E)).^2).^(-0.5);
i = find(abs(XI)<pi/2);
I1(i,1) = (5*10^-9+AI(i).*cos(XI(i))).*F(i);
i = find(abs(XI)>=pi/2);
I1(i,1) = (5*10^-9)*F(i);

figure(2)
for j = 1:32
    i = find(sat_gps(:,1) == j);
    time = 3600*sat_gps(i,5)+60*sat_gps(i,6)+sat_gps(i,7);
    I = c*I1(i);
    plot(time,I,'.')
    hold on
end
hold off
grid on 
title('Klobuchar Ionophere model')
ylabel('Error_(_m_)')
xlabel('time_(_s_)')

%%
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


figure(3)
for j = 1:32
    i = find(sat_gps(:,1) == j);
    time = 3600*sat_gps(i,5)+60*sat_gps(i,6)+sat_gps(i,7);
    I = Tr(i);
    plot(time,I,'.')
    grid on
    hold on
end
hold off
title('Collins Troposphere model')
ylabel('Error_(_m_)')
xlabel('time_(_s_)')
x_sat = sat_gps(:,21);
y_sat = sat_gps(:,22);
z_sat = sat_gps(:,23);
ro0 = sqrt((pos(1)-x_sat).^2+(pos(2)-y_sat).^2+(pos(3)-z_sat).^2);
A = [(pos(1)-x_sat)./ro0 (pos(2)-y_sat)./ro0 (pos(3)-z_sat)./ro0];
P = 400*(A'*A)^-1;



j = [1 -1 0 0 0 0;0 0 1 -1 0 0;0 0 0 0 1 -1];
d = j*[P zeros(3);zeros(3) P]*j';
rmse = sqrt(trace(d))*1000

%% Sustamiinen
hp = hp/1000;
g = 1 - 0.0026*cosd(2*latp)-0.00000028*hp;
% pr = 1013.25;
% Tr = 18+273;
% hr = 0;
% Hr = 50/100;
% p = pr*(1-0.0000226*(hp-hr)).^5.225;
% T = Tr - 0.0065*(hp-hr);
% H = Hr*e.^(-0.0006396*(hp-hr));
[T,P,e,B,l] = linear_interpolation(deg2rad(latp),D);
Tr_dw = 0.002277./(g).*(P+(1255./T+0.05).*e);
Tr_s = (Tr_dw).*M;

figure(4)
for j = 1:32
    i = find(sat_gps(:,1) == j);
    time = 3600*sat_gps(i,5)+60*sat_gps(i,6)+sat_gps(i,7);
    I = Tr_s(i);
    plot(time,I,'.')
    grid on
    hold on
end
title('Sustamiinen Troposphere model')
ylabel('Error_(_m_)')
xlabel('time_(_s_)')
hold off


%% Hopfield

ZHD = (0.62291/T+0.0023081)*P;
ZWD = (555.7+1.792*10^-4*exp((T-273)/22.9))*e/T^2;
Tr_h = 0.5*(ZHD+ZWD).*M;

figure(5)
for j = 1:32
    i = find(sat_gps(:,1) == j);
    time = 3600*sat_gps(i,5)+60*sat_gps(i,6)+sat_gps(i,7);
    I = Tr_h(i);
    plot(time,I,'.')
    grid on
    hold on
end
hold off
title('Hopfield Troposphere model')
ylabel('Error_(_m_)')
xlabel('time_(_s_)')
% m = (1+(a/(1+b/(1+c))))./((sin(E)+a./(sin(E)+b./(sin(E)+c))));
% 
% .
% 
% M_neil = M_w+M_d;