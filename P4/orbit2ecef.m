function coordinate_ecef = orbit2ecef( M0, delta_n, e, sqrtA, OMEGA, i0, omega, OMEGA_DOT, ...
    i_dot, cuc, cus, crc, crs, cic, cis, toe,...
    year,month,day,hour,min,sec)
JD = juliand(year,month,day,hour,min,sec) ;
JD0 = juliand(1980,1,6,0,0,0) ;
gps_week = floor( (JD - JD0) / 7 );
gps_seconds=round(((((JD-2444244.5)/7)-gps_week)*24*60*60*7)/0.5)*0.5; %% GPS secend
t = gps_seconds - toe;
GM = 3.986004418e+14;
we = 7.29211514671*10^-5;
a = sqrtA^2;
M = M0 + (sqrt(GM/a^3) + delta_n)*t;
if or((-pi < M < 0) , (M > pi))
    E = M - e;
else
    E = M + e;
end
E_new=1000;
while abs(E_new - E) > 10^-15
    E_new = E + ((M - E + e .* sin(E))./(1 - e .* cos(E)));
    
    E = E_new;
end
v = atan2((sqrt(1-e^2)*sin(E)),(cos(E) - e));
u = omega +v+ cuc*cos(2*(omega + v)) + cus*sin(2*(omega + v));
r = a*(1-e*cos(E)) + crc*cos(2*(omega + v)) + crs*sin(2*(omega + v));
i = i0 + i_dot*t + cic*cos(2*(omega + v)) + cis*sin(2*(omega + v));
omg = OMEGA + (OMEGA_DOT - we)*t - we*toe;

v_pf = [r; 0; 0];
R_ecef =[(cos(omg)*cos(u))-(sin(omg)*sin(u)*cos(i)) , (-cos(omg)*sin(u))-(sin(omg)*cos(u)*cos(i)) , (sin(omg)*sin(i));
    (sin(omg)*cos(u))+(cos(omg)*sin(u)*cos(i)) , (-sin(omg)*sin(u))+(cos(omg)*cos(u)*cos(i)) , (-cos(omg)*sin(i));
    sin(u)*sin(i)  ,  cos(u)*sin(i)   ,  cos(i)];
coordinate_ecef = R_ecef*v_pf;




