function sat_gps = gps_time_delay_coordiante(data,sat_gps,pos,c,we)
for t1 = 1:numel(data(:,1))/8
    
    prn1(t1,1) = data(-7+8*t1,1);
    crs1(t1,1) = data(-6+8*t1,2);
    delta_n1(t1,1) = data(-6+8*t1,3);
    M01(t1,1) = data(-6+8*t1,4);
    cuc1(t1,1) =  data(-5+8*t1,1);
    e1(t1,1) = data(-5+8*t1,2);
    cus1(t1,1) = data(-5+8*t1,3);
    sqrtA1(t1,1) = data(-5+8*t1,4);
    toe1(t1,1) = data(-4+8*t1,1);
    cic1(t1,1) = data(-4+8*t1,2);
    OMEGA1(t1,1) =data(-4+8*t1,3);
    cis1(t1,1) = data(-4+8*t1,4);
    i01(t1,1) = data(-3+8*t1,1);
    
    omega1(t1,1) = data(-3+8*t1,3);
    OMEGA_DOT1(t1,1) = data(-3+8*t1,4);
    i_dot1(t1,1) = data(-2+8*t1,1) ;
    crc1(t1,1) = data(-3+8*t1,2);
    
    year1(t1,1) = data(-7+8*t1,2) ;
    month1(t1,1) = data(-7+8*t1,3) ;
    day1(t1,1) = data(-7+8*t1,4) ;
    hour1(t1,1) = data(-7+8*t1,5) ;
    minute1(t1,1) = data(-7+8*t1,6) ;
    sec1(t1,1) = data(-7+8*t1,7) ;
    a01(t1,1) = data(-7+8*t1,8) ;
    a11(t1,1) =  data(-7+8*t1,9) ;
    a21(t1,1) =  data(-7+8*t1,10) ;
    TGD1(t1,1) = data(-1+8*t1,3);
end
[k1,k2] = size(sat_gps);
for k = 1 : size(sat_gps,1)
    prn = sat_gps(k,1);
    j = find(prn1 == prn);
    crs2 = crs1(j,1);
    delta_n2 = delta_n1(j,1);
    M02 = M01(j,1);
    cuc2 = cuc1(j,1);
    e2 = e1(j,1);
    cus2 = cus1(j,1);
    sqrtA2 = sqrtA1(j,1);
    toe2 = toe1(j,1);
    cic2 = cic1(j,1);
    OMEGA2 = OMEGA1(j,1);
    cis2 = cis1(j,1);
    i02 = i01(j,1);
    
    omega2 = omega1(j,1);
    OMEGA_DOT2 = OMEGA_DOT1(j,1);
    i_dot2 = i_dot1(j,1);
    crc2 = crc1(j,1);
    
    year2 = year1(j,1);
    month2 = month1(j,1);
    day2 = day1(j,1);
    hour2 = hour1(j,1);
    minute2 = minute1(j,1);
    sec2 = sec1(j,1);
    a02 = a01(j,1);
    a12 = a11(j,1);
    a22 = a21(j,1);
    TGD22 = TGD1(j,1);
    t1 = sat_gps(k,5)*3600 + sat_gps(k,6)*60 + sat_gps(k,7);
    
    
    time = hour2*3600+minute2*60+sec2;
    if time(end) == 0
        time(end) = 3600*24;
    end
    dt = -time+t1;
    i3 = find(dt<0);
    dt(i3) =  [];
    i2 = find(dt == min(dt));
    if isempty(i2)
        i2 = 1;
    end
    crs = crs2(i2,1);
    delta_n = delta_n2(i2,1);
    M0 = M02(i2,1);
    cuc = cuc2(i2,1);
    e = e2(i2,1);
    cus = cus2(i2,1);
    sqrtA = sqrtA2(i2,1);
    toe = toe2(i2,1);
    cic = cic2(i2,1);
    OMEGA = OMEGA2(i2,1);
    cis = cis2(i2,1);
    i0 = i02(i2,1);
    TGD = TGD22(i2,1);
    
    omega = omega2(i2,1);
    OMEGA_DOT = OMEGA_DOT2(i2,1);
    i_dot = i_dot2(i2,1);
    crc = crc2(i2,1);
    
    year = year2(i2,1);
    month = month2(i2,1);
    day = day2(i2,1);
    hour = hour2(i2,1);
    minute = minute2(i2,1);
    sec = sec2(i2,1);
    h = sat_gps(k,5);
    m = sat_gps(k,6);
    s = sat_gps(k,7);
    a0 = a02(i2,1);
    a1 = a12(i2,1);
    a2 = a22(i2,1);
    coordinate_ecef = orbit2ecef( M0, delta_n, e, sqrtA, OMEGA, i0, omega, OMEGA_DOT, ...
        i_dot, cuc, cus, crc, crs, cic, cis, toe,...
        year,month,day,h,m,s);
    N = 1;
    R = norm(-pos+coordinate_ecef);
    delta_t = R/c;
    delta_t_total = delta_t;
    while true
        N = N+1;
        coordinate_ecef1 = orbit2ecef( M0, delta_n, e, sqrtA, OMEGA, i0, omega, OMEGA_DOT, ...
            i_dot, cuc, cus, crc, crs, cic, cis, toe + delta_t,...
            year,month,day,sat_gps(k,5),sat_gps(k,6),s);
        R = norm(-pos+coordinate_ecef1);
        delta_t = R/c;
        new_R = norm(coordinate_ecef-coordinate_ecef1);
        if abs(new_R)<10^-6
            teta = -we*delta_t;
            coordinate_ecef1 = [cos(teta), sin(teta), 0; -sin(teta), cos(teta), 0; 0, 0, 1]*coordinate_ecef1;
            sat_gps(k,k2+1:k2+3) = coordinate_ecef1;
            sat_gps(k,k2+4) = h*3600+m*60+sec-delta_t;
            JD = juliand(year,month,day,h,m,s) ;
            JD0 = juliand(1980,1,6,0,0,0) ;
            gps_week = floor( (JD - JD0) / 7 );
            gps_seconds=round(((((JD-2444244.5)/7)-gps_week)*24*60*60*7)/0.5)*0.5; %% GPS secend
            t_gps_sec = gps_seconds-delta_t;
            sat_gps (k,k2+5) = a0+a1*(t_gps_sec-toe)+a2*(t_gps_sec-toe)^2;
            sat_gps (k,k2+6) = (a0+a1*(t_gps_sec-toe)+a2*(t_gps_sec-toe)^2)*c;
            sat_gps (k,k2+7) = -delta_t;
            sat_gps (k,k2+8) = sqrt(sum((coordinate_ecef1 - pos).^2));
            sat_gps (k,k2+9) = TGD;
            break
        end
        delta_t_total = delta_t_total+delta_t;
        coordinate_ecef = coordinate_ecef1;
    end
end