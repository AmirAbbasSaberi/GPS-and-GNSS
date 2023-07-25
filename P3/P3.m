clc
clear
close all

% filename1 = 'igs19822.sp3';
% file1 = Extract_sp3_data(filename1);
tic
filename = 'brdc0770.15n';
[data] = Extract_Navigation_data (filename);


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
end

sat_num = unique(prn1);

for in = 1:numel(sat_num)
    prn = sat_num(in);
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
    for ii = 0:24*3600/30-1
        t1 = ii*30;
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
        
        
        JD = juliand(year,month,day,hour,minute,sec) ;
        JD0 = juliand(1980,1,6,0,0,0) ;
        gps_week = floor( (JD - JD0) / 7 );

        gps_seconds=round(((((JD-2444244.5)/7)-gps_week)*24*60*60*7)/0.5)*0.5; %% GPS secend

        t = gps_seconds + t1 - toe - (hour*3600+minute*60+sec);
        
        if t<0
            t = t+24*60*60*7;
        end

        GM = 3.986004418e+14;
        we = 7.2921151467e-5;
        a = sqrtA^2;


        M = M0 + (sqrt(GM/a^3) + delta_n)*t;

        if or((-pi < M )&&(M< 0) , (M > pi))
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
        
        ecef = R_ecef*v_pf;
        x = ecef(1,1);           
        y = ecef(2,1);
        z = ecef(3,1);
        coordinate_ecef_navigation(ii+1,1:11) = [prn,x,y,z,t1,year,month,day,hour,minute,sec]; %% Navigation ECEF Coordinate 

        R_eci= [cos(OMEGA) sin(OMEGA) 0;-sin(OMEGA) cos(OMEGA) 0;0 0 1]...
                  *[1 0 0;0 cos(i0) sin(i0);0 -sin(i0) cos(i0)]...
                  *[cos(omega) sin(omega) 0;-sin(omega) cos(omega) 0;0 0 1];
        
        eci = R_eci*[r*cos(v) ;r*sin(v); 0];
        x = eci(1,1);           
        y = eci(2,1);
        z = eci(3,1); 
        coordinate_eci_navigation(ii+1,1:11) = [prn,x,y,z,year,t1,month,day,hour,minute,sec]; %% Navigation ECEF Coordinate 

        
        
    end
    clear d
    xyz(sat_num(in),1:2) = {coordinate_ecef_navigation,coordinate_eci_navigation};
end

filename1 = 'igs18363.sp3';
file1 = Extract_sp3_data(filename1);
date1 = zeros(1,6);
 for i = 1:numel(file1(:,1))/33
     date1(i,:) = file1(32*(i-1)+1,1:6);
     coordinate_ecef_sp3(31*i-30:31*i,1:4) = file1(32*i-30:32*i,1:4);
 end
 for i = 1:numel(file1(:,1))/33
     d = date1(i,:);
     
     coordinate_ecef_sp3(31*i-30:31*i,5:10)= [d(1)*ones(31,1) d(2)*ones(31,1) d(3)*ones(31,1) d(4)*ones(31,1) d(5)*ones(31,1) d(6)*ones(31,1)];
 end


[coordinate_eci_navigation1,coordinate_ecef_navigation1 ]= orbit2ecef(data);
sat_num = unique(coordinate_ecef_navigation1(:,1));

figure(1)
for i = 1:numel(sat_num)
    j = sat_num(i);
    x_eci = xyz{j,2};
    if isempty(x_eci) == 0
    plot3(x_eci(:,2),x_eci(:,3),x_eci(:,4))
    axis equal
    hold on
    end
end
title('ECI Coordiante of all GPS Sats')
xlabel('x_(_m_)')
ylabel('y_(_m_)')
zlabel('z_(_m_)')
figure(2)
for i = 1:numel(sat_num)
    j = sat_num(i);
    x_ecef = xyz{j,1};
    if isempty(x_ecef) == 0
    plot3(x_ecef(:,2),x_ecef(:,3),x_ecef(:,4))
    axis equal
    hold on
    end
end
title('ECEF Coordiante of all GPS Sats')
xlabel('x_(_m_)')
ylabel('y_(_m_)')
xlabel('z_(_m_)')

coordinate_ecef_sp3(:,11) = coordinate_ecef_sp3(:,8)*3600 + coordinate_ecef_sp3(:,9)*60 +coordinate_ecef_sp3(:,10);
i = find(coordinate_ecef_sp3(:,1) == 12);
x_sp3_1 = coordinate_ecef_sp3(i,:);
x_ecef1 = xyz{12,1};
[~,k] = intersect(x_ecef1(:,5),x_sp3_1(:,11));
x_ecef1 = x_ecef1(k,:);
dx = x_ecef1(:,2:4) - x_sp3_1(:,2:4)*1000;
norm_dx = sqrt(dx(:,1).^2 + dx(:,2).^2 + dx(:,3).^2);
time = x_sp3_1(:,11);

figure(3)
plot(time,norm_dx,'linewidth',3)
title('Difference between sp3 and calculated ECEF by navigation GPS12' )
xlabel('Time_(_s_)')
ylabel('dx_(_m_)')



