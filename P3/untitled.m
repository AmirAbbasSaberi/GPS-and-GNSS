clc
clear
close all
format long g
% filename1 = 'igs19822.sp3';
% file1 = Extract_sp3_data(filename1);
filename = 'brdc0770.15n';
[prn1, M01, delta_n1, e1, sqrtA1, OMEGA1, i01, omega1, OMEGA_DOT1, ...
    i_dot1, cuc1, cus1, crc1, crs1, cic1, cis1, toe1,...
    year1,month1,day1,hour1,minute1,sec1,data] = Extract_Navigation_data (filename);
xyz = {};

for prn_s = 1:1
    for j = 0:3600*24/30-1
        i1 = find(prn1 == prn_s);
        if isempty(i1) == 0
            i1 = find(prn1 == prn_s);
            prn = prn1(i1);
            hour = hour1(i1);
            minute = minute1(i1);
            sec = sec1(i1);
            time = 3600*hour+60*minute+sec;
            t1 = 30*j;
            M0 = M01(i1);
            M0 = M0(1);
            delta_n = spline(time,10^9*delta_n1(i1),t1)/10^9;
            e = spline(time,10^5*e1(i1),t1)/10^5;
            sqrtA = spline(time,sqrtA1(i1),t1);
            OMEGA = spline(time,OMEGA1(i1),t1);
            i0 = spline(time,i01(i1),t1);
            omega = spline(time,omega1(i1),t1);
            OMEGA_DOT = spline(time,10^9*OMEGA_DOT1(i1),t1)/10^9;
            i_dot = spline(time,10^9*i_dot1(i1),t1)/10^9;
            cuc = spline(time,10^9*cuc1(i1),t1)/10^9;
            cus = spline(time,10^9*cus1(i1),t1)/10^9;
            crc = spline(time,10^9*crc1(i1),t1)/10^9;
            crs = spline(time,10^9*crs1(i1),t1)/10^9;
            cic = spline(time,10^9*cic1(i1),t1)/10^9;
            cis = spline(time,10^9*cis1(i1),t1)/10^9;
            toe = toe1(i1);
            
            
            i_t = 1;
            
            t = t1 - toe(i_t);
            if t<0
                t = t+24*60*60*7;
            end
            
            GM = 3.986004418e+14;
            we = 7.2921151467e-5;
            a = sqrtA(i_t)^2;
            M = M0(i_t) + (sqrt(GM/a^3) + delta_n(i_t))*t;
            if or(((-pi < M) & (M < 0)) , (M > pi))
                E = M - e(i_t);
            else
                E = M + e(i_t);
            end
            E_new=1000;
            while abs(E_new - E) > 10^-15
                E_new = E + ((M - E + e(i_t) .* sin(E))./(1 - e(i_t) .* cos(E)));
                
                E = E_new;
            end
            
            v = atan2((sqrt(1-e(i_t)^2)*sin(E)),(cos(E) - e(i_t)));
            
            u = omega(i_t) +v+ cuc(i_t)*cos(2*(omega(i_t) + v)) + cus(i_t)*sin(2*(omega(i_t) + v));
            
            r = a*(1-e(i_t)*cos(E)) + crc(i_t)*cos(2*(omega(i_t) + v)) + crs(i_t)*sin(2*(omega(i_t) + v));
            
            i = i0(i_t) + i_dot(i_t)*t + cic(i_t)*cos(2*(omega(i_t) + v)) + cis(i_t)*sin(2*(omega(i_t) + v));
            
            omg = OMEGA(i_t) + (OMEGA_DOT(i_t) - we)*t - we*toe(i_t);
            
            v_pf = [r; 0; 0];
            
            R_ecef =[(cos(omg)*cos(u))-(sin(omg)*sin(u)*cos(i)) , (-cos(omg)*sin(u))-(sin(omg)*cos(u)*cos(i)) , (sin(omg)*sin(i));
                (sin(omg)*cos(u))+(cos(omg)*sin(u)*cos(i)) , (-sin(omg)*sin(u))+(cos(omg)*cos(u)*cos(i)) , (-cos(omg)*sin(i));
                sin(u)*sin(i)  ,  cos(u)*sin(i)   ,  cos(i)];
            
            ecef = R_ecef*v_pf;
            x = ecef(1,1);
            y = ecef(2,1);
            z = ecef(3,1);
            coordinate_ecef_navigation = [prn(i_t),x,y,z,t1]; %% Navigation ECEF Coordinate
            
            R_eci= [cos(OMEGA(i_t)) sin(OMEGA(i_t)) 0;-sin(OMEGA(i_t)) cos(OMEGA(i_t)) 0;0 0 1]...
                *[1 0 0;0 cos(i0(i_t)) sin(i0(i_t));0 -sin(i0(i_t)) cos(i0(i_t))]...
                *[cos(omega(i_t)) sin(omega(i_t)) 0;-sin(omega(i_t)) cos(omega(i_t)) 0;0 0 1];
            
            eci = R_eci*[r*cos(v) ;r*sin(v); 0];
            x = eci(1,1);
            y = eci(2,1);
            z = eci(3,1);
            coordinate_eci_navigation = [prn(i_t),x,y,z,t1]; %% Navigation ECEF Coordinate
            x_eci(1+j,1:5) = coordinate_eci_navigation;
            x_ecef(1+j,1:5) = coordinate_ecef_navigation;
        else
            x_eci = [];
            x_ecef = [];
        end
    end
    xyz(prn_s,1:3) = {prn_s,x_eci,x_ecef};
end
figure(1)
for i = 1:1
    x_eci = xyz{i,2};
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
for i = 1:1
    x_ecef = xyz{i,3};
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
