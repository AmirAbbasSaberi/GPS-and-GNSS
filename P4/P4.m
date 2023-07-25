clc
clear
close all
format long g
% git remote add origin git@github.com:Arsenic2000/GPS_matlab_master.git
% git branch -M main
% git push -u origin main
%%
file = fopen('mizu0770.15o');%mizu3260.16o  mizu0770.15o
sat = {};
while true
    
    line = fgetl(file);
    if contains(line,'APPROX POSITION XYZ')
        i = find(ismember(line,'APPROXPOSITIONXYZ'));
        j = min(i);
        pos = str2num(line(1:j-1));
        
    end
    if contains(line,'TYPES OF OBSERV')
        
        sat(1,1:16) = {'PRN','year','month','day','hour','minute','second',line(11:12),line(17:18),line(23:24),line(29:30) ...
            line(35:36) , line(41:42) ,line(47:48) , line(53:54) , line(59:60)};
        line = fgetl(file);
        sat(1,17:19) = {line(11:12),line(17:18),'Time'};
        break
    end
    
    
end
pos = pos';
while true
    
    line = fgetl(file);
    
    if contains(line,'END OF HEADER')
        break
    end
    
    
end

k = 1 ;
p = 2;
while feof(file) ~=1
    line = fgetl(file);
    year = str2double(line(2:3));
    month = str2double(line(5:6));
    day = str2double(line(8:9));
    hour = str2double(line(11:12));
    minute = str2double(line(14:15));
    sec = str2double(line(17:26));
    sat_num = str2double(line(31:32));
    t = hour*3600+minute*60+sec;
    
    
    
    if sat_num <= 12
        for i = 33:3:number
            number = length(line);
            sat(p,1) = {line(i:i+2)};
            p = p+1;
            
        end
    elseif and(sat_num<=24 , sat_num>12)
        number = 66;
        number1 = sat_num-12;
        for i = 33:3:number
            number = length(line);
            sat(p,1) = {line(i:i+2)};
            p = p+1;
            
        end
        line = fgetl(file);
        for i = 33:3:number1*3-3+33
            number = length(line);
            sat(p,1) = {line(i:i+2)};
            p = p+1;
            
        end
        s = 1;
    elseif and(sat_num>24 , sat_num<=36)
        number = 66;
        number1 = 12;
        number2 = sat_num-24;
        for i = 33:3:number
            number = length(line);
            sat(p,1) = {line(i:i+2)};
            p = p+1;
            
        end
        line = fgetl(file);
        for i = 33:3:number1*3-3+33
            number = length(line);
            sat(p,1) = {line(i:i+2)};
            p = p+1;
            
        end
        line = fgetl(file);
        for i = 33:3:number2*3-3+33
            number = length(line);
            sat(p,1) = {line(i:i+2)};
            p = p+1;
            
        end
        
    end
    
    for i = 1:sat_num
        line = fgetl(file);
        l = numel(line);
        if l == 0
            C1 = NaN;
            C2 = NaN;
            C5 = NaN;
            L1 = NaN;
            L2 = NaN;
        elseif and(l<=16,l>0)
            C1 = str2double(line(1:16));
            C2 = NaN;
            C5 = NaN;
            L1 = NaN;
            L2 = NaN;
        elseif and(l<=32,l>16)
            C1 = str2double(line(1:16));
            C2 = str2double(line(17:32));
            C5 = NaN;
            L1 = NaN;
            L2 = NaN;
        elseif and(l<=48,l>32)
            C1 = str2double(line(1:16));
            C2 = str2double(line(17:32));
            C5 = str2double(line(33:48));
            L1 = NaN;
            L2 = NaN;
        elseif and(l<=64,l>48)
            C1 = str2double(line(1:16));
            C2 = str2double(line(17:32));
            C5 = str2double(line(33:48));
            L1 = str2double(line(49:64));
            L2 = NaN;
        elseif and(l<=80,l>64)
            C1 = str2double(line(1:16));
            C2 = str2double(line(17:32));
            C5 = str2double(line(33:48));
            L1 = str2double(line(49:64));
            L2 = str2double(line(65:80));
        end
        line = fgetl(file);
        l = numel(line);
        if l == 0
            L5 = NaN;
            P1 = NaN;
            P2 = NaN;
            S1 = NaN;
            S2 = NaN;
        elseif and(l<=16,l>0)
            L5 = str2double(line(1:16));
            P1 = NaN;
            P2 = NaN;
            S1 = NaN;
            S2 = NaN;
        elseif and(l<=32,l>16)
            L5 = str2double(line(1:16));
            P1 = str2double(line(17:32));
            P2 = NaN;
            S1 = NaN;
            S2 = NaN;
        elseif and(l<=48,l>32)
            L5 = str2double(line(1:16));
            P1 = str2double(line(17:32));
            P2 = str2double(line(33:48));
            S1 = NaN;
            S2 = NaN;
        elseif and(l<=64,l>48)
            L5 = str2double(line(1:16));
            P1 = str2double(line(17:32));
            P2 = str2double(line(33:48));
            S1 = str2double(line(49:64));
            S2 = NaN;
        elseif and(l<=80,l>64)
            L5 = str2double(line(1:16));
            P1 = str2double(line(17:32));
            P2 = str2double(line(33:48));
            S1 = str2double(line(49:64));
            S2 = str2double(line(65:80));
        end
        line = fgetl(file);
        l = numel(line);
        if l == 0
            S5 = NaN;
        elseif and(l<=16,l>0)
            S5 = str2double(line(1:16));
        end
        
        sat(k+i,2:19) = {year,month,day,hour,minute,sec,C1,C2,C5,L1,L2,L5,P1,P2,S1,S2,S5,t};
    end
    k = size(sat,1);
end
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
        if abs(new_R)<10^-4
            teta = -we*delta_t;
            coordinate_ecef1 = [cos(teta), sin(teta), 0; -sin(teta), cos(teta), 0; 0, 0, 1]*coordinate_ecef1;
            sat_gps(k,k2+1:k2+3) = coordinate_ecef1;
            sat_gps(k,k2+4) = h*3600+m*60+sec-delta_t;
            JD = juliand(year,month,day,h,m,s) ;
            JD0 = juliand(1980,1,6,0,0,0) ;
            gps_week = floor( (JD - JD0) / 7 );
            gps_seconds=round(((((JD-2444244.5)/7)-gps_week)*24*60*60*7)/0.5)*0.5; %% GPS secend
            t_gps_sec = gps_seconds-delta_t;
            sat_gps(k,k2+5) = a0+a1*(t_gps_sec-toe)+a2*(t_gps_sec-toe)^2;
            sat_gps (k,k2+6) = (a0+a1*(t_gps_sec-toe)+a2*(t_gps_sec-toe)^2)*c;
            sat_gps (k,k2+7) = -delta_t;
            break
        end
        delta_t_total = delta_t_total+delta_t;
        coordinate_ecef = coordinate_ecef1;
    end
end

sat_gps (:,k2+8) = sat_gps (:,24)+sat_gps (:,25);

i = find(sat_gps(:,1) == 5);
deltaT_sat = sat_gps(i,end-2);
time = sat_gps(i,5)*3600 + sat_gps(i,6)*60 + sat_gps(i,7);

figure(1)
plot(time,deltaT_sat,'.')
title('Satellite clock dist offset c\delta_s_a_t PRN = GPS5')
xlabel('time_(_s_)')
ylabel('offset dist_(_m_)')
grid on
figure(2)
plot3(sat_gps(1:11,k2+1),sat_gps(1:11,k2+2),sat_gps(1:11,k2+3),'^','markerfacecolor','red')
hold on
for i=1:11
    
    plot3([pos(1) sat_gps(i,k2+1)],[pos(2) sat_gps(i,k2+2)],[pos(3) sat_gps(i,k2+3)],'--','markerfacecolor','red') 
    hold on
end

plot3(pos(1),pos(2),pos(3),'o','markerfacecolor','b','markersize' , 8)
grid on
title('Satellite Coordiante PRN = GPS5')
hold on 
[X,Y,Z] = sphere(20);
hold on
r = 6371000;
X = X * r;
Y = Y * r;
Z = Z * r;
surf(X,Y,Z)
axis equal
colormap summer
xlabel('X_E_C_E_F')
ylabel('Y_E_C_E_F')
zlabel('Z_E_C_E_F')
%%
figure(3)
for i = 1:32
    j = find(sat_gps(:,1) == i);
    if isempty(j)       
    else
        d = sqrt(sum((sat_gps(j,k2+1:k2+3)-[pos(1)*ones(numel(j),1) pos(2)*ones(numel(j),1) pos(3)*ones(numel(j),1)]).^2,2));
        time = sat_gps(j,5)*3600 + sat_gps(j,6)*60 + sat_gps(j,7);        
        plot(time,d,'.')
        hold on
    end    
end
title('Satellite Dist')
xlabel('time_(_s_)')
ylabel('Dist(_m_e_t_e_r_)')
grid on

figure(4) 
plot(sat_gps(1:11,20) - sat_gps(1:11,28),'.','markersize',8,'color','red')
xlabel('number of obs (1:11)')
ylabel('dt_(_s_)')
title('Diference time between geomtric and pseudorange method')
grid on






