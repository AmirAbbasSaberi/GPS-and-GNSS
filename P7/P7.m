clc
clear 
close all 


path = 'mizu0770.15o';
[pos,sat] = mizuRinexExtraction(path);



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

time = sat_gps(:,5)*3600 + sat_gps(:,6)*60 + sat_gps(:,7);
time = time/30;
i = find(time == 0);
sat_gps1 = sat_gps(i,:);


[data] = Extract_Navigation_data (filename);

dx1 = [10 10 10];
Cdt_rcv = 0;
while norm(dx1) > 10^-4

    sat_gps1 = sat_gps(i,:);
    sat_gps1 = gps_time_delay_coordiante(data,sat_gps1,pos,c,we);
    sat_gps1 (:,30) = sat_gps1 (:,24)+sat_gps1 (:,25);


    x = sat_gps1(:,21);
    y = sat_gps1(:,22);
    z = sat_gps1(:,23);

    wgs84 = wgs84Ellipsoid('kilometer');
    [lat,lon,h] = ecef2geodetic(wgs84,x,y,z);
    [I1,E,~] =  Klobuchar(sat_gps1,pos,filename,c);
    Tropsphere = Sustamiinen(pos,E);
    Ionosphere = I1*c; 
    TGD = sat_gps1(:,end-1);
    Ecludian_dist = sqrt((x-pos(1)).^2 + (y-pos(2)).^2 + (z-pos(3)).^2);
    Psudorange = sat_gps1(:,14);

    A = [(x-pos(1))./Ecludian_dist (y-pos(2))./Ecludian_dist (z-pos(3))./Ecludian_dist];
    A(:,end+1) = 1;

    Cdt_sat = sat_gps1(:,26);
    dl = Psudorange - Ecludian_dist + Cdt_sat + TGD*c - Ionosphere - Tropsphere - Cdt_rcv;

    dx = (A'*A)^-1*A'*dl;
    pos = pos - dx(1:3);
    dx1 = dx(1:3);
    Cdt_rcv = Cdt_rcv + dx(4);

end

Qx__pos_tartibi_without = 0.04*(A'*A)^-1;

dX = dx;
for i = 1:max(time)

    j = find(time == i);
    sat_gps1 = sat_gps(j,:);
    
    sat_gps1 = gps_time_delay_coordiante(data,sat_gps1,pos,c,we);
    sat_gps1 (:,30) = sat_gps1 (:,24)+sat_gps1 (:,25);


    x = sat_gps1(:,21);
    y = sat_gps1(:,22);
    z = sat_gps1(:,23);

    wgs84 = wgs84Ellipsoid('kilometer');
    [lat,lon,h] = ecef2geodetic(wgs84,x,y,z);
    [I1,E,~] =  Klobuchar(sat_gps1,pos,filename,c);
    Tropsphere = Sustamiinen(pos,E);
    Ionosphere = I1*c; 
    TGD = sat_gps1(:,end-1);
    Ecludian_dist = sqrt((x-pos(1)).^2 + (y-pos(2)).^2 + (z-pos(3)).^2);
    Psudorange = sat_gps1(:,14);


    Ai = [-(x-pos(1))./Ecludian_dist -(y-pos(2))./Ecludian_dist -(z-pos(3))./Ecludian_dist];
    Ai(:,end+1) = 1;


    Cdt_sat = sat_gps1(:,26);
    dli =  Psudorange - Ecludian_dist + Cdt_sat + TGD*c - Ionosphere - Tropsphere - Cdt_rcv;
    
    k = find(abs(dli) > 100);
    dli(k) = [];
    Ai(k,:) = [];

    Mi = Qx__pos_tartibi_without*Ai';
    Pi = (0.04*eye(size(Ai,1))+Ai*Mi)^-1;
    Si = dli - (Ai*dX);
    pso_1 = pos;
    dx = Mi*Pi*Si;
    
    pos = pos + dx(1:3);
    trend_tartibi_without(:,i) = pos;
    dds(:,i) = dx;
    dX = dx;
    Cdt_rcv = Cdt_rcv + dX(4);
    Qx__pos_tartibi_without = Qx__pos_tartibi_without - Mi*Pi*Mi';

end


figure
plot(trend_tartibi_without(1,:)-trend_tartibi_without(1,end),'r.')
xlabel('Step')
ylabel('x_(_m_)')
hold on 
n = numel(trend_tartibi_without(1,:));
x = ones(n,1);
plot(0.3*x,'b')
xlabel('Step')
ylabel('x_(_m_)')
hold on
plot(0.1*x,'k')
xlabel('Step')
ylabel('x_(_m_)')

figure
plot(trend_tartibi_without(2,:)-trend_tartibi_without(2,end),'g.')
xlabel('Step')
ylabel('y_(_m_)')

figure
plot(trend_tartibi_without(3,:)-trend_tartibi_without(3,end),'b.')
xlabel('Step')
ylabel('z_(_m_)')