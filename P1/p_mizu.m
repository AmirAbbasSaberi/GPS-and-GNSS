clear
clc
close all
format long g
%%
file = fopen('mizu0770.15o');%mizu3260.16o  mizu0770.15o
sat = {};
while true
    
    line = fgetl(file);
    
    if contains(line,'TYPES OF OBSERV')
            
        sat(1,1:16) = {'PRN','year','month','day','hour','minute','second',line(11:12),line(17:18),line(23:24),line(29:30) ...
            line(35:36) , line(41:42) ,line(47:48) , line(53:54) , line(59:60)};
        line = fgetl(file);
        sat(1,17:19) = {line(11:12),line(17:18),'Time'};
        break
    end
    
    
end

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


%% Data
index = find(strcmpi(char(sat(:,1)),"G18") ); 

f1 = 1575.42; 
f2 = 1227.6;
iL1 = find(strcmpi(sat(1,:) , 'L1'));
iL2 = find(strcmpi(sat(1,:) , 'L2'));
iP1 = find(strcmpi(sat(1,:) , 'P1'));
iP2 = find(strcmpi(sat(1,:) , 'P2'));
itime= find(strcmpi(sat(1,:) , 'Time'));
time = cell2mat(sat(index,itime));
L1 = cell2mat(sat(index,iL1));
L2 = cell2mat(sat(index,iL2));
P1 = cell2mat(sat(index,iP1));
P2 = cell2mat(sat(index,iP2));
il1nan = find(isnan(L1));
ip1nan = find(isnan(P1));
il2nan = find(isnan(L2));
ip2nan = find(isnan(P2));
L1([il1nan;ip1nan;il2nan;ip2nan]) = [];
L2([il1nan;ip1nan;il2nan;ip2nan]) = [];
P1([il1nan;ip1nan;il2nan;ip2nan]) = [];
P2([il1nan;ip1nan;il2nan;ip2nan]) = [];
time([il1nan;ip1nan;il2nan;ip2nan]) = [];
%% Normal Hatch
N = 600;

R1 = P1(1);
for k = 2: numel(P1)
    if k<N
        n = k;
    else 
        n = N;
    end
    R1(k,1) = (1/n)*P1(k) + ((n-1)/n)*(R1(k-1)+L1(k)-L1(k-1));
    
end

R2 = P2(1);
for k = 2: numel(P2)
    if k<N
        n = k;
    else 
        n = N;
    end
    R2(k,1) = (1/n)*P2(k) + ((n-1)/n)*(R2(k-1)+L2(k)-L2(k-1));
    
end


figure(1) 
subplot(2,1,1)
hold on
plot(time,R1,'.')
plot(time,L1,'.')
xlabel('times')
ylabel('Meters')
legend('R_1','L_1','Location','best')
title('MIZU Normal Hatch Algorithm first frequency')
hold off
grid on
subplot(2,1,2)
hold on
plot(time,R1,'.')
plot(time,P1,'.')
xlabel('times')
ylabel('Meters')
legend('R_1','P_1','Location','best')
hold off
grid on

figure(2) 
subplot(2,1,1)
hold on
plot(time,R2,'.')
plot(time,L2,'.')
xlabel('times')
ylabel('Meters')
legend('R_2','L_2','Location','best')
title('MIZU Normal Hatch Algorithm second frequency')
hold off
grid on
subplot(2,1,2)
hold on
plot(time,R2,'.')
plot(time,P2,'.')
xlabel('times')
ylabel('Meters')
legend('R_2','P_2','Location','best')
hold off
grid on


%% Hatch Ionosphere Free

Lc=(f1^2*L1-f2^2*L2)/(f1^2-f2^2);
Pc=(f1^2*P1-f2^2*P2)/(f1^2-f2^2);

Rc = Pc(1);
for k = 2: numel(Pc)
    if k<N
        n = k;
    else 
        n = N;
    end
    Rc(k,1) = (1/n)*Pc(k) + ((n-1)/n)*(Rc(k-1)+Lc(k)-Lc(k-1));
    
end


figure(3) 
subplot(2,1,1)
hold on
plot(time,Rc,'.')
plot(time,Lc,'.')
xlabel('times')
ylabel('Meters')
legend('R_c','L_c','Location','best')
title('Hatch Ion-free Algorithm frequency')
hold off
grid on
subplot(2,1,2)
hold on
plot(time,Rc,'.')
plot(time,Pc,'.')
xlabel('times')
ylabel('Meters')
legend('R_c','P_c','Location','best')
hold off
grid on


%% Hatch Divergence Free

alfa = ((40.3*(f1^2-f2^2))/(f1*f2))*10^16;
Ldf1 = L1+2*alfa*(L1-L2);
Ldf2 = L2+2*alfa*(L2-L1);

Rdf1 = P1(1);
for k = 2: numel(Pc)
    if k<N
        n = k;
    else 
        n = N;
    end
    Rdf1(k,1) = (1/n)*P1(k) + ((n-1)/n)*(Rdf1(k-1)+Ldf1(k)-Ldf1(k-1));
    
end

Rdf2 = P2(1);
for k = 2: numel(Pc)
    if k<N
        n = k;
    else 
        n = N;
    end
    Rdf2(k,1) = (1/n)*P2(k) + ((n-1)/n)*(Rdf2(k-1)+Ldf2(k)-Ldf2(k-1));
    
end


figure(4) 
subplot(2,1,1)
hold on
plot(time,Rdf1,'.')
plot(time,Ldf1,'.')
xlabel('times')
ylabel('Meters')
legend('Rdf_1','Ldf_1','Location','best')
title('Hatch Div-free Algorithm first frequency')
hold off
grid on
subplot(2,1,2)
hold on
plot(time,Rdf1,'.')
plot(time,P1,'.')
xlabel('times')
ylabel('Meters')
legend('Rdf_1','P_1','Location','best')
hold off
grid on


figure(5) 
subplot(2,1,1)
hold on
plot(time,Rdf2,'.')
plot(time,Ldf2,'.')
xlabel('times')
ylabel('Meters')
legend('Rdf_2','Ldf_2','Location','best')
title('Hatch Div-free Algorithm second frequency')
hold off
grid on
subplot(2,1,2)
hold on
plot(time,Rdf1,'.')
plot(time,P1,'.')
xlabel('times')
ylabel('Meters')
legend('Rdf_2','P_2','Location','best')
hold off
grid on

close all
Li = L1-L2;
figure(1)
plot(time,Li,'.')
i1 = find (time == 6300);
i2 = find (time == 6360);
k1 = Li(i2-1) - Li(i1)
k2 = Li(i2) - Li(i1)
n = 1;
p = 1;
spli = {};
p = 1;
k2 = 1;
k = 1;
while true
    
    for i = k : numel(Li)-1
        d = time(i+1) - time(i);
        k1 = i;
        if d > 15000
            spli(p,1) = {time(k2:k1)};
            k = k1+1;
            k2 = k1+1;
            p = p+1;
            break
        elseif k1 == numel(Li)-1
            spli(p,1) = {time(k2:k1)};
            p = p+1;
            break
        end
        
    end
    if k1 == numel(Li)-1
        break
    end
    
    
end


for j = 1: numel(spli)
    t = spli{j,1};
    i = find(ismember(time , t));
    l = Li(i);
    tt = min(t):30:max(t);
    li = transpose(spline(t,l,tt));
    spli(j,2:3) = {tt,li};
    figure()
    plot(tt,li)
    hold on
    plot(t,l,'ro')
end


clear h
for i = 1:size(spli,1)
    li = spli{i,3};
    tt = spli{i,2};
    p =1;
    n = 1;

    while true
        
        while true
            x = tt(n:n+2);
            y = li(n:n+2);
            an = polyfit(x,y,2);
            P = polyval(an,tt(n+3));
            %         T0 = 10;
            %         a0 = 5.50033939790375;
            %         deltaT = tt(n+3)-tt(n+2);
            %         treshold = a0 - a0*exp(-deltaT/T0)/2;
            if abs(li(n+3)-P)> 400
                n = n+3;
                h(p,1:3) = [tt(n),tt(n+3),i];
                p = p+1;
                break
            elseif  n == numel(li)-4
                n = n+3;
                break
                %         elseif time(n) == 5700
                %             break
            end
            n = n+1;
        end
        if n >= numel(li)-4
            break
            %     elseif time(n) == 5700
            %             break
        end
    end
end

point = mean2(h(:,1:2));

i = find(abs(time-point) == min(abs(time-point)));
near_time = time(max(i));










