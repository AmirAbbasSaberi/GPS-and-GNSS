clear
clc
close all
format long g
%%
file = fopen('mtka3260.16o');%mtka3260.16o mtka0770.15o
sat = {};
sat(1,1:15) = {'PRN','year','month','day','hour','minute','second','L1','L2','P1','P2','C1','S1','S2','Time'};

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
        number = length(line);
        for i = 33:3:number
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
            C1 = str2double(line(1:15));
            C2 = NaN;
            C5 = NaN;
            L1 = NaN;
            L2 = NaN;
        elseif and(l<=32,l>16)
            C1 = str2double(line(1:15));
            C2 = str2double(line(17:30));
            C5 = NaN;
            L1 = NaN;
            L2 = NaN;
        elseif and(l<=48,l>32)
            C1 = str2double(line(1:15));
            C2 = str2double(line(17:30));
            C5 = str2double(line(33:46));
            L1 = NaN;
            L2 = NaN;
        elseif and(l<=62,l>48)
            C1 = str2double(line(1:15));
            C2 = str2double(line(17:30));
            C5 = str2double(line(33:48));
            L1 = str2double(line(49:62));
            L2 = NaN;
        elseif and(l<=78,l>62)
            C1 = str2double(line(1:15));
            C2 = str2double(line(17:30));
            C5 = str2double(line(33:46));
            L1 = str2double(line(49:64));
            L2 = str2double(line(65:78));
        end
        line = fgetl(file);
        l = numel(line);
        
        
        if l == 0
            S1 = NaN;
            S2 = NaN;
        elseif and(l<=16,l>0)
            S1 = str2double(line(1:14));
            S2 = NaN;
        elseif and(l<=31,l>17)
            S1 = str2double(line(1:14));
            S2 = str2double(line(17:30));
        end
        
        sat(k+i,2:15) = {year,month,day,hour,minute,sec,C1,C2,C5,L1,L2,S1,S2,t};
    end
    k = size(sat,1);
end
clear P1 P2 L1 L2 L5 n p  k i hour minute sec month line s S1 S2 S5 year t C1 C2 C5 number number1 number2 sat_num day end_of_header


close all


%% Data
index = find(strcmpi(char(sat(:,1)),"G28") );

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

close all
Li = L1-L2;
figure(1)
plot(time,Li,'.')

% this code works just for specific PRN obsrevation
% this algorithm splits all arcline of Rinex obs for preparing for finding cycle-slip time 
% the most important point is the Obs being free of NaN data so Plz eleminate them first
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

% the "for loop" is interping gaps in arclines of Rinex and illustrating them in figure
for j = 1: numel(spli)
    t = spli{j,1};
    i = find(ismember(time , t));
    l = Li(i);
    tt = min(t):30:max(t);
    li = transpose(spline(t,l,tt));
    spli(j,2:9) = {tt,li,t,l,L1(i),L2(i),P1(i),P2(i)};
    figure()
    plot(tt,li)
    hold on
    plot(t,l,'ro')
end



% Sec-order polynomial algorithm for finding cycle-slip
clear cycle_slip_poly2 
for i = 1:size(spli,1)
    li = spli{i,3}; %Ln
    tt = spli{i,2}; %Time
    p =1;
    n = 1;

    while true
        
        while true
            x = tt(n:n+2);
            y = li(n:n+2);
            an = polyfit(x,y,2);
            P = polyval(an,tt(n+3));
            if abs(li(n+3)-P)> 3.5*std(y)
                cycle_slip_poly2(p,1:3) = [tt(n),tt(n+3),i];
                p = p+1;
                n = n+3;
                break
            elseif  n >= numel(li)-4
                n = n+3;
                break
            end
            n = n+1;
        end
        if n >= numel(li)-4
            break
        end
    end
end
i = cycle_slip_poly2(1,3);
t = spli{i,2};
j = find(spli{i,2} == cycle_slip_poly2(1,1));
li = spli{i,3};
figure
plot(spli{i,2},spli{i,3},'.','markersize',3)
hold on
plot(t(j),li(j),'*','markersize')




clear cycle_slip_spline 
p = 1;
for j = 1:numel(spli(:,1))
x = spli{j,5}; %Li
t = spli{j,4}/30-min(spli{j,4}/30)+1; %Time
k = 5;
while true
    for i = k:numel(x)-6
        y1 = x(k:i+4);
        x1 = t(k:i+4);
        ii = i;
        P = transpose(spline(x1,y1,i+5));
        if abs(x(i+5)-P) > 6*std(x(i-3:i+3))
            cycle_slip_spline(p,1:3) = [x(i+5),j,i];
            p = p+1;
            k = i+5;
            break           
        end        
    end
    if ii >= numel(x)-10
        break
    end
        
        
end

end

clear cycle_slip_BW 
p =1;
for j = 1:size(spli,1)
    L1 = spli{j,6};
    L2 = spli{j,7};
    
    P1 = spli{j,8};
    P2 = spli{j,9};
    LW = (f1*L1-f2*L2)/(f1-f2);
    PN = (f1*P1+f2*P2)/(f1+f2);
    BW = LW-PN;
    landa_W=299792458/(f1-f2);
    m_BW = BW(1);
    s_BW = landa_W/10;

    
    n =1;

    while true
        
        while true
            if n >= numel(BW)-25
                break
            end
            m_BW  = mean(BW(n:n+7));
            s_BW = std(BW(n:n+7));
            if abs(BW(n+8)-m_BW)> 4*s_BW  
                cycle_slip_BW(p,1:2) = [n+7,j];
                p = p+1;
                n = n+8;
                break
            end
            n = n+1;
        end
        if n >= numel(BW)-25
            break
        end
    end
    
end




