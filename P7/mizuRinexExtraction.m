function [pos,sat] = mizuRinexExtraction(path)
file = fopen(path);%mizu3260.16o  mizu0770.15o
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
end