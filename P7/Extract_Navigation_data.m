function[data] = Extract_Navigation_data (filename)
%______________________________________________________________________________________________________%
%%%  this function perduced by AmirAbbas Saberi %%%
%----- input = file navgation  from  <<http://www.gnsscalendar.com/>>
%----- par nav data
%-----|   PRN   |          satellite PRN         |-----%%%-----|Cuc| cosine term, arg. of latitude
%-----|    M0   | mean anomaly at reference time |-----%%%-----|Cus|sine term, arg. of latitude
%-----| delta_n |     mean motion difference     |-----%%%-----|Crc|cosine term, radius
%-----|    e    |          eccentricity          |-----%%%-----|Crs|sine term, radius
%-----|  sqrtA  |   where A is semimajor axis    |-----%%%-----|Cic|cosine term, inclination
%-----|  OMEGA  |      LoAN at weekly epoch      |-----%%%-----|toe|
%-----|OMEGA_dot|    rate of right ascension     |-----%%%-----|IODE|
%-----|  i_dot  |   rate of inclination angle    |-----%%%-----|GPS_wk|
%-----|   toe   |       time of ephemeris        |-----%%%-----|af0| ....
%-----|   af2   |              .....             |-----%%%-----|TGD|......
%-----|    Y    |               year             |-----%%%-----|M|month
%-----|    D    |              day               |-----%%%-----|H|hour
%-----|   sec   |            second              |-----%%%-----|omega  |argument of perigee
%______________________________________________________________________________________________________%
fid = fopen(filename);
if fid == -1
    errordlg(['The file ''' filename ''' does not exist.']);
    return;
end
end_of_header = 0;
while end_of_header == 0
    current_line = fgetl(fid);
    if contains(current_line,'END OF HEADER')
        end_of_header=1;
    end
end
j = 0;
p = cell(1,1);
while feof(fid) ~= 1
    j = j+1;
    l = fgetl(fid);
    p{j,1} = char(l);
end
fileID = fopen('celldata.dat','w');
[nrows,~] = size(p);
for row = 1:nrows
    fprintf(fileID,'%s\n' ,p{row,:});
end
fclose(fileID);
clear ;

data = dlmread('celldata.dat');
end