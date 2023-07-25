function file1 = Extract_sp3_data(filename1)
%______________________________________________________________________________________________________%
%%%  this function perduced by AmirAbbas Saberi %%%
%----- input = file sp3  from  <<http://www.gnsscalendar.com/>> 
%----- part of GPS Final Orbits (IGS) 
%----- output [X Y Z sigmax sigmay sigmaz .....] n*p marix 
%----- output CalenderTime ----------------------------------
%______________________________________________________________________________________________________%
fid = fopen(filename1);
if fid == -1
    errordlg(['The file ''' filename1 ''' does not exist.']);
    return;
end
end_of_header = 0;
while end_of_header == 0
    current_line = fgetl(fid);
    if contains(current_line,'PCV') %% find PVC in file 
        end_of_header=1;
    end
end
j = 0;
p = cell(1,1);
while feof(fid) ~= 1
    j = j+1;
    l = fgetl(fid);
    p{j,1} = char(l(3:end)); %% clear the header and strings in file
end

fileID = fopen('celldata1.dat','w'); 
[nrows,~] = size(p);
for row = 1:nrows-1
fprintf(fileID,'%s\n' ,p{row,:});
end
fclose(fileID);
clear ;
file1 = dlmread('celldata1.dat');  %% save file in dat file 
end