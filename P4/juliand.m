function J_Date = juliand(Y,M,D,H,min,sec) 
%______________________________________________________________________________________________________%
%%%  this function perduced by AmirAbbas Saberi %%%
%-----  year Y , month M , Day D , minute min , second sec 
%----- output J_Date time 
%______________________________________________________________________________________________________%
const = 365.25;
if and(Y >= 80 , Y <= 99)
    Y = 1900 + Y;
end
if and(Y >= 0 , Y <= 79)
    Y = 2000 + Y;
end
if (M <= 2)
    y = Y - 1;
    m = M + 12;
end
if (M > 2)
    y = Y;
    m = M;
end
J_Date = fix( ( const* y) ) + floor( (30.6001 * (m+1)) ) + D + ( (H + min / 60 + sec / 3600) / 24 ) + 1720981.5;
end