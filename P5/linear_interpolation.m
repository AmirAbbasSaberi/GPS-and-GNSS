function [T,P,e,beta,lambda] = linear_interpolation(phi,D)

phi = abs(rad2deg(phi));

Lat = [15, 30, 45, 60, 75]';

T0 = [299.65, 294.15, 283.15, 272.15, 263.65]';
P0 = [1013.25, 1017.25, 1015.75, 1011.75, 1013.00]';
e0 = [26.31, 21.79, 11.66, 6.78, 4.11]';
beta0 = [6.3, 6.05, 5.58, 5.39, 4.53]'*10^(-3);
lambda0 = [2.77, 3.15, 2.57, 1.81, 1.55]';

dT = [0 , 7, 11, 15, 14.5]';
dP = [0 , -3.75, -2.25, -1.75, -0.5]';
de = [0 , 8.85, 7.24, 5.36, 3.39]';
dbeta = [0 , 0.25, 0.32, 0.81, 0.62]'*10^(-3);
dlambda = [0 , 0.33, 0.46, 0.74, 0.30]';


state = phi/15;
index = floor(phi/15);
if state-index==0 || index==0 || index==5 || index>5
    if index==0 
        state=1;
    elseif index==5 || index>5
        state=5;
    end
    intrp_T0 = T0(state);
    intrp_P0 = P0(state);
    intrp_e0 = e0(state);
    intrp_beta0 = beta0(state);
    intrp_lambda0 = lambda0(state);
    
    intrp_dT = dT(state);
    intrp_dP = dP(state);
    intrp_de = de(state);
    intrp_dbeta = dbeta(state);
    intrp_dlambda = dlambda(state);
    
else
    
   
    x = [Lat(index), Lat(index+1)];
    
    y_T0 = [T0(index), T0(index+1)];
    intrp_T0 = interp1(x,y_T0,phi);
    y_P0 = [P0(index), P0(index+1)];
    intrp_P0 = interp1(x,y_P0,phi);
    y_e0 = [e0(index), e0(index+1)];
    intrp_e0 = interp1(x,y_e0,phi);
    y_beta0 = [beta0(index), beta0(index+1)];
    intrp_beta0 = interp1(x,y_beta0,phi);
    y_lambda0 = [lambda0(index), lambda0(index+1)];
    intrp_lambda0 = interp1(x,y_lambda0,phi);
    
    y_dT = [dT(index), dT(index+1)];
    intrp_dT = interp1(x,y_dT,phi);
    y_dP = [dP(index), dP(index+1)];
    intrp_dP = interp1(x,y_dP,phi);
    y_de = [de(index), de(index+1)];
    intrp_de = interp1(x,y_de,phi);
    y_dbeta = [dbeta(index), dbeta(index+1)];
    intrp_dbeta = interp1(x,y_dbeta,phi);
    y_dlambda = [dlambda(index), dlambda(index+1)];
    intrp_dlambda = interp1(x,y_dlambda,phi);

end


if phi<0
    D_min = 211;            
else
    D_min =28;              
end

T = intrp_T0 - intrp_dT*cos(2*pi*(D-D_min)/365.25);
P = intrp_P0 - intrp_dP*cos(2*pi*(D-D_min)/365.25);
e = intrp_e0- intrp_de*cos(2*pi*(D-D_min)/365.25);
beta = intrp_beta0- intrp_dbeta*cos(2*pi*(D-D_min)/365.25);
lambda = intrp_lambda0- intrp_dlambda*cos(2*pi*(D-D_min)/365.25);

end